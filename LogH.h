#ifndef __LOGH_H__
#define __LOGH_H__

#include <sstream>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>

/* This is LogH only */

/* use this one if run on x86_64 */
#if 1
#define _SSE_
#endif

#if 0
#ifdef _SSE_
#error "_SSE_ is not supported with _SYMPCORR_"
#endif
#define _SYMPCORR_
#endif

/* use this option to compute forces in N^2/s steps rather than N^2 steps */
#if 1
#define _N2FAST_
#endif

/* use this option to use rational functoin interpolation instead of polynomial */
#if 1
#define _RATIONAL_FUNCTIONS_
#endif

#include "particle.h"
#include "mytimer.h"

#define SQR(x) ((x)*(x))

#define SMALLM 1.0e-64

struct Force
{
  typedef std::vector<Force, __gnu_cxx::malloc_allocator<Force, 64> > Vector;
  vec3 acc, iPad;
  Force() {}
  Force(const vec3 &_acc) : acc(_acc) {}
};

#ifdef _SSE_
struct ForceSIMD
{
  typedef std::vector<ForceSIMD, __gnu_cxx::malloc_allocator<ForceSIMD, 64> > Vector;
  v2df accx, accy, accz, iPad;
  ForceSIMD() {}
  ForceSIMD(const real v) 
  {
    accx = accy = accz = iPad = (v2df){v, v};
  }
};
#endif

struct Nbody
{
  unsigned long long iteration;
  real time, h;

  int ntry, ntry1;

  Particle::Vector ptcl;
  Force   ::Vector force;

#ifdef _SSE_
  ParticleSIMD::Vector  ptclSIMD;
  ForceSIMD   ::Vector forceSIMD;
#endif

  unsigned long long flops;
  double tbeg;
  double dt_force;
  double dt_step, dt_multi, dt_extra, dt_err;

  real dt;
  real negE0;

  void reset_counters()
  {
    ntry = ntry1 = 0;
    flops = 0.0;
    dt_force = 0.0;
    dt_step = dt_multi = dt_extra = dt_err = 0.0;
    tbeg = mytimer::get_wtime();
  }

  real get_gflops() const
  {
    return flops/(mytimer::get_wtime() - tbeg)/1e9;
  }

  real get_gflops_force() const
  {
    assert(dt_force > 0.0);
    return flops / dt_force / 1e9;
  }


  Nbody(const unsigned long long i, const real tepoch, const real _h, const Particle::Vector &_ptcl) : 
    iteration(i), time(tepoch), h(_h), ptcl(_ptcl) 
  {
    reset_counters();
    force.resize(ptcl.size());
    vec3 cm_pos(0.0), cm_vel(0.0);
    real Mtot = 0.0;
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      it->mass = std::max(it->mass, SMALLM);
      Mtot   += it->mass;
      cm_pos += it->mass*it->pos;
      cm_vel += it->momentum();
    }

    cm_pos *= 1.0/Mtot;
    cm_vel *= 1.0/Mtot;

    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      it->pos -= cm_pos;
      it->vel -= cm_vel;
    }

#ifdef _SSE_
    assert((ptcl.size() & 1) == 0);
    ptclSIMD.resize(ptcl.size());
    forceSIMD.resize(ptclSIMD.size());
#endif

    negE0 = -Etot();
  }

  real Ekin(const Particle::Vector &ptcl) const
  {
    real Ekin = 0.0;
    vec3 cmom(0.0);
    for (Particle::Vector::const_iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      Ekin += it->Ekin();
      cmom += it->momentum();
    }
    return Ekin;
  }

  real Epot() const { return Epot(ptcl); }
  real Ekin() const { return Ekin(ptcl); }
  real Epot(const Particle::Vector &ptcl) const
  {
    const int nbody = ptcl.size();
    real gpot = 0.0;

    for (int i = 0; i < nbody-1; i++)
      for (int j = i+1; j < nbody; j++)
      {
        const real ds2 = (ptcl[i].pos - ptcl[j].pos).norm2();
        assert(ds2 > 0.0);
        gpot -= ptcl[i].mass*ptcl[j].mass/std::sqrt(ds2);
      }

    return gpot;
  }

  real Etot() const {return Ekin(ptcl) + Epot(ptcl);}
  real E0  () const {return -negE0;}

  std::string print_orbit(const int i) const
  {
    const Particle &p0 = ptcl[0];
    const Particle &p1 = ptcl[i];
    const vec3 R = p1.pos - p0.pos;
    const vec3 V = p1.vel - p0.vel;
    const real Mtot = p0.mass + p1.mass;

    const vec3 L = R%V;
    const real h2 = L.norm2();
    const real h  = std::sqrt(h2); 

    const real inc = std::acos(L.z/h);
    const real fac = std::sqrt(L.x*L.x + L.y*L.y)/h;
    const real LSMALL = 1.0e-10;


    double capom, u;
    if (fac < LSMALL)
    {
      capom = 0.0;
      u     = std::atan2(R.y, R.x);
      if (std::abs(inc - M_PI) < 10*LSMALL) 
        u = -u;
    } 
    else
    {
      capom = std::atan2(L.x, -L.y);
      u = std::atan2( R.z/std::sin(inc) , R.x*std::cos(capom) + R.y*std::sin(capom));
    }

    if (capom < 0.0) capom += 2.0*M_PI;
    if (u     < 0.0) u     += 2.0*M_PI;

    const double r  = R.abs();
    const double v2 = V.norm2();
    const double vdotr = R*V;
    const double  energy = 0.5*v2 -  Mtot/r; 
    assert(energy < 0.0);

    double e, f, a, omega, capm;
    {
      
      a = -0.5*Mtot/energy;
      const double fac = 1.0 - h2/(Mtot*a);
      assert(a > 0.0);

    
      double cape; 
      if (fac > LSMALL)
      {
        e = std::sqrt(fac);
        const double face = (a-r)/(a*e);
        cape = face > 0.9999999 ? 0.0 : (face > -0.9999999 ? std::acos(face) : M_PI);

        if (vdotr < 0.0) cape = 2.0*M_PI - cape;

        const double cf = (cos(cape) - e)/(1.0 - e*cos(cape));
        const double sf = std::sqrt(1-e*e)*std::sin(cape)/(1.0 - e*std::cos(cape));
        f = std::atan2(sf, cf);
        if (f < 0.0) f += 2.0*M_PI;
      }
      else
      {
        e = 0.0;
        f = u;
        cape = u;
      }

      capm = cape - e*std::sin(cape);
      omega = u - f;
      if (omega < 0) omega += 2.0*M_PI;
      omega = omega - int(omega/(2.0*M_PI))*2.0*M_PI;  /* longitude of pericentre */
    }

    double wp = capom + omega;
    double lambda = capm + wp;
    lambda = lambda - int(lambda/(2.0*M_PI))*2.0*M_PI;

    std::stringstream oss;
    oss <<  "#" <<  i << ": I= " << inc << " a= " << a << " e= " <<  e << 
      " o= " << wp << " l= " << lambda;
    return oss.str();
  }

#ifndef _SSE_

  real compute_force(const Particle::Vector &ptcl, Force::Vector &force)
  {
    const int PP_FLOP = 38;
    const int n = ptcl.size();
    const double t0 = mytimer::get_wtime();
    real Utot = 0.0;

#ifdef _N2FAST_

    for (int i = 0; i < n; i++)
      force[i] = Force(0.0);

    for (int i = 0; i < n-1; i++)
      for (int j = i+1; j < n; j++)
      {
        const vec3 dr   = ptcl[j].pos - ptcl[i].pos;
        const real r2   = dr*dr;
        const real rinv = 1.0/std::sqrt(r2);
        const real rinv1 = rinv * ptcl[i].mass*ptcl[j].mass;
        const real rinv2 = rinv*rinv;
        const real rinv3 = rinv1*rinv2;

        const real pot = rinv1;
        const vec3 acc = rinv3 * dr;
        force[i].acc += acc;
        force[j].acc -= acc;
        Utot         += pot;
      }

    for (int i = 0; i < n; i++)
      force[i].acc *= 1.0/ptcl[i].mass;

#else /* _N2FAST_ */

    for (int i = 0; i < n; i++)
    {
      const Particle &pi = ptcl[i];
      vec3 iacc(0.0);
      real ipot(0.0);
      for (int j = 0; j < n; j++)
      {
        const Particle &pj = ptcl[j];
        if (i == j) continue;

        const vec3   dr = pj.pos - pi.pos;    
        const real   r2 = dr*dr;             
        const real rinv2 = 1.0/r2;          
        const real rinv1 = pj.mass * std::sqrt(rinv2);  
        const real rinv3 = rinv1*rinv2;                

        const real pot = rinv1; 
        const vec3 acc = rinv3 * dr;      

        iacc += acc;                     
        ipot += pot;                    
      }                                
      force[i] = Force(iacc);
      Utot    += 0.5*pi.mass*ipot;
    }

#endif /* _N2FAST_ */

    Utot = -Utot;
    dt_force += mytimer::get_wtime() - t0;
    flops += PP_FLOP*n*n;
    return Utot;
  }

#else /* _SSE_ */

  v2df compute_force(const ParticleSIMD::Vector &ptcl, ForceSIMD::Vector &force) 
  {
    const int PP_FLOP = 38;

#ifdef _N2FAST_

    asm("#SSE-checkpoint1"); 
    const double t0 = mytimer::get_wtime();
    const int n  = this->ptcl.size();
    assert((n&1) == 0);

    const int nh = n >> 1;

    for (int i = 0; i < nh; i++)
      force[i] = ForceSIMD(0.0);

    v2df simdU = {0.0, 0.0};
    for (int i = 0; i < nh; i++)
      for (int j = i; j < nh; j++)
      {
        asm("#SSE-checkpoint3a"); 

        const ParticleSIMD pi = ptcl[i];
        const ParticleSIMD p1(  ptcl[j], true );
        const ParticleSIMD p2(  ptcl[j], false);

        const v2df dx_1 = p1.posx - pi.posx;
        const v2df dy_1 = p1.posy - pi.posy;
        const v2df dz_1 = p1.posz - pi.posz;
        const v2df r2_1 = dx_1*dx_1 + dy_1*dy_1 + dz_1*dz_1;   

        const v2df dx_2 = p2.posx - pi.posx;
        const v2df dy_2 = p2.posy - pi.posy;
        const v2df dz_2 = p2.posz - pi.posz;
        const v2df r2_2 = dx_2*dx_2 + dy_2*dy_2 + dz_2*dz_2;   

        v2df rinv_1 = r2_1;
        v2df rinv_2 = r2_2;
        __rsqrtpd(rinv_1, rinv_2);

        const long long Imask = -(i!=j);
        const v2di mask = (v2di){Imask, Imask};

        const v2df rinv1_1 = __builtin_ia32_andpd((v2df)mask, pi.mass*p1.mass * rinv_1);
        const v2df rinv1_2 =                                  pi.mass*p2.mass * rinv_2 ;
        const v2df rinv3_1 = rinv_1*rinv_1 * rinv1_1;
        const v2df rinv3_2 = rinv_2*rinv_2 * rinv1_2;

        simdU += rinv1_1 + rinv1_2;

        const v2df fx_1233 = rinv3_1*dx_1;
        const v2df fx_1244 = rinv3_2*dx_2;
        const v2df fy_1233 = rinv3_1*dy_1;
        const v2df fy_1244 = rinv3_2*dy_2;
        const v2df fz_1233 = rinv3_1*dz_1;
        const v2df fz_1244 = rinv3_2*dz_2;

        force[i].accx += fx_1233 + fx_1244;
        force[i].accy += fy_1233 + fy_1244;
        force[i].accz += fz_1233 + fz_1244;

        const v2df fx_3411 = __builtin_ia32_unpcklpd(fx_1233, fx_1244);
        const v2df fx_3422 = __builtin_ia32_unpckhpd(fx_1233, fx_1244);
        const v2df fy_3411 = __builtin_ia32_unpcklpd(fy_1233, fy_1244);
        const v2df fy_3422 = __builtin_ia32_unpckhpd(fy_1233, fy_1244);
        const v2df fz_3411 = __builtin_ia32_unpcklpd(fz_1233, fz_1244);
        const v2df fz_3422 = __builtin_ia32_unpckhpd(fz_1233, fz_1244);

        force[j].accx -= fx_3411 + fx_3422;
        force[j].accy -= fy_3411 + fy_3422;
        force[j].accz -= fz_3411 + fz_3422;

        asm("#SSE-checkpoint3b"); 
      }

    for (int i = 0; i < nh; i++)
    {
      const v2df invM = (v2df){1.0, 1.0}/ptcl[i].mass;
      force[i].accx *= invM;
      force[i].accy *= invM;
      force[i].accz *= invM;
    }
    asm("#SSE-checkpoint4"); 
    const real Utot = __builtin_ia32_vec_ext_v2df(__builtin_ia32_haddpd(simdU, simdU), 0);

#else /* _N2FAST_ */

    asm("#SSE-checkpoint1"); 
    const double t0 = mytimer::get_wtime();
    const int n  = this->ptcl.size();
    assert((n&1) == 0);

    const int nh = n >> 1;

    for (int i = 0; i < nh; i++)
      force[i] = ForceSIMD(0.0);

    v2df simdU = {0.0, 0.0};
    for (int i = 0; i < nh; i++)
    {
      asm("#SSE-checkpoint2"); 
      const ParticleSIMD pi = ptcl[i];

      v2df accx = {0.0, 0.0};
      v2df accy = {0.0, 0.0};
      v2df accz = {0.0, 0.0};
      for (int j = 0; j < nh; j++)
      {
        asm("#SSE-checkpoint3a"); 

        const ParticleSIMD p1(ptcl[j], true);
        const ParticleSIMD p2(ptcl[j], false);

        const v2df dx_1 = p1.posx - pi.posx;
        const v2df dy_1 = p1.posy - pi.posy;
        const v2df dz_1 = p1.posz - pi.posz;
        const v2df r2_1 = dx_1*dx_1 + dy_1*dy_1 + dz_1*dz_1;   

        const v2df dx_2 = p2.posx - pi.posx;
        const v2df dy_2 = p2.posy - pi.posy;
        const v2df dz_2 = p2.posz - pi.posz;
        const v2df r2_2 = dx_2*dx_2 + dy_2*dy_2 + dz_2*dz_2;   

        v2df rinv_1 = r2_1;
        v2df rinv_2 = r2_2;
        __rsqrtpd(rinv_1, rinv_2);

        const v2df rinv1_1 = p1.mass * rinv_1;
        const v2df rinv2_1 = rinv_1*rinv_1;
        const v2df rinv3_1 = rinv2_1 * rinv1_1;
        const v2df rinv1_2 = p2.mass * rinv_2;
        const v2df rinv2_2 = rinv_2*rinv_2;
        const v2df rinv3_2 = rinv2_2 * rinv1_2;

        simdU += (v2df){0.5, 0.5}*pi.mass*(rinv1_1 + rinv1_2);
        accx += rinv3_1 * dx_1 + rinv3_2*dx_2;
        accy += rinv3_1 * dy_1 + rinv3_2*dy_2;
        accz += rinv3_1 * dz_1 + rinv3_2*dz_2;

        asm("#SSE-checkpoint3b"); 
      }                     
      force[i].accx = accx;
      force[i].accy = accy;
      force[i].accz = accz;
    }
    asm("#SSE-checkpoint4"); 
    const real Utot = __builtin_ia32_vec_ext_v2df(__builtin_ia32_haddpd(simdU, simdU), 0);
#endif  /* _N2FAST */

    dt_force += mytimer::get_wtime() - t0;
    flops += PP_FLOP*n*n;
    return (v2df){Utot, Utot};
  }
#endif  /* _SSE_ */

#if 0
  void iterate(const real EPS = 1.0e-13, const real atol = 1.0e-32)
  {
    const double t0 = mytimer::get_wtime();
    const real EP[] = {.4e-1, .16e-2, .64e-4, .256e-5};
    static real D[7];

    const int n = ptcl.size();
    const int N = n*6 + 1;

    const int NMAX = 1024*6 + 1;
    assert(N <= NMAX);
    static real DT[NMAX][7], YR[NMAX], YS[NMAX], Y[NMAX], S[NMAX];


    Particle::Vector y(ptcl);
    for (int i = 0; i < n; i++)
    {
      const int i6 = i*6;
      Y[i6 + 0] = y[i].pos[0];
      Y[i6 + 1] = y[i].pos[1];
      Y[i6 + 2] = y[i].pos[2];
      Y[i6 + 3] = y[i].vel[0];
      Y[i6 + 4] = y[i].vel[1];
      Y[i6 + 5] = y[i].vel[2];
    }
    Y[N-1] = y[0].time;

    for (int i = 0; i < N; i++)
      S[i] = std::abs(YS[i]);

    bool failed = true;
    real H = h;

    int JTI = 0;
    real FY = 1.0;
    real redu = 0.8;
    real Odot7 = 0.7;

    const int jmax = 10;

    ntry1++;
    while (failed)
    {
      bool BO = false;

      int M  = 1;
      int JR = 2;
      int JS = 3;

      for (int j = 0; j < jmax; j++)
      {
        ntry++;

        for (int i = 0; i < N; i++)
        {
          YS[i] = Y[i];
          S[i]  = std::max(std::abs(YS[i]), std::abs(S[i]));
        }

        if (BO)
        {
          D[1]=1.777777777777778e0;
          D[3]=7.111111111111111e0;
          D[5]=2.844444444444444e1;
        }
        else
        {
          D[1] = 2.25e0;
          D[3] = 9.e0;
          D[5] = 36.0e0;
        }

        int L;
        if (j+1 > 7)
        {
          L = 7-1;
          D[6] = 6.4e1;
        } 
        else
        {
          L = j;
          D[L] = M*M;
        }

        bool KONV = L+1 > 3;
        Multistep(M, H, ptcl, y);
        for (int i = 0; i < n; i++)
        {
          const int i6 = i*6;
          YS[i6 + 0] = y[i].pos[0];
          YS[i6 + 1] = y[i].pos[1];
          YS[i6 + 2] = y[i].pos[2];
          YS[i6 + 3] = y[i].vel[0];
          YS[i6 + 4] = y[i].vel[1];
          YS[i6 + 5] = y[i].vel[2];
        }
        YS[N-1] = y[0].time;

        const bool KL = L+1 > 2;
        const bool GR = L+1 > 5;

        real FS = 0.0;

        for (int i = 0; i < N; i++)
        {
          real V = DT[i][0];
          real C = YS[i];
          DT[i][0] = C;
          real TA = C;
          if (KL)
          {
            real W;
            for (int k = 1; k < L; k++)
            {
              real B1 = D[k]*V;
              real B  = B1 - C;
              W  = C -  V;
              real U  = V;
              if (B != 0.0)
              {
                B = W/B;
                U = C*B;
                C = B1*B;
              }
              V = DT[i][k];
              DT[i][k] = U;
              TA = U + TA;
            }
            real SI = std::max(S[i], std::abs(TA));
            if (std::abs(YR[i] - TA) > atol + SI*EPS) KONV = false;
            if (!(GR || SI == 0.0))
            {
              real FV = std::abs(W)/SI;
              if (FS < FV) FS = FV;
            }
          }
          YR[i] = TA;
        }

        if (FS != 0.0) 
        {
          real FA = FY;
          int  K  = L - 1;
          assert(K >= 0);
          FY = std::pow(EP[K]/FS, 1.0/(L + K));
          //          fprintf(stderr, " -- M= %d  FY= %g  j= %d \n", M, FY, j);
          FY = std::min(FY, 1.4);
          if (!((L+1 != 2 && FY<Odot7*FA)|| FY > Odot7))
          {
            H = H*FY;
            JTI = JTI+1;
            if (JTI+1 > 25)
            {
              H = 0.0;
              assert(0);
              return;
            }
            failed = true;
            break;
          }
        }

        if (KONV)
        {
          for (int i = 0; i < n; i++)
          {
            const int i6 = i*6;
            y[i].pos[0] = YR[i6 + 0];
            y[i].pos[1] = YR[i6 + 1];
            y[i].pos[2] = YR[i6 + 2];
            y[i].vel[0] = YR[i6 + 3];
            y[i].vel[1] = YR[i6 + 4];
            y[i].vel[2] = YR[i6 + 5];
          }
          y[0].time = YR[N-1];
          //          fprintf(stderr, "M= %d  h_old= %g h_new= %g  FY= %g \n", M, h, H*FY, FY);
          h    = H*FY;

          ptcl = y;

          dt = ptcl[0].time - time;
          time = ptcl[0].time;
          iteration++;
#if 0
          assert(0);
#endif
          dt_step += mytimer::get_wtime() - t0;
          return;
        }

        D[2] = 4.0e0;
        D[4] = 1.6e1;

        BO = !BO;
        M  = JR;
        JR = JS;
        JS = M + M;
      }

      redu = redu*redu + 0.001;
      H = H*redu;
      failed = true;
    }

  }
#else
  void iterate(const real rtol = 1.0e-13, const real atol = 1.0e-32, const bool dense = false)
  {
    assert(rtol >= 1.0e-15);
    const double t0 = mytimer::get_wtime();
    const  int  KMAXX = 7;
    const  int  IMAXX = KMAXX + 1;

    static int  nseq [IMAXX];
    static int  cost [IMAXX];
    static real coeff[IMAXX][IMAXX];
    static real errfac[2*IMAXX+2];
    static int  ipoint[IMAXX+1];
    static Particle::Vector table[KMAXX];

    static real hnext = -1;
    static int k_targ = -1;
    if (k_targ == -1)
    {
      if (dense)
        for (int i = 0; i < IMAXX; i++)
          nseq[i] = 4*i + 2;
      else
      {
#if 0  /* sequences suggested by Numerical Recipes, 3rd edition */
        for (int i = 0; i < IMAXX; i++)
          nseq[i] = 2*(i+1);
#endif

#if 0  /* Harmonic sequence */
        for (int i = 0; i < IMAXX; i++)
          nseq[i] = i + 1;
#endif

#if 1   /* Bulirsch sequence works the best, especially with Rational function extrapolation */
        nseq[0] = 2;
        nseq[1] = 3;
        for (int i = 2; i < IMAXX; i++)
          nseq[i] = nseq[i-2] << 1;
#endif
      }
      cost[0] = nseq[0] + 1;
      for (int k = 0; k < KMAXX; k++)
        cost[k+1] = cost[k] + nseq[k+1];

      const real logfact = -log(std::max(1.0e-16, rtol))/log(10.0)*0.6 + 0.5;
      k_targ = std::max(1, std::min(KMAXX-1, (int)logfact));

      for (int k = 0; k < IMAXX; k++)
        for (int l = 0; l < k; l++)
        {
          const real ratio = real(nseq[k])/(real)nseq[l];
#ifndef _RATIONAL_FUNCTIONS_
          coeff[k][l] = 1.0/(ratio*ratio - 1.0);
#else
          coeff[k][l] = ratio*ratio;
#endif
        }

      for (int i = 0; i < 2*IMAXX + 1; i++)
      {
        const int ip5 = i + 5;
        errfac[i] = 1.0/(ip5*ip5);
        const real e = 0.5*std::sqrt(real(i+1.0)/ip5);
        for (int j = 0; j <= i; j++)
          errfac[i] *= e/(j+1);
      }
      ipoint[0] = 0;
      for (int i = 1; i <= IMAXX; i++)
      {
        int njadd = 4*i - 2;
        if (nseq[i-1] > njadd) 
          njadd++;
        ipoint[i] = ipoint[i-1] + njadd;
      }
    }

    const real STEPFAC1 = 0.65;
    const real STEPFAC2 = 0.94;
    const real STEPFAC3 = 0.02;
    const real STEPFAC4 = 4.0;
    const real    KFAC1 = 0.8;
    const real    KFAC2 = 0.9;

    static real hopt[IMAXX], work[IMAXX];

    static bool first_step = true, last_step = false; 
    static bool forward, reject = false, prev_reject = false;

    const real htry = h;

    work[0] = 0;
    h = htry;
    forward = h > 0;
    assert(!dense);
    assert(forward);

    Particle::Vector y;
    const Particle::Vector &ysav = ptcl;

    if (h != hnext && !first_step)
      last_step = true;

    if (reject)
    {
      prev_reject = true;
      last_step   = false;
    }

    reject = false;
    bool firstk = true;

    real hnew = std::abs(h);

    int k;
    ntry1++;
    while (firstk || reject)
    {
      h = forward ? hnew : -hnew;

      firstk = false;
      reject = false;

      //      int ipt = -1;
      for (k = 0; k <= k_targ + 1; k++)
      {
        ntry++;

        if (k == 0)  Multistep(nseq[k], h, ptcl, y         );
        else         Multistep(nseq[k], h, ptcl, table[k-1]);

        if (k != 0)
        {
          Extrapolate<KMAXX, IMAXX>(k, coeff, table, y);
          const real err    = Error(y, table[0], ysav, rtol, atol);
          const real expo   = 1.0/(2*k + 1);
          const real facmin = std::pow(STEPFAC3, expo);
          real fac;
          if (err == 0.0) 
            fac = 1.0/facmin;
          else
          {
            fac = STEPFAC2/std::pow(err/STEPFAC1, expo);
            fac = std::max(facmin/STEPFAC4, std::min(1.0/facmin, fac));
          }
          hopt[k] = std::abs(h*fac);
          work[k] = cost[k]/hopt[k];

          if ((first_step || last_step) && err <= 1.0)
            break;

          if (k == k_targ - 1 && !prev_reject && !first_step && !last_step)
          {
            if (err <= 1.0)
              break;
            else if (err > SQR(nseq[k_targ]*nseq[k_targ+1]/SQR(nseq[0])))
            {
              reject = true;
              k_targ = k;
              if (k_targ > 1 && work[k-1] < KFAC1*work[k])
                k_targ--;
              hnew = hopt[k_targ];
              break;
            }
          }

          if (k == k_targ)
          {
            if (err <= 1.0)
              break;
            else if (err > SQR(nseq[k+1]/nseq[0]))
            {
              reject = true;
              if (k_targ > 1 && work[k-1] < KFAC1*work[k])
                k_targ--;
              hnew = hopt[k_targ];
              break;
            }
          }

          if (k == k_targ+1)
          {
            if (err > 1.0)
            {
              reject = true;
              if (k_targ > 1 && work[k_targ - 1] < KFAC1*work[k_targ])
                k_targ--;
              hnew = hopt[k_targ];
            }
            break;
          }

        }
      }
      if (reject)
        prev_reject = true;
    }

    ptcl = y;

    first_step = false;
    int kopt;

    if (k == 1)
      kopt = 2;
    else if (k <= k_targ)
    {
      kopt = k;
      if (work[k-1] < KFAC1*work[k])
        kopt = k-1;
      else if (work[k] < KFAC2*work[k-1])
        kopt = std::min(k+1, KMAXX-1);
    }
    else
    {
      kopt = k-1;
      if (k > 2 && work[k-2] < KFAC1*work[k-1])
        kopt = k-2;
      if (work[k] < KFAC2*work[kopt])
        kopt = std::min(k, KMAXX-1);
    }

    if (prev_reject)
    {
      k_targ = std::min(kopt, k);
      hnew = std::min(std::abs(h), hopt[k_targ]);
      prev_reject = false;
    }
    else
    {
      if (kopt <= k)
        hnew = hopt[kopt];
      else
      {
        if (k < k_targ && work[k] < KFAC2*work[k-1])
          hnew = hopt[k] * cost[kopt+1]/cost[k];
        else
          hnew = hopt[k] * cost[kopt  ]/cost[k];
      }
      k_targ = kopt;
    }
    if (forward)
      hnext =  hnew;
    else
      hnext = -hnew;

    h  = hnext;

    dt   = ptcl[0].time - time;
    time = ptcl[0].time;
    iteration++;
    dt_step += mytimer::get_wtime() - t0;
  }
#endif

  void SymplecticCorrector(const real h, const Particle::Vector &ptcl, Particle::Vector &dPQ, real &dT)
  {
    const int n = ptcl.size();

    std::vector< std::pair<vec3, vec3> > force(n);

    real U = 0.0;
    real T = 0.0;
    real W = 0.0;
    const real P = negE0; 
    for (int i = 0; i < n; i++)
    {
      const Particle &pi = ptcl[i];
      T += pi.Ekin();
      const real mih = 0.5*pi.mass;
      vec3 acc = 0.0;
      vec3 jrk = 0.0;
      for (int j = 0; j < n; j++)
      {
        if (i == j) continue;
        const Particle &pj = ptcl[j];
        const vec3 dr = pj.pos - pi.pos;
        const vec3 dv = pj.vel - pi.vel;
        const real r2 = dr.norm2();
        const real rv = dr*dv;
 
        const real rinv  = 1.0/std::sqrt(r2);
        const real rinv2 = rinv    * rinv;
        const real rinv1 = pj.mass * rinv;
        const real rinv3 = rinv1   * rinv2;
        const real alpha = rv*rinv2;

        const vec3 Aij = rinv3*dr;
        const vec3 Jij = rinv3*dv - (3.0*alpha)*Aij;

        acc += Aij;
        jrk += Jij;
        U   += mih*rinv1;
      }
      force[i] = std::make_pair(acc, jrk);
      W       += ptcl[i].mass * (ptcl[i].vel * acc);
    }
    const real invU = 1.0/U;
    const real invP = 1.0/(T + P);
    const real invF = h*h/24.0*invU*invP;

    W *= -invF;

    for (int i = 0; i < n; i++)
    {
      dPQ[i].pos = -1*((W*invP)* ptcl[i].vel   + invF*force[i].first );
      dPQ[i].vel =  1*((W*invU)*force[i].first + invF*force[i].second);
    }

    dT = -W*invP;
  }

  void cvt2symp(const real h, Particle::Vector &ptcl, real &T)
  {
    Particle::Vector ptcl0(ptcl), dPQ(ptcl);
    real T0(T);
    real dT;

    const int n = ptcl.size();

    int iter = 100;
    while (1)
    {
      SymplecticCorrector(h, ptcl, dPQ, dT);
      real err = 0.0;

      for (int i = 0; i < n; i++)
      {
        vec3 dpos = ptcl[i].pos;
        vec3 dvel = ptcl[i].vel;

        ptcl[i].pos = ptcl0[i].pos - dPQ[i].pos;
        ptcl[i].vel = ptcl0[i].vel - dPQ[i].vel;
        
        dpos -= ptcl[i].pos;
        dvel -= ptcl[i].vel;

        err = std::max(err, dpos.norm2()/ptcl[i].pos.norm2());
        err = std::max(err, dvel.norm2()/ptcl[i].vel.norm2());
      }
      real del = T;
      T = T0 - dT;
      del -= T;
      err = std::max(err, SQR(del/T));


      err = std::sqrt(err);

      if (err < 1.0e-15) break;

      iter--;
      assert(iter > 0);
    }
  }
  
  void cvt2phys(const real h, Particle::Vector &ptcl, real &T)
  {
    Particle::Vector dPQ(ptcl);
    real dT;

    const int n = ptcl.size();

    SymplecticCorrector(h, ptcl, dPQ, dT);

    for (int i = 0; i < n; i++)
    {
      ptcl[i].pos += dPQ[i].pos;
      ptcl[i].vel += dPQ[i].vel;

    }
    T += dT;
    assert(T > 0.0);
  }

#ifdef _SYMPCORR_ /* with symplectic corrector */
  void Multistep(const int nsteps, const real tau, const Particle::Vector &ptcl_in, Particle::Vector &ptcl)
  {
    const double t0 = mytimer::get_wtime();
    const int nbody = ptcl_in.size();
    ptcl = ptcl_in;
    real time = this->time;
    real h = tau/(real)nsteps;

#if 1
    cvt2symp(h, ptcl, time);
#endif

    real T = negE0 + Ekin(ptcl);

    real dth = 0.5*h/T;
    for (int i = 0; i < nbody; i++)
      if (ptcl[i].mass > SMALLM)
        ptcl[i].update_pos(dth);
    time += dth;
    for (int k = 0; k < nsteps; k++)
    {
      const real U = compute_force(ptcl, force);
      const real h_over_U = h/(-U);
      T = 0;
      for (int i = 0; i < nbody; i++)
      {
        if (ptcl[i].mass > SMALLM)
          ptcl[i].update_vel(h_over_U, force[i].acc);
        T += ptcl[i].mass*ptcl[i].vel.norm2();
      }
      T = 0.5*T + negE0;

      real dth = h/T;
      if (k == nsteps-1) dth *= 0.5;
      for (int i = 0; i < nbody; i++)
        if (ptcl[i].mass > SMALLM)
          ptcl[i].update_pos(dth);
      time += dth;
    }

#if 1
    cvt2phys(h, ptcl, time);
#endif

    ptcl[0].time = time;


    dt_multi += mytimer::get_wtime() - t0;
  }

#else /* with symplectic corrector */

  void Multistep(const int nsteps, const real tau, const Particle::Vector &ptcl_in, Particle::Vector &ptcl)
  {
    const double t0 = mytimer::get_wtime();
    const int nbody = ptcl_in.size();
    ptcl.resize(nbody);

#ifndef _SSE_

    real h = tau/(real)nsteps;
    real T = negE0 + Ekin(ptcl_in);
    real time = this->time;

    real dth = 0.5*h/T;
    for (int i = 0; i < nbody; i++)
    {
      ptcl[i] = ptcl_in[i];
      if (ptcl[i].mass > SMALLM)
        ptcl[i].update_pos(dth);
    }
    time += dth;
    for (int k = 0; k < nsteps; k++)
    {
      const real U = compute_force(ptcl, force);
      const real h_over_U = h/(-U);
      T = 0;
      for (int i = 0; i < nbody; i++)
      {
        if (ptcl[i].mass > SMALLM)
          ptcl[i].update_vel(h_over_U, force[i].acc);
        T += ptcl[i].mass*ptcl[i].vel.norm2();
      }
      T = 0.5*T + negE0;

      real dth = h/T;
      if (k == nsteps-1) dth *= 0.5;
      for (int i = 0; i < nbody; i++)
        if (ptcl[i].mass > SMALLM)
          ptcl[i].update_pos(dth);
      time += dth;
    }

    ptcl[0].time = time;

#else /* _SSE_ */

    for (int i = 0; i < nbody; i += 2)
      ptclSIMD[i>>1] = ParticleSIMD(&ptcl_in[i]);

    const int n = nbody >> 1;

    const v2df h = {tau/(real)nsteps, tau/(real)nsteps}; 
    const v2df negE0v = {negE0, negE0};
    v2df  time = {this->time, this->time};

    const v2df zero = {0.0, 0.0};
    const v2df half = {0.5, 0.5};

    v2df T = zero;
    for (int i = 0; i < n; i++)
    {
      const ParticleSIMD &p = ptclSIMD[i];
      T += p.mass * (p.velx*p.velx + p.vely*p.vely + p.velz*p.velz);
    }
    T = half*reduce(T) + negE0v;

    v2df dth = half*h/T;
    for (int i = 0; i < n; i++)
    {
      ParticleSIMD &p = ptclSIMD[i];
      const v2df mask = __builtin_ia32_cmpgtpd(p.mass, (v2df){SMALLM, SMALLM});
      const v2df dt1  = __builtin_ia32_andpd(mask, dth);
      p.posx += p.velx * dt1;
      p.posy += p.vely * dt1;
      p.posz += p.velz * dt1;
    }
    time += dth;

    for (int k = 0; k < nsteps; k++)
    {
      const v2df U = compute_force(ptclSIMD, forceSIMD);
      v2df dt = h/U;
      T = zero;
      for (int i = 0; i < n; i++)
      {
        ParticleSIMD &p = ptclSIMD[i];
        const v2df mask = __builtin_ia32_cmpgtpd(p.mass, (v2df){SMALLM, SMALLM});
        const v2df dt1  = __builtin_ia32_andpd(mask, dt);
        p.velx += forceSIMD[i].accx * dt1;
        p.vely += forceSIMD[i].accy * dt1;
        p.velz += forceSIMD[i].accz * dt1;
        T += p.mass * (p.velx*p.velx + p.vely*p.vely + p.velz*p.velz);
      }
      T = half*reduce(T) + negE0v;

      dt = h/T;
      if (k == nsteps-1) 
        dt *= half;
      for (int i = 0; i < n; i++)
      {
        ParticleSIMD &p = ptclSIMD[i];
        const v2df mask = __builtin_ia32_cmpgtpd(p.mass, (v2df){SMALLM, SMALLM});
        const v2df dt1  = __builtin_ia32_andpd(mask, dt);
        p.posx += p.velx * dt1;
        p.posy += p.vely * dt1;
        p.posz += p.velz * dt1;
      }
      time += dt;
    }


    for (int i = 0; i < nbody; i += 2)
    {
      ptcl[i  ] = ptclSIMD[i>>1].scalar(0);
      ptcl[i+1] = ptclSIMD[i>>1].scalar(1);
    }

    ptcl[0].time = __builtin_ia32_vec_ext_v2df(time, 0);

#endif /* _SSE_ */


    dt_multi += mytimer::get_wtime() - t0;
  }
#endif


#ifndef _RATIONAL_FUNCTIONS_

  template<int KMAXX, int IMAXX>
    void Extrapolate(const int k, const real coeff[IMAXX][IMAXX], Particle::Vector table[KMAXX], Particle::Vector &last)
    {
      const double t0 = mytimer::get_wtime();
      const int n = last.size();
      for (int j = k-1; j > 0; j--)
      {
        for (int i = 0; i < n; i++)
        {
          table[j-1][i].pos = table[j][i].pos + coeff[k][j]*(table[j][i].pos - table[j-1][i].pos);
          table[j-1][i].vel = table[j][i].vel + coeff[k][j]*(table[j][i].vel - table[j-1][i].vel);
        }
        table[j-1][0].time = table[j][0].time + coeff[k][j]*(table[j][0].time - table[j-1][0].time);
      }
      for (int i = 0; i < n; i++)
      {
        last[i].pos = table[0][i].pos + coeff[k][0]*(table[0][i].pos - last[i].pos);
        last[i].vel = table[0][i].vel + coeff[k][0]*(table[0][i].vel - last[i].vel);
      }
      last[0].time = table[0][0].time + coeff[k][0]*(table[0][0].time - last[0].time);
      dt_extra += mytimer::get_wtime() - t0;
    }

#else

  template<int KMAXX, int IMAXX>
    void Extrapolate(const int k, const real coeff[IMAXX][IMAXX], Particle::Vector table[KMAXX], Particle::Vector &last)
    {
      const double t0 = mytimer::get_wtime();
      const int n = last.size();
      assert(k > 0);

      const real SMALL = 1.0e-30;
     
      static Particle::Vector tmp(n);

      for (int j = k-1; j > 0; j--)
      {
        for (int i = 0; i < n; i++)
        {
          const vec3 dr0 = table[j][i].pos - table[j-1][i].pos;
          const vec3 dv0 = table[j][i].vel - table[j-1][i].vel;
          const vec3 dr1 = table[j][i].pos - tmp       [i].pos;
          const vec3 dv1 = table[j][i].vel - tmp       [i].vel;

          tmp[i].pos = table[j-1][i].pos;
          tmp[i].vel = table[j-1][i].vel;


          for (int l = 0; l < 3; l++)
          {
            asm("#-DIV-0");
            const real dq = coeff[k][j]*(dr1[l] - dr0[l]) - dr1[l] + SMALL;
            const real dp = coeff[k][j]*(dv1[l] - dv0[l]) - dv1[l] + SMALL;
#ifndef _SSE_
            const real idq = 1.0/dq;
            const real idp = 1.0/dp;
#else
            real idq = dq, idp = dp;
            __divpd(idq, idp);
#endif
            table[j-1][i].pos[l] = table[j][i].pos[l] + dr0[l]*dr1[l] * idq;
            table[j-1][i].vel[l] = table[j][i].vel[l] + dv0[l]*dv1[l] * idp;
            asm("#-DIV-1");
          }
        }
        const real dt0 = table[j][0].time - table[j-1][0].time;
        const real dt1 = table[j][0].time - tmp       [0].time;
        const real dt_ = coeff[k][j]*(dt1 - dt0) - dt1 + SMALL;
        table[j-1][0].time = table[j][0].time + dt0*dt1 / dt_;
      }

      for (int i = 0; i < n; i++)
      {
        const vec3 dr0 = table[0][i].pos - last[i].pos;
        const vec3 dv0 = table[0][i].vel - last[i].vel;
        const vec3 dr1 = table[0][i].pos - tmp [i].pos;
        const vec3 dv1 = table[0][i].vel - tmp [i].vel;

        tmp[i].pos = 0.0;
        tmp[i].vel = 0.0;

        for (int l = 0; l < 3; l++)
        {
          const real dq = coeff[k][0]*(dr1[l] - dr0[l]) - dr1[l] + SMALL;
          const real dp = coeff[k][0]*(dv1[l] - dv0[l]) - dv1[l] + SMALL;
#ifndef _SSE_
          const real idq = 1.0/dq;
          const real idp = 1.0/dp;
#else
          real idq = dq, idp = dp;
          __divpd(idq, idp);
#endif
          last[i].pos[l] = table[0][i].pos[l] + dr0[l]*dr1[l] * idq;
          last[i].vel[l] = table[0][i].vel[l] + dv0[l]*dv1[l] * idp;
        }
      }
      const real dt0 = table[0][0].time - last[0].time;
      const real dt1 = table[0][0].time - tmp [0].time;
      const real dt_ = coeff[k][0]*(dt1 - dt0) - dt1 + SMALL;
      tmp[0].time = 0.0;
      last[0].time = table[0][0].time + dt0*dt1 / dt_;

      dt_extra += mytimer::get_wtime() - t0;
    }
#endif

  real Error(const Particle::Vector &y, const Particle::Vector &y1, const Particle::Vector &ysav,
      const real rtol, const real atol)
  {
    const double t0 = mytimer::get_wtime();
    assert(y.size() == ysav.size());
    const int n = y.size();

    float err = 0.0;
    for (int i = 0; i < n; i++)
      for (int k = 0; k < 3; k++)
      {
        const float scale_pos = atol + rtol*std::max(std::abs(ysav[i].pos[k]), std::abs(y[i].pos[k]));
        const float scale_vel = atol + rtol*std::max(std::abs(ysav[i].vel[k]), std::abs(y[i].vel[k]));
        err = std::max(err, std::abs(float(y[i].pos[k] - y1[i].pos[k])/scale_pos));
        err = std::max(err, std::abs(float(y[i].vel[k] - y1[i].vel[k])/scale_vel));
      }
#if 1
    const float scale_time = atol + rtol*std::max(std::abs(ysav[0].time), std::abs(y[0].time));
    err = std::max(err, std::abs(float(y[0].time - y1[0].time)/scale_time));
#endif

    dt_err += mytimer::get_wtime() - t0;
    return err;
  }

};


#endif /* __LOGH_H__ */
