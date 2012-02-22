#ifndef __KEPLER_H__
#define __KEPLER_H__

#include <iostream>
#include <cassert>
#include "vector3.h"

#if 0
struct Kepler
{
  vec3 pos, vel;
  
  Kepler(const vec3 &_pos, const vec3 &_vel, const double Mo, const double dt = 0.0) : 
    pos(_pos), vel(_vel)
  {
    if (dt == 0.0) return;
    
    assert(step(Mo, dt));
  }

  bool step(const double Mo, const double dt)
  {
    const double    R = pos.abs();
    assert(R > 0.0);
    const double invR = 1.0/R;
    const double invM = 1.0/Mo;
    const vec3 Lvec = pos%vel;
    const vec3 Evec = (vel%Lvec)*invM - pos*invR; 
    const vec3 Bvec = Lvec%Evec;

    const double hs   = (vel*vel)*0.5 - Mo*invR;
    const double k    = -2.0*hs*invM;
    const double j    = std::sqrt(2.0*std::abs(hs));
    const double n    = j*j*j*invM;

    const double e = Evec.abs();
    if (e == 0.0)
    {
      assert(hs < 0.0);
      assert(false);
    }
    else if (e < 1.0)
    {
      const double inve = 1.0/e;
      assert(j > 0.0);
      const double invj = 1.0/j;
      assert(k > 0.0);
      const double invk = 1.0/k;

      const double x = inve * (1.0 - n*R*invj);
      const double y = n*invk*invM*inve * (pos*vel);
      const double E0 = std::atan2(y,x);
      const double nP = -(E0 - e*sin(E0));

      double M = n*dt - nP;
      while (M < 0.0)     M += 2*M_PI;
      while (M >= 2*M_PI) M -= 2*M_PI;

      const double E = EccentricAnomaly(e, M);
      const double cv = std::sin(E);
      const double cr = std::cos(E);

     
#if 0 
      assert((pos%vel).abs() > 0.0);
      
      fprintf(stderr, "pos= %g %g %g  vel= %g %g %g \n",
          pos.x, pos.y, pos.z,
          vel.x, vel.y, vel.z);
      
      pos = invk*(cr*inve - 1.0)*Evec + cv*inve*invj*Bvec;
      
      const double invr = 1.0/pos.abs();
      
      vel = inve*invr*(-Mo*cv*invj*Evec + cr*Bvec); 

      fprintf(stderr, "M= %g cr= %g  cv= %g  pos= %g %g %g  vel= %g %g %g \n",
          M,
          cr, cv,
          pos.x, pos.y, pos.z,
          vel.x, vel.y, vel.z);
      {
        const vec3 cr = Evec%Bvec;
        fprintf(stderr," Evec= %g %g %g \n", Evec.x, Evec.y, Evec.z);
        fprintf(stderr," Bvec= %g %g %g \n", Bvec.x, Bvec.y, Bvec.z);
        fprintf(stderr," cr= %g %g %g \n", cr.x, cr.y, cr.z);
      }
      assert((Evec%Bvec).abs() > 0.0);
      assert((pos%vel).abs() > 0.0);
#else
      const double invr = k/(1.0 - e*cr);
      pos = invk*(cr*inve - 1.0)*Evec + cv*inve*invj*Bvec;
      vel = inve*invr*(-Mo*cv*invj*Evec + cr*Bvec); 
#endif
    }
    else if (e == 1.0)
    {
      assert(0);
    }
    else  /* e > 1.0 */
    {
      assert(0);
    }

    return true;
  }

  static double sign(const double x) {return x < 0.0 ? -1.0 : (x > 0.0 ? 1.0 : 0.0);}

  static double EccentricAnomaly(const double e, const double M)
  {
    int niter = 100;
    const double tol = 1.0e-10;

    double Ei = M + 0.85*e*sign(std::sin(M));
    double E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));

#if 0

    while (std::abs(E - Ei) > tol*std::abs((E + Ei)*0.5))
    {
      if (niter < 10)
      {
        fprintf(stderr, "Ei= %g  E= %g  diff= %g \n",
            Ei, E, Ei-E);
      }
      assert(niter >= 0);
      Ei = E;
      E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));
      niter--;
    }
#else

    while (std::abs(E - Ei) > tol*E)
    {
      if (niter < 10)
      {
        fprintf(stderr, "Ei= %g  E= %g  diff= %g \n",
            Ei, E, (Ei-E)/E);
      }
      assert(niter >= 0);
      const double sinE = std::sin(E);
      const double cosE = std::cos(E);
      const double f3 = e*cosE;
      const double f2 = e*sinE;
      const double f1 = 1.0 - f3;
      const double f0 = E - e*sinE - M;
      const double d1 = -f0/ f1;
      const double d1h = 0.5*d1;
      const double d2 = -f0/(f1 + d1h* f2);
      const double d3 = -f0/(f1 + d1h*(f2 + (1.0/3.0)*d2*d2*f3));

      Ei = E;
      E  = Ei + d3;

      niter--;
    }
#endif

    return E;
  }

};

#else


struct Kepler
{
  vec3 pos, vel;

  double G_func2(double q) 
  {
    int l = 3;
    int d = 15;
    int n = 0;
    double A, B, G;

    if(q==0.0) return 1.0;	/* this isn't necessary when first
                               Newt-Raph iteration is done by hand */

    A = B = G = 1.0;

    while (fabs(B/G)>1e-15) 
    {
      l += 2;
      d += 4*l;
      n += 10*l;

      A = d/(d-n*A*q);
      B *= A-1.0;
      G += B;

      l += 2;
      d += 4*l;
      n -= 8*l;

      A = d/(d-n*A*q);
      B *= A-1.0;
      G += B;
    }

    return G;
  };

  Kepler(const vec3 &_pos, const vec3 &_vel, const double Mo, const double dt = 0.0) : 
    pos(_pos), vel(_vel)
  {
    if (dt == 0.0) return;

    if (!step(Mo, dt))
    {
      assert(step(Mo, dt/2.0));
      assert(step(Mo, dt/2.0));
    }
  }

  bool step(const double Mo, const double dt)
  {
    double r0mag, v0mag2;
    double r0v0;	/* r dot v */
    double rcalc, dtcalc, terr;
    double u;	/* (?) universal variable */
    double beta;	/* (?) vis-a-vis integral */
    double P;	/* period (for elliptic orbits only) */
    double dU;
    int n;
    double q;
    double U0w2, U1w2;
    double U, U0, U1, U2, U3;
    double f, g, F, G;
    int no_iter;

    double du1, du2, du3, dqdu, d2qdu2, drdu, d2rdu2, fn, fnp, fnpp, fnppp;

    r0mag  = pos.abs();   // sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
    v0mag2 = vel.norm2(); // v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2];
    r0v0   = pos*vel;     // r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2];
    beta   = 2*Mo/r0mag - v0mag2;

    if (beta > 0) 
    {	
      /* elliptic orbit */
      P = 2*M_PI*Mo/sqrt(beta*beta*beta);
      n = floor((dt + P/2.0 -2*r0v0/beta)/P);
      dU = 2*n*M_PI/sqrt(beta*beta*beta*beta*beta);
    } 
    else 
    {
      dU = 0.0;
    }

    u = 0;	/* a "better" guess is possible, see footnote at Battin p.219 */
    //	u = dt/(4*r0mag); 				  /* N-R step by hand */
    //	u = dt/(4*r0mag-2*dt*r0v0/r0mag);
    //	u = init_u(dt, mu, r0mag, r0v0, beta);

    no_iter = 0;
    do 
    {
      q = beta*u*u/(1+beta*u*u);
      if (q > 0.5 || no_iter > 12) 
      {
        return false;
        assert(0);
      }
      //        return DRIFT_FAIL;
      dqdu = 2*beta*u/(1+beta*u*u)/(1+beta*u*u);
      d2qdu2 = 2*beta/(1+beta*u*u) 
        - 8*beta*beta*u*u / (1+beta*u*u)/(1+beta*u*u)/(1+beta*u*u);
      U0w2 = 1 - 2*q;
      U1w2 = 2*(1-q)*u;
      U = 16.0/15 * U1w2*U1w2*U1w2*U1w2*U1w2 * G_func2(q) + dU;
      U0 = 2*U0w2*U0w2 - 1;
      U1 = 2*U0w2*U1w2;
      U2 = 2*U1w2*U1w2;
      U3 = beta*U + U1*U2/3.0;
      rcalc = r0mag*U0 + r0v0*U1 + Mo*U2;
      drdu   = 4*(1-q)*(r0v0*U0 + (Mo-beta*r0mag)*U1);
      d2rdu2 = -4*dqdu*(r0v0*U0 + (Mo-beta*r0mag)*U1)
        + (4*(1-q)*4*(1-q))*(-beta*r0v0*U1 + (Mo-beta*r0mag)*U0);
      dtcalc = r0mag*U1 + r0v0*U2 + Mo*U3;

      fn    = dtcalc-dt;
      fnp   = 4*(1-q)*rcalc;
      fnpp  = 4*(drdu*(1-q) - rcalc*dqdu);
      fnppp = -8*drdu*dqdu - 4*rcalc*d2qdu2 + 4*(1-q)*d2rdu2;

      du1  = -fn/fnp;
      du2  = -fn/(fnp + du1*fnpp/2);
      du3  = -fn/(fnp + du2*fnpp/2 + du2*du2*fnppp/6);

      u += du3;
      no_iter++;

      terr = fabs((dt-dtcalc)/dt);
    } 
    while (terr > 1e-15);

    f = 1 - (Mo/r0mag)*U2;
    g = r0mag*U1 + r0v0*U2;
    F = -Mo*U1/(rcalc*r0mag);
    G = 1 - (Mo/rcalc)*U2;

    const vec3 r1 = f*pos + g*vel;
    const vec3 v1 = F*pos + G*vel;

    pos = r1;
    vel = v1;
    return true;
  }


  static double EccentricAnomaly(const double e, const double M)
  {
    const int niter_max = 100;

    double Ei = M/2.0;
    double E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));

    const double tol = 1.0e-13;
    int iter = 0;
    while (std::abs(E - Ei) > tol*std::abs((E + Ei)*0.5))
    {
      Ei = E;
      E  = Ei + (M + e*sin(Ei) - Ei)/(1 - e*cos(Ei));
      assert(iter < niter_max);
      iter++;
    }

    return E;
  }

};
#endif


#endif // __KEPLER_H__
