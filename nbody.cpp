#include <cstdlib>
#include <fstream>

#if 0
#include "LogH.h"
#else
#include "LogH+TTL.h"
#endif
#include "kepler.h"
#include "mytimer.h"

#ifndef __MACOSX_
#define __LINUX__
#endif

#ifdef __MACOSX__
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
inline void fpe_catch() {
  _mm_setcsr( _MM_MASK_MASK &~
              (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
}
#elif defined __LINUX__
#include <fenv.h>
void fpe_catch(void) {
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#else
crap
void fpe_catch(void) {}
#endif

typedef double real;

/* position are in AU, G = 1, time is in years */
/* Period @ R = 1 AU is 1 year */

Particle::Vector read_xyz(const int nbody)
{
  Particle::Vector ptcl;
  assert(nbody > 1);

  real mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);

  for (int i = 0; i < nbody; i++)
  {
    int idummy;
    real mass;
    vec3 pos, vel;
    std::cin >> idummy >> mass >> pos.x >> pos.y >> pos.z >> vel.x >> vel.y >> vel.z;
    mass *= mass_scale;
    pos  *= pos_scale;
    vel  *= vel_scale;
#if 0
    fprintf(stderr, " i =%d : mass= %g  pos= %g %g %g (r= %g)   vel= %g %g %g  (v= %g)   cos(r,v)= %g\n",
        i, mass,
        pos.x, pos.y, pos.z, pos.abs(),
        vel.x, vel.y, vel.z, vel.abs(), 
        pos*vel/(pos.abs()*vel.abs()));
#endif

    ptcl.push_back(Particle(mass, pos, vel));
  }
  
  fprintf(stderr, "read= %d  particles \n", nbody);
  return ptcl;
}

Particle::Vector read_aei(const int nbody)
{
  Particle::Vector ptcl;
  assert(nbody > 0);
  
  real mass_scale, pos_scale, vel_scale;
  std::cin >> mass_scale >> pos_scale >> vel_scale;
  fprintf(stderr, " scaling mass by %g \n", mass_scale);
  fprintf(stderr, " scaling position by %g \n", pos_scale);
  fprintf(stderr, " scaling velocity by %g \n", vel_scale);

  real Mcentre;
  int index;
  std::cin >> index;
  std::cin >> Mcentre;

  fprintf(stderr, " Mcentre= %g \n", Mcentre);
  ptcl.push_back(Particle(Mcentre, 0.0, 0.0));

  for (int i = 0; i < nbody-1; i++)
  {
    real mass;
    real a; // semi-major axis
    real e; // eccentricity
    real I; // inclination in degrees;
    real w; // argument of the pericentre in degrees
    real O; // longitude of the ascending node in degrees
    real M; // mean anomaly in degrees;
    std::cin >> index >> 
      mass >> 
      a    >> 
      e    >> 
      I    >> 
      w    >> 
      O    >> 
      M;
    fprintf(stderr, " i =%d : index= %d mass= %g  a= %g  e= %g  I= %g  w= %g  O= %g  M= %g :: r_peri= %g  r_apo= %g\n",
        i, index, mass,
        a, e, I, w, O, M,
        a*(1.0 -e), a*(1.0 +e));

    w *= M_PI/180.0;
    I *= M_PI/180.0;
    M *= M_PI/180.0;
    O *= M_PI/180.0;

#if 0
    const real x  = sma[i]*(1+ecc[i]);
    const real y  = 0.0;
    const real vx = 0.0;
    const real vy = sqrt(1.0/sma[i]*(1-ecc[i])/(1+ecc[i]));
#else
    const real Mt   = Mcentre + mass;
    const real E    = Kepler::EccentricAnomaly(e, M);
    const real Edot = sqrt(Mt/(a*a*a))/(1 - e*cos(E));

    const real x    =  a*(cos(E) - e);
    const real y    =  a*sqrt(1.0 - e*e)*sin(E);
    const real vx   = -a*sin(E)*Edot;
    const real vy   =  a*sqrt(1.0 - e*e)*cos(E)*Edot;
#endif

    real r, f;

    r = sqrt(x*x+y*y);
    f = atan2(y,x);

    const vec3 pos(
        r*(cos(O)*cos(w+f) - sin(O)*sin(w+f)*cos(I)),
        r*(sin(O)*cos(w+f) + cos(O)*sin(w+f)*cos(I)),
        r*(sin(w+f)*sin(I)));

    r = sqrt(vx*vx+vy*vy);
    f = atan2(vy,vx);

    const vec3 vel(
        r*(cos(O)*cos(w+f) - sin(O)*sin(w+f)*cos(I)),
        r*(sin(O)*cos(w+f) + cos(O)*sin(w+f)*cos(I)),
        r*(sin(w+f)*sin(I)));

    fprintf(stderr, " m= %g  pos= %g %g %g  vel= %g %g %g \n",
        mass,
        pos.x, pos.y, pos.z,
        vel.x, vel.y, vel.z);

    ptcl.push_back(Particle(mass, pos, vel));
  }

  // positions are heliocentric
  // momenta   are  barycentric

  fprintf(stderr, "read= %d  particles \n", nbody-1);
  assert(nbody == (int)ptcl.size());

  return ptcl;
}

int main(int argc, char * argv[])
{
#ifndef _SSE_
  fpe_catch();
#endif

  const real h = 1.0e-5;
  
  unsigned long long iteration = 0;
  Particle::Vector ptcl;
  
  int nbody, npl;
  std::cin >> nbody >> npl;
  assert(nbody != 0);
  assert(npl >= 0);

  real tepoch, Tscale;
  std::cin >> tepoch >> Tscale;

  int dI;
  real Tend, dt_log, dt_out, dt_snap;
  std::cin >> Tend >> dt_out >> dt_log >> dt_snap >> dI;
  assert(dt_out > 0.0);
  assert(dt_log > 0.0);
  assert(dt_snap > 0.0);
  assert(Tscale > 0.0);
  const unsigned long long di_iter_max = 1LLU << 62;
  const unsigned long long di_iter = dI < 0 ? di_iter_max : dI;

  fprintf(stderr, " Tepoch= %g   Tscale= %g \n", tepoch, Tscale);
  fprintf(stderr, " Tend= %g \n", Tend);
  fprintf(stderr, " dTout= %g \n", dt_out);
  fprintf(stderr, " dTlog= %g \n", dt_log);
  fprintf(stderr, " dTsnap= %g \n", dt_snap);

  real tolerance;
  std::cin >> tolerance;
  fprintf(stderr, " Tolerance= %g \n", tolerance);


  if (nbody < 0)
  {
    fprintf(stderr, " Reading aei format ... \n");
    ptcl = read_aei(-nbody);
  }
  else
  {
    fprintf(stderr, " Reading xyz format ... \n");
    ptcl = read_xyz(nbody);
  }

  Nbody s(iteration, tepoch, h, ptcl);

  const real E0 = s.Etot();
  fprintf(stderr, " tepoch= %g   Etot= %g \n", tepoch, E0);

  real Ep = E0;

  real t_log  = s.time/Tscale;
  real t_out  = s.time/Tscale;
  real t_snap = s.time/Tscale;


  real de_max = 0.0;

  double t0 = mytimer::get_wtime();
  const double t_start = t0;
  while (s.time/Tscale < Tend)
  {
    s.iterate(tolerance);
    const double t1 = mytimer::get_wtime() + 1.0/HUGE;
    if (s.time/Tscale > t_log || s.time/Tscale >= Tend || s.iteration % di_iter == 0)
    {
      const real E1 = s.Epot() + s.Ekin();
      de_max = std::max(de_max, std::abs((E1 - E0)/E0));
      fprintf(stderr, "iter= %llu :: t= %g dt= %4.3g  dE= %4.3g ddE= %4.3g dEmax= %4.3g  <n>= %3.2g  Twall= %4.3g hr | <T>= %4.3g [%4.3g %4.3g %4.3g %4.3g %4.3g] sec  <P>= %4.3g [ %4.3g ] GFLOP/s] \n",
          s.iteration,
          s.time/Tscale, s.dt/Tscale,
          (E1 - E0)/std::abs(E0),
          (E1 - Ep)/std::abs(Ep),
          de_max, s.ntry1 > 0 ? (real)s.ntry/(real)s.ntry1 : 1, 
          (t1 - t_start) / 3600.0,
          t1 - t0,
          s.dt_step,
          s.dt_multi,
          s.dt_force,
          s.dt_extra,
          s.dt_err,
          s.get_gflops(), s.get_gflops_force());
      t0 = t1;
      s.reset_counters();
      t_log += dt_log;
      Ep = E1;
    }

    if (s.time/Tscale > t_out || s.time/Tscale >= Tend)
    {
      fprintf(stdout, "%g ", s.time/Tscale);
      for (int ipl = 1; ipl <= npl; ipl++)
        fprintf(stdout,"%s ", s.print_orbit(ipl).c_str());
      fprintf(stdout, "\n");
      fflush(stdout);
      t_out += dt_out;
    }
   
    if (s.time/Tscale >= t_snap || s.time/Tscale >= Tend)
    {
      t_snap += dt_snap;
      /*************/

      char fn[256];
      const char path[] = "output";
      sprintf(fn, "%s/snap_%.8d.out", path, int(t_snap/dt_snap));
      std::ofstream fout(fn);
     
     
      char o1[256]; 
      sprintf(o1, "%lld \n", di_iter == di_iter_max ? -1 : di_iter);

      fout.precision(15);
      fout << s.ptcl.size() << " " << npl << std::endl;
      fout << s.time << " " << Tscale << std::endl;
      fout << Tend << " " << dt_out << " " << dt_log  << " " << dt_snap << " " << o1;
      fout << tolerance << std::endl;
      fout << "1.0 1.0 1.0 \n";
      fout << s.print_output();

      fout.close();

      /*************/
    }
  }

  return 0;

}
