#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <iostream>
#include <vector>
#include <cassert>
#include "vector3.h"
#include "memalign_allocator.h"

struct Particle
{
  typedef std::vector<Particle, __gnu_cxx::malloc_allocator<Particle, 64> > Vector;
  real mass;        // 1
  vec3 pos;         // 4
  vec3 vel;         // 7
  real time;        // 8

  Particle() {}
  Particle(const real c) : mass(c), pos(c), vel(c), time(c) {}
  Particle(const real m, const vec3 &_pos, const vec3 _vel) : mass(m), pos(_pos), vel(_vel) {}

  real Ekin() const {return mass*vel.norm2()*0.5;}
  vec3 momentum() const {return mass*vel;}
	
  void update_vel(const real h_over_U, const vec3 &acc)
	{
    vel += h_over_U * acc;
	}
	void update_pos(const real h_half_over_TPt)
	{
		pos += h_half_over_TPt * vel;
	}
};


#ifdef _SSE_
typedef double v2df __attribute__((vector_size(16)));
typedef float  v4sf __attribute__((vector_size(16)));
typedef long long v2di __attribute__((vector_size(16)));
struct ParticleSIMD
{
  typedef std::vector<ParticleSIMD, __gnu_cxx::malloc_allocator<ParticleSIMD, 128> > Vector;
  v2df mass;
  v2df posx;
  v2df posy;
  v2df posz;
  v2df velx;
  v2df vely;
  v2df velz;
  v2df time;
  ParticleSIMD() {}
  ParticleSIMD(const Particle &ptcl) 
  { 
    mass = (v2df){ptcl.mass,  ptcl.mass };
    posx = (v2df){ptcl.pos.x, ptcl.pos.x};
    posy = (v2df){ptcl.pos.y, ptcl.pos.y};
    posz = (v2df){ptcl.pos.z, ptcl.pos.z};
    velx = (v2df){ptcl.vel.x, ptcl.vel.x};
    vely = (v2df){ptcl.vel.y, ptcl.vel.y};
    velz = (v2df){ptcl.vel.z, ptcl.vel.z};
  }
  ParticleSIMD(const Particle *p)
  { 
    mass = (v2df){p[0].mass,  p[1].mass };
    posx = (v2df){p[0].pos.x, p[1].pos.x};
    posy = (v2df){p[0].pos.y, p[1].pos.y};
    posz = (v2df){p[0].pos.z, p[1].pos.z};
    velx = (v2df){p[0].vel.x, p[1].vel.x};
    vely = (v2df){p[0].vel.y, p[1].vel.y};
    velz = (v2df){p[0].vel.z, p[1].vel.z};
  }

  ParticleSIMD(const ParticleSIMD &p, const bool flag)
  {
    if (flag)
    {
      mass = __builtin_ia32_unpcklpd(p.mass, p.mass);
      posx = __builtin_ia32_unpcklpd(p.posx, p.posx);
      posy = __builtin_ia32_unpcklpd(p.posy, p.posy);
      posz = __builtin_ia32_unpcklpd(p.posz, p.posz);
    }
    else
    {
      mass = __builtin_ia32_unpckhpd(p.mass, p.mass);
      posx = __builtin_ia32_unpckhpd(p.posx, p.posx);
      posy = __builtin_ia32_unpckhpd(p.posy, p.posy);
      posz = __builtin_ia32_unpckhpd(p.posz, p.posz);
    }
  }

  Particle scalar(const int i) const
  {
    Particle p;
    p.mass = __builtin_ia32_vec_ext_v2df(mass, i);
    p.time = __builtin_ia32_vec_ext_v2df(time, i);
    p.pos = vec3(
        __builtin_ia32_vec_ext_v2df(posx, i),
        __builtin_ia32_vec_ext_v2df(posy, i),
        __builtin_ia32_vec_ext_v2df(posz, i));
    p.vel = vec3(
        __builtin_ia32_vec_ext_v2df(velx, i),
        __builtin_ia32_vec_ext_v2df(vely, i),
        __builtin_ia32_vec_ext_v2df(velz, i));
    return p;
  }

};

inline void __divpd(double &a, double &b)
{
  const v2df one = {1.0, 1.0};
#if 1
  const v2df v = one/(v2df){a, b};
#else
  const v2df A = {a, b};
  const v4sf x = __builtin_ia32_rcpps((v4sf){a, b, a, b});
  v2df v = __builtin_ia32_cvtps2pd(x);
  const v2df h  = one - A*v;
  const v2df h2 = h*h;
//  v *= (one + h) * (one + h2);
//  v *= (one + (one + h2)*(h + h2));
  v *= (one + h)*(one + h2)*(one + h2*h2);
#endif

  a = __builtin_ia32_vec_ext_v2df(v, 0);
  b = __builtin_ia32_vec_ext_v2df(v, 1);
}


inline void __rsqrtpd(v2df &r2_1, v2df &r2_2)
{
#if 0
  const v2df flag_1  = __builtin_ia32_cmpgtpd(r2_1, (v2df){0.0, 0.0});
  const v2df flag_2  = __builtin_ia32_cmpgtpd(r2_2, (v2df){0.0, 0.0});
  r2_1 = __builtin_ia32_andpd(flag_1, (v2df){1.0, 1.0}/__builtin_ia32_sqrtpd(r2_1));
  r2_2 = __builtin_ia32_andpd(flag_2, (v2df){1.0, 1.0}/__builtin_ia32_sqrtpd(r2_2));
#endif

#if 1
  const v4sf x = __builtin_ia32_movlhps(__builtin_ia32_cvtpd2ps(r2_1),__builtin_ia32_cvtpd2ps(r2_2));
  const v4sf ysp = __builtin_ia32_rsqrtps(x);
  v2df y1 = __builtin_ia32_cvtps2pd(ysp);
  v2df y2 = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(ysp, ysp));
  const v2df c1 = {-0.5, -0.5};
  const v2df c2 = {-3.0, -3.0};
  y1 = (c1 * y1) * (r2_1*y1*y1 + c2);
  y2 = (c1 * y2) * (r2_2*y2*y2 + c2);
  y1 = (c1 * y1) * (r2_1*y1*y1 + c2);
  y2 = (c1 * y2) * (r2_2*y2*y2 + c2);
  const v2df flag_1  = __builtin_ia32_cmpgtpd(r2_1, (v2df){0.0, 0.0});
  const v2df flag_2  = __builtin_ia32_cmpgtpd(r2_2, (v2df){0.0, 0.0});
  r2_1  = __builtin_ia32_andpd(flag_1, y1);
  r2_2  = __builtin_ia32_andpd(flag_2, y2);
#endif

#if 0
  const v4sf x = __builtin_ia32_movlhps(__builtin_ia32_cvtpd2ps(r2_1),__builtin_ia32_cvtpd2ps(r2_2));
  const v4sf ysp = __builtin_ia32_rsqrtps(x);
  const v2df y1 = __builtin_ia32_cvtps2pd(ysp);
  const v2df y2 = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(ysp, ysp));
  const v2df c1 = {1.0, 1.0};
  const v2df c2 = {0.5, 0.5};
  const v2df c3 = {0.375, 0.375};
  const v2df c4 = {0.3125, 0.3125};
  const v2df z1 = c1 - r2_1*y1*y1;
  const v2df z2 = c1 - r2_2*y2*y2;
  const v2df flag_1  = __builtin_ia32_cmpgtpd(r2_1, (v2df){0.0, 0.0});
  const v2df flag_2  = __builtin_ia32_cmpgtpd(r2_2, (v2df){0.0, 0.0});
  r2_1  = __builtin_ia32_andpd(flag_1, y1*(c1 + z1*(c2 + z1*(c3 + z1*c4))));
  r2_2  = __builtin_ia32_andpd(flag_2, y2*(c1 + z2*(c2 + z2*(c3 + z2*c4))));
#endif
}

inline v2df reduce(const v2df val)
{
  const double v = __builtin_ia32_vec_ext_v2df(__builtin_ia32_haddpd(val, val), 0);
  return (v2df){v,v};
}


#endif  /* _SSE_ */


#endif //  __PARTICLE_H__
