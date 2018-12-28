#ifndef PARTICLE_INCLUDE
#define PARTICLE_INCLUDE

#include "TVector3.h"

class Particle {
public:
  Particle(TVector3, TVector3, double);
  TVector3 pos;
  TVector3 dir;
  double beta;
  double theta();
  double phi();
  double dist();
};

#endif