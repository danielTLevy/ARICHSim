#ifndef PHOTON_INCLUDE
#define PHOTON_INCLUDE

#include "particle.h"

class Photon : public Particle {
public:
  Photon(TVector3, TVector3, double);
  double wav;
  double getWavelength();
  int numScatters = 0;
};

#endif