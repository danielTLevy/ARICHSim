#include "particle.h"

Particle::Particle(TVector3 pos, TVector3 dir) {
  this->pos = pos;
  this->dir = dir;
}

double Particle::theta() {
  return atan(sqrt(dir[0]*dir[0] + dir[1]*dir[1]) / dir[2]);
}

double Particle::phi() {
  return atan(dir[0] / dir[1]);
}

double Particle::dist(double z) {
  // distance to plane in direction of travel
  return (z - pos[2]) / cos(theta());
}