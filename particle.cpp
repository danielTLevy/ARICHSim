#include "particle.h"

Particle::Particle(TVector3 pos, TVector3 dir, double beta) {
  this->pos = pos;
  this->dir = dir;
  this->beta = beta;
}

double Particle::theta() {
	return atan(sqrt(dir[0]*dir[0] + dir[1]*dir[1]) / dir[2]);
}

double Particle::phi() {
	return  atan(dir[0] / dir[1]);
}

double Particle::dist() {
    return pos[2] / cos(theta());
}