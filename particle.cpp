#include "particle.h"

Particle::Particle(TVector3 pos, TVector3 dir, double beta) {
  this->pos = pos;
  this->dir = dir;
  this->beta = beta;
}