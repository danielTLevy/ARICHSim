#include "photon.h"

Photon::Photon(TVector3 pos, TVector3 dir, double wavelength) : Particle(pos, dir) {
  this->pos0 = pos;
  this->dir0 = dir;
  this->pos = pos;
  this->dir = dir;
  this->wav = wavelength;
}

double Photon::getWavelength() {
  return wav;
}
