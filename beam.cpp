#include "beam.h"

Beam::Beam(TVector3 pos0, TVector3 dir0, double beta, double errX, double errY, double errDirX, double errDirY) {
  this->pos0 = pos0;
  this->dir0 = dir0;
  this->beta = beta;
  this->errX = errX;
  this->errY = errY;
  this->errDirX = errDirX;
  this->errDirY = errDirY;
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  //particles.reserve(10000);
}

double Beam::getBeta() {
  return beta;
}

Particle Beam::generateParticle() {
  double paX = pos0(0) + randomGenerate->Gaus(0.0, errX);
  double paY = pos0(1) + randomGenerate->Gaus(0.0, errY);
  TVector3 pos = TVector3(paX, paY, 0);
  double paDirX = pos0(0) + randomGenerate->Gaus(0.0, errDirX);
  double paDirY = pos0(1) + randomGenerate->Gaus(0.0, errDirY);
  double paDirZ =  sqrt(1 - paDirX*paDirX - paDirY*paDirY);
  TVector3 dir = TVector3(paDirX, paDirY, paDirZ);
  double paTheta =  atan(sqrt(paDirX*paDirX + paDirY*paDirY) / paDirZ);
  double paPhi = atan(paDirX / paDirY);

  // optional: Project beam trajectory onto detector
  double z = pos[2];
  double r = z / cos(paTheta);
  //beamHist->Fill(paX+r*paDirX, paY+r*paDirY);

  return Particle(pos, dir, beta);
}

void Beam::makeParticles(int N) {

}

void Beam::plotParticles() {

};