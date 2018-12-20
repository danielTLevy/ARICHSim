#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
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
  particles.reserve(10000);
}

double Beam::getBeta() {
  return beta;
}

Particle* Beam::generateParticle() {
  double paX = pos0(0) + randomGenerate->Gaus(0.0, errX);
  double paY = pos0(1) + randomGenerate->Gaus(0.0, errY);
  TVector3 pos = TVector3(paX, paY, pos0(2));
  double paDirX = dir0(0) + randomGenerate->Gaus(0.0, errDirX);
  double paDirY = dir0(1) + randomGenerate->Gaus(0.0, errDirY);
  double paDirZ =  sqrt(1 - paDirX*paDirX - paDirY*paDirY);
  TVector3 dir = TVector3(paDirX, paDirY, paDirZ);

  return new Particle(pos, dir, beta);
}

void Beam::makeParticles(int N) {
  for (int i = 0; i < N; i++) {
    particles.push_back(generateParticle());
  }
  n_particles = N;
}

TH2D* Beam::plotParticles() {
  TH2D *beamHist = new TH2D("beamHist2","beamHist2",200,-15,15,200,-15,15);
  for (int i = 0; i < n_particles; i++) {
    Particle p = *particles[i];
    double paTheta =  atan(sqrt(p.dir[0]*p.dir[0] + p.dir[1]*p.dir[1]) / p.dir[2]);
    double z = p.pos[2];
    double r = z / cos(paTheta);
    beamHist->Fill(p.pos[0]+r*p.dir[0], p.pos[1]+r*p.dir[1]);
  }
  return beamHist;
};