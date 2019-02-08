#ifndef BEAM_INCLUDE
#define BEAM_INCLUDE

#include "TRandom3.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "particle.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"

class Beam {
private:
  TVector3 pos0;
  TVector3 dir0;
  double beta;
  // Todo: Better errs
  double errX;
  double errY;
  double errDirX;
  double errDirY;
  std::vector<Particle*> particles;
  std::shared_ptr<TRandom3> randomGenerate;
  int n_particles = 0;



public:
  Beam(TVector3, TVector3, double, double, double, double, double);
  double getBeta();
  Particle* generateParticle();
  Particle* getParticle(int i);
  TH2D* plotParticles(double);
  void makeParticles(int N);
};

#endif