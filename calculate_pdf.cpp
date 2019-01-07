#include "stdlib.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TEllipse.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVector3.h"
#include "beam.h"
#include "aerogel.h"
#include "particle.h"
#include "photon.h"


int main(int argc, char *argv[]) {
  // Beam parameters:
  double beta = 0.99; // velocity of particle
  double n = 1.035; // index of refraction
  double aeroPos[2] = {0., 2.0}; // positions of aerogel planes
  double dist = 21.0; // dist to detector plane
  double thickness = 2.0; // thickness of aerogel layer
  double x_0 = -3.0; // beam initial position
  double y_0 = -2.0;
  TVector3 pos_0 = TVector3(x_0, y_0, 0);
  double dirX_0 = 0.05; // beam initial direction
  double dirY_0 = -0.1;
  double dirZ_0 = sqrt(1 - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 dir_0 = TVector3(dirX_0, dirY_0, dirZ_0);
  double errX = 0.0; // beam position error
  double errY = 0.0;
  double errDirX = 0.01; // beam direction error
  double errDirY = 0.01;
  // Generate beam
  int nIter = 10000; // number of particles simulated for beam
  Beam *beam = new Beam(pos_0, dir_0, beta, errX, errY, errDirX, errDirY);
  beam->makeParticles(nIter);

  // Make Aerogel layer
  Aerogel* aerogel = new Aerogel(thickness, n, dist-aeroPos[0], beta);
  // Get plots ready
  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  TH2D *beamHist = beam->plotParticles(dist);
  TH2D *photonHist = new TH2D("photonHist","photonHist",200,-15,15,200,-15,15);

  Photon* ph;

  for (int i = 0; i < nIter; i++) {
    // Generate particles

    Particle *pa = beam->getParticle(i);
    std::vector<Photon*> photons = aerogel->generatePhotons(pa);

    for (int j = 0; j < photons.size(); j++) {
      ph = photons[j];
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;
      double phDist = ph->dist(dist);
      photonHist->Fill(phPos[0] + phDist*phDir[0]/phDir[2], phPos[1] + phDist*phDir[1]/phDir[2]);

    }
  }
    c1->cd();
    beamHist->Draw("colz");
    photonHist->Draw("samecolz");
    c1->SaveAs("./output/beamAndPhoton.pdf");

};
