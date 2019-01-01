#include "stdlib.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TH2D.h"
#include "TFile.h"
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
  double beta = 0.99;
  double n = 1.06;
  double dist[2] = {21.0,19.0};
  double thickness = dist[0] - dist[1];
  double x_0 = -3.0;
  double y_0 = -2.0;
  TVector3 pos_0 = TVector3(x_0, y_0, dist[0]);
  double dirX_0 = 0.05;
  double dirY_0 = -0.1;
  double dirZ_0 = sqrt(1 - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 dir_0 = TVector3(dirX_0, dirY_0, dirZ_0);
  double errX = 0.0;
  double errY = 0.0;
  double errDirX = 0.01;
  double errDirY = 0.01;
  // Generate beam
  int nIter = 10000;
  Beam *beam = new Beam(pos_0, dir_0, beta, errX, errY, errDirX, errDirY);
  beam->makeParticles(nIter);

  // Make Aerogel layer
  Aerogel* aerogel = new Aerogel(2.0, n, dist[0], beam, beta);
  // Get plots ready
  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  TH2D *beamHist = beam->plotParticles();
  TH2D *photonHist = new TH2D("photonHist","photonHist",200,-15,15,200,-15,15);

  //double xs[100]; double ys[100]; double ps[100]; double ts[100];

  for (int i = 0; i < nIter; i++) {
    // Generate particles

    Particle *pa = beam->getParticle(i);

    std::vector<Photon*> photons = aerogel->generatePhotons(pa);

    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;

      double phDist = ph->dist();
      photonHist->Fill(phPos[0] + phDist*phDir[0]/phDir[2], phPos[1] + phDist*phDir[1]/phDir[2]);

    }

  }

    c1->cd();
    beamHist->Draw("colz");
    photonHist->Draw("samecolz");
    c1->SaveAs("beamAndPhoton.pdf");

};
