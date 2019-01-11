#include "stdlib.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
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
#include <chrono>
using namespace std;
using namespace std::chrono;

struct phStructruct {
  double dirx;
  double diry;
  double dirz;
  double posx;
  double posy;
  double posz;
  int paId;
};

struct particleStruct {
  double dirx;
  double diry;
  double dirz;
  double posx;
  double posy;
  double posz;
  int id;
};

int main(int argc, char *argv[]) {
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  // Beam parameters:
  double beta = 0.99; // velocity of particle
  double n = 1.035; // index of refraction
  double aeroPos[2] = {0., 2.0}; // positions of aerogel planes
  double dist = 21.0; // dist to detector plane
  double thickness = 2.0; // thickness of aerogel layer
  double x_0 = -3.0; // beam initial position
  double y_0 = -2.0;
  TVector3 pos_0 = TVector3(x_0, y_0, 0);
  double dirX_0 = 0.1; // beam initial direction
  double dirY_0 = -0.05;
  double dirZ_0 = sqrt(1 - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 dir_0 = TVector3(dirX_0, dirY_0, dirZ_0);
  double errX = 0.0; // beam position error
  double errY = 0.0;
  double errDirX = 0.01; // beam direction error
  double errDirY = 0.01;
  int nIter = 10000; // number of particles simulated for beam

  // Generate beam
  Beam *beam = new Beam(pos_0, dir_0, beta, errX, errY, errDirX, errDirY);
  beam->makeParticles(nIter);

  // Make Aerogel layer
  Aerogel* aerogel = new Aerogel(thickness, n, dist-aeroPos[0], beta);

  // Get plots ready
  TH2D *beamHist = beam->plotParticles(dist);
  TH2D *photonHist = new TH2D("photonHist","photonHist",200,-15,15,200,-15,15);
  TH1D *scatterHist = new TH1D("scatterHist", "scatterHist", 101, 0, 100);
  TH1D *distHist = new TH1D("distHist", "distHist", 200, 0., 5.);


  // Make a tree to save the photons

  phStructruct phStruct;
  particleStruct paStruct;
  TFile *f = new TFile("./output/photons.root","RECREATE");
  TTree *tree = new TTree("T","Output photon data");
  TBranch *phBranch = tree->Branch("photons",&phStruct.dirx,"dirx/D:diry:dirz:posx:posy:posz:parentid/i");
  TBranch *paBranch = tree->Branch("particles",&paStruct.dirx,"dirx/D:diry:dirz:posx:posy:posz:id/i");

  for (int i = 0; i < nIter; i++) {
    // Generate particles
    Particle *pa = beam->getParticle(i);
    paStruct.dirx = pa->dir[0];
    paStruct.diry = pa->dir[1];
    paStruct.dirz = pa->dir[2];
    paStruct.posx = pa->pos[0];
    paStruct.posy = pa->pos[1];
    paStruct.posz = pa->pos[2];
    paStruct.id = i;
    paBranch->Fill();

    std::vector<Photon*> photons = aerogel->generatePhotons(pa);
    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;
      phStruct.dirx = phDir[0];
      phStruct.diry = phDir[1];
      phStruct.dirz = phDir[2];
      phStruct.posx = phPos[0];
      phStruct.posy = phPos[1];
      phStruct.posz = phPos[2];
      phStruct.paId = i;
      phBranch->Fill();
      double phDist = ph->dist(dist);
      photonHist->Fill(phPos[0] + phDist*phDir[0]/phDir[2], phPos[1] + phDist*phDir[1]/phDir[2]);
      scatterHist->Fill(ph->numScatters);
      distHist->Fill(aerogel->getDistInGel(ph));
    }
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>( t2 - t1 ).count();
  cout << "Time taken: " << duration / 1000000. << endl;

  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  c1->cd();
  //beamHist->Draw("colz");
  photonHist->Draw("samecolz");
  c1->SaveAs("./output/beamAndPhotonScatter.pdf");

  TCanvas *c2 = new TCanvas("c2","c2",600,500);
  c2->cd();
  c2->SetLogy();
  scatterHist->Draw();
  c2->SaveAs("./output/scatterCount.pdf");


  TCanvas *c3 = new TCanvas("c3","c3",600,500);
  c3->cd();
  //c2->SetLogy();
  distHist->Draw();
  c3->SaveAs("./output/distHist.pdf");

  f->Write();
  tree->Print();
  return 0;
};
