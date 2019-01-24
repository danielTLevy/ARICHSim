#include "stdlib.h"
#include <iostream>
#include <chrono>
#include "TMath.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVector3.h"
#include "beam.h"
#include "aerogel.h"
#include "particle.h"
#include "photon.h"
#include "detector.h"
using namespace std;
using namespace std::chrono;

struct phStructruct {
  // Initial direction and position
  double dirxi;
  double diryi;
  double dirzi;
  double posxi;
  double posyi;
  double poszi;
  // Direction and position upon exit of aerogel
  double dirxe;
  double dirye;
  double dirze;
  double posxe;
  double posye;
  double posze; 
  // Parent particle
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
  double n1 = 1.035; // index of refraction
  double n2 = 1.045; // index of refraction

  double aeroPos[2] = {0., 2.0}; // positions of aerogel planes
  double dist = 21.0; // dist to detector plane
  double thickness = 2.0; // thickness of aerogel layer
  double x_0 = -0.0; // beam initial position
  double y_0 = -0.0;
  TVector3 pos_0 = TVector3(x_0, y_0, 0);
  double dirX_0 = 0.2; // beam initial direction
  double dirY_0 = -0.00;
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
  Aerogel* aerogel1 = new Aerogel(thickness, n1, aeroPos[0], beta);
  Aerogel* aerogel2 = new Aerogel(thickness, n2, aeroPos[1], beta);

  // Get plots ready
  TH2D *beamHist = beam->plotParticles(dist);
  TH2D *photonHist = new TH2D("photonHist","photonHist",200,-15,15,200,-15,15);
  TH1D *scatterHist = new TH1D("scatterHist", "scatterHist", 101, 0, 100);
  TH1D *distHist = new TH1D("distHist", "distHist", 200, 0., 5.);
  TH1D *numPhotonHist = new TH1D("numPhotonHist", "numPhotonHist", 200, 1.8, 2.1);


  // Make a tree to save the photons

  phStructruct phStruct;
  particleStruct paStruct;
  TFile *f = new TFile("./output/photons.root","RECREATE");
  TTree *tree = new TTree("T","Output photon data");
  TBranch *phBranch = tree->Branch("photons",&phStruct.dirxe,"dirx/D:diry:dirz:posx:posy:posz:parentid/i");
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
    numPhotonHist->Fill(aerogel1->getDistInGel(pa));

    // Make photons in first aerogel and scatter them
    std::vector<Photon*> photons1 = aerogel1->generatePhotons(pa);
    aerogel1->applyPhotonScatters(photons1);
    // Advance particle forward to next aerogel and generate photons
    pa->travelZDist(aeroPos[1] - aeroPos[0]);
    std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa);

    // Combine photons from both aerogels
    std::vector<Photon*> photons;
    photons.reserve(photons1.size() + photons2.size());
    photons.insert(photons.end(), photons1.begin(), photons1.end());
    photons.insert(photons.end(), photons2.begin(), photons2.end());
    // Include scattering in second aerogel
    aerogel2->applyPhotonScatters(photons);

    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;
      phStruct.dirxe = phDir[0];
      phStruct.dirye = phDir[1];
      phStruct.dirze = phDir[2];
      phStruct.posxe = phPos[0];
      phStruct.posye = phPos[1];
      phStruct.posze = phPos[2];
      phStruct.paId = i;
      phBranch->Fill();
      double phDist = ph->dist(dist);
      photonHist->Fill(phPos[0] + phDist*phDir[0], phPos[1] + phDist*phDir[1]);

      scatterHist->Fill(ph->numScatters);
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
  distHist->Draw();
  c3->SaveAs("./output/distHist.pdf");

  TCanvas *c4 = new TCanvas("c4","c4",600,500);
  c4->cd();
  numPhotonHist->Draw();
  c4->SaveAs("./output/numPhotonHist.pdf");

  f->Write();
  tree->Print();
  return 0;
};
