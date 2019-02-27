#include "stdlib.h"
#include <iostream>
#include <chrono>
#include <TROOT.h>
#include <TStyle.h>
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
#include "TEllipse.h"
#include "TCutG.h"
#include "beam.h"
#include "aerogel.h"
#include "particle.h"
#include "photon.h"
#include "detector.h"
using namespace std;
using namespace std::chrono;

const int nEvents = 1000; // number of particles simulated for beam
const double aeroPos[2] = {0., 2.0}; // positions of aerogel planes
const double thickness = 2.0; // thickness of aerogel layer
const double width = 10.; // x width of aerogel
const double height = 10.; // y height of aerogel
const double n1 = 1.035; // outer index of refraction
const double n2 = 1.045; // inner index of refraction
const double dist = 21.0; // dist to detector plane
const double errDirX = 0.000; // beam direction error
const double errDirY = 0.000;
const double errX = 0.001; // beam position error
const double errY = 0.001;
const double errBeta = 0.001;

const char* particleNames[3] = {"Pi", "Proton", "Kaon"};
const double particleMasses[3] = {0.1395701, 0.938272, 0.493677};

struct photonStruct {
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
  // Wavelength
  double wav;
  // Parent particle
  int paid;
  // Number of scatters
  int numscat;
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

double calcBeta(int particlei, double mom) {
  double M = particleMasses[particlei];
  return sqrt(1/(1 + M*M/(mom*mom)));
}

double integrateAndDrawEllipse(TVector3 pos0, TVector3 dir0, double beta, TH2D* photonHist, TCanvas* canvas, Aerogel* aerogel) {
  // Define Ellipse and integrate over this ring
  double dirX_0 = dir0[0];
  double dirY_0 = dir0[1];
  double dirZ_0 = dir0[2];

  // Naively assume that all photons are generated in middle of second aerogel layer.
  double middlePoint = thickness + thickness / 2.;
  double newX_0 = pos0[0] + middlePoint*dirX_0/dirZ_0;
  double newY_0 = pos0[1] + middlePoint*dirY_0/dirZ_0;
  double newDist = dist - middlePoint;
  double dirTheta =  atan(sqrt(dirX_0*dirX_0 +dirY_0*dirY_0) / dirZ_0);
  double chAngle = aerogel->getChAngle();
  // Calculate ellipse parameters
  double radiusA = 0.5*newDist*TMath::Abs(tan(dirTheta+chAngle) - tan(dirTheta-chAngle));
  double radiusB = 0.5*newDist*TMath::Abs(tan(chAngle) - tan(-chAngle));
  // Distance travelled on detector plane with respect to travel direction
  double deltaR = newDist*tan(dirTheta-chAngle) + radiusA;
  // Multiply by X and Y components of direction to get final x and y position
  double ringX, ringY;
  if (dirX_0 != 0 || dirY_0 != 0) {
    ringX = newX_0 + deltaR*dirX_0/sqrt(dirX_0*dirX_0 +dirY_0*dirY_0);
    ringY = newY_0 + deltaR*dirY_0/sqrt(dirX_0*dirX_0 +dirY_0*dirY_0);
  } else {
    ringX = newX_0;
    ringY = newY_0;
  }
  // Rotate the ellipse by the particle's phi direction
  double dirPhiDeg;
  if (dirX_0 != 0) {
    dirPhiDeg = atan(dirY_0/dirX_0)*180./TMath::Pi();
  } else {
    dirPhiDeg = 0.;
  }
  // Create two ellipses to encapsulate the photons
  double ringOuterA = radiusA+1.5;
  double ringOuterB = radiusB+1.5;
  TEllipse *elOuter = new TEllipse(ringX,ringY,ringOuterA,ringOuterB, 0, 360, dirPhiDeg);
  double ringInnerA = radiusA-1.5;
  double ringInnerB = radiusB-1.5;
  TEllipse *elInner = new TEllipse(ringX,ringY,ringInnerA,ringInnerB, 0, 360, dirPhiDeg);
  TCutG *outerCut = createCutFromEllipse(elOuter);
  TCutG *innerCut = createCutFromEllipse(elInner);

  canvas->cd();
  outerCut->Draw("same");
  innerCut->Draw("same");

  return outerCut->IntegralHist(photonHist) - innerCut->IntegralHist(photonHist);
}


TH2D* generateEvent(TVector3 pos0, TVector3 dir0, double beta) {
  // Generate a single particle event
  Beam *beam = new Beam(pos0, dir0, beta, errX, errY, errDirX, errDirY);
  Particle *pa = beam->generateParticle();

  // Make Aerogel layer
  Aerogel* aerogel1 = new Aerogel(thickness, n1, aeroPos[0], beta);
  Aerogel* aerogel2 = new Aerogel(thickness, n2, aeroPos[1], beta);
  // Make them aware of each other for refraction purposes
  aerogel1->setDownIndex(n2);
  aerogel2->setUpIndex(n1);

  // Make the detector
  Detector* detector = new Detector(dist);

  // Make photons in first aerogel
  std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);
  // Advance particle forward to next aerogel and generate photons
  pa->travelZDist(aeroPos[1] - aeroPos[0]);
  std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);

  aerogel1->applyPhotonScatters(photons);
  aerogel1->applyPhotonScatters(photons);

  // Combine photons from both aerogels
  photons.insert(photons.end(), photons2.begin(), photons2.end());

  // Include scattering in second aerogel
  aerogel2->applyPhotonScatters(photons);


  // Throw out photons based off fill factor
  int numPhotonsDetected = (int) (detector->getFillFactor() * photons.size());
  photons.resize(numPhotonsDetected);

  // Project photons onto detector and plot distribution
  TH2D *photonHist = new TH2D("generatedEvent","generatedEvent",48,-15,15,48,-15,15);
  detector->projectPhotons(photonHist, photons);

  // Draw out photon histogram and ellipse outline
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd();
  TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
  center_pad->Draw();
  TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
  right_pad->Draw();
  TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
  top_pad->Draw();
  gStyle->SetOptStat(0);

  center_pad->cd();
  photonHist->Draw("colz");

  right_pad->cd();
  photonHist->ProjectionY()->Draw("hbar");

  top_pad->cd();
  photonHist->ProjectionX()->Draw("bar");

  photonHist->Draw("samecolz");

  double nPhotons = integrateAndDrawEllipse(pos0, dir0, beta, photonHist, c1, aerogel2);
  cout << "SINGLE EXAMPLE EVENT: Integrated number of photons in ring: " << nPhotons << endl;
  photonHist->SaveAs("./output/generatedEvent.root");
  c1->SaveAs("./output/generatedEvent.pdf");

  return photonHist;
}

TH2D* calculatePdf(TVector3 pos0, TVector3 dir0, double beta) {
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  // Generate beam
  Beam *beam = new Beam(pos0, dir0, beta, errX, errY, errDirX, errDirY);
  beam->makeParticles(nEvents);

  // Make Aerogel layer
  Aerogel* aerogel1 = new Aerogel(thickness, n1, aeroPos[0], beta);
  Aerogel* aerogel2 = new Aerogel(thickness, n2, aeroPos[1], beta);
  // Make them aware of each other for refraction purposes
  aerogel1->setDownIndex(n2);
  aerogel2->setUpIndex(n1);

  // Make the detector
  Detector* detector = new Detector(dist);
  // Get plots ready
  TH2D *photonHist = new TH2D("photonHist","photonHist",48,-15,15,48,-15,15);
  TH1D *scatterHist = new TH1D("scatterHist", "scatterHist", 101, 0, 100);

  TH1D *numPhotonHist = new TH1D("numPhotonHist", "numPhotonHist", 400, 0, 400);
  TH1D *wavHist = new TH1D("wavHist", "wavHist", 200, 250E-9, 700E-9);

  // Make a tree to save the photons
  photonStruct phStruct;
  particleStruct paStruct;
  TFile *f = new TFile("./output/photons.root","RECREATE");
  TTree *tree = new TTree("T","Output photon data");
  TBranch *phBranch = tree->Branch("photons",&phStruct.dirxi,
    "dirxi/D:diryi:dirzi:posxi:posyi:poszi:dirxe:dirye:dirze:posxe:posye:posze:wav:paid/i:numscat");
  TBranch *paBranch = tree->Branch("particles",&paStruct.dirx,"dirx/D:diry:dirz:posx:posy:posz:id/i");

  for (int i = 0; i < nEvents; i++) {
    // Generate particles
    Particle *pa = beam->getParticle(i);
    paStruct.dirx = pa->dir0[0];
    paStruct.diry = pa->dir0[1];
    paStruct.dirz = pa->dir0[2];
    paStruct.posx = pa->pos0[0];
    paStruct.posy = pa->pos0[1];
    paStruct.posz = pa->pos0[2];
    paStruct.id = i;

    // Make photons in first aerogel
    std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);

    // Advance particle forward to next aerogel and generate photons
    pa->travelZDist(aeroPos[1] - aeroPos[0]);
    std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);

    // Scatter photons in first aerogel, move them forwards out of aerogel
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);

    // Combine photons from both aerogels
    photons.insert(photons.end(), photons2.begin(), photons2.end());

    // Include scattering in second aerogel
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);

    // Do some more scattering, for the heck of it
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);

    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      // Save parent particle ID
      phStruct.paid = i;
      // Save initial direction and positions
      TVector3 phPos0 = ph->pos0;
      TVector3 phDir0 = ph->dir0;
      phStruct.dirxi = phDir0[0];
      phStruct.diryi = phDir0[1];
      phStruct.dirzi = phDir0[2];
      phStruct.posxi = phPos0[0];
      phStruct.posyi = phPos0[1];
      phStruct.poszi = phPos0[2];
      // Save exit direction and positions
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;
      phStruct.dirxe = phDir[0];
      phStruct.dirye = phDir[1];
      phStruct.dirze = phDir[2];
      phStruct.posxe = phPos[0];
      phStruct.posye = phPos[1];
      phStruct.posze = phPos[2];

      phStruct.wav = ph->wav;

      // Save number of scattering events
      phStruct.numscat = ph->numScatters;

      tree->Fill();
      scatterHist->Fill(ph->numScatters);
    }
    numPhotonHist->Fill(photons.size());

    // Project photons onto detector and plot distribution
    detector->projectPhotons(photonHist, photons);
  }

  // Scale photon histogram to the number of iterations
  //photonHist->Scale(1. / nEvents);
  // Scale photon to fill factor of detector
  // photonHist->Scale(detector->getFillFactor());


  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>( t2 - t1 ).count();
  cout << "Time taken: " << duration / 1000000. << endl;


  // Draw out photon histogram and ellipse outline
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd();
  TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
  center_pad->Draw();
  TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
  right_pad->Draw();
  TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
  top_pad->Draw();
  gStyle->SetOptStat(0);

  center_pad->cd();
  photonHist->Draw("colz");

  right_pad->cd();
  photonHist->ProjectionY()->Draw("hbar");

  top_pad->cd();
  photonHist->ProjectionX()->Draw("bar");

  photonHist->Draw("samecolz");

  double nPhotons = integrateAndDrawEllipse(pos0, dir0, beta, photonHist, c1, aerogel2);
  cout << "PHOTON DISTRIBUTION: Integrated number of photons in ring: " << nPhotons << endl;

  photonHist->SaveAs("./output/photonHist.root");
  c1->SaveAs("./output/photonHist.pdf");


  TCanvas *c4 = new TCanvas("c4","c4",600,500);
  c4->cd();
  numPhotonHist->Draw();
  c4->SaveAs("./output/numPhotonHist.pdf");

  f->Write();
  return photonHist;
}


int mainPID(int argc, char *argv[]) {
  // Given particle, momentum, simulate example event
  // Assume centered for now

  int particlei = atoi(argv[1]);
  double particleMom = atof(argv[2]);
  double beta = calcBeta(particlei, particleMom);
  TVector3 pos0 = TVector3(0.,0.,0.);
  TVector3 dir0 = TVector3(0.,0.,1.);
  cout << "Particle: " << particleNames[particlei] << endl;
  cout << "Momentum: " << particleMom << " GeV" << endl;
  cout << "Beta: " << beta << endl;
  TH2D* generatedEvent = generateEvent(pos0, dir0, beta);
  for (int i = 0; i < 3; i++) {
    double beta0 = calcBeta(i, particleMom);
    // calculate within 2 standard deviations of each particle hypothesis
    for (int i = -2; i < 3; i++) {
      TH2D *calculatedPdf = calculatePdf(pos0, dir0, beta + i*errBeta);
    }
  }

}

int mainWithBeamParameters(int argc, char *argv[]) {
  // Default beam parameters:
  double beta = 0.999; // velocity of particle
  double dirX_0 = 0.0; // beam initial direction
  double dirY_0 = 0.0;
  double x_0 = -0.0; // beam initial position
  double y_0 = -0.0;
  // Optional changes
  if (argc >= 2) {
    beta = atof(argv[1]);
  }
  if (argc >= 4) {
    dirX_0 = atof(argv[2]);
    dirY_0 = atof(argv[3]);
  }
  if (argc >= 6) {
    x_0 = atof(argv[4]);
    y_0 = atof(argv[5]);
  }
  cout << "Beta: " << beta << endl;
  cout << "X Dir: " << dirX_0 << endl;
  cout << "Y Dir: " << dirY_0 << endl;
  cout << "X Pos: " << x_0 << endl;
  cout << "Y Pos: " << y_0 << endl;
  double dirZ_0 = sqrt(1. - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 pos0 = TVector3(x_0, y_0, 0);
  TVector3 dir0 = TVector3(dirX_0, dirY_0, dirZ_0).Unit();

  // Generate a single candidate event to look at
  //TH2D* eventHist = generateEvent(pos0, dir0, beta);
  // Generate photon PDF of beta hypothesis
  TH2D* pdfHist = calculatePdf(pos0, dir0, beta);

  return 0;
};

int main(int argc, char *argv[]) {
  // Given particle and momentum, simulate centered beam
  int particlei = atoi(argv[1]);
  double particleMom = atof(argv[2]);
  double beta = calcBeta(particlei, particleMom);
  cout << "Particle: " << particleNames[particlei] << endl;
  cout << "Momentum: " << particleMom << " GeV" << endl;
  cout << "Beta: " << beta << endl;
  calculatePdf(TVector3(0.,0.,0.), TVector3(0.,0.,1.), beta);
  return 0;
}
