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

int generateEvent(int argc, char *argv[]) {
  double aeroPos[2] = {0., 2.0}; // positions of aerogel planes
  double thickness = 2.0; // thickness of aerogel layer
  double n1 = 1.035; // outer index of refraction
  double n2 = 1.045; // inner index of refraction
  double dist = 21.0; // dist to detector plane

  // Beam parameters:
  double beta = 0.999; // velocity of particle

  double dirX_0 = 0.2; // beam initial direction
  double dirY_0 = -0.00;
  double x_0 = -0.0; // beam initial position
  double y_0 = -0.0;
  double errDirX = 0.0001; // beam direction error
  double errDirY = 0.0001;
  double errX = 0.0; // beam position error
  double errY = 0.0;

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
  if (argc >= 10) {
    errDirX =atof(argv[6]);
    errDirY = atof(argv[7]);
    errX = atof(argv[8]);
    errY = atof(argv[9]);
  }
  cout << "Beta: " << beta << endl;
  cout << "X Dir: " << dirX_0 << endl;
  cout << "Y Dir: " << dirY_0 << endl;
  cout << "X Pos: " << x_0 << endl;
  cout << "Y Pos: " << y_0 << endl;
  cout << "X Dir Err: " << errDirX << endl;
  cout << "Y Dir Err: " << errDirY << endl;
  cout << "X Pos Err: " << errX << endl;
  cout << "Y Pos Err: " << errY << endl;

  double dirZ_0 = sqrt(1 - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 pos_0 = TVector3(x_0, y_0, 0);
  TVector3 dir_0 = TVector3(dirX_0, dirY_0, dirZ_0);

  // Generate beam
  Beam *beam = new Beam(pos_0, dir_0, beta, errX, errY, errDirX, errDirY);
  Particle *pa = beam->generateParticle();

  // Make Aerogel layer
  Aerogel* aerogel1 = new Aerogel(thickness, n1, aeroPos[0], beta);
  Aerogel* aerogel2 = new Aerogel(thickness, n2, aeroPos[1], beta);

  // Make the detector
  Detector* detector = new Detector(dist);

  // Make photons in first aerogel and scatter them
  std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);
  aerogel1->applyPhotonScatters(photons);
  // Advance particle forward to next aerogel and generate photons
  pa->travelZDist(aeroPos[1] - aeroPos[0]);
  std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);

  // Combine photons from both aerogels
  photons.insert(photons.end(), photons2.begin(), photons2.end());

  // Include scattering in second aerogel
  aerogel2->applyPhotonScatters(photons);

  // Throw out photons based off fill factor

  int numPhotonsDetected = (int) (detector->getFillFactor() * photons.size());
  photons.resize(numPhotonsDetected); 

  // Project photons onto detector and plot distribution
  // Get plot ready
  TH2D *photonHist = new TH2D("generatedEvent","generatedEvent",48,-15,15,48,-15,15);
  detector->projectPhotons(photonHist, photons);


  //gStyle->SetOptStat(0);
  // Define Ellipse and integrate over this ring
  // Naively assume that all photons are generated in middle of second aerogel layer.
  double middlePoint = thickness + thickness / 2.;
  double newX_0 = x_0 + middlePoint*dirX_0/dirZ_0;
  double newY_0 = y_0 + middlePoint*dirY_0/dirZ_0;
  double newDist = dist - middlePoint;
  double dirTheta =  atan(sqrt(dirX_0*dirX_0 +dirY_0*dirY_0) / dirZ_0);
  double chAngle = aerogel2->getChAngle();
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

  // Draw out photon histogram and ellipse outline
  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  c1->cd();
  photonHist->Draw("samecolz");
  outerCut->Draw("same");
  innerCut->Draw("same");
  double nPhotons = outerCut->IntegralHist(photonHist) - innerCut->IntegralHist(photonHist);
  cout << "Integrated number of photons in ring: " << nPhotons << endl;
  photonHist->SaveAs("./output/generatedEvent.root");
  c1->SaveAs("./output/generatedEvent.pdf");


  return 0;
};
