#include "stdlib.h"
#include <iostream>
#include <chrono>
#include <TROOT.h>
#include <TStyle.h>
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
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

const int nEvents = 10000; // number of particles simulated for beam
const double aeroPos[2] = {0., 2.0}; // positions of aerogel planes
const double thickness = 2.0; // thickness of aerogel layer
const double width = 10.; // x width of aerogel
const double height = 10.; // y height of aerogel
const double n1 = 1.035; // outer index of refraction
const double n2 = 1.045; // inner index of refraction
const double detectorDist = 21.0; // dist to detector plane
const double errDirX = 0.000; // beam direction error
const double errDirY = 0.000;
const double errX = 0.001; // beam position error
const double errY = 0.001;

const char* particleNames[3] = {"Pion", "Kaon", "Proton"};
const int NUMPARTICLES = 3;
const double particleMasses[3] = {0.1395701, 0.493677, 0.938272};
const double PLANCKSCONST = 4.135667662E-15;
const double SPEEDOFLIGHT = 2.99792458E8;
const double DETECTORFILL = 0.8*0.87;

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

double integrateAndDrawEllipse(TVector3 pos0, TVector3 dir0, double beta, TH2D* photonHist, TPad* pad, Aerogel* aerogel) {
  // Define Ellipse and integrate over this ring
  double dirX_0 = dir0[0];
  double dirY_0 = dir0[1];
  double dirZ_0 = dir0[2];

  // Naively assume that all photons are generated in middle of second aerogel layer.
  double middlePoint = thickness + thickness / 2.;
  double newX_0 = pos0[0] + middlePoint*dirX_0/dirZ_0;
  double newY_0 = pos0[1] + middlePoint*dirY_0/dirZ_0;
  double newDist = detectorDist - middlePoint;
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

  pad->cd();
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
  Detector* detector = new Detector(detectorDist);
  // Make photons in first aerogel
  std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);
  aerogel1->applyPhotonScatters(photons);
  aerogel1->exitAerogel(photons, true);
  // Advance particle forward to next aerogel and generate photons
  pa->travelZDist(aeroPos[1] - aeroPos[0]);
  std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);
  // Combine photons from both aerogels
  photons.insert(photons.end(), photons2.begin(), photons2.end());
  // Include scattering in second aerogel
  aerogel2->applyPhotonScatters(photons);
  aerogel2->exitAerogel(photons, true);
  aerogel1->applyPhotonScatters(photons);
  aerogel1->exitAerogel(photons, true);
  aerogel2->applyPhotonScatters(photons);
  aerogel2->exitAerogel(photons, true);
  // Throw out photons based off fill factor
  int numPhotonsDetected = (int) (detector->getFillFactor() * photons.size());
  photons.resize(numPhotonsDetected);
  // Project photons onto detector and plot distribution
  TH2D *photonHist = new TH2D("generatedEvent","generatedEvent",48,-15,15,48,-15,15);
  detector->projectPhotons(photonHist, photons);
  // Draw out photon histogram and ellipse outline
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  double nPhotons = integrateAndDrawEllipse(pos0, dir0, beta, photonHist, c1, aerogel2);
  cout << "SINGLE EXAMPLE EVENT: Integrated number of photons in ring: " << nPhotons << endl;
  photonHist->SaveAs("./output/generatedEvent.root");
  photonHist->Draw("colz");
  c1->SaveAs("./output/generatedEvent.pdf");
  return photonHist;
}

TH2D* geant4Pdf(TFile* g4File, int particlei) {
  const char* particle = particleNames[particlei];
  char* filename = Form("%sG4Pdf", particle);
  // Prepare values to update in our loop
  TTreeReader reader("h1000", g4File);
  TTreeReaderValue<Int_t> raPid(reader, "Pid");
  TTreeReaderArray<Double_t> raMom(reader, "Mom");
  TTreeReaderArray<Double_t> raPos(reader, "Pos");
  TTreeReaderArray<Double_t> raDir(reader, "Dir");
  Detector* detector = new Detector(0);
  TH2D *allEventHist = new TH2D(filename, filename,48,-15,15,48,-15,15);

  while (reader.Next()) {
      // Get only forward-exiting optical photons
      if (*raPid != 0 || raPos[2] < 829.99) {
          continue;
      }
      // Momentum in GeV/c:
      Double_t mom = sqrt(raMom[0]*raMom[0] + raMom[1]*raMom[1] + raMom[2]*raMom[2]);
      // Wavelength in m (multiply by c and planck's constant in eV*s)
      Double_t wav = SPEEDOFLIGHT * PLANCKSCONST / (1E9 * mom);
      Double_t efficiency = DETECTORFILL * detector->evalQEff(wav);
        // Project onto detector, and convert from mm to cm
      Double_t xFinal = 0.1*(raPos[0] + (1000. - raPos[2])*raDir[0]/raDir[2]);
      Double_t yFinal = 0.1*(raPos[1] + (1000. - raPos[2])*raDir[1]/raDir[2]);
      allEventHist->Fill(xFinal, yFinal, efficiency/10000.);
  }

  allEventHist->SaveAs(Form("./output/geant4pdfs/%s.root", filename));
}

TH2D* calculatePdf(TVector3 pos0, TVector3 dir0, double beta, bool save = false) {
  // Make beam
  Beam *beam = new Beam(pos0, dir0, beta, errX, errY, errDirX, errDirY);
  // Make Aerogel layer
  Aerogel* aerogel1 = new Aerogel(thickness, n1, aeroPos[0], beta);
  Aerogel* aerogel2 = new Aerogel(thickness, n2, aeroPos[1], beta);
  // Make them aware of each other for refraction purposes
  aerogel1->setDownIndex(n2);
  aerogel2->setUpIndex(n1);
  // Make the detector
  Detector* detector = new Detector(detectorDist);


  // Get plots ready
  TH2D *photonHist = new TH2D("photonHist","photonHist",48,-15,15,48,-15,15);
  TH1D *rHist = new TH1D("rHist", "rHist", 500, 0., 10.);

  // Make a tree to save the photons
  photonStruct phStruct;
  particleStruct paStruct;
  TFile *f = new TFile("./output/photons.root","RECREATE");
  TTree *tree = new TTree("T","Output photon data");
  TBranch *phBranch = tree->Branch("photons",&phStruct.dirxi,
    "dirxi/D:diryi:dirzi:posxi:posyi:poszi:dirxe:dirye:dirze:posxe:posye:posze:wav:paid/i:numscat");
  TBranch *paBranch = tree->Branch("particles",&paStruct.dirx,"dirx/D:diry:dirz:posx:posy:posz:id/i");

  // Make events and loop over them
  for (int i = 0; i < nEvents; i++) {
    Particle *pa = beam->generateParticle();
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
    // Do some more scattering
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);

    // Project photons onto detector and plot distribution
    detector->projectPhotons(photonHist, rHist, photons);

    // Save photons to TTree and delete
    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      phStruct.paid = i;
      TVector3 phPos0 = ph->pos0;
      TVector3 phDir0 = ph->dir0;
      phStruct.dirxi = phDir0[0];
      phStruct.diryi = phDir0[1];
      phStruct.dirzi = phDir0[2];
      phStruct.posxi = phPos0[0];
      phStruct.posyi = phPos0[1];
      phStruct.poszi = phPos0[2];
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;
      phStruct.dirxe = phDir[0];
      phStruct.dirye = phDir[1];
      phStruct.dirze = phDir[2];
      phStruct.posxe = phPos[0];
      phStruct.posye = phPos[1];
      phStruct.posze = phPos[2];
      phStruct.wav = ph->wav;
      tree->Fill();
      delete ph;
    }
    delete pa;
  }

  // Scale photon histogram to the number of iterations
  photonHist->Scale(1. / nEvents);
  rHist->Scale(1. / nEvents);

  // Scale photon to fill factor of detector
  photonHist->Scale(detector->getFillFactor());
  rHist->Scale(detector->getFillFactor());

  // Draw out photon histogram and ellipse outline
  if (save) {
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1","c1",900,900);
    c1->cd();
    TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.55,0.55);
    center_pad->Draw();
    TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
    right_pad->Draw();
    TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
    top_pad->Draw();

    TPad *corner_pad = new TPad("corner_pad", "corner_pad",0.55,0.55,1.0,1.0);
    corner_pad->Draw();

    center_pad->cd();
    center_pad->SetGrid(1);
    photonHist->SetTitle("");
    photonHist->SetXTitle("x [cm]");
    photonHist->SetYTitle("y [cm]");
    photonHist->SetZTitle("Mean Photon Count");
    photonHist->Draw("colz");

    double nPhotons = integrateAndDrawEllipse(pos0, dir0, beta, photonHist, center_pad, aerogel2);
    cout << "PHOTON DISTRIBUTION: Integrated number of photons in ring: " << nPhotons << endl;

    right_pad->cd();
    right_pad->SetGrid(1);
    TH1D* yProj = (TH1D*) photonHist->ProjectionY()->Clone("yProj");
    yProj->SetYTitle("y [cm]");
    yProj->Draw("hbar");

    top_pad->cd();
    top_pad->SetGrid(1);
    TH1D* xProj = (TH1D*) photonHist->ProjectionY()->Clone("xProj");
    xProj->SetYTitle("x [cm]");

    xProj->Draw("bar");

    corner_pad->cd();
    corner_pad->SetGrid(1);
    rHist->SetTitle("");
    rHist->SetXTitle("Distance from origin [cm]");
    rHist->Draw("bar");


    c1->SaveAs("./output/photonHistWithProjections.root");
    c1->SaveAs("./output/photonHistWithProjections.pdf");
    rHist->SaveAs("./output/rHist.root");
    photonHist->SaveAs("./output/photonHist.root");

    f->Write();
  }
  delete aerogel1;
  delete aerogel2;
  delete detector;

  return photonHist;
}

double computeLogLikelihood(TH2D* event, TH2D* distribution) {
  int nBins = event->GetSize();
  if (nBins != distribution->GetSize()) {
      throw "Error: Bin Mismatch";
  }
  double logLikelihood = 0.;
  for (int i = 0; i < nBins; i++) {
      double lambda = distribution->GetBinContent(i);
      bool pixelHit = event->GetBinContent(i) > 0;
      if (pixelHit) {
        logLikelihood += log(1 - exp(-lambda));
      } else {
        logLikelihood += log(exp(-lambda));
      }
  }
  return -2*logLikelihood;
}


void calculateSeparation(double particleMom, TVector3 pos0, TVector3 dir0, char* analysisDir) {
  /*
  Compare each particle geant4 root file to a different simulated particle PDF
  Requires an analysisDir in the out directory
    - This contains a directory "g4" which contains all the ttrees for the geant4
  Outputs loglikelihood ratio for each particle pair, and .pdfs of the PDFs
  */

  // First generate 3 particle hypothesis
  TH2D* particlePdfs[NUMPARTICLES];
  for (int i = 0; i < NUMPARTICLES; i++) {
    char* namei = (char*) particleNames[i];
    double massi = particleMasses[i];
    double betai = calcBeta(i, particleMom);
    TH2D* particleiPdf = calculatePdf(pos0, dir0, betai);
    particleiPdf->SetName(Form("%sPdf", namei));
    particleiPdf->SetTitle(Form("%s PDF", namei));
    particlePdfs[i] = particleiPdf;
    TCanvas* pdfCanvas = new TCanvas();
    particleiPdf->Draw("colz");
    pdfCanvas->SaveAs(Form("./output/%s/%sPdfHist.pdf", analysisDir, namei));
    delete pdfCanvas;
  }
  // Next, compare each of the 3 geant4 outputs to these 3 particle hypotheses
  TH1D* likelihoodRatios[NUMPARTICLES][NUMPARTICLES];
  for (int j = 0; j < NUMPARTICLES; j++) {
    for (int i = 0; i < NUMPARTICLES; i++) {
      if (i != j) {
        char* ratioHistName = Form("%s%sRatio%s", (char*)particleNames[max(i,j)],
                                                  (char*)particleNames[min(i,j)],
                                                  (char*)particleNames[j]);
        // Ratio of particle likelihoods for max(i,j) to min(i,j) for a given particle j
        // max and min are used so that the ratio is consistent for both particle types looked at
        likelihoodRatios[i][j] = new TH1D(ratioHistName, ratioHistName, 300, -150, 150);
      }
    }
  }
  TRandom3 randomGen = TRandom3();
  Detector* detector = new Detector(0);
  TH2D *g4EventHist = new TH2D("g4EventHist","g4EventHist",48,-15,15,48,-15,15);
  for (int j = 0; j < NUMPARTICLES; j++) {
    // For each particle type, run through the whole geant4-generated ttree.
    g4EventHist->Reset();
    char* namej = (char*) particleNames[j];
    cout << "Checking likelihoods for " << namej << endl;
    double massj = particleMasses[j];
    double betaj = calcBeta(j, particleMom);
    TFile* g4File = TFile::Open(Form("./output/%s/g4/%s.root", analysisDir, namej));
    // Prepare values to update in our loop
    TTreeReader reader("h1000", g4File);
    TTreeReaderValue<Int_t> rvEvent(reader, "EventNumber");
    TTreeReaderValue<Int_t> rvPid(reader, "Pid");
    TTreeReaderArray<Double_t> raMom(reader, "Mom");
    TTreeReaderArray<Double_t> raPos(reader, "Pos");
    TTreeReaderArray<Double_t> raDir(reader, "Dir");
    int currEventId = 0;
    while (reader.Next()) {
      // Get only forward-exiting optical photons
      if (*rvPid != 0 || raPos[2] < 829.99) {
          continue;
      }
      // After we have all our photons in an event, fill loglikelihood histograms for each particle type
      if (*rvEvent != currEventId) {
        double likelihoods[NUMPARTICLES];
        for (int i = 0; i < NUMPARTICLES; i++) {
          likelihoods[i] = computeLogLikelihood(g4EventHist, particlePdfs[i]);
        }
        for (int i = 0; i < NUMPARTICLES; i++) {
          if (j != i) {
            likelihoodRatios[i][j]->Fill(likelihoods[max(i,j)] - likelihoods[min(i,j)]);
          }
        }
        g4EventHist->Reset();
        currEventId = *rvEvent;
      }
      // Fill our event Histogram
      // Momentum in GeV/c:
      double momMag = sqrt(raMom[0]*raMom[0] + raMom[1]*raMom[1] + raMom[2]*raMom[2]);
      // Wavelength in m
      double wav = SPEEDOFLIGHT * PLANCKSCONST / (1E9 * momMag);
      double efficiency = DETECTORFILL * detector->evalQEff(wav);
      if ( randomGen.Uniform() < efficiency) {
        // Project onto detector, and convert from mm to cm
        double xFinal = 0.1*(raPos[0] + (1000. - raPos[2])*raDir[0]/raDir[2]);
        double yFinal = 0.1*(raPos[1] + (1000. - raPos[2])*raDir[1]/raDir[2]);
        g4EventHist->Fill(xFinal, yFinal);
      }
    }
  }
  // Save the output for each of the ratios.
  TFile *likelihoodFile = new TFile(Form("./output/%s/likelihoods.root", analysisDir), "RECREATE"); 
  likelihoodFile->cd();
  for (int j = 0; j < NUMPARTICLES; j++) {
    for (int i = 0; i < NUMPARTICLES; i++) {
      if (i != j) {
        likelihoodRatios[i][j]->Write();
      }
    }
  }
}


void identifyParticle(int particlei, double particleMom, TVector3 pos0, TVector3 dir0, double errMom = 0.5) {
  // Given particle, momentum, simulate example event
  double realBeta = calcBeta(particlei, particleMom);
  cout << "Real Beta: " << realBeta << endl;
  TH2D* generatedEvent = generateEvent(pos0, dir0, realBeta);
  cout << endl;

  vector<double> betas;
  vector<double> loglikes;
  for (int particleId = 0; particleId < 3; particleId++) {
    // calculate within 2 standard deviations of each particle hypothesis
    cout << "Guess: " << particleNames[particleId] << endl;
    for (int i = -2; i < 3; i++) {
      double betaGuess = calcBeta(particleId, particleMom + i*errMom);
      cout << "Beta: " << betaGuess << endl;
      TH2D *calculatedPdf = calculatePdf(pos0, dir0, betaGuess);
      double logLikelihood = computeLogLikelihood(generatedEvent, calculatedPdf);
      cout << "logLikelihood: " << logLikelihood << endl << endl;
      betas.push_back(betaGuess);
      loglikes.push_back(logLikelihood);
    }
  }
  TGraph* betalikelihoods = new TGraph(betas.size(), &betas[0], &loglikes[0]);
  TCanvas* betagraph = new TCanvas("betacanvas", "betacanvas", 900, 900);
  betalikelihoods->Draw();
  betagraph->SaveAs("./output/BETA.pdf");
}


int main(int argc, char *argv[]) {
  if (argc == 1) {
    cerr << "Usage: " << argv[0] << endl
         << "Make pdf given particle and momentum:  -g <pid> <mom> [xdir ydir xpos ypos]" << endl
         << "Make pdf given beta:                   -b [beta xdir ydir xpos ypos]" << endl
         << "Run particle identification:           -p <pid> <mom> [xdir ydir xpos ypos]" << endl
         << "Make PDF given Geant4 TTree:           -gpdf <g4filename> <pid>" << endl
         << "Check particle separation:             -s <analysisdir> <mom> [xdir ydir xpos ypos]" << endl;
    return -1;
  }

  string mode = string(argv[1]);

  // Set particle parameters
  double beta = 0.999; // velocity of particle
  double dirX_0, dirY_0, x_0, y_0 ;
  dirX_0 = dirY_0 = x_0 = y_0 = 0;
  int particlei = 0;
  double particleMom = 0;
  TFile *g4File = nullptr;
  char* analysisDir;
  bool g4 = false;

  int argi = 2;

  if (mode == "-gpdf") {
    g4File = TFile::Open(argv[argi]);
    argi = argi + 1;
  }
  if (mode == "-s") {
    analysisDir = argv[argi];
    argi = argi + 1;
  }
  if (mode ==  "-g" || mode == "-p" || mode == "-gpdf") {
    particlei = atoi(argv[argi]);
    cout << "Particle: " << particleNames[particlei] << endl;
    argi = argi + 1;
  }
  if (mode ==  "-g" || mode == "-p" || mode == "-s" || mode == "-gpdf" || mode == "-sd") {
    particleMom = atof(argv[argi]);
    cout << "Momentum: " << particleMom << " GeV" << endl;
    argi = argi + 1;
  }
  if (mode == "-b" && argc > 2) {
    beta = atof(argv[argi]);
    argi = argi + 1;
  }
  if (argc >= argi + 1) {
    dirX_0 = atof(argv[argi]);
    dirY_0 = atof(argv[argi + 1]);
    argi = argi + 2;
  }
  if (argc >= argi + 1) {
    x_0 = atof(argv[argi]);
    y_0 = atof(argv[argi + 1]);
    argi = argi + 2;
  }
  cout << "X Dir: " << dirX_0 << endl;
  cout << "Y Dir: " << dirY_0 << endl;
  cout << "X Pos: " << x_0 << endl;
  cout << "Y Pos: " << y_0 << endl;
  double dirZ_0 = sqrt(1. - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 pos0 = TVector3(x_0, y_0, 0);
  TVector3 dir0 = TVector3(dirX_0, dirY_0, dirZ_0).Unit();

  // Do the thing
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  if (mode == "-g") {
    // Given particle and momentum, simulate centered beam à la Geant4
    beta = calcBeta(particlei, particleMom);
    cout << "Beta: " << beta << endl;
    calculatePdf(pos0, dir0, beta, true);
  }
  if (mode == "-b") {
    calculatePdf(pos0, dir0, beta, true);
  }

  if (mode == "-p") {
    identifyParticle(particlei, particleMom, pos0,  dir0);
  }

  if (mode == "-s") {
    calculateSeparation(particleMom, pos0, dir0, analysisDir);
  }

  if (mode == "-gpdf") {
    geant4Pdf(g4File, particlei);
  }

  auto duration = duration_cast<microseconds>( high_resolution_clock::now() - t1 ).count();
  cout << "Time taken: " << duration / 1000000. << endl;
  return 0;
}
