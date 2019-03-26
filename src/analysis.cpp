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
#include "THStack.h"
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
#include "detector.h"
#include "arich.h"

using namespace std;
using namespace std::chrono;

const char* pNames[3] = {"Pion", "Kaon", "Proton"};
const int NUMPARTICLES = 3;
const double particleMasses[3] = {0.1395701, 0.493677, 0.938272};
const double PLANCKSCONST = 4.135667662E-15;
const double SPEEDOFLIGHT = 2.99792458E8;
const double DETECTORFILL = 0.8*0.87;


double calcBeta(int particlei, double mom) {
  double M = particleMasses[particlei];
  return sqrt(1/(1 + M*M/(mom*mom)));
}


TH2D* geant4Pdf(TFile* g4File, int particlei) {
  const char* particle = pNames[particlei];
  char* filename = Form("%sG4Pdf", particle);
  // Prepare values to update in our loop
  TTreeReader reader("h1000", g4File);
  TTreeReaderValue<Int_t> raPid(reader, "Pid");
  TTreeReaderArray<Double_t> raMom(reader, "Mom");
  TTreeReaderArray<Double_t> raPos(reader, "Pos");
  TTreeReaderArray<Double_t> raDir(reader, "Dir");
  Detector* detector = new Detector(0);
  TH2D *allEventHist = new TH2D(filename, filename,200,-20,20,200,-20,20);

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
  return allEventHist;
}



double computeLogLikelihood(TH2D* event, TH2D* distribution) {
  /*
  Compare every bin of an event histogram to a bin in some probability distribution,
  to give the likelihood of that event under that distribution
  */
  int nBins = event->GetSize();
  if (nBins != distribution->GetSize()) {
      cerr << "ERROR: Bin Mismatch" << endl;
      throw "Bin Mismatch Error";
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

void calculateAllLoglikes(double particleMom, TVector3 pos0, TVector3 dir0, char* analysisDir) {
  /*
  Compare each particle geant4 root file to a different simulated particle PDF
  Requires an analysisDir, containing g4 ttrees for each particle
  Outputs each loglikelihood for each particle
  */

  double mom = particleMom;
  double xdir = dir0[0];
  double ydir = dir0[1];
  double xpos = pos0[0];
  double ypos = pos0[1];
  int particle;
  double piloglike;
  double kloglike;
  double ploglike;

  TString filename = Form("%sloglikes.root", analysisDir);
  TFile *hfile = 0;
  TTree *tree = new TTree("T","loglikes");
  tree->Branch("mom",&mom,"mom/D");
  tree->Branch("xdir",&xdir,"xdir/D");
  tree->Branch("ydir",&ydir,"ydir/D");
  tree->Branch("xpos",&xpos,"xpos/D");
  tree->Branch("ypos",&ypos,"ypos/D");
  tree->Branch("particle", &particle, "particle/I");
  tree->Branch("piloglike",&piloglike,"piloglike/D");
  tree->Branch("kloglike",&kloglike,"kloglike/D");
  tree->Branch("ploglike",&ploglike,"ploglike/D");



  // First generate 3 particle hypothesis
  TH2D* particlePdfs[NUMPARTICLES];
  for (int i = 0; i < NUMPARTICLES; i++) {
    char* namei = (char*) pNames[i];
    double betai = calcBeta(i, particleMom);
    TH2D* particleiPdf = Arich::calculatePdf(pos0, dir0, betai);
    particleiPdf->SetName(Form("%sPdf", namei));
    particleiPdf->SetTitle(Form("%s PDF", namei));
    particlePdfs[i] = particleiPdf;
    particleiPdf->SaveAs(Form("%s%sPdfHist.root", analysisDir, namei));
  }

  // Next, compare each of the 3 geant4 outputs to these 3 particle hypotheses
  TRandom3 randomGen = TRandom3();
  Detector* detector = new Detector(0);

  TH2D *g4EventHist = new TH2D("g4EventHist","g4EventHist",48,-15,15,48,-15,15);
  for (int j = 0; j < NUMPARTICLES; j++) {
    particle = j;

    // For each particle type, run through the whole geant4-generated ttree.
    g4EventHist->Reset();
    char* namej = (char*) pNames[j];
    TFile* g4File = TFile::Open(Form("%s%s.root", analysisDir, namej));
    char* g4PdfName = Form("g4%sPdf", namej);
    TH2D *g4Pdf = new TH2D(g4PdfName,g4PdfName,48,-15,15,48,-15,15);
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
        piloglike  = computeLogLikelihood(g4EventHist, particlePdfs[0]);
        kloglike  = computeLogLikelihood(g4EventHist, particlePdfs[1]);
        ploglike  = computeLogLikelihood(g4EventHist, particlePdfs[2]);
        tree->Fill();
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
        g4Pdf->Fill(xFinal, yFinal);
      }
    }
    g4Pdf->SaveAs(Form("%s%s.root", analysisDir, g4PdfName));
  }
  // Save the likelihoods
  hfile = TFile::Open(filename,"RECREATE");
  tree->Write();
}



void calculateSeparationHists(double particleMom, TVector3 pos0, TVector3 dir0, char* analysisDir) {
  /*
  Compare each particle geant4 root file to a different simulated particle PDF
  //Requires an analysisDir in the out directory
    - This contains a directory "g4" which contains all the ttrees for the geant4
  Outputs loglikelihood ratio for each particle pair, and .pdfs of the PDFs
  */

  // First generate 3 particle hypothesis
  TH2D* particlePdfs[NUMPARTICLES];
  for (int i = 0; i < NUMPARTICLES; i++) {
    char* namei = (char*) pNames[i];
    double betai = calcBeta(i, particleMom);
    TH2D* particleiPdf = Arich::calculatePdf(pos0, dir0, betai);
    particleiPdf->SetName(Form("%sPdf", namei));
    particleiPdf->SetTitle(Form("%s PDF", namei));
    particlePdfs[i] = particleiPdf;
    particleiPdf->SaveAs(Form("%s%sPdfHist.root", analysisDir, namei));
  }

  // Next, compare each of the 3 geant4 outputs to these 3 particle hypotheses
  TH1D* likelihoodRatios[NUMPARTICLES][NUMPARTICLES];
  for (int j = 0; j < NUMPARTICLES; j++) {
    for (int i = 0; i < NUMPARTICLES; i++) {
      if (i != j) {
        // Ratio of particle likelihoods for max(i,j) to min(i,j) for a given particle j
        // max and min are used so that the ratio is consistent for both particle types looked at
        char* histName = Form("%s%sRatio%s",
                                (char*)pNames[max(i,j)],
                                (char*)pNames[min(i,j)],
                                (char*) pNames[j]);
        char* ratioName = Form("%s/%s Ratio", (char*)pNames[max(i,j)],
                                              (char*)pNames[min(i,j)]);
        char* histTitle = Form("%s, %ss", ratioName, (char*) pNames[j]);
        likelihoodRatios[i][j] = new TH1D(histName, histTitle, 300, -150, 150);
        likelihoodRatios[i][j]->SetXTitle(ratioName);
        likelihoodRatios[i][j]->SetYTitle("Count");
      }
    }
  }
  TRandom3 randomGen = TRandom3();
  Detector* detector = new Detector(0);
  TH2D *g4EventHist = new TH2D("g4EventHist","g4EventHist",48,-15,15,48,-15,15);
  for (int j = 0; j < NUMPARTICLES; j++) {
    // For each particle type, run through the whole geant4-generated ttree.
    g4EventHist->Reset();
    char* namej = (char*) pNames[j];
    //cout << "Checking likelihoods for " << namej << endl;
    TFile* g4File = TFile::Open(Form("%s%s.root", analysisDir, namej));
    char* g4PdfName = Form("g4%sPdf", namej);
    TH2D *g4Pdf = new TH2D(g4PdfName,g4PdfName,48,-15,15,48,-15,15);
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
        g4Pdf->Fill(xFinal, yFinal);
      }
    }
    g4Pdf->SaveAs(Form("%s%s.root", analysisDir, g4PdfName));
  }
  // Save the output for each of the ratios.
  TFile *likelihoodFile = new TFile(Form("%s/likelihoods.root", analysisDir), "RECREATE"); 
  likelihoodFile->cd();
  for (int j = 0; j < NUMPARTICLES; j++) {
    for (int i = 0; i < NUMPARTICLES; i++) {
      if (i != j) {
        likelihoodRatios[i][j]->Write();
        if (i > j) {
          cout << pNames[i] << pNames[j] << ": " << endl;
          double deltaMean = likelihoodRatios[i][j]->GetMean() - likelihoodRatios[j][i]->GetMean();
          double width = sqrt(TMath::Power(likelihoodRatios[i][j]->GetRMS(),2) +
                            TMath::Power(likelihoodRatios[j][i]->GetRMS(),2));
          cout << "Separation: " << deltaMean / width << endl;
          double pctjMisidentified = 100.*likelihoodRatios[i][j]->Integral(0,   149)/likelihoodRatios[i][j]->GetSum();
          double pctiMisidentified = 100.*likelihoodRatios[j][i]->Integral(150, 300)/likelihoodRatios[i][j]->GetSum();
          cout << pNames[i] << "ErrPct: " << pctiMisidentified << endl;
          cout << pNames[j] << "ErrPct: " << pctjMisidentified << endl;
        }
      }
    }
  }
}


void identifyParticle(TH2D* eventHist, double particleMom, TVector3 pos0, TVector3 dir0, double errMom = 0.5) {
  vector<double> betas;
  vector<double> loglikes;
  for (int particleId = 0; particleId < NUMPARTICLES; particleId++) {
    // calculate within 2 standard deviations of each particle hypothesis
    cout << "Guess: " << pNames[particleId] << endl;
    for (int i = -2; i < 3; i++) {
      double betaGuess = calcBeta(particleId, particleMom + i*errMom);
      cout << "Beta: " << betaGuess << endl;
      TH2D *calculatedPdf = Arich::calculatePdf(pos0, dir0, betaGuess);
      double logLikelihood = computeLogLikelihood(eventHist, calculatedPdf);
      delete calculatedPdf;
      cout << "logLikelihood: " << logLikelihood << endl << endl;
      betas.push_back(betaGuess);
      loglikes.push_back(logLikelihood);
    }
  }
  TGraph(betas.size(), &betas[0], &loglikes[0]).SaveAs("./output/Beta.root");
}


void testIdentifyParticle(int particlei, double particleMom, TVector3 pos0, TVector3 dir0, double errMom = 0.5) {
  // Given particle, momentum, simulate example event, and see if we can identify it
  double realBeta = calcBeta(particlei, particleMom);
  cout << "Real Beta: " << realBeta << endl;
  TH2D* generatedEvent = Arich::generateEvent(pos0, dir0, realBeta);
  identifyParticle(generatedEvent, particleMom, pos0, dir0);
}


void identifyMultiParticle(TH2D* eventHist, int nParticles, vector<int> particleis, vector<double> particleMoms,
                           vector<TVector3> pos0s, vector<TVector3> dir0s, double errMom=0.5) {
  // Run multidimensional particle identification 
  THStack *hs = new THStack("pdfStack","");
  for (int i = 0; i < nParticles; i++) {
    double realBeta = calcBeta(particleis[i], particleMoms[i]);
    hs->Add(Arich::calculatePdf(pos0s[i], dir0s[i], realBeta));
  }
  hs->GetStack()->Last()->SaveAs("./output/stackedPdfs.root");

}

void testIdentifyMultiParticle(int nParticles) {
  // Simulate multiparticle event, and test our ability to identify it
  THStack histStack("eventStack","");
  vector<int> particles;
  vector<TVector3> pos0s;
  vector<TVector3> dir0s;
  vector<double> moms;
  TRandom3 randomGen = TRandom3();
  for (int i = 0; i < nParticles; i++) {
    int particlei = randomGen.Integer(NUMPARTICLES);
    double momentumi = randomGen.Gaus(10, 2);
    double betai = calcBeta(particlei, momentumi);
    TVector3 pos0i;
    TVector3 dir0i;
    for (int d = 0; d < 2; d++) {
      pos0i[d] = min(5.,max(-5.,(randomGen.Gaus(0., 2.))));
      dir0i[d] = min(0.6,max(-0.6, randomGen.Gaus(0,0.3)));
    }
    dir0i[2] = sqrt(1-dir0i[0]*dir0i[0]-dir0i[1]*dir0i[0]);
    pos0i[2] = 0;
    histStack.Add(Arich::generateEvent(pos0i, dir0i, betai, false));
    particles.push_back(particlei);
    pos0s.push_back(pos0i);
    dir0s.push_back(dir0i);
    moms.push_back(momentumi);
  }
  TH2D* eventHist = (TH2D*) histStack.GetStack()->Last();
  eventHist->SaveAs("./output/stackedEvents.root");
  identifyMultiParticle(eventHist, nParticles, particles, moms, pos0s, dir0s);

}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    cerr << "Usage: " << argv[0] << endl
         << "Make pdf given particle and momentum:  -g <pid> <mom> [xdir ydir xpos ypos]" << endl
         << "Make pdf given beta:                   -b [beta xdir ydir xpos ypos]" << endl
         << "Run particle identification:           -p <pid> <mom> [xdir ydir xpos ypos]" << endl
         << "Run multi-particle identification:     -mp <nparticles>" << endl
         << "Make PDF given Geant4 TTree:           -gpdf <g4filename> <pid>" << endl
         << "Check loglikes for each particle:      -s <analysisdir> <mom> [xdir ydir xpos ypos]" << endl;
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
  TVector3 pos0;
  TVector3 dir0;
  int nParticles = 0;
  int argi = 2;

  if (mode == "-mp") {
    nParticles = atoi(argv[argi]);
    argi = argi + 1;
  }
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
    cout << "Particle: " << pNames[particlei] << endl;
    argi = argi + 1;
  }
  if (mode ==  "-g" || mode == "-p" || mode == "-s" || mode == "-sd") {
    particleMom = atof(argv[argi]);
    cout << "Mom: " << particleMom << endl;
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
  if (mode != "-mp") {
    cout << "XDir: " << dirX_0 << endl;
    cout << "YDir: " << dirY_0 << endl;
    cout << "XPos: " << x_0 << endl;
    cout << "YPos: " << y_0 << endl;
    double dirZ_0 = sqrt(1. - dirX_0*dirX_0 - dirY_0*dirY_0);
    pos0 = TVector3(x_0, y_0, 0);
    dir0 = TVector3(dirX_0, dirY_0, dirZ_0).Unit();    
  }

  // Do the thing
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  if (mode == "-mp") {
    testIdentifyMultiParticle(nParticles);
  }
  if (mode == "-g" || mode == "-b") {
    // Given particle and momentum, simulate centered beam à la Geant4
    if (mode == "-g") {
      beta = calcBeta(particlei, particleMom);
      cout << "Beta: " << beta << endl;
    }
    Arich::simulateBeam(pos0, dir0, beta);
  }

  if (mode == "-p") {
    testIdentifyParticle(particlei, particleMom, pos0,  dir0);
  }

  if (mode == "-s") {
    calculateAllLoglikes(particleMom, pos0, dir0, analysisDir);
  }

  if (mode == "-gpdf") {
    geant4Pdf(g4File, particlei);
  }

  auto duration = duration_cast<microseconds>( high_resolution_clock::now() - t1 ).count();
  cout << "Time taken: " << duration / 1000000. << endl;
  return 0;
}
