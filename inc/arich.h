#ifndef ARICH_INCLUDE
#define ARICH_INCLUDE
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

class Arich {

private:
  const int nEvents = 10000;
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

  double integrateAndDrawEllipse(TVector3 pos0, TVector3 dir0, double beta, TH2D* photonHist, TPad* pad, Aerogel* aerogel);

public:
  Arich();
  TH2D* calculatePdf(TVector3 pos0, TVector3 dir0, double beta);
  TH2D* generateEvent(TVector3 pos0, TVector3 dir0, double beta);
  TH2D* simulateBeam(TVector3 pos0, TVector3 dir0, double beta);
};

#endif