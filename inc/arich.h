#ifndef ARICH_INCLUDE
#define ARICH_INCLUDE
#include "stdlib.h"
#include <iostream>
#include <chrono>
#include <TROOT.h>
#include <TStyle.h>
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
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

  struct particleInfoStruct {
    // Input parameters to a simple fast simulation
    TVector3 pos;
    TVector3 dir;
    double beta;
  };

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
    double beta;
    int id;
  };

class Arich {

private:
  static constexpr double aeroPos1 = 0.; // positions of upstream edges of aerogel layers
  static constexpr double aeroPos2 = 2.;
  static constexpr double thickness1 = 2.0; // thickness of aerogel layers
  static constexpr double thickness2 = 2.0;
  static constexpr double width = 10.; // x width of aerogel
  static constexpr double height = 10.; // y height of aerogel
  static constexpr double n1 = 1.0352; // outer index of refraction
  static constexpr double n2 = 1.0452; // inner index of refraction
  static constexpr double detectorDist = 21.0; // dist to detector plane

  Aerogel* aerogel1;
  Aerogel* aerogel2;
  Detector* detector;

  double integrateAndDrawEllipse(particleInfoStruct params, TH2D* photonHist, TPad* pad);

public:
  Arich(bool mirror = false);
  TH2D* calculatePdf(particleInfoStruct params, char* histName="photonHist");
  TH2D* generateEvent(particleInfoStruct params, bool save=true, char* histName="generatedEvent");
  TH2D* simulateBeam(particleInfoStruct params, char* outputDir="simulatedBeam");
};

#endif