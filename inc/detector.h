#ifndef DETECTOR_INCLUDE
#define DETECTOR_INCLUDE

//Root
#include "TMath.h"
#include "TGraph.h"
#include "utility.h"
#include "beam.h"
#include "photon.h"
#include "particle.h"

class Detector {
private:
  double zPos;
  TGraph *quantumEff;
  double fillFactor;
  static TGraph *createQEff();
  bool mirror;
  std::shared_ptr<TRandom3> randomGenerate;

public:
  double xmin = -15.;
  double xmax = 15.;
  double ymin = -15.;
  double ymax = 15.;
  double xpixels = 48;
  double ypixels = 48;

  Detector(double zPos, bool mirror = false);
  TH2D* makeDetectorHist(char* name, char* title);
  double evalQEff(double);
  double getFillFactor();
  void projectPhotons(TH2D*,  std::vector<Photon*>, TH1D* rHist = nullptr);

};

#endif