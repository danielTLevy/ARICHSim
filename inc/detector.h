#ifndef DETECTOR_INCLUDE
#define DETECTOR_INCLUDE

//Root
#include "TMath.h"
#include "TGraph.h"
#include "beam.h"
#include "photon.h"
#include "particle.h"

class Detector {
private:
  double zPos;
  TGraph *quantumEff;
  double fillFactor;
  static TGraph *createQEff();
  std::shared_ptr<TRandom3> randomGenerate;

public:
  Detector(double);
  double evalQEff(double);
  double getFillFactor();
  void projectPhotons(TH2D*,  std::vector<Photon*> );
  void projectPhotons(TH2D*, TH1D*, std::vector<Photon*> );

};

#endif