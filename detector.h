#ifndef DETECTOR_INCLUDE
#define DETECTOR_INCLUDE

//Root
#include "TGraph.h"
#include "beam.h"
#include "photon.h"
#include "particle.h"

class Detector {
private:
  TGraph *quantumEff;
  static TGraph *createQEff();

public:
  Detector();
  double evalQEff(double);
};

#endif