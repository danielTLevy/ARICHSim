#ifndef AEROGEL_INCLUDE
#define AEROGEL_INCLUDE

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>
#include <memory>

//Root
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "beam.h"
#include "photon.h"
#include "particle.h"
#include "utility.h"

class Aerogel {
private:
  double refractiveIndex=0.0;
  double thickness=0.0;
  double dist=0.0;
  Beam* beam;
  double chAngle = 0.0;
  TF1 *wavPdf;
  std::shared_ptr<TRandom3> randomGenerate;
  static double calcChAngle(double, double);
  std::vector<double> getRandomTheta(int);
  std::vector<double> getRandomEnergy(int, double, double);
  static TF1* calcWavPdf(double, double);
  double getRandomWav();

public:
  Aerogel(double, double, double, Beam*, double);
  double getRefractiveIndex();
  double getThickness();
  double getDistance();
  std::vector<Photon*> generatePhotons(Particle* pa);
};

#endif