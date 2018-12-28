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
#include "TRandom3.h"
#include "TSpline.h"
#include "beam.h"
#include "photon.h"
#include "utility.h"

class Aerogel {
private:
  double refractiveIndex=0.0;
  double thickness=0.0;
  double dist=0.0;
  Beam* beam;
  double chAngle = 0.0;
  std::shared_ptr<TRandom3> randomGenerate;

  std::vector<double> getRandomTheta(int);
  std::vector<double> getRandomEnergy(int, double, double);

public:
  Aerogel(double, double, double, Beam*, double);
  double getRefractiveIndex();
  double getThickness();
  double getDistance();
  std::vector<Particle*> generatePhotons(Particle* pa);
};

#endif