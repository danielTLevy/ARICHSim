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

class Aerogel {
private:
  double refractiveIndex=0.0;
  double thickness=0.0;
  double dist=0.0;
  std::shared_ptr<TRandom3> randomGenerate;

  std::vector<double> getRandomTheta(int);
  std::vector<double> getRandomEnergy(int, double, double);
  double getCherenkovAngle(double);

public:
  Aerogel(double, double, double);
  double getRefractiveIndex();
  double getThickness();
  double getDistance();
};
