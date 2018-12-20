#include "aerogel.h"

Aerogel::Aerogel(double thickness, double refractiveIndex, double dist) {
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  this->thickness = thickness;
  this->refractiveIndex = refractiveIndex;
  this->dist = dist;
}


double Aerogel::getThickness() {
  return thickness;
}

double Aerogel::getRefractiveIndex() {
  return refractiveIndex;
}


double Aerogel::getDistance() {
  return dist;
}

double Aerogel::getCherenkovAngle(double beta) {
  return acos(1.0 / (refractiveIndex * beta));
}
