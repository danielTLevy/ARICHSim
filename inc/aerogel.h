#ifndef AEROGEL_INCLUDE
#define AEROGEL_INCLUDE

#include <fstream>

//Root
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "beam.h"
#include "photon.h"
#include "particle.h"
#include "detector.h"
#include "utility.h"

class Aerogel {
private:
  double refractiveIndex=1.0;
  double thickness=0.0;
  double height=10.0;
  double width=10.0;
  double zPos=0.0;
  double chAngle = 0.0;
  double upIndex = 1.0; //Upstream refractive index
  double downIndex = 1.0; //Downstream refractive
  TF1 *wavPdf;
  TF1 *scatAngleFunc;
  std::shared_ptr<TRandom3> randomGenerate;
  double dNdX;
  std::vector<double> interactionLengths;

  static std::vector<double> readInteractionLength(double);
  static double calcChAngle(double, double);
  std::vector<double> getRandomTheta(int);
  std::vector<double> getRandomEnergy(int, double, double);
  static double calcdNdX(double, double);
  double getRandomWav();
  double getRandomScatAngle();
  void refractPhoton(Photon*);
  void applyPhotonScatter(Photon*);
  double getIntLengthForWav(double);
  double getRandomIntDistance(double);

public:
  Aerogel(double, double, double, double);
  int calcNumPhotons(double);
  double getChAngle();
  double getRefractiveIndex();
  double getThickness();
  double getZPos();
  void setUpIndex(double);
  void setDownIndex(double);
  void exitAerogel(Photon*,  bool refract = true);
  void exitAerogel(std::vector<Photon*>, bool refract = true);
  bool isInAerogel(TVector3);
  void applyPhotonScatters(std::vector<Photon*>);
  std::vector<Photon*> generatePhotons(Particle*, Detector*);
  double getDistInGel(Particle*);

};

#endif