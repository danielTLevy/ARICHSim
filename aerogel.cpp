#include "aerogel.h"

Aerogel::Aerogel(double thickness, double refractiveIndex, double dist, Beam* beam, double beta) {
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  this->thickness = thickness;
  this->refractiveIndex = refractiveIndex;
  this->dist = dist;
  this->beam = beam;
  this->chAngle = acos(1.0 / (refractiveIndex * beta));
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

std::vector<Particle*> Aerogel::generatePhotons(Particle* pa) {
  TMatrixD rot = makeRotationMatrix(pa->dir);
  double paTheta = pa->theta();
  double paPhi = pa->phi();
  double r = dist / cos(paTheta);
  double paDist = thickness / cos(paTheta);

  std::vector<Particle*> photons;
  int nPhotons = 10;
  photons.reserve(nPhotons);

  for (int i = 0; i < nPhotons; i++) {

    double phIntDist = randomGenerate->Uniform(paDist);
    double phIntR = r - phIntDist;

    double phX = pa->pos[0] + phIntDist*sin(paTheta)*cos(paPhi);
    double phY = pa->pos[1] + phIntDist*sin(paTheta)*sin(paPhi);
    TVector3 phPos;
    phPos[0] = pa->pos[0] + phIntDist*sin(paTheta)*cos(paPhi);
    phPos[1] = pa->pos[1] + phIntDist*sin(paTheta)*sin(paPhi);
    phPos[2] = phIntR * cos(paTheta);
    // Photons
    double phPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());

    TVector3 dirC;
    TVector3 dirCR;
    // Beam frame angle of photon
    dirC.SetX(cos(phPhi)*sin(chAngle));
    dirC.SetY(sin(phPhi)*sin(chAngle));
    dirC.SetZ(cos(chAngle));

    // photon direction rotated onto particle direction
    //TVector3 dirCR = rot * dirC;
    dirCR[0] = rot[0][0]*dirC[0]+rot[0][1]*dirC[1]+rot[0][2]*dirC[2];
    dirCR[1] = rot[1][0]*dirC[0]+rot[1][1]*dirC[1]+rot[1][2]*dirC[2];
    dirCR[2] = rot[2][0]*dirC[0]+rot[2][1]*dirC[1]+rot[2][2]*dirC[2];

    Particle* photon = new Particle(phPos, dirCR, 1);

    photons.push_back(photon);
  }
  return photons;

}