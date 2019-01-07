#include "aerogel.h"

Aerogel::Aerogel(double thickness, double refractiveIndex, double dist, double beta) {
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  this->thickness = thickness;
  this->refractiveIndex = refractiveIndex;
  this->dist = dist;
  this->chAngle = calcChAngle(refractiveIndex, beta);
  this->wavPdf = calcWavPdf(refractiveIndex, beta);
  this->dNdX = calcdNdX(refractiveIndex, beta);
}


double Aerogel::calcChAngle(double n, double beta) {
  return acos(1.0 / (n * beta));
}

TF1* Aerogel::calcWavPdf(double n, double beta) {
  TF1 *pdf = new TF1("pdf", "(1. - 1./([0]*[1]*[0]*[1]))/(x*x)",
                      300E-9, 700E-9);
  pdf->SetParameter(0, n); pdf->SetParName(0, "n");
  pdf->SetParameter(1, beta); pdf->SetParName(1, "beta");
  return pdf;
}

double Aerogel::calcdNdX(double n, double beta) {
  double alpha = 1./137;
  double lowWav = 300E-9;
  double highWav = 700E-9; // TODO: make these constants!
  return 2*TMath::Pi()*alpha*(1. - 1./(n*beta*n*beta))*(1./lowWav - 1./highWav);
}

double Aerogel::getRandomWav() {
  return wavPdf->GetRandom();
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

double Aerogel::getDistInGel(Particle* pa) {
  // calculate how far a particle has remaining in the gel
  return (thickness - pa->pos[2]) / cos(pa->theta());

}

int Aerogel::calcNumPhotons(double particleDist) {
  // N ~= dN/dX * X (in meters)
  return particleDist*0.01*dNdX;
}

double Aerogel::getRandomIntDistance(double wav) {
  // TODO: use transmission distance to calculate random interaction distance
  // by sampling from exponential function
  return 3;
}

double Aerogel::getRandomScatAngle(double wav) {
  // Rayleigh scattering is proportional to the inverse 4th power of the wavelength
  // and 1 + cos^2(theta)
  TF1 *scatPdf = new TF1("scatPdf", "(1 + pow(cos(x), 2))*pow([0], -4)",
                     0, 2*TMath::Pi());
  scatPdf->SetParameter(0, wav);
  return scatPdf->GetRandom();
}

void Aerogel::applyPhotonScatter(Photon* photon) {
  // continuously scatter photon until it exits gel
  double intLength = getRandomIntDistance(photon->wav);
  double dist = getDistInGel(photon);
  while (intLength < dist) {
    double scatAngle = getRandomScatAngle(photon->wav);
    photon->dir = photon->dir;     // TODO: apply new theta
    double intLength = getRandomIntDistance(photon->wav);
    double dist = getDistInGel(photon);
  }
}

std::vector<Photon*> Aerogel::generatePhotons(Particle* pa) {
  // Create Cherenkov photons as the particle passes through the gel
  TMatrixD rotMatrix = makeRotationMatrix(pa->dir);
  double paDist = getDistInGel(pa);

  std::vector<Photon*> photons;
  int nPhotons = calcNumPhotons(paDist);
  photons.reserve(nPhotons);
  for (int i = 0; i < nPhotons; i++) {
    double phIntDist = randomGenerate->Uniform(paDist);
    TVector3 phPos = pa->pos + phIntDist*pa->dir;
    double phPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());
    TVector3 dirCR = rotateVector(rotMatrix, chAngle, phPhi);
    double wav = getRandomWav();
    Photon* photon = new Photon(phPos, dirCR, wav);
    photons.push_back(photon);
  }

  return photons;

}