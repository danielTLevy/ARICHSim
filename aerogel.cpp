#include "aerogel.h"


Aerogel::Aerogel(double thickness, double refractiveIndex, double dist, double beta) {
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  this->thickness = thickness;
  this->refractiveIndex = refractiveIndex;
  this->dist = dist;
  this->chAngle = calcChAngle(refractiveIndex, beta);
  this->wavPdf = calcWavPdf(refractiveIndex, beta);
  this->scatAngleFunc = new TF1("scatPdf", "(1 + pow(cos(x), 2))", 0, 2*TMath::Pi());
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
  return abs((thickness - pa->pos[2]) / cos(pa->theta()));

}

int Aerogel::calcNumPhotons(double particleDist) {
  // N ~= dN/dX * X (in meters)
  return particleDist*0.01*dNdX;
}

double Aerogel::getRandomIntDistance(double wav) {
  // TODO: use transmission distance to calculate random interaction distance
  // by sampling from exponential function
  return randomGenerate->Uniform(0.,4.);
}

double Aerogel::getRandomScatAngle() {
  // Rayleigh scattering is proportional to  1 + cos^2(theta)
  return scatAngleFunc->GetRandom();
}

void Aerogel::applyPhotonScatter(Photon* photon) {
  // Continuously scatter photon while the distance travelled before interacting
  // is less than the distance to exit the gel
  double intDist = getRandomIntDistance(photon->wav);
  double gelDist = getDistInGel(photon);

  while ((intDist < gelDist) && (photon->numScatters <= 100)) {
    // update position
    photon->pos = photon->pos + intDist*photon->dir;
    // update direction
    TMatrixD rotMatrix = makeRotationMatrix(photon->dir);
    double scatTheta = getRandomScatAngle();
    double scatPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());
    photon->dir = rotateVector(rotMatrix, scatTheta, scatPhi);

    photon->numScatters += 1;
    intDist = getRandomIntDistance(photon->wav);
    gelDist = getDistInGel(photon);
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
    // Get the point where the photon was emitted
    double phIntDist = randomGenerate->Uniform(paDist);
    TVector3 phPos = pa->pos + phIntDist*pa->dir;
    // Get the direection of the new photon
    double phPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());
    TVector3 dirCR = rotateVector(rotMatrix, chAngle, phPhi);
    // Get the wavelength of the photon
    double wav = getRandomWav();

    Photon* photon = new Photon(phPos, dirCR, wav);

    applyPhotonScatter(photon);
    photons.push_back(photon);
  }

  return photons;

}