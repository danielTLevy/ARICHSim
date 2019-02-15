#include "aerogel.h"

const double lowWav = 290E-9;
const double highWav = 700E-9;
const double fineStructConst = 1./137;

Aerogel::Aerogel(double thickness, double refractiveIndex, double zPos, double beta) {
  randomGenerate=std::make_shared<TRandom3>();
  randomGenerate->SetSeed(0);
  this->thickness = thickness;
  this->refractiveIndex = refractiveIndex;
  this->zPos = zPos;
  this->chAngle = calcChAngle(refractiveIndex, beta);
  this->wavPdf = calcWavPdf(refractiveIndex, beta);
  this->scatAngleFunc = new TF1("scatPdf", "(1 + pow(cos(x), 2))", 0, 2*TMath::Pi());
  this->dNdX = calcdNdX(refractiveIndex, beta);
  this->interactionLengths = readInteractionLength(refractiveIndex);
}

double Aerogel::getChAngle() {
  return this->chAngle;
}

double Aerogel::calcChAngle(double n, double beta) {
  return acos(1.0 / (n * beta));
}

TF1* Aerogel::calcWavPdf(double n, double beta) {
  TF1 *pdf = new TF1("pdf", "(1. - 1./([0]*[1]*[0]*[1]))/(x*x)", lowWav, highWav);
  pdf->SetParameter(0, n); pdf->SetParName(0, "n");
  pdf->SetParameter(1, beta); pdf->SetParName(1, "beta");
  return pdf;
}

double Aerogel::calcdNdX(double n, double beta) {
  return 2*TMath::Pi()*fineStructConst*(1. - 1./(n*beta*n*beta))*(1./lowWav - 1./highWav);
}

std::vector<double> Aerogel::readInteractionLength(double n) {
  std::string files[5] = {"leps1-1b","btr4-1a","hds2-3b","leps2-1a","leps6-1a"};
  float ns[5] = {1.0505,1.0452,1.0401,1.0352,1.0297};
  // Get nearest index of refraction
  float minDiff = 1.;
  int minI = 0;
  for (int i = 0; i < 5; i++) {
    if (abs(n - ns[i]) < minDiff) {
      minDiff = abs(n - ns[i]);
      minI = i;
    }
  }
  // Get the corresponding file of interaction lengths
  std::string fileName = "./data/" + files[minI] + "IntLength.csv";
  std::cout << "Reading aerogel data from: " << fileName << std::endl;

  std::ifstream intLengthFile(fileName);

  std::vector<double> wavs;
  std::vector<double> intLengths;
  std::string line;

  while(std::getline(intLengthFile,line, '\n')) {
    std::string wavString = line.substr(0, line.find(' '));
    std::string intLengthString = line.substr(line.find(' ')+1, -1);

    wavs.push_back(atof(wavString.c_str()));
    intLengths.push_back(atof(intLengthString.c_str()));
  }

  return intLengths;
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

double Aerogel::getZPos() {
  return zPos;
}

bool Aerogel::isInAerogel(TVector3 pos) {
  return (TMath::Abs(pos[0]) <= width/2.)
          && (TMath::Abs(pos[1]) <= height/2.)
          && ((pos[2] - zPos) <= thickness)
          && ((pos[2] - zPos) > 0);
}

double Aerogel::getDistInGel(Particle* pa) {
  // calculate how far a particle has remaining in the gel
  return abs((zPos + thickness - pa->pos[2]) / pa->dir[2]);
}

int Aerogel::calcNumPhotons(double particleDist) {
  // N ~= dN/dX * X (in meters)
  return (int) particleDist*0.01*dNdX;
}

double Aerogel::getIntLengthForWav(double wav) {
  // Get index of nearest wavelengths
  double maxWav = 800.;
  double deltaWav = 0.5;
  wav = wav*1E9; // convert to nm
  int index = floor((maxWav - wav)/deltaWav);

  return interactionLengths[index];
}

double Aerogel::getRandomIntDistance(double wav) {
  double intLength = getIntLengthForWav(wav);
  // Our CDF of whether photon has interacted is 1 - exp(- x / interactionL)
  // Sample randomly for value of CDF
  // Invert the equation to get the x value that would get this random value
  double randY = randomGenerate->Uniform(0,1);
  return - intLength * log(1. - randY);

}

double Aerogel::getRandomScatAngle() {
  // Rayleigh scattering is proportional to  1 + cos^2(theta)
  return scatAngleFunc->GetRandom();
}

void Aerogel::applyPhotonScatter(Photon* photon) {
  // Continuously scatter photon while the distance travelled before interacting
  // is less than the distance to exit the gel - update position and direction
  double intDist = getRandomIntDistance(photon->wav);
  double gelDist = getDistInGel(photon);
  TVector3 newPos = photon->pos + intDist*photon->dir;
  TVector3 newDir = photon->dir;
  while (isInAerogel(newPos) && (photon->numScatters<=100)) {
    // update position
    photon->pos = newPos;
    // update direction
    TMatrixD rotMatrix = makeRotationMatrix(photon->dir);
    double scatTheta = getRandomScatAngle();
    double scatPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());
    photon->dir = rotateVector(rotMatrix, scatTheta, scatPhi);

    // Predict where it might scatter next
    photon->numScatters += 1;
    intDist = getRandomIntDistance(photon->wav);
    newPos = photon->pos + intDist*photon->dir;
  }
}

void Aerogel::applyPhotonScatters(std::vector<Photon*> photons) {
  for(int i = 0; i < photons.size(); i++) {
    applyPhotonScatter(photons[i]);
  }
}

std::vector<Photon*> Aerogel::generatePhotons(Particle* pa, Detector* detector) {
  // Create Cherenkov photons as the particle passes through the gel
  // Hacky, but requires detector
  TMatrixD rotMatrix = makeRotationMatrix(pa->dir);
  double paDist = getDistInGel(pa);

  std::vector<Photon*> photons;
  int nPhotons = calcNumPhotons(paDist);
  photons.reserve(nPhotons);
  for (int i = 0; i < nPhotons; i++) {
    // First decide if we want to throw it out yet, based of quantum efficiency
    // Can do this now because the photon's wavelength won't change
    double wav = getRandomWav();
    if (randomGenerate->Uniform(0,1) > detector->evalQEff(wav)) {
      continue;
    } else {
      // Get the point where the photon was emitted
      double phIntDist = randomGenerate->Uniform(paDist);
      TVector3 phPos = pa->pos + phIntDist*pa->dir;
      // Get the direection of the new photon
      double phPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());
      TVector3 dirCR = rotateVector(rotMatrix, chAngle, phPhi);
      // Get the wavelength of the photon

      Photon* photon = new Photon(phPos, dirCR, wav);
      photons.push_back(photon);
    }
  }
  return photons;
}