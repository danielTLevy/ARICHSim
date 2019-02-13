#include "detector.h"

Detector::Detector(double zPos) {
  this->zPos = zPos;
  this->quantumEff = createQEff();
  this->randomGenerate=std::make_shared<TRandom3>();
  this->randomGenerate->SetSeed(0);
  this->fillFactor = 0.8*0.87;
};

TGraph* Detector::createQEff() {
  TGraph *qe = new TGraph();
  qe->SetPoint(0,267.60597413356646, 13.28179513808478);
  qe->SetPoint(1,275.5511575459106, 18.57988946896094);
  qe->SetPoint(2,284.8598517973029, 24.01437582097361);
  qe->SetPoint(3,298.2015464629263, 29.830659978841663);
  qe->SetPoint(4,315.59925438394623, 32.90657208607398);
  qe->SetPoint(5,337.04146913977996, 33.53499361680753);
  qe->SetPoint(6,357.1546923183136, 32.85256215587456);
  qe->SetPoint(7,373.2464215031987, 32.18913466589778);
  qe->SetPoint(8,398.7331430938463, 30.305446608103644);
  qe->SetPoint(9,420.20412390113677, 27.97804240426051);
  qe->SetPoint(10,436.31311271689594, 25.834853963434693);
  qe->SetPoint(11,455.1203571592949, 22.480051555019603);
  qe->SetPoint(12,472.5813503935195, 19.952411696799835);
  qe->SetPoint(13,492.7290928338012, 17.360557858087155);
  qe->SetPoint(14,507.54936254429964, 13.41903446952999);
  qe->SetPoint(15,515.6613890550927, 10.58224261870728);
  qe->SetPoint(16,526.4544115616513, 8.34426589284665);
  qe->SetPoint(17,541.2689280618584, 6.578539947918845);
  qe->SetPoint(18,556.0834445620656, 5.1864583897625645);
  qe->SetPoint(19,568.2112118562158, 4.171018298153135);
  qe->SetPoint(20,581.7024899894141, 3.099239370381881);
  qe->SetPoint(21,593.8590233350209, 2.257913856121514);
  qe->SetPoint(22,601.9825562663966, 1.7115761852025484);
  qe->SetPoint(23,611.4638468265291, 1.2226741190629429);
  qe->SetPoint(24,627.7396787407373, 0.6364576389506086);
  qe->SetPoint(25,644.0730427578588, 0.27188821897218746);
  qe->SetPoint(26,654.975376259953, 0.1472681410200624);
  qe->SetPoint(27,667.2814930731348, 0.06417745922393044);
  qe->SetPoint(28,675.5028305794633, 0.03476537833796253);
  qe->SetPoint(29,686.4857090256363, 0.01427881998143405);
  return qe;
}

double Detector::getFillFactor() {
  return fillFactor;
}

double Detector::evalQEff(double wav) {
  wav = wav*1E9;
  if ((wav >= 267.) && (wav <= 687.)) {
    return quantumEff->Eval(wav) / 100.;
  } else {
    return 0.;
  }
}

void Detector::projectPhotons(TH2D* photonHist, std::vector<Photon*> photons) {
  for (int j = 0; j < photons.size(); j++) {
    Photon* ph = photons[j];
    double phDist = ph->dist(zPos);
    photonHist->Fill(ph->pos[0] + phDist*ph->dir[0], ph->pos[1] + phDist*ph->dir[1]);
  }
}


