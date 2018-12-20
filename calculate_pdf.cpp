#include "stdlib.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TEllipse.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVector3.h"

TMatrixD makeRotationMatrix(TVector3 dir) {
  // https://math.stackexchange.com/questions/180418
  // https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation

  TMatrixD rot = TMatrixD(3, 3);
  // a = (0, 0, 1)
  // b = (dirx, diry, dirz)
  // v = a x b = (-diry, dirx, 0)
  TVector3 v(-dir.Y(), dir.X(), 0);
  double s = v.Mag(); // sin(a,b) = ||v|| 
  double c = dir.Z(); // cos(a,b) = a . b = (0,0,1) . (dirx, diry, dirz) = dirz

  // rot =   I   + v    + v^2                    * (1-c)/s^2
  TMatrixD V = TMatrixD(3, 3);
  /*
  V(0,0) = 0;
  V(1,0) = v(2);
  V(2,0) = -v(1);
  V(0,1) = -v(2);
  v(1,1) = 0;
  v(2,1) = v(0);
  v(0,2) = v(1);
  v(1,2) = -v(0);
  v(2,2) = 0;
  */
  double k = 1./(1.+c);
  rot(0,0) = 1.0        - (v[2]*v[2]+v[1]*v[1])*k;
  rot(0,1) =     - v[2] - (v[0]*v[1])          *k;
  rot(0,2) =     + v[1] + (v[0]*v[2])          *k;
  rot(1,0) =     + v[2] + (v[0]*v[1])          *k;
  rot(1,1) = 1.0        - (v[2]*v[2]+v[0]*v[0])*k; 
  rot(1,2) =     - v[0] - (v[2]*v[1])          *k;
  rot(2,0) =     - v[1] - (v[0]*v[2])          *k;
  rot(2,1) =     + v[0] + (v[1]*v[2])          *k;
  rot(2,2) = 1.0        - (v[1]*v[1]+v[0]*v[0])*k;
  return rot;
}

int main(int argc, char *argv[]) {
//void calculate_pdf() {
  double beta = 0.9999;
  double n = 1.06;
  double dist[2] = {21.0,19.0};
  double x_0 = -2.0;
  double y_0 = -2.0;
  TVector3 pos_0 = TVector3(x_0, y_0, dist[0]);
  double dirX_0 = 0.01;
  double dirY_0 = -0.1;
  double dirZ_0 = sqrt(1 - dirX_0*dirX_0 - dirY_0*dirY_0);
  TVector3 dir_0 = TVector3(dirX_0, dirY_0, dirZ_0);

  double errX = 0.0;
  double errY = 0.0;
  double errDirX = 0.01;
  double errDirY = 0.01;



  double thetaCh = 0.;
  if(n*beta>=1.0){ 
    thetaCh =  acos(1.0/(n*beta));
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  TH2D *beamHist = new TH2D("beamHist","beamHist",200,-15,15,200,-15,15);

  //TCanvas *c2 = new TCanvas("c2","c2",600,500);
  //TH2D *photonHist = new TH2D("photonHist","photonHist",200,-15,15,200,-15,15);

  std::shared_ptr<TRandom3> randomGenerate(std::make_shared<TRandom3>());
  randomGenerate->SetSeed(1);

  int nIter = 10000;
  //double xs[100]; double ys[100]; double ps[100]; double ts[100];
  for (int i = 0; i < nIter; i++) {
    // Generate particles

    double paX = x_0 + randomGenerate->Gaus(0.0, errX);
    double paY = y_0 + randomGenerate->Gaus(0.0, errY);
    double paDirX = dirX_0 + randomGenerate->Gaus(0.0, errDirX);
    double paDirY = dirY_0 + randomGenerate->Gaus(0.0, errDirY);
    double paDirZ =  sqrt(1 - paDirX*paDirX - paDirY*paDirY);
    double paTheta =  atan(sqrt(paDirX*paDirX + paDirY*paDirY) / paDirZ);
    double paPhi = atan(paDirX / paDirY);

    // optional: Project beam trajectory onto detector
    double z = dist[0];
    double r = z / cos(paTheta);
    beamHist->Fill(paX+r*paDirX, paY+r*paDirY);


    // Make the rotation matrix: 
    TVector3 paDir(paDirX, paDirY, paDirZ);
    TMatrixD rot = makeRotationMatrix(paDir);

    // dist travelled through the aerogel by particle
    double paDist = (dist[0] - dist[1]) / cos(paTheta);

    // Generate photons
    int nPhotons = 10;
    for (int j = 0; j < nPhotons; j++) {
      // Distance into aerogel that photon is emitted
      double phIntDist = randomGenerate->Uniform(paDist);
      double phIntR = r - phIntDist;

      double phX = paX+phIntDist*sin(paTheta)*cos(paPhi);
      double phY = paY+phIntDist*sin(paTheta)*sin(paPhi);

      // Photons
      double phPhi = randomGenerate->Uniform(0., 2.*TMath::Pi());

      TVector3 dirC;
      TVector3 dirCR;
      dirC.SetX(cos(phPhi)*sin(thetaCh)); 
      dirC.SetY(sin(phPhi)*sin(thetaCh));  
      dirC.SetZ(cos(thetaCh));
      // x y z of cherenkov rotated onto particle
      //TVector3 dirCR = rot * dirC;
      dirCR[0] = rot[0][0]*dirC[0]+rot[0][1]*dirC[1]+rot[0][2]*dirC[2]; 
      dirCR[1] = rot[1][0]*dirC[0]+rot[1][1]*dirC[1]+rot[1][2]*dirC[2];
      dirCR[2] = rot[2][0]*dirC[0]+rot[2][1]*dirC[1]+rot[2][2]*dirC[2]; 
      phIntDist = randomGenerate->Rndm()*(dist[1]-dist[0])+dist[0];
      beamHist->Fill(phX + phIntDist*dirCR[0]/dirCR[2], phY + phIntDist*dirCR[1]/dirCR[2]);
    }


  }
  //beamHist->SetStats(false);
  c1->cd();
  beamHist->Draw("colz");
  c1->SaveAs("beamhist.pdf");
  //c2->cd();
  //photonHist->Draw("colz");


};


