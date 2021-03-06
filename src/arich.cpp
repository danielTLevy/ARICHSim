#include "arich.h"
using namespace std;

Arich::Arich(bool mirror) {
  this->aerogel1 = new Aerogel(n1, thickness1, aeroPos1);
  this->aerogel2 = new Aerogel(n2, thickness2, aeroPos2);
  this->detector = new Detector(detectorDist, mirror);
  aerogel1->setDownIndex(n2);
  aerogel2->setUpIndex(n1);
}


double Arich::integrateAndDrawEllipse(particleInfoStruct params, TH2D* photonHist, TPad* pad) {
  /*
    Draw an ellipse over a photonhist showing where particles are expected to go
    Return number of photon hits found in this ring
  */
  // Define Ellipse and integrate over this ring
  double dirX_0 = params.dir[0];
  double dirY_0 = params.dir[1];
  double dirZ_0 = params.dir[2];
  // Naively assume that all photons are generated in middle of second aerogel layer.
  double middlePoint = thickness1 + thickness2 / 2.;
  double newX_0 = params.pos[0] + middlePoint*dirX_0/dirZ_0;
  double newY_0 = params.pos[1] + middlePoint*dirY_0/dirZ_0;
  double newDist = detectorDist - middlePoint;
  double dirTheta =  atan(sqrt(dirX_0*dirX_0 +dirY_0*dirY_0) / dirZ_0);
  double chAngle = aerogel2->calcChAngle(params.beta);
  // Calculate ellipse parameters
  double radiusA = 0.5*newDist*TMath::Abs(tan(dirTheta+chAngle) - tan(dirTheta-chAngle));
  double radiusB = 0.5*newDist*TMath::Abs(tan(chAngle) - tan(-chAngle));
  // Distance travelled on detector plane with respect to travel direction
  double deltaR = newDist*tan(dirTheta-chAngle) + radiusA;
  // Multiply by X and Y components of direction to get final x and y position
  double ringX, ringY;
  if (dirX_0 != 0 || dirY_0 != 0) {
    ringX = newX_0 + deltaR*dirX_0/sqrt(dirX_0*dirX_0 +dirY_0*dirY_0);
    ringY = newY_0 + deltaR*dirY_0/sqrt(dirX_0*dirX_0 +dirY_0*dirY_0);
  } else {
    ringX = newX_0;
    ringY = newY_0;
  }
  // Rotate the ellipse by the particle's phi direction
  double dirPhiDeg;
  if (dirX_0 != 0) {
    dirPhiDeg = atan(dirY_0/dirX_0)*180./TMath::Pi();
  } else {
    dirPhiDeg = 0.;
  }
  // Create two ellipses to encapsulate the photons
  double ringOuterA = radiusA+1.5;
  double ringOuterB = radiusB+1.5;
  TEllipse *elOuter = new TEllipse(ringX,ringY,ringOuterA,ringOuterB, 0, 360, dirPhiDeg);
  double ringInnerA = radiusA-1.5;
  double ringInnerB = radiusB-1.5;
  TEllipse *elInner = new TEllipse(ringX,ringY,ringInnerA,ringInnerB, 0, 360, dirPhiDeg);
  TCutG *outerCut = createCutFromEllipse(elOuter);
  TCutG *innerCut = createCutFromEllipse(elInner);

  pad->cd();
  outerCut->Draw("same");
  innerCut->Draw("same");

  return outerCut->IntegralHist(photonHist) - innerCut->IntegralHist(photonHist);
}



TH2D* Arich::calculatePdf(particleInfoStruct params, char* histName) {
  /*
  Calculate mean distribution of photons over detector for given beta
  */
  // Number of events to simulate
  int nEvents = 10000;
  Beam* beam = new Beam(params.pos, params.dir, params.beta);
  // Get hist ready
  TH2D *photonHist = detector->makeDetectorHist(histName, histName);
  // Make events and loop over them
  for (int i = 0; i < nEvents; i++) {
    Particle *pa = beam->generateParticle();
    // Make photons in first aerogel
    std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);
    // Advance particle forward to next aerogel and generate photons
    pa->travelZDist(aeroPos2 - aeroPos1);
    std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);
    // Scatter photons in first aerogel, move them forwards out of aerogel
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);
    // Combine photons from both aerogels
    photons.insert(photons.end(), photons2.begin(), photons2.end());
    // Include scattering in second aerogel
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);
    // Do some more scattering
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, true);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, true);
    // Project photons onto 00 and plot distribution
    detector->projectPhotons(photonHist, photons);
    // delete
    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      delete ph;
    }
    delete pa;
  }
  // Scale photon histogram to the number of iterations
  photonHist->Scale(1. / nEvents);
  // Scale photon to fill factor of detector
  photonHist->Scale(detector->getFillFactor());
  return photonHist;
}

TH2D* Arich::generateEvent(particleInfoStruct params, bool save, char* histName, char* outputDir) {
  /*
  Simulate a photon distribution resulting from a single particle
  */
  Beam *beam = new Beam(params.pos, params.dir, params.beta);
  Particle *pa = beam->generateParticle();
  // Make photons in first aerogel
  std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);
  bool refract = true;
  aerogel1->applyPhotonScatters(photons);
  aerogel1->exitAerogel(photons, refract);
  // Advance particle forward to next aerogel and generate photons
  pa->travelZDist(aeroPos2 - aeroPos1);
  std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);
  // Combine photons from both aerogels
  photons.insert(photons.end(), photons2.begin(), photons2.end());
  // Include scattering in second aerogel
  aerogel2->applyPhotonScatters(photons);
  aerogel2->exitAerogel(photons, refract);
  aerogel1->applyPhotonScatters(photons);
  aerogel1->exitAerogel(photons, refract);
  aerogel2->applyPhotonScatters(photons);
  aerogel2->exitAerogel(photons, refract);
  // Throw out photons based off fill factor
  int numPhotonsDetected = (int) (detector->getFillFactor() * photons.size());
  photons.resize(numPhotonsDetected);
  // Project photons onto detector and plot distribution
  TH2D *photonHist = detector->makeDetectorHist(histName,histName);
  detector->projectPhotons(photonHist, photons);
  if (save) {
    // Draw out photon histogram and ellipse outline
    TCanvas *c1 = new TCanvas("c1","c1",900,900);
    gStyle->SetOptStat(0);
    double nPhotons = Arich::integrateAndDrawEllipse(params, photonHist, c1);
    cout << "SINGLE EXAMPLE EVENT: Integrated number of photons in ring: " << nPhotons << endl;
    photonHist->SaveAs(Form("%s/%s.root", outputDir, histName));
    photonHist->Draw("colz");
    c1->SaveAs(Form("%s/%s.pdf", outputDir, histName));
  }

  return photonHist;
}

TH2D* Arich::simulateBeam(particleInfoStruct params, char* outputDir) {
  /*
  Plot and save info into some output directory
  Including TTree holding photon and particle info, and final photon distribution
  */
  double nEvents = 10000;
  // Make beam
  Beam *beam = new Beam(params.pos, params.dir, params.beta);

  // Get plots ready
  TFile *f = new TFile(Form("%s.root", outputDir), "RECREATE");
  TH2D *photonHist = detector->makeDetectorHist("photonHist","photonHist");
  photonHist->SetXTitle("x [cm]");
  photonHist->SetYTitle("y [cm]");
  TH1D *rHist = new TH1D("rHist", "rHist", 500, 0., 10.);

  // Make a tree to save the photons
  photonStruct phStruct;
  particleStruct paStruct;
  TTree *tree = new TTree("T","Output photon data");
  TBranch *phBranch = tree->Branch("photons",&phStruct.dirxi,
    "dirxi/D:diryi:dirzi:posxi:posyi:poszi:dirxe:dirye:dirze:posxe:posye:posze:wav:paid/i:numscat");
  TBranch *paBranch = tree->Branch("particles",&paStruct.dirx,"dirx/D:diry:dirz:posx:posy:posz:id/i");

  // Make events and loop over them
  for (int i = 0; i < nEvents; i++) {
    Particle *pa = beam->generateParticle();
    paStruct.dirx = pa->dir0[0];
    paStruct.diry = pa->dir0[1];
    paStruct.dirz = pa->dir0[2];
    paStruct.posx = pa->pos0[0];
    paStruct.posy = pa->pos0[1];
    paStruct.posz = pa->pos0[2];
    paStruct.id = i;

    // Make photons in first aerogel 
    std::vector<Photon*> photons = aerogel1->generatePhotons(pa, detector);
    // Advance particle forward to next aerogel and generate photons
    pa->travelZDist(aeroPos2 - aeroPos1);
    std::vector<Photon*> photons2 = aerogel2->generatePhotons(pa, detector);
    // Scatter photons in first aerogel, move them forwards out of aerogel
    bool refract = true;
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, refract);
    // Combine photons from both aerogels
    photons.insert(photons.end(), photons2.begin(), photons2.end());
    // Include scattering in second aerogel
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, refract);
    // Do some more scattering

    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, refract);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, refract);
    aerogel1->applyPhotonScatters(photons);
    aerogel1->exitAerogel(photons, refract);
    aerogel2->applyPhotonScatters(photons);
    aerogel2->exitAerogel(photons, refract);

    // Project photons onto detector and plot distribution
    detector->projectPhotons(photonHist, photons, rHist);

    // Save photons to TTree and delete
    for (int j = 0; j < photons.size(); j++) {
      Photon* ph = photons[j];
      phStruct.paid = i;
      TVector3 phPos0 = ph->pos0;
      TVector3 phDir0 = ph->dir0;
      phStruct.dirxi = phDir0[0];
      phStruct.diryi = phDir0[1];
      phStruct.dirzi = phDir0[2];
      phStruct.posxi = phPos0[0];
      phStruct.posyi = phPos0[1];
      phStruct.poszi = phPos0[2];
      TVector3 phPos = ph->pos;
      TVector3 phDir = ph->dir;
      phStruct.dirxe = phDir[0];
      phStruct.dirye = phDir[1];
      phStruct.dirze = phDir[2];
      phStruct.posxe = phPos[0];
      phStruct.posye = phPos[1];
      phStruct.posze = phPos[2];
      phStruct.wav = ph->wav;
      phStruct.numscat = ph->numScatters;
      tree->Fill();
      delete ph;
    }
    delete pa;
  }

  // Scale photon histogram to the number of iterations
  photonHist->Scale(1. / nEvents);
  rHist->Scale(1. / nEvents);

  // Scale photon to fill factor of detector
  photonHist->Scale(detector->getFillFactor());
  rHist->Scale(detector->getFillFactor());

  // Draw out photon histogram and ellipse outline
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd();
  TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.55,0.55);
  center_pad->Draw();
  TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
  right_pad->Draw();
  TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
  top_pad->Draw();

  TPad *corner_pad = new TPad("corner_pad", "corner_pad",0.55,0.55,1.0,1.0);
  corner_pad->Draw();

  center_pad->cd();
  center_pad->SetGrid(1);
  photonHist->SetTitle("");
  photonHist->SetZTitle("Mean Photon Count");
  photonHist->Draw("colz");

  double nPhotons = Arich::integrateAndDrawEllipse(params, photonHist, center_pad);
  cout << "PHOTON DISTRIBUTION: Integrated number of photons in ring: " << nPhotons << endl;

  right_pad->cd();
  right_pad->SetGrid(1);
  TH1D* yProj = (TH1D*) photonHist->ProjectionY()->Clone("yProj");
  yProj->SetYTitle("y [cm]");
  yProj->Draw("hbar");

  top_pad->cd();
  top_pad->SetGrid(1);
  TH1D* xProj = (TH1D*) photonHist->ProjectionY()->Clone("xProj");
  xProj->SetYTitle("x [cm]");

  xProj->Draw("bar");

  corner_pad->cd();
  corner_pad->SetGrid(1);
  rHist->SetTitle("");
  rHist->SetXTitle("Distance from origin [cm]");
  rHist->Draw("bar");


  c1->Write();
  rHist->Write();
  photonHist->Write();

  f->Write();

  return photonHist;
}