#include "stdlib.h"
#include <iostream>
#include "TFile.h"
#include "TMath.h"
using namespace std;
using namespace std::chrono;


void calcInteractionDist() {
  TString file;
  float n;
  float nerr;
  float thickness;
  TTree y("data", "data");
  y.ReadFile("info.csv", "file[30]/C:n/F:nerr:thickness");
  y.SetBranchAddress("file",&file);
  y.SetBranchAddress("n",&n);
  y.SetBranchAddress("nerr",&nerr);
  y.SetBranchAddress("thickness",&thickness);

  for(int i = 0; i < y.GetEntries(); i++) {
    y.GetEntry(i);
  }
}