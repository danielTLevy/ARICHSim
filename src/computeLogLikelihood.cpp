#include "stdlib.h"
#include <iostream>
#include <chrono>
#include <TROOT.h>
#include <TStyle.h>
#include "TMath.h"
#include "TH2D.h"
#include "TH1D.h"

double computeLogLikelihood(TH2D* event, TH2D* distribution) {
    int nBins = event->GetSize();
    if (nBins != distribution->GetSize()) {
        throw "Error: Bin Mismatch";
    }
    double logLikelihood = 0.;
    for (int i = 0; i < nBins; i++) { 
        double lambda = distribution->GetBinContent(i);
        bool pixelHit = event->GetBinContent(i) > 0;
        if (pixelHit) {
            logLikelihood += log(exp(-lambda));
        } else {
            logLikelihood += log(1 - exp(-lambda));
        }
    }
    return logLikelihood;
}

/*
TFile *f = new TFile("photonHist.root")
TH2D* photonHist = (TH2D*) f->Get("photonHist")
TFile *g = new TFile("generatedEvent.root")
TH2D* generatedEvent =(TH2D*) g->Get("generatedEvent")
computeLogLikelihood(generatedEvent, photonHist)
*/