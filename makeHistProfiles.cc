#include "TGraph.h"
#include "TTreeReader.h"

struct EventEntry{
  Int_t eventNumber;
  Int_t trackNumber;
  Int_t pid;
  Double_t pos[3]; 
  Double_t mom[3];
  Double_t dir[3];
  Int_t stateID;
  Int_t prodCode;
  string material1;
  string material2;
  string volume1;
  string volume2;
  Double_t fourMomTransfer;
};

Double_t fillFactor = 0.8*0.87;

void makeHistProfiles() {
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

	TH2D *photonHist = new TH2D("photonHist","photonHist",48,-150,150,48,-150,150);
	TH1D *rHist = new TH1D("rHist", "rHist", 500, 0., 100.);

	TFile *runFile = TFile::Open("./EMPHATICSim-build/test_exitonly.root");
	TTreeReader reader("h1000", runFile);
	TTreeReaderValue<Int_t> raPid(reader, "Pid");
	TTreeReaderArray<Double_t> raMom(reader, "Mom");
	TTreeReaderArray<Double_t> raPos(reader, "Pos");
	TTreeReaderArray<Double_t> raDir(reader, "Dir");
	TRandom3 randomGen = TRandom3();
	while (reader.Next()) {
		if (*raPid != 0 || raPos[2] < 829.99) {
			continue;
		}
		// Momentum in GeV/c:
		Double_t mom = sqrt(raMom[0]*raMom[0] + raMom[1]*raMom[1] + raMom[2]*raMom[2]);
		// Wavelength in nm (multiply by c and planck's constant in eV*s)
		Double_t wav = 1E9 *299792458 * 4.135667662E-15 / (1E9 * mom);
		Double_t qEff = max(qe->Eval(wav) / 100., 0.);
		Double_t xFinal = raPos[0] + (1000. - raPos[2])*raDir[0]/raDir[2];
		Double_t yFinal = raPos[1] + (1000. - raPos[2])*raDir[1]/raDir[2];
		
		photonHist->Fill(xFinal, yFinal, fillFactor*qEff/1000.);
		rHist->Fill(sqrt(xFinal*xFinal + yFinal*yFinal), fillFactor*qEff/1000.);
		
	}
	TCanvas *c1 = new TCanvas("c1", "c1",900,900);
	TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
	center_pad->Draw();
	TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
	right_pad->Draw();
	TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
	top_pad->Draw();
 	TPad *corner_pad = new TPad("corner_pad", "corner_pad",0.55,0.55,1.0,1.0);
	corner_pad->Draw();

	center_pad->cd();
	center_pad->SetGrid(1);
	photonHist->Draw("colz");

	right_pad->cd();
	right_pad->SetGrid(1);
  	photonHist->ProjectionY()->Draw("hbar");

	top_pad->cd();
	top_pad->SetGrid(1);
	photonHist->ProjectionX()->Draw("bar");


	corner_pad->cd();
	corner_pad->SetGrid(1);
	rHist->Draw("bar");

	int leftbin = rHist->GetXaxis()->FindBin(48.);
	int rightbin = rHist->GetXaxis()->FindBin(60.);
	cout << "Number of photons: " << rHist->Integral(leftbin, rightbin);
	//gStyle->SetOptStat(0);
	
	c1->SaveAs("histWithProjections.root");
	photonHist->SaveAs("photonHist.root");
	rHist->SaveAs("rHist.root");
}
