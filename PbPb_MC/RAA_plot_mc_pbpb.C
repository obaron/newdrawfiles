// Owen Baron
// August 29

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

static const Int_t nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
const int radius = 3;

void RAA_plot_mc_pbpb(){
//TFile *f = TFile::Open("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/PbPb_MC_histograms_FromForest_akPu3_anabin_20_eta_20.root");
TFile *f = new TFile(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/PbPb_MC_histograms_FromForest_akPu%d_anabin_20_eta_20.root",radius),"r");

TCanvas *EtaCan55[nbins_cent];
TCanvas *EtaCan65[nbins_cent];
TCanvas *EtaCan80[nbins_cent];
TCanvas *EtaPrecutCan55[nbins_cent];
TCanvas *EtaPrecutCan65[nbins_cent];
TCanvas *EtaPrecutCan80[nbins_cent];

TCanvas *pTCan[nbins_cent];
TCanvas *pTPrecutCan[nbins_cent];
TCanvas *pTCan80[nbins_cent];
TCanvas *pTPrecutCan80[nbins_cent];
TCanvas *pTCan65[nbins_cent];
TCanvas *pTPrecutCan65[nbins_cent];
TCanvas *pTCan55[nbins_cent];
TCanvas *pTPrecutCan55[nbins_cent];

TCanvas *EtaRatCan[nbins_cent];

TH1F *hEta80[nbins_cent];
TH1F *hPhi80[nbins_cent];
TH1F *hPhieta80[nbins_cent];
TH1F *hEta65[nbins_cent];
TH1F *hPhi65[nbins_cent];
TH1F *hPhieta65[nbins_cent];
TH1F *hEta55[nbins_cent];
TH1F *hPhi55[nbins_cent];
TH1F *hPhieta55[nbins_cent];

TH1F *hEtaPrecut80[nbins_cent];
TH1F *hPhiPrecut80[nbins_cent];
TH1F *hPhiEtaPrecut80[nbins_cent];
TH1F *hEtaPrecut65[nbins_cent];
TH1F *hPhiPrecut65[nbins_cent];
TH1F *hPhiEtaPrecut65[nbins_cent];
TH1F *hEtaPrecut55[nbins_cent];
TH1F *hPhiPrecut55[nbins_cent];
TH1F *hPhiEtaPrecut55[nbins_cent];

TH1F *pTGen[nbins_cent];
TH1F *pTPrecutGen[nbins_cent];
TH1F *pTGen80[nbins_cent];
TH1F *pTPrecutGen80[nbins_cent];
TH1F *pTGen65[nbins_cent];
TH1F *pTPrecutGen65[nbins_cent];
TH1F *pTGen55[nbins_cent];
TH1F *pTPrecutGen55[nbins_cent];
TH1F *pTReco[nbins_cent];
TH1F *pTPrecutReco[nbins_cent];
TH1F *pTReco80[nbins_cent];
TH1F *pTPrecutReco80[nbins_cent];
TH1F *pTReco65[nbins_cent];
TH1F *pTPrecutReco65[nbins_cent];
TH1F *pTReco55[nbins_cent];
TH1F *pTPrecutReco55[nbins_cent];

TH1F *EtaRatio80[nbins_cent];
TH1F *EtaRatio65[nbins_cent];
TH1F *EtaRatio55[nbins_cent];



for(int i=0;i<nbins_cent;i++){
//reduce margins if possible
	
	//cout<<"start loop"<<endl;
	EtaCan80[i] = new TCanvas(Form("EtaCan80_%d",i),Form("Jet80 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan80[i]->Divide(3,1);
	EtaCan80[i]->cd(1);
	hEta80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_etacent_R%d_cent%d",radius,i));
	hEta80[i]->GetXaxis()->SetTitle("Eta");
	hEta80[i]->GetXaxis()->CenterTitle();
	hEta80[i]->GetYaxis()->SetTitle("Multiplicity");
	hEta80[i]->GetYaxis()->CenterTitle();
	hEta80[i]->Draw();
	EtaCan80[i]->cd(2);
	hPhi80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_phicent_R%d_cent%d",radius,i));
	hPhi80[i]->GetXaxis()->SetTitle("Phi");
	hPhi80[i]->GetXaxis()->CenterTitle();
	hPhi80[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhi80[i]->GetYaxis()->CenterTitle();
	hPhi80[i]->Draw();
	EtaCan80[i]->cd(3);
	hPhieta80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_phieta_R%d_cent%d",radius,i));
	hPhieta80[i]->GetXaxis()->SetTitle("Eta");
	hPhieta80[i]->GetXaxis()->CenterTitle();
	hPhieta80[i]->GetYaxis()->SetTitle("Phi");
	hPhieta80[i]->GetYaxis()->CenterTitle();
	hPhieta80[i]->Draw("colz,logz");

	EtaCan65[i] = new TCanvas(Form("EtaCan65_%d",i),Form("Jet65 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan65[i]->Divide(3,1);
	EtaCan65[i]->cd(1);
	hEta65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_etacent_R%d_cent%d",radius,i));
	hEta65[i]->GetXaxis()->SetTitle("Eta");
	hEta65[i]->GetXaxis()->CenterTitle();
	hEta65[i]->GetYaxis()->SetTitle("Multiplicity");
	hEta65[i]->GetYaxis()->CenterTitle();
	hEta65[i]->Draw();
	EtaCan65[i]->cd(2);
	hPhi65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_phicent_R%d_cent%d",radius,i));
	hPhi65[i]->GetXaxis()->SetTitle("Phi");
	hPhi65[i]->GetXaxis()->CenterTitle();
	hPhi65[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhi65[i]->GetYaxis()->CenterTitle();
	hPhi65[i]->Draw();
	EtaCan65[i]->cd(3);
	hPhieta65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_phieta_R%d_cent%d",radius,i));
	hPhieta65[i]->GetXaxis()->SetTitle("Eta");
	hPhieta65[i]->GetXaxis()->CenterTitle();
	hPhieta65[i]->GetYaxis()->SetTitle("Phi");
	hPhieta65[i]->GetYaxis()->CenterTitle();
	hPhieta65[i]->Draw("colz,logz");
	
	EtaCan55[i] = new TCanvas(Form("EtaCan55_%d",i),Form("Jet55 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan55[i]->Divide(3,1);
	EtaCan55[i]->cd(1);
	hEta55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_etacent_R%d_cent%d",radius,i));
	hEta55[i]->GetXaxis()->SetTitle("Eta");
	hEta55[i]->GetXaxis()->CenterTitle();
	hEta55[i]->GetYaxis()->SetTitle("Multiplicity");
	hEta55[i]->GetYaxis()->CenterTitle();
	hEta55[i]->Draw();
	EtaCan55[i]->cd(2);
	hPhi55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_phicent_R%d_cent%d",radius,i));
	hPhi55[i]->GetXaxis()->SetTitle("Phi");
	hPhi55[i]->GetXaxis()->CenterTitle();
	hPhi55[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhi55[i]->GetYaxis()->CenterTitle();
	hPhi55[i]->Draw();
	EtaCan55[i]->cd(3);
	hPhieta55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_phieta_R%d_cent%d",radius,i));
	hPhieta55[i]->GetXaxis()->SetTitle("Eta");
	hPhieta55[i]->GetXaxis()->CenterTitle();
	hPhieta55[i]->GetYaxis()->SetTitle("Phi");
	hPhieta55[i]->GetYaxis()->CenterTitle();
	hPhieta55[i]->Draw("colz,logz");
	
	pTCan[i] = new TCanvas(Form("pTCan%d",i),Form("pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTCan[i]->Divide(2,1);
	pTCan[i]->cd(1);
	pTGen[i] = (TH1F*)f->Get(Form("hpbpb_gen_R%d_20_eta_20_cent%d",radius,i));
	pTGen[i]->GetXaxis()->SetTitle("Gen pT");
	pTGen[i]->GetXaxis()->CenterTitle();
	pTGen[i]->GetXaxis()->SetRange(0,300);
	pTGen[i]->Draw();
	pTCan[i]->cd(2);
	pTReco[i] = (TH1F*)f->Get(Form("hpbpb_reco_R%d_20_eta_20_cent%d",radius,i));
	pTReco[i]->GetXaxis()->SetTitle("Reco pT");
	pTReco[i]->GetXaxis()->CenterTitle();
	pTReco[i]->GetXaxis()->SetRange(0,300);
	pTReco[i]->Draw("logy");
	
	pTCan80[i] = new TCanvas(Form("pTCan80_%d",i),Form("Jet 80 pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTCan80[i]->Divide(2,1);
	pTCan80[i]->cd(1);
	pTGen80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_gen_R%d_20_eta_20_cent%d",radius,i));
	pTGen80[i]->GetXaxis()->SetTitle("Jet 80 Gen pT");
	pTGen80[i]->GetXaxis()->CenterTitle();
	pTGen80[i]->GetXaxis()->SetRange(0,300);
	pTGen80[i]->Draw("logy");
	pTCan80[i]->cd(2);
	pTReco80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_reco_R%d_20_eta_20_cent%d",radius,i));
	pTReco80[i]->GetXaxis()->SetTitle("Jet 80 Reco pT");
	pTReco80[i]->GetXaxis()->CenterTitle();
	pTReco80[i]->GetXaxis()->SetRange(0,300);
	pTReco80[i]->Draw("logy");
	
	pTCan65[i] = new TCanvas(Form("pTCan65_%d",i),Form("Jet 65 pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTCan65[i]->Divide(2,1);
	pTCan65[i]->cd(1);
	pTGen65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_gen_R%d_20_eta_20_cent%d",radius,i));
	pTGen65[i]->GetXaxis()->SetTitle("Jet 65 Gen pT");
	pTGen65[i]->GetXaxis()->CenterTitle();
	pTGen65[i]->GetXaxis()->SetRange(0,300);
	pTGen65[i]->Draw("logy");
	pTCan65[i]->cd(2);
	pTReco65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_reco_R%d_20_eta_20_cent%d",radius,i));
	pTReco65[i]->GetXaxis()->SetTitle("Jet 65 Reco pT");
	pTReco65[i]->GetXaxis()->CenterTitle();
	pTReco65[i]->GetXaxis()->SetRange(0,300);
	pTReco65[i]->Draw("logy");
	
	pTCan55[i] = new TCanvas(Form("pTCan55_%d",i),Form("Jet 55 pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTCan55[i]->Divide(2,1);
	pTCan55[i]->cd(1);
	pTGen55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_gen_R%d_20_eta_20_cent%d",radius,i));
	pTGen55[i]->GetXaxis()->SetTitle("Jet 55 Gen pT");
	pTGen55[i]->GetXaxis()->CenterTitle();
	pTGen55[i]->GetXaxis()->SetRange(0,300);
	pTGen55[i]->Draw("logy");
	pTCan55[i]->cd(2);
	pTReco55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_reco_R%d_20_eta_20_cent%d",radius,i));
	pTReco55[i]->GetXaxis()->SetTitle("Jet 55 Reco pT");
	pTReco55[i]->GetXaxis()->CenterTitle();
	pTReco55[i]->GetXaxis()->SetRange(0,300);
	pTReco55[i]->Draw("logy");
	
	
	//PRECUT LOOP
	
	EtaPrecutCan80[i] = new TCanvas(Form("EtaPrecutCan80_%d",i),Form("Jet80 Precut R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaPrecutCan80[i]->Divide(3,1);
	EtaPrecutCan80[i]->cd(1);
	hEtaPrecut80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_etaprecut_cent_R%d_cent%d",radius,i));
	hEtaPrecut80[i]->GetXaxis()->SetTitle("Eta");
	hEtaPrecut80[i]->GetXaxis()->CenterTitle();
	hEtaPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	hEtaPrecut80[i]->GetYaxis()->CenterTitle();
	hEtaPrecut80[i]->Draw();
	EtaPrecutCan80[i]->cd(2);
	hPhiPrecut80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_phiprecut_cent_R%d_cent%d",radius,i));
	hPhiPrecut80[i]->GetXaxis()->SetTitle("Phi");
	hPhiPrecut80[i]->GetXaxis()->CenterTitle();
	hPhiPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhiPrecut80[i]->GetYaxis()->CenterTitle();
	hPhiPrecut80[i]->Draw();
	EtaPrecutCan80[i]->cd(3);
	hPhiEtaPrecut80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_phieta_precut_R%d_cent%d",radius,i));
	hPhiEtaPrecut80[i]->GetXaxis()->SetTitle("Eta");
	hPhiEtaPrecut80[i]->GetXaxis()->CenterTitle();
	hPhiEtaPrecut80[i]->GetYaxis()->SetTitle("Phi");
	hPhiEtaPrecut80[i]->GetYaxis()->CenterTitle();
	hPhiEtaPrecut80[i]->Draw("colz,logz");

	EtaPrecutCan65[i] = new TCanvas(Form("EtaPrecutCan65_%d",i),Form("Jet65 Precut R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaPrecutCan65[i]->Divide(3,1);
	EtaPrecutCan65[i]->cd(1);
	hEtaPrecut65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_etaprecut_cent_R%d_cent%d",radius,i));
	hEtaPrecut65[i]->GetXaxis()->SetTitle("Eta");
	hEtaPrecut65[i]->GetXaxis()->CenterTitle();
	hEtaPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	hEtaPrecut65[i]->GetYaxis()->CenterTitle();
	hEtaPrecut65[i]->Draw();
	EtaPrecutCan65[i]->cd(2);
	hPhiPrecut65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_phiprecut_cent_R%d_cent%d",radius,i));
	hPhiPrecut65[i]->GetXaxis()->SetTitle("Phi");
	hPhiPrecut65[i]->GetXaxis()->CenterTitle();
	hPhiPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhiPrecut65[i]->GetYaxis()->CenterTitle();
	hPhiPrecut65[i]->Draw();
	EtaPrecutCan65[i]->cd(3);
	hPhiEtaPrecut65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_phieta_precut_R%d_cent%d",radius,i));
	hPhiEtaPrecut65[i]->GetXaxis()->SetTitle("Eta");
	hPhiEtaPrecut65[i]->GetXaxis()->CenterTitle();
	hPhiEtaPrecut65[i]->GetYaxis()->SetTitle("Phi");
	hPhiEtaPrecut65[i]->GetYaxis()->CenterTitle();
	hPhiEtaPrecut65[i]->Draw("colz,logz");
	
	EtaPrecutCan55[i] = new TCanvas(Form("EtaPrecutCan55_%d",i),Form("Jet55 Precut R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaPrecutCan55[i]->Divide(3,1);
	EtaPrecutCan55[i]->cd(1);
	hEtaPrecut55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_etaprecut_cent_R%d_cent%d",radius,i));
	hEtaPrecut55[i]->GetXaxis()->SetTitle("Eta");
	hEtaPrecut55[i]->GetXaxis()->CenterTitle();
	hEtaPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	hEtaPrecut55[i]->GetYaxis()->CenterTitle();
	hEtaPrecut55[i]->Draw();
	EtaPrecutCan55[i]->cd(2);
	hPhiPrecut55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_phiprecut_cent_R%d_cent%d",radius,i));
	hPhiPrecut55[i]->GetXaxis()->SetTitle("Phi");
	hPhiPrecut55[i]->GetXaxis()->CenterTitle();
	hPhiPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhiPrecut55[i]->GetYaxis()->CenterTitle();
	hPhiPrecut55[i]->Draw();
	EtaPrecutCan55[i]->cd(3);
	hPhiEtaPrecut55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_phieta_precut_R%d_cent%d",radius,i));
	hPhiEtaPrecut55[i]->GetXaxis()->SetTitle("Eta");
	hPhiEtaPrecut55[i]->GetXaxis()->CenterTitle();
	hPhiEtaPrecut55[i]->GetYaxis()->SetTitle("Phi");
	hPhiEtaPrecut55[i]->GetYaxis()->CenterTitle();
	hPhiEtaPrecut55[i]->Draw("colz,logz");
	
	pTPrecutCan[i] = new TCanvas(Form("pTPrecutCan_%d",i),Form("pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTPrecutCan[i]->Divide(2,1);
	pTPrecutCan[i]->cd(1);
	pTPrecutGen[i] = (TH1F*)f->Get(Form("hpbpb_gen_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutGen[i]->GetXaxis()->SetTitle("Precut Gen pT");
	pTPrecutGen[i]->GetXaxis()->CenterTitle();
	pTPrecutGen[i]->GetXaxis()->SetRange(0,300);
	pTPrecutGen[i]->Draw("logy");
	pTPrecutCan[i]->cd(2);
	pTPrecutReco[i] = (TH1F*)f->Get(Form("hpbpb_reco_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutReco[i]->GetXaxis()->SetTitle("Precut Reco pT");
	pTPrecutReco[i]->GetXaxis()->CenterTitle();
	pTPrecutReco[i]->GetXaxis()->SetRange(0,300);
	pTPrecutReco[i]->Draw("logy");
	
	pTPrecutCan80[i] = new TCanvas(Form("pTPrecutCan80_%d",i),Form("Jet 80 pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTPrecutCan80[i]->Divide(2,1);
	pTPrecutCan80[i]->cd(1);
	pTPrecutGen80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_gen_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutGen80[i]->GetXaxis()->SetTitle("Precut Jet 80 Gen pT");
	pTPrecutGen80[i]->GetXaxis()->CenterTitle();
	pTPrecutGen80[i]->GetXaxis()->SetRange(0,300);
	pTPrecutGen80[i]->Draw("logy");
	pTPrecutCan80[i]->cd(2);
	pTPrecutReco80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_reco_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutReco80[i]->GetXaxis()->SetTitle("Precut Jet 80 Reco pT");
	pTPrecutReco80[i]->GetXaxis()->CenterTitle();
	pTPrecutReco80[i]->GetXaxis()->SetRange(0,300);
	pTPrecutReco80[i]->Draw("logy");
		
	pTPrecutCan65[i] = new TCanvas(Form("pTPrecutCan65_%d",i),Form("Jet 65 pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTPrecutCan65[i]->Divide(2,1);
	pTPrecutCan65[i]->cd(1);
	pTPrecutGen65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_gen_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutGen65[i]->GetXaxis()->SetTitle("Precut Jet 65 Gen pT");
	pTPrecutGen65[i]->GetXaxis()->CenterTitle();
	pTPrecutGen65[i]->GetXaxis()->SetRange(0,300);
	pTPrecutGen65[i]->Draw("logy");
	pTPrecutCan65[i]->cd(2);
	pTPrecutReco65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_reco_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutReco65[i]->GetXaxis()->SetTitle("Precut Jet 65 Reco pT");
	pTPrecutReco65[i]->GetXaxis()->CenterTitle();
	pTPrecutReco65[i]->GetXaxis()->SetRange(0,300);
	pTPrecutReco65[i]->Draw("logy");
	
	pTPrecutCan55[i] = new TCanvas(Form("pTPrecutCan55_%d",i),Form("Jet 55 pT plots R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1167,500);
	pTPrecutCan55[i]->Divide(2,1);
	pTPrecutCan55[i]->cd(1);
	pTPrecutGen55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_gen_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutGen55[i]->GetXaxis()->SetTitle("Precut Jet 55 Gen pT");
	pTPrecutGen55[i]->GetXaxis()->CenterTitle();
	pTPrecutGen55[i]->GetXaxis()->SetRange(0,300);
	pTPrecutGen55[i]->Draw("logy");
	pTPrecutCan55[i]->cd(2);
	pTPrecutReco55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_reco_precut_R%d_20_eta_20_cent%d",radius,i));
	pTPrecutReco55[i]->GetXaxis()->SetTitle("Precut Jet 55 Reco pT");
	pTPrecutReco55[i]->GetXaxis()->CenterTitle();
	pTPrecutReco55[i]->GetXaxis()->SetRange(0,300);
	pTPrecutReco55[i]->Draw("logy");

	EtaRatCan[i] = new TCanvas(Form("EtaRatios_%d",i),Form("Eta Cut/Precut Ratios R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaRatCan[i]->Divide(3,1);
	EtaRatCan[i]->cd(1);
	EtaRatio80[i] = new TH1F(Form("E80%d",i),Form("Eta 80 Cut/Precut"),30,-2.5,2.5);
	EtaRatio80[i]->Divide(hEta80[i],hEtaPrecut80[i]);
	EtaRatio80[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio80[i]->GetXaxis()->CenterTitle();
	EtaRatio80[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio80[i]->Draw();
	EtaRatCan[i]->cd(2);
	EtaRatio65[i] = new TH1F(Form("E65%d",i),Form("Eta 65 Cut/Precut"),30,-2.5,2.5);
	EtaRatio65[i]->Divide(hEta65[i],hEtaPrecut65[i]);
	EtaRatio65[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio65[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio65[i]->GetXaxis()->CenterTitle();
	EtaRatio65[i]->Draw();
	EtaRatCan[i]->cd(3);
	EtaRatio55[i] = new TH1F(Form("E55%d",i),Form("Eta 55 Cut/Precut"),30,-2.5,2.5);
	EtaRatio55[i]->Divide(hEta55[i],hEtaPrecut55[i]);
	EtaRatio55[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio55[i]->GetXaxis()->CenterTitle();
	EtaRatio55[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio55[i]->Draw();
	
	EtaCan80[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet80_R%d_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaCan65[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet65_R%d_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaCan55[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet55_R%d_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaPrecutCan80[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet80_R%d_phietaprecut_cent%d.pdf",radius,i),"RECREATE");
	EtaPrecutCan65[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet65_R%d_phietaprecut_cent%d.pdf",radius,i),"RECREATE");
	EtaPrecutCan55[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet55_R%d_phietaprecut_cent%d.pdf",radius,i),"RECREATE");
	pTCan[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_gen_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTCan80[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet80_gen_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTCan65[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet65_gen_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTCan55[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet55_gen_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTPrecutCan[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_pTprecut_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTPrecutCan80[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet80_pTprecut_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTPrecutCan65[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet65_pTprecut_R%d_cent%d.pdf",radius,i),"RECREATE");
	pTPrecutCan55[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_Jet55_pTprecut_R%d_cent%d.pdf",radius,i),"RECREATE");
	EtaRatCan[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/output/MC/hpbpb_etaratio_R%d_cent%d.pdf",radius,i),"RECREATE");

	
	cout<<"end loop"<<endl;	
}

//pbpb plots here:
//TCanvas *t[radius][nbins_eta]

}
