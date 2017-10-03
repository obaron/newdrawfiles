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
const int radius = 2;

//leg->SetTextSize(0.04);
//leg->Draw();


void pbpb_Plot_RAA(){
TFile *fdata = new TFile(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/PbPb_Data_histograms_FromForest_30GeV_akPu%d_anabin_20_eta_20.root",radius),"r");
TFile *fmc = new TFile(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/PbPb_MC_histograms_FromForest_30GeV_akPu%d_anabin_20_eta_20.root",radius),"r");

TCanvas *EtaCan55;
TCanvas *EtaCan65;
TCanvas *EtaCan80;
TCanvas *EtaPrecutCan55;
TCanvas *EtaPrecutCan65;
TCanvas *EtaPrecutCan80;
TCanvas *EtaRatCan55;
TCanvas *EtaRatCan65;
TCanvas *EtaRatCan80;

TCanvas *PhiCan55;
TCanvas *PhiCan65;
TCanvas *PhiCan80;
TCanvas *PhiPrecutCan55;
TCanvas *PhiPrecutCan65;
TCanvas *PhiPrecutCan80;

TCanvas *PhiEtaCan55;
TCanvas *PhiEtaCan65;
TCanvas *PhiEtaCan80;
TCanvas *PhiEtaPrecutCan55;
TCanvas *PhiEtaPrecutCan65;
TCanvas *PhiEtaPrecutCan80;


TH1F *h_mc_Eta80[nbins_cent];
TH1F *h_mc_Phi80[nbins_cent];
TH1F *h_mc_PhiEta80[nbins_cent];
TH1F *h_mc_Eta65[nbins_cent];
TH1F *h_mc_Phi65[nbins_cent];
TH1F *h_mc_PhiEta65[nbins_cent];
TH1F *h_mc_Eta55[nbins_cent];
TH1F *h_mc_Phi55[nbins_cent];
TH1F *h_mc_PhiEta55[nbins_cent];

TH1F *h_mc_EtaPrecut80[nbins_cent];
TH1F *h_mc_PhiPrecut80[nbins_cent];
TH1F *h_mc_PhiEtaPrecut80[nbins_cent];
TH1F *h_mc_EtaPrecut65[nbins_cent];
TH1F *h_mc_PhiPrecut65[nbins_cent];
TH1F *h_mc_PhiEtaPrecut65[nbins_cent];
TH1F *h_mc_EtaPrecut55[nbins_cent];
TH1F *h_mc_PhiPrecut55[nbins_cent];
TH1F *h_mc_PhiEtaPrecut55[nbins_cent];

TH1F *h_mc_EtaRatio80[nbins_cent];
TH1F *h_mc_EtaRatio65[nbins_cent];
TH1F *h_mc_EtaRatio55[nbins_cent];

TH1F *h_data_Eta80[nbins_cent];
TH1F *h_data_Phi80[nbins_cent];
TH1F *h_data_PhiEta80[nbins_cent];
TH1F *h_data_Eta65[nbins_cent];
TH1F *h_data_Phi65[nbins_cent];
TH1F *h_data_PhiEta65[nbins_cent];
TH1F *h_data_Eta55[nbins_cent];
TH1F *h_data_Phi55[nbins_cent];
TH1F *h_data_PhiEta55[nbins_cent];

TH1F *h_data_EtaPrecut80[nbins_cent];
TH1F *h_data_PhiPrecut80[nbins_cent];
TH1F *h_data_PhiEtaPrecut80[nbins_cent];
TH1F *h_data_EtaPrecut65[nbins_cent];
TH1F *h_data_PhiPrecut65[nbins_cent];
TH1F *h_data_PhiEtaPrecut65[nbins_cent];
TH1F *h_data_EtaPrecut55[nbins_cent];
TH1F *h_data_PhiPrecut55[nbins_cent];
TH1F *h_data_PhiEtaPrecut55[nbins_cent];

TH1F *h_data_EtaRatio80[nbins_cent];
TH1F *h_data_EtaRatio65[nbins_cent];
TH1F *h_data_EtaRatio55[nbins_cent];

//TLegend *leg[nbins_cent] = myLegend(0.15,0.15,0.25,0.25);    

EtaCan80 = new TCanvas(Form("EtaCan80"),Form("Jet80 R%d Eta plots for PbPb MC and Data",radius),1750,1000);
EtaCan80->Divide(3,2);
EtaCan65 = new TCanvas(Form("EtaCan65"),Form("Jet65 R%d Eta plots for PbPb MC and Data",radius),1750,1000);
EtaCan65->Divide(3,2);
EtaCan55 = new TCanvas(Form("EtaCan55"),Form("Jet55 R%d Eta plots for PbPb MC and Data",radius),1750,1000);
EtaCan55->Divide(3,2);
EtaPrecutCan80 = new TCanvas(Form("EtaPrecutCan80"),Form("Jet80 R%d Eta plots without CutA for PbPb MC and Data",radius),1750,1000);
EtaPrecutCan80->Divide(3,2);
EtaPrecutCan65 = new TCanvas(Form("EtaPrecutCan65"),Form("Jet65 R%d Eta plots without CutA for PbPb MC and Data",radius),1750,1000);
EtaPrecutCan65->Divide(3,2);
EtaPrecutCan55 = new TCanvas(Form("EtaPrecutCan55"),Form("Jet55 R%d Eta plots without CutA for PbPb MC and Data",radius),1750,1000);
EtaPrecutCan55->Divide(3,2);

PhiCan80 = new TCanvas(Form("PhiCan80"),Form("Jet80 R%d Phi plots for PbPb MC and Data",radius),1750,1000);
PhiCan80->Divide(3,2);
PhiCan65 = new TCanvas(Form("PhiCan65"),Form("Jet65 R%d Phi plots for PbPb MC and Data",radius),1750,1000);
PhiCan65->Divide(3,2);
PhiCan55 = new TCanvas(Form("PhiCan55"),Form("Jet55 R%d Phi plots for PbPb MC and Data",radius),1750,1000);
PhiCan55->Divide(3,2);
PhiPrecutCan80 = new TCanvas(Form("PhiPrecutCan80"),Form("Jet80 R%d Phi plots without CutA for PbPb MC and Data",radius),1750,1000);
PhiPrecutCan80->Divide(3,2);
PhiPrecutCan65 = new TCanvas(Form("PhiPrecutCan65"),Form("Jet65 R%d Phi plots without CutA for PbPb MC and Data",radius),1750,1000);
PhiPrecutCan65->Divide(3,2);
PhiPrecutCan55 = new TCanvas(Form("PhiPrecutCan55"),Form("Jet55 R%d Phi plots without CutA for PbPb MC and Data",radius),1750,1000);
PhiPrecutCan55->Divide(3,2);

PhiEtaCan80 = new TCanvas(Form("PhiEtaCan80"),Form("Jet80 R%d Phi vs Eta plots for PbPb MC and Data",radius),1750,1000);
PhiEtaCan80->Divide(3,2);
PhiEtaCan65 = new TCanvas(Form("PhiEtaCan65"),Form("Jet65 R%d Phi vs Eta plots for PbPb MC and Data",radius),1750,1000);
PhiEtaCan65->Divide(3,2);
PhiEtaCan55 = new TCanvas(Form("PhiEtaCan55"),Form("Jet55 R%d Phi vs Eta plots for PbPb MC and Data",radius),1750,1000);
PhiEtaCan55->Divide(3,2);
PhiEtaPrecutCan80 = new TCanvas(Form("PhiEtaPrecutCan80"),Form("Jet80 R%d Phi vs Eta plots without CutA for PbPb MC and Data",radius),1750,1000);
PhiEtaPrecutCan80->Divide(3,2);
PhiEtaPrecutCan65 = new TCanvas(Form("PhiEtaPrecutCan65"),Form("Jet65 R%d Phi vs Eta plots without CutA for PbPb MC and Data",radius),1750,1000);
PhiEtaPrecutCan65->Divide(3,2);
PhiEtaPrecutCan55 = new TCanvas(Form("PhiEtaPrecutCan55"),Form("Jet55 R%d Phi vs Eta plots without CutA for PbPb MC and Data",radius),1750,1000);
PhiEtaPrecutCan55->Divide(3,2);

EtaRatCan80 = new TCanvas(Form("EtaRatCan80"),Form("Jet80 R%d Eta ratio plots for PbPb MC and Data",radius),1750,1000);
EtaRatCan80->Divide(3,2);
EtaRatCan65 = new TCanvas(Form("EtaRatCan65"),Form("Jet65 R%d Eta ratio plots for PbPb MC and Data",radius),1750,1000);
EtaRatCan65->Divide(3,2);
EtaRatCan55 = new TCanvas(Form("EtaRatCan55"),Form("Jet55 R%d Eta ratio plots for PbPb MC and Data",radius),1750,1000);
EtaRatCan55->Divide(3,2);

//TLegend leg[nbins_cent] = TLegend(0.15,0.15,0.25,0.25);
TLegend *leg[nbins_cent];

for(int i=0;i<nbins_cent;i++){

	EtaCan80->cd(i+1);
	h_mc_Eta80[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet80_etacent_R%d_cent%d",radius,i));
	h_mc_Eta80[i]->GetXaxis()->SetTitle("Eta");
	h_mc_Eta80[i]->GetXaxis()->CenterTitle();
	h_mc_Eta80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_Eta80[i]->GetYaxis()->CenterTitle();
	h_data_Eta80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_etacent_R%d_cent%d",radius,i));
	h_data_Eta80[i]->SetLineColor(2);
	h_mc_Eta80[i]->DrawNormalized();
	h_data_Eta80[i]->DrawNormalized("same");
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_Eta80[0],"MC","pl");
    leg[i]->AddEntry(h_data_Eta80[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
    
	
	EtaCan65->cd(i+1);
	h_mc_Eta65[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet65_etacent_R%d_cent%d",radius,i));
	h_mc_Eta65[i]->GetXaxis()->SetTitle("Eta");
	h_mc_Eta65[i]->GetXaxis()->CenterTitle();
	h_mc_Eta65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_Eta65[i]->GetYaxis()->CenterTitle();
	h_mc_Eta65[i]->DrawNormalized();
	h_data_Eta65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_etacent_R%d_cent%d",radius,i));
	h_data_Eta65[i]->SetLineColor(2);
	h_data_Eta65[i]->DrawNormalized("same");	
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_Eta65[0],"MC","pl");
    leg[i]->AddEntry(h_data_Eta65[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");	
	
	EtaCan55->cd(i+1);
	h_mc_Eta55[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet55_etacent_R%d_cent%d",radius,i));
	h_mc_Eta55[i]->GetXaxis()->SetTitle("Eta");
	h_mc_Eta55[i]->GetXaxis()->CenterTitle();
	h_mc_Eta55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_Eta55[i]->GetYaxis()->CenterTitle();
	h_mc_Eta55[i]->DrawNormalized();
	h_data_Eta55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_etacent_R%d_cent%d",radius,i));
	h_data_Eta55[i]->SetLineColor(2);
	h_data_Eta55[i]->DrawNormalized("same");
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_Eta55[0],"MC","pl");
    leg[i]->AddEntry(h_data_Eta55[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	PhiCan80->cd(i+1);
	h_mc_Phi80[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet80_phicent_R%d_cent%d",radius,i));
	h_mc_Phi80[i]->GetXaxis()->SetTitle("Phi");
	h_mc_Phi80[i]->GetXaxis()->CenterTitle();
	h_mc_Phi80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_Phi80[i]->GetYaxis()->CenterTitle();
	h_data_Phi80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phicent_R%d_cent%d",radius,i));
	h_data_Phi80[i]->SetLineColor(2);
	h_mc_Phi80[i]->DrawNormalized();
	h_data_Phi80[i]->DrawNormalized("same");
	leg[i] = new TLegend(.78,.80,0.88,0.88);
	leg[i]->AddEntry(h_mc_Phi80[0],"MC","pl");
    leg[i]->AddEntry(h_data_Phi80[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	PhiCan65->cd(i+1);
	h_mc_Phi65[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet65_phicent_R%d_cent%d",radius,i));
	h_mc_Phi65[i]->GetXaxis()->SetTitle("Phi");
	h_mc_Phi65[i]->GetXaxis()->CenterTitle();
	h_mc_Phi65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_Phi65[i]->GetYaxis()->CenterTitle();
	h_mc_Phi65[i]->DrawNormalized();
	h_data_Phi65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phicent_R%d_cent%d",radius,i));
	h_data_Phi65[i]->SetLineColor(2);
	h_data_Phi65[i]->DrawNormalized("same");
	leg[i] = new TLegend(.78,.80,0.88,0.88);
	leg[i]->AddEntry(h_mc_Phi65[0],"MC","pl");
    leg[i]->AddEntry(h_data_Phi65[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	PhiCan55->cd(i+1);
	h_mc_Phi55[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet55_phicent_R%d_cent%d",radius,i));
	h_mc_Phi55[i]->GetXaxis()->SetTitle("Phi");
	h_mc_Phi55[i]->GetXaxis()->CenterTitle();
	h_mc_Phi55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_Phi55[i]->GetYaxis()->CenterTitle();
	h_mc_Phi55[i]->DrawNormalized();
	h_data_Phi55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phicent_R%d_cent%d",radius,i));
	h_data_Phi55[i]->SetLineColor(2);
	h_data_Phi55[i]->DrawNormalized("same");
	leg[i] = new TLegend(.78,.80,0.88,0.88);
	leg[i]->AddEntry(h_mc_Phi55[0],"MC","pl");
    leg[i]->AddEntry(h_data_Phi55[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	/*
	PhiEtaCan80->cd(i+1);
	h_mc_PhiEta80[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet80_phieta_R%d_cent%d",radius,i));
	h_mc_PhiEta80[i]->GetXaxis()->SetTitle("PhiEta");
	h_mc_PhiEta80[i]->GetXaxis()->CenterTitle();
	h_mc_PhiEta80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiEta80[i]->GetYaxis()->CenterTitle();
	h_data_PhiEta80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phieta_R%d_cent%d",radius,i));
	h_data_PhiEta80[i]->SetLineColor(2);
	h_mc_PhiEta80[i]->DrawNormalized();
	h_data_PhiEta80[i]->DrawNormalized("same");
		
	PhiEtaCan65->cd(i+1);
	h_mc_PhiEta65[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet65_phieta_R%d_cent%d",radius,i));
	h_mc_PhiEta65[i]->GetXaxis()->SetTitle("PhiEta");
	h_mc_PhiEta65[i]->GetXaxis()->CenterTitle();
	h_mc_PhiEta65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiEta65[i]->GetYaxis()->CenterTitle();
	h_mc_PhiEta65[i]->DrawNormalized();
	h_data_PhiEta65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phieta_R%d_cent%d",radius,i));
	h_data_PhiEta65[i]->SetLineColor(2);
	h_data_PhiEta65[i]->DrawNormalized("same");	
	
	PhiEtaCan55->cd(i+1);
	h_mc_PhiEta55[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet55_phieta_R%d_cent%d",radius,i));
	h_mc_PhiEta55[i]->GetXaxis()->SetTitle("PhiEta");
	h_mc_PhiEta55[i]->GetXaxis()->CenterTitle();
	h_mc_PhiEta55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiEta55[i]->GetYaxis()->CenterTitle();
	h_mc_PhiEta55[i]->DrawNormalized();
	h_data_PhiEta55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phieta_R%d_cent%d",radius,i));
	h_data_PhiEta55[i]->SetLineColor(2);
	h_data_PhiEta55[i]->DrawNormalized("same");	
	*/
	//precut
	
	
	EtaPrecutCan80->cd(i+1);
	h_mc_EtaPrecut80[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet80_etaprecut_cent_R%d_cent%d",radius,i));
	h_mc_EtaPrecut80[i]->GetXaxis()->SetTitle("Eta");
	h_mc_EtaPrecut80[i]->GetXaxis()->CenterTitle();
	h_mc_EtaPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_EtaPrecut80[i]->GetYaxis()->CenterTitle();
	h_data_EtaPrecut80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_etaprecut_cent_R%d_cent%d",radius,i));
	h_data_EtaPrecut80[i]->SetLineColor(2);
	h_mc_EtaPrecut80[i]->DrawNormalized();
	h_data_EtaPrecut80[i]->DrawNormalized("same");
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_EtaPrecut80[0],"MC","pl");
    leg[i]->AddEntry(h_data_EtaPrecut80[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	EtaPrecutCan65->cd(i+1);
	h_mc_EtaPrecut65[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet65_etaprecut_cent_R%d_cent%d",radius,i));
	h_mc_EtaPrecut65[i]->GetXaxis()->SetTitle("Eta");
	h_mc_EtaPrecut65[i]->GetXaxis()->CenterTitle();
	h_mc_EtaPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_EtaPrecut65[i]->GetYaxis()->CenterTitle();
	h_mc_EtaPrecut65[i]->DrawNormalized();
	h_data_EtaPrecut65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_etaprecut_cent_R%d_cent%d",radius,i));
	h_data_EtaPrecut65[i]->SetLineColor(2);
	h_data_EtaPrecut65[i]->DrawNormalized("same");	
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_EtaPrecut65[0],"MC","pl");
    leg[i]->AddEntry(h_data_EtaPrecut65[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	EtaPrecutCan55->cd(i+1);
	h_mc_EtaPrecut55[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet55_etaprecut_cent_R%d_cent%d",radius,i));
	h_mc_EtaPrecut55[i]->GetXaxis()->SetTitle("Eta");
	h_mc_EtaPrecut55[i]->GetXaxis()->CenterTitle();
	h_mc_EtaPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_EtaPrecut55[i]->GetYaxis()->CenterTitle();
	h_mc_EtaPrecut55[i]->DrawNormalized();
	h_data_EtaPrecut55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_etaprecut_cent_R%d_cent%d",radius,i));
	h_data_EtaPrecut55[i]->SetLineColor(2);
	h_data_EtaPrecut55[i]->DrawNormalized("same");
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_EtaPrecut55[0],"MC","pl");
    leg[i]->AddEntry(h_data_EtaPrecut55[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	PhiPrecutCan80->cd(i+1);
	h_mc_PhiPrecut80[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet80_phiprecut_cent_R%d_cent%d",radius,i));
	h_mc_PhiPrecut80[i]->GetXaxis()->SetTitle("Phi");
	h_mc_PhiPrecut80[i]->GetXaxis()->CenterTitle();
	h_mc_PhiPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiPrecut80[i]->GetYaxis()->CenterTitle();
	h_data_PhiPrecut80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phiprecut_cent_R%d_cent%d",radius,i));
	h_data_PhiPrecut80[i]->SetLineColor(2);
	h_mc_PhiPrecut80[i]->DrawNormalized();
	h_data_PhiPrecut80[i]->DrawNormalized("same");
	leg[i] = new TLegend(.78,.80,0.88,0.88);
	leg[i]->AddEntry(h_mc_PhiPrecut80[0],"MC","pl");
    leg[i]->AddEntry(h_data_PhiPrecut80[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
		
	PhiPrecutCan65->cd(i+1);
	h_mc_PhiPrecut65[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet65_phiprecut_cent_R%d_cent%d",radius,i));
	h_mc_PhiPrecut65[i]->GetXaxis()->SetTitle("Phi");
	h_mc_PhiPrecut65[i]->GetXaxis()->CenterTitle();
	h_mc_PhiPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiPrecut65[i]->GetYaxis()->CenterTitle();
	h_mc_PhiPrecut65[i]->DrawNormalized();
	h_data_PhiPrecut65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phiprecut_cent_R%d_cent%d",radius,i));
	h_data_PhiPrecut65[i]->SetLineColor(2);
	h_data_PhiPrecut65[i]->DrawNormalized("same");	
	leg[i] = new TLegend(.78,.80,0.88,0.88);
	leg[i]->AddEntry(h_mc_PhiPrecut65[0],"MC","pl");
    leg[i]->AddEntry(h_data_PhiPrecut65[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	PhiPrecutCan55->cd(i+1);
	h_mc_PhiPrecut55[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet55_phiprecut_cent_R%d_cent%d",radius,i));
	h_mc_PhiPrecut55[i]->GetXaxis()->SetTitle("Phi");
	h_mc_PhiPrecut55[i]->GetXaxis()->CenterTitle();
	h_mc_PhiPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiPrecut55[i]->GetYaxis()->CenterTitle();
	h_mc_PhiPrecut55[i]->DrawNormalized();
	h_data_PhiPrecut55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phiprecut_cent_R%d_cent%d",radius,i));
	h_data_PhiPrecut55[i]->SetLineColor(2);
	h_data_PhiPrecut55[i]->DrawNormalized("same");
	leg[i] = new TLegend(.78,.80,0.88,0.88);
	leg[i]->AddEntry(h_mc_PhiPrecut55[0],"MC","pl");
    leg[i]->AddEntry(h_data_PhiPrecut55[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");

/*	
	PhiEtaPrecutCan80->cd(i+1);
	h_mc_PhiEtaPrecut80[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet80_phieta_precut_R%d_cent%d",radius,i));
	h_mc_PhiEtaPrecut80[i]->GetXaxis()->SetTitle("PhiEta");
	h_mc_PhiEtaPrecut80[i]->GetXaxis()->CenterTitle();
	h_mc_PhiEtaPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiEtaPrecut80[i]->GetYaxis()->CenterTitle();
	h_data_PhiEtaPrecut80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phieta_precut_R%d_cent%d",radius,i));
	h_data_PhiEtaPrecut80[i]->SetLineColor(2);
	h_mc_PhiEtaPrecut80[i]->DrawNormalized();
	h_data_PhiEtaPrecut80[i]->DrawNormalized("same");
		
	PhiEtaPrecutCan65->cd(i+1);
	h_mc_PhiEtaPrecut65[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet65_phieta_precut_R%d_cent%d",radius,i));
	h_mc_PhiEtaPrecut65[i]->GetXaxis()->SetTitle("PhiEta");
	h_mc_PhiEtaPrecut65[i]->GetXaxis()->CenterTitle();
	h_mc_PhiEtaPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiEtaPrecut65[i]->GetYaxis()->CenterTitle();
	h_mc_PhiEtaPrecut65[i]->DrawNormalized();
	h_data_PhiEtaPrecut65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phieta_precut_R%d_cent%d",radius,i));
	h_data_PhiEtaPrecut65[i]->SetLineColor(2);
	h_data_PhiEtaPrecut65[i]->DrawNormalized("same");	
	
	PhiEtaPrecutCan55->cd(i+1);
	h_mc_PhiEtaPrecut55[i] = (TH1F*)fmc->Get(Form("hpbpb_Jet55_phieta_precut_R%d_cent%d",radius,i));
	h_mc_PhiEtaPrecut55[i]->GetXaxis()->SetTitle("PhiEta");
	h_mc_PhiEtaPrecut55[i]->GetXaxis()->CenterTitle();
	h_mc_PhiEtaPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_mc_PhiEtaPrecut55[i]->GetYaxis()->CenterTitle();
	h_mc_PhiEtaPrecut55[i]->DrawNormalized();
	h_data_PhiEtaPrecut55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phieta_precut_R%d_cent%d",radius,i));
	h_data_PhiEtaPrecut55[i]->SetLineColor(2);
	h_data_PhiEtaPrecut55[i]->DrawNormalized("same");	
*/
//ratio plots

	EtaRatCan80->cd(i+1);
	h_mc_EtaRatio80[i] = new TH1F(Form("E80_mc_%d",i),Form("Eta 80 MC Cut/Precut"),30,-2.5,2.5);
	h_mc_EtaRatio80[i]->Divide(h_mc_Eta80[i],h_mc_EtaPrecut80[i]);	
	h_mc_EtaRatio80[i]->GetXaxis()->SetTitle("Eta");
	h_mc_EtaRatio80[i]->GetXaxis()->CenterTitle();
	h_mc_EtaRatio80[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	h_mc_EtaRatio80[i]->Draw();	
	h_data_EtaRatio80[i] = new TH1F(Form("E80_data_%d",i),Form("Eta 80 Data Cut/Precut"),30,-2.5,2.5);
	h_data_EtaRatio80[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	h_data_EtaRatio80[i]->Divide(h_data_Eta80[i],h_data_EtaPrecut80[i]);
	h_data_EtaRatio80[i]->SetLineColor(2);
	h_data_EtaRatio80[i]->Draw("same");	
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_EtaRatio80[0],"MC","pl");
    leg[i]->AddEntry(h_data_EtaRatio80[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	EtaRatCan65->cd(i+1);
	h_mc_EtaRatio65[i] = new TH1F(Form("E65_mc_%d",i),Form("Eta 65 MC Cut/Precut"),30,-2.5,2.5);
	h_mc_EtaRatio65[i]->Divide(h_mc_Eta65[i],h_mc_EtaPrecut65[i]);
	h_mc_EtaRatio65[i]->GetXaxis()->SetTitle("Eta");
	h_mc_EtaRatio65[i]->GetXaxis()->CenterTitle();
	h_mc_EtaRatio65[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	h_mc_EtaRatio65[i]->Draw();	
	h_data_EtaRatio65[i] = new TH1F(Form("E65_data_%d",i),Form("Eta 65 Data Cut/Precut"),30,-2.5,2.5);
	h_data_EtaRatio65[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	h_data_EtaRatio65[i]->Divide(h_data_Eta65[i],h_data_EtaPrecut65[i]);
	h_data_EtaRatio65[i]->SetLineColor(2);
	h_data_EtaRatio65[i]->Draw("same");	
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_EtaRatio65[0],"MC","pl");
    leg[i]->AddEntry(h_data_EtaRatio65[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
	EtaRatCan55->cd(i+1);
	h_mc_EtaRatio55[i] = new TH1F(Form("E55_mc_%d",i),Form("Eta 55 MC Cut/Precut"),30,-2.5,2.5);
	h_mc_EtaRatio55[i]->Divide(h_mc_Eta55[i],h_mc_EtaPrecut55[i]);
	h_mc_EtaRatio55[i]->GetXaxis()->SetTitle("Eta");
	h_mc_EtaRatio55[i]->GetXaxis()->CenterTitle();
	h_mc_EtaRatio55[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	h_mc_EtaRatio55[i]->Draw();	
	h_data_EtaRatio55[i] = new TH1F(Form("E55_data_%d",i),Form("Eta 55 Data Cut/Precut"),30,-2.5,2.5);
	h_data_EtaRatio55[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	h_data_EtaRatio55[i]->Divide(h_data_Eta55[i],h_data_EtaPrecut55[i]);
	h_data_EtaRatio55[i]->SetLineColor(2);
	h_data_EtaRatio55[i]->Draw("same");	
	leg[i] = new TLegend(0.15,0.15,0.25,0.25);
	leg[i]->AddEntry(h_mc_EtaRatio55[0],"MC","pl");
    leg[i]->AddEntry(h_data_EtaRatio55[0],"Data","pl");
    leg[i]->SetTextSize(0.04);
    leg[i]->Draw("same");
	
}
EtaCan80->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet80_30GeV_R%d_eta.pdf",radius),"RECREATE");
EtaCan65->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet65_30GeV_R%d_eta.pdf",radius),"RECREATE");
EtaCan55->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet55_30GeV_R%d_eta.pdf",radius),"RECREATE");
EtaPrecutCan80->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet80_30GeV_R%d_eta_precut.pdf",radius),"RECREATE");
EtaPrecutCan65->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet65_30GeV_R%d_eta_precut.pdf",radius),"RECREATE");
EtaPrecutCan55->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet55_30GeV_R%d_eta_precut.pdf",radius),"RECREATE");
EtaRatCan80->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet80_30GeV_R%d_etaratio.pdf",radius),"RECREATE");
EtaRatCan65->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet65_30GeV_R%d_etaratio.pdf",radius),"RECREATE");
EtaRatCan55->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet55_30GeV_R%d_etaratio.pdf",radius),"RECREATE");

PhiCan80->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet80_30GeV_R%d_phi.pdf",radius),"RECREATE");
PhiCan65->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet65_30GeV_R%d_phi.pdf",radius),"RECREATE");
PhiCan55->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet55_30GeV_R%d_phi.pdf",radius),"RECREATE");
PhiPrecutCan80->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet80_30GeV_R%d_phi_precut.pdf",radius),"RECREATE");
PhiPrecutCan65->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet65_30GeV_R%d_phi_precut.pdf",radius),"RECREATE");
PhiPrecutCan55->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/output/hpbpb_Jet55_30GeV_R%d_phi_precut.pdf",radius),"RECREATE");

}