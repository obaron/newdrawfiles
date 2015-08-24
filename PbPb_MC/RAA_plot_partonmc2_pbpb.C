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
TFile *f = new TFile("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/PbPb_MC_partonhistograms_FromForest_akPu3_20_eta_20.root","r");
TCanvas *EtaCan55[nbins_cent];
TCanvas *EtaCan65[nbins_cent];
TCanvas *EtaCan80[nbins_cent];

TH1F *hEta80[nbins_cent];
TH1F *hPhi80[nbins_cent];
TH1F *hPhieta80[nbins_cent];

TH1F *hEta65[nbins_cent];
TH1F *hPhi65[nbins_cent];
TH1F *hPhieta65[nbins_cent];

TH1F *hEta55[nbins_cent];
TH1F *hPhi55[nbins_cent];
TH1F *hPhieta55[nbins_cent];

TCanvas *ChargedMS80[nbins_cent];
TCanvas *ChargedMS65[nbins_cent];
TCanvas *ChargedMS55[nbins_cent];

/*TCanvas *ChargedSum80[nbins_cent];
TCanvas *ChargedSum65[nbins_cent];
TCanvas *ChargedSum55[nbins_cent];*/

TH1F *hChargedmax_q_80[nbins_cent];
TH1F *hChargedmax_g_80[nbins_cent];
TH1F *hChargedsum_q_80[nbins_cent];
TH1F *hChargedsum_g_80[nbins_cent];

TH1F *hChargedmax_q_65[nbins_cent];
TH1F *hChargedmax_g_65[nbins_cent];
TH1F *hChargedsum_q_65[nbins_cent];
TH1F *hChargedsum_g_65[nbins_cent];

TH1F *hChargedmax_q_55[nbins_cent];
TH1F *hChargedmax_g_55[nbins_cent];
TH1F *hChargedsum_q_55[nbins_cent];
TH1F *hChargedsum_g_55[nbins_cent];

for(int i=0;i<nbins_cent;i++){
//reduce margins if possible
	
	
	ChargedMS80[i] = new TCanvas(Form("ChargedMS80%d",i),Form("Jet80 Q and G ChargedSum and ChargedMax plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1500,500);
	ChargedMS80[i]->Divide(2,1);
	ChargedMS80[i]->cd(1);
	hChargedmax_q_80[i] = (TH1F*)f->Get(Form("hpbpb_g_chargedsum_R3_20_eta_20_cent%d",i));
	hChargedmax_g_80[i] = (TH1F*)f->Get(Form("hpbpb_g_chargedsum_R3_20_eta_20_cent%d",i));
	hChargedmax_q_80[i]->Draw();
	hChargedmax_g_80[i]->SetMarkerColor(2);
	hChargedmax_g_80[i]->Draw("same");
	ChargedMS80[i]->cd(2);
	hChargedsum_q_80[i] = (TH1F*)f->Get(Form("hpbpb_g_chargedsum_R3_20_eta_20_cent%d",i));
	hChargedsum_g_80[i] = (TH1F*)f->Get(Form("hpbpb_g_chargedsum_R3_20_eta_20_cent%d",i));
	hChargedsum_q_80[i]->Draw();
	hChargedsum_g_80[i]->SetMarkerColor(2);
	hChargedsum_g_80[i]->Draw("same");

	
	
	/* 
	//cout<<"start loop"<<endl;
	EtaCan80[i] = new TCanvas(Form("EtaCan80%d",i),Form("Jet80 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan80[i]->Divide(3,1);
	EtaCan80[i]->cd(1);
	hEta80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_etacent_R3_cent%d",i));
	hEta80[i]->GetXaxis()->SetTitle("Eta");
	hEta80[i]->GetXaxis()->CenterTitle();
	hEta80[i]->GetYaxis()->SetTitle("Multiplicity");
	hEta80[i]->GetYaxis()->CenterTitle();
	hEta80[i]->Draw();
	EtaCan80[i]->cd(2);
	hPhi80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_phicent_R3_cent%d",i));
	hPhi80[i]->GetXaxis()->SetTitle("Phi");
	hPhi80[i]->GetXaxis()->CenterTitle();
	hPhi80[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhi80[i]->GetYaxis()->CenterTitle();
	hPhi80[i]->Draw();
	EtaCan80[i]->cd(3);
	hPhieta80[i] = (TH1F*)f->Get(Form("hpbpb_Jet80_phieta_R3_cent%d",i));
	hPhieta80[i]->GetXaxis()->SetTitle("Eta");
	hPhieta80[i]->GetXaxis()->CenterTitle();
	hPhieta80[i]->GetYaxis()->SetTitle("Phi");
	hPhieta80[i]->GetYaxis()->CenterTitle();
	hPhieta80[i]->Draw("colz");

	EtaCan65[i] = new TCanvas(Form("EtaCan65%d",i),Form("Jet65 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan65[i]->Divide(3,1);
	EtaCan65[i]->cd(1);
	hEta65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_etacent_R3_cent%d",i));
	hEta65[i]->GetXaxis()->SetTitle("Eta");
	hEta65[i]->GetXaxis()->CenterTitle();
	hEta65[i]->GetYaxis()->SetTitle("Multiplicity");
	hEta65[i]->GetYaxis()->CenterTitle();
	hEta65[i]->Draw();
	EtaCan65[i]->cd(2);
	hPhi65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_phicent_R3_cent%d",i));
	hPhi65[i]->GetXaxis()->SetTitle("Phi");
	hPhi65[i]->GetXaxis()->CenterTitle();
	hPhi65[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhi65[i]->GetYaxis()->CenterTitle();
	hPhi65[i]->Draw();
	EtaCan65[i]->cd(3);
	hPhieta65[i] = (TH1F*)f->Get(Form("hpbpb_Jet65_phieta_R3_cent%d",i));
	hPhieta65[i]->GetXaxis()->SetTitle("Eta");
	hPhieta65[i]->GetXaxis()->CenterTitle();
	hPhieta65[i]->GetYaxis()->SetTitle("Phi");
	hPhieta65[i]->GetYaxis()->CenterTitle();
	hPhieta65[i]->Draw("colz");
	
	EtaCan55[i] = new TCanvas(Form("EtaCan55%d",i),Form("Jet55 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan55[i]->Divide(3,1);
	EtaCan55[i]->cd(1);
	hEta55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_etacent_R3_cent%d",i));
	hEta55[i]->GetXaxis()->SetTitle("Eta");
	hEta55[i]->GetXaxis()->CenterTitle();
	hEta55[i]->GetYaxis()->SetTitle("Multiplicity");
	hEta55[i]->GetYaxis()->CenterTitle();
	hEta55[i]->Draw();
	EtaCan55[i]->cd(2);
	hPhi55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_phicent_R3_cent%d",i));
	hPhi55[i]->GetXaxis()->SetTitle("Phi");
	hPhi55[i]->GetXaxis()->CenterTitle();
	hPhi55[i]->GetYaxis()->SetTitle("Multiplicity");
	hPhi55[i]->GetYaxis()->CenterTitle();
	hPhi55[i]->Draw();
	EtaCan55[i]->cd(3);
	hPhieta55[i] = (TH1F*)f->Get(Form("hpbpb_Jet55_phieta_R3_cent%d",i));
	hPhieta55[i]->GetXaxis()->SetTitle("Eta");
	hPhieta55[i]->GetXaxis()->CenterTitle();
	hPhieta55[i]->GetYaxis()->SetTitle("Phi");
	hPhieta55[i]->GetYaxis()->CenterTitle();
	hPhieta55[i]->Draw("colz");
	
	EtaCan80[i]->SaveAs(Form("hpbpb_Jet80_R%d_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaCan65[i]->SaveAs(Form("hpbpb_Jet65_R%d_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaCan55[i]->SaveAs(Form("hpbpb_Jet55_R%d_phieta_cent%d.pdf",radius,i),"RECREATE");
	
	//cout<<"end loop"<<endl;	
	*/
	
	
	
}

//pbpb plots here:
//TCanvas *t[radius][nbins_eta]

}
