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
const int radius = 4;

void RAA_plot_data_pbpb(){
TFile *f = new TFile(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/PbPb_Data_histograms_FromForest_akPu%d_anabin_20_eta_20.root",radius),"r");

TCanvas *EtaCan55[nbins_cent];
TCanvas *EtaCan65[nbins_cent];
TCanvas *EtaCan80[nbins_cent];
TCanvas *EtaPrecutCan55[nbins_cent];
TCanvas *EtaPrecutCan65[nbins_cent];
TCanvas *EtaPrecutCan80[nbins_cent];
TCanvas *EtaRatCan[nbins_cent];
TCanvas *pTCan[nbins_cent];


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

TH1F *EtaRatio80[nbins_cent];
TH1F *EtaRatio65[nbins_cent];
TH1F *EtaRatio55[nbins_cent];

TH1F *simplepT[nbins_cent];

for(int i=0;i<nbins_cent;i++){
//reduce margins if possible
	cout<<"make canvases"<<endl;
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
	hEta55[i]->GetYaxis()->SetRange(-0.5,2.5);
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
	EtaRatio65[i]->GetXaxis()->CenterTitle();
	EtaRatio65[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio65[i]->Draw();
	EtaRatCan[i]->cd(3);
	EtaRatio55[i] = new TH1F(Form("E55%d",i),Form("Eta 55 Cut/Precut"),30,-2.5,2.5);
	EtaRatio55[i]->Divide(hEta55[i],hEtaPrecut55[i]);
	EtaRatio55[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio55[i]->GetXaxis()->CenterTitle();
	EtaRatio55[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio55[i]->Draw();
	
	EtaCan80[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_Jet80_R%d_data_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaCan65[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_Jet65_R%d_data_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaCan55[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_Jet55_R%d_data_phieta_cent%d.pdf",radius,i),"RECREATE");
	EtaPrecutCan80[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_Jet80_R%d_data_phietaprecut_cent%d.pdf",radius,i),"RECREATE");
	EtaPrecutCan65[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_Jet65_R%d_data_phietaprecut_cent%d.pdf",radius,i),"RECREATE");
	EtaPrecutCan55[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_Jet55_R%d_data_phietaprecut_cent%d.pdf",radius,i),"RECREATE");
	EtaRatCan[i]->SaveAs(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/output/hpbpb_etaratio_R%d_data_cent%d.pdf",radius,i),"RECREATE");

	
	cout<<"end loop"<<endl;	
}

//pbpb plots here:
//TCanvas *t[radius][nbins_eta]

}
