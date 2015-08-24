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

void pbpb_Plot_RAA(){
TFile *fdata = new TFile(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_Data/PbPb_Data_histograms_FromForest_akPu%d_anabin_20_eta_20.root",radius),"r");
TFile *fmc = new TFile(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_20/newdrawfiles/PbPb_MC/PbPb_MC_histograms_FromForest_akPu%d_anabin_20_eta_20.root",radius),"r");

TCanvas *EtaRatCan55;
TCanvas *EtaRatCan65;
TCanvas *EtaRatCan80;

TH1F *h_data_EtaRatio80[nbins_cent];
TH1F *h_data_EtaRatio65[nbins_cent];
TH1F *h_data_EtaRatio55[nbins_cent];

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



}