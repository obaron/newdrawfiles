// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN
// for questions or comments: raghav.k.e at CERN dot CH

// 
// read all the MC files for PbPb and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
// 

// July 19 - all pp histograms will have 2D arrays with [radius][eta_bin]. the PbPb histograms will be defined by 3D arrays with [radius][eta_bin][centrality]. 
// July 20 - the loop structure(s) are defined as follows for the several histograms 
//            Radius    : iteration variable: k;       number of iterations: no_radius;                                 values for the radii: list_radius
//            Eta bins  : iteration variable: j;       number of iterations: nbins_eta;                                 values of the bins  : boundaries_eta
//            Centrality: iteration variable: i;       number of iterations: nbins_cent +1;                             values of the bins  : boundaries_cent (till nbins_cent) + 0-200 (the whole range) for the final iteration. 
//            p_T Hats  : iteration variable: h;       number of iterations: nbins_pthat (PbPb) and nbinsPP_pthat (pp); values of the bins  : boundaries_pthat (PbPb) and boundariesPP_pthat (pp)  
//            jets      : iteration variable: g;       number of iterations: no of jets in Data[k][h];
//            p_T       : defined just below as nbins_pt with 39 bins. to match our NLO and jet RpA analysis bins. 

// Oct 23 - removed the cuts from the MC -> like the noisefilter etc... 

// Nov 4th - added the supernova event cut rejection based on the no of hits in the pixel. 

// Dec 9th - going to PU for the Jet RAA. 

// Dec 17th - changing the file list to smaller 50k files on which JEC were derived to check for PF electron problems, requested by Marguerite.

// Jan 13th 2015 - adding in the official pp mc (from Dragos) 
//               - this is going to be a bit tricky since each file is split up into 4 smaller files. so each pthat will have a TChain!


// Feb 12th - cleaned up the macro to make it usable (hopefuly) by others.

// Jun 22th - going back to the HiForest with the trees from Pawan's for the event selection cuts and PF electron cuts. 

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

#define pi 3.14159265

static const int nbins_pt = 39;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300, 
  330, 362, 395, 430,
  468, 507, 548, 592,
  638, 686, 1000 
};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 3; 
//static const int list_radius[no_radius] = {2,3,4};

static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 
static const Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
static const Float_t effecPrescl = 2.047507;
static const char * etaWidth = (char*)"20_eta_20";

int findBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  if(bin<10)ibin=0; //! 0-5%
  else if(bin>=10  && bin<20 )ibin=1; //! 5-10%
  else if(bin>=20  && bin<60 )ibin=2;  //! 10-30%
  else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
  else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
  else if(bin>=140 && bin<180)ibin=5;  //! 70-90%
  else if(bin>=180 && bin<200)ibin=6;  //! 90-100%
  return ibin;
}

static const double cent_HF_bound[] = {0, 0, 4.18352, 6.93443, 7.81838, 8.54289, 9.23664, 9.88781, 10.5473, 11.1902, 11.8706, 12.5891, 13.319, 14.0674, 14.8253, 15.6046, 16.4299, 17.2737, 18.1496, 19.0479, 20.0015, 20.912, 21.9531, 23.0394, 24.1431, 25.328, 26.5151, 27.8044, 29.0709, 30.4214, 31.8381, 33.283, 34.6493, 36.0984, 37.6395, 39.2199, 40.9588, 42.6406, 44.569, 46.3474, 48.2599, 50.2071, 52.3226, 54.5107, 56.9493, 59.3141, 61.7028, 64.244, 66.8831, 69.5266, 72.3771, 75.4711, 78.4157, 81.5896, 84.769, 88.2238, 91.8252, 95.4077, 98.98, 102.798, 106.782, 110.906, 115.16, 119.625, 124.341, 129.065, 133.87, 138.692, 143.708, 149.042, 154.287, 159.825, 165.462, 171.257, 176.947, 183.167, 189.585, 195.883, 202.593, 209.261, 216.241, 223.538, 231.183, 238.531, 246.204, 254.287, 262.356, 270.439, 278.961, 287.651, 297.26, 306.368, 315.948, 325.339, 335.032, 345.099, 355.337, 365.822, 376.405, 387.344, 399.021, 411.008, 422.029, 433.915, 446.391, 458.551, 470.541, 483.077, 495.431, 508.121, 521.027, 534.43, 548.443, 562.924, 576.951, 592.656, 608.19, 623.398, 639.439, 654.628, 670.442, 686.71, 702.497, 719.216, 736.875, 753.611, 771.997, 790, 809.284, 827.785, 846.916, 866.011, 885.882, 904.369, 924.095, 944.936, 966.012, 987.122, 1009.07, 1031.27, 1054.05, 1076.25, 1099.91, 1122.04, 1144.78, 1167.66, 1191.41, 1216.09, 1241.13, 1267.38, 1293.67, 1321.55, 1350.66, 1378.46, 1407.54, 1434.55, 1462.55, 1490.14, 1520.59, 1550.78, 1581.24, 1613.24, 1645.85, 1678.99, 1711.05, 1745.42, 1778.74, 1812.54, 1846.99, 1884.43, 1922.16, 1958.71, 1997.48, 2037.29, 2076.29, 2114.1, 2155.37, 2196.62, 2237.9, 2279.55, 2324.23, 2371.58, 2415.58, 2465.63, 2515.93, 2562.36, 2612.32, 2663.82, 2719.3, 2770.5, 2829.21, 2890.39, 2948.11, 3011.59, 3073.29, 3135.95, 3202.43, 3270.81, 3340.31, 3429.82, 5805.99};
static const int nbins_HF_bound = sizeof(cent_HF_bound)/sizeof(double) -1;

int findHFbin(float HF_energy){
  int ibin = -1;
  for(int k = 0; k<nbins_HF_bound; ++k){
    if(HF_energy > cent_HF_bound[k])
      ibin = 200 - k;
  }
  return ibin;
}


// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

struct Jet{
  int id;
  float pt;
};
bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}

using namespace std;

void RAA_read_data_pbpb(int startfile = 0,
			int endfile = 1,
			int radius = 4,
			std::string kFoname="test_output.root",
			std::string kFonametxt1 = "test1.txt",
			std::string kFonametxt2 = "test2.txt"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;
  std::string infile_Select;

  infile_Forest = "jetRAA_PbPb_data_forest.txt";
  infile_Select = Form("jetRAA_PbPb_data_akPu%d_select.txt",radius);
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  std::ifstream instr_Select(infile_Select.c_str(),std::ifstream::in);
  std::string filename_Select;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
    instr_Select>>filename_Select;
  }

  const int N = 5; //6
  float pfptmin=30;

  TChain * jetpbpb[N];
  TChain * evt_select, * jet_select; 

  Float_t ScaleFactor_55_NeqScalePerCent[nbins_cent] = {60.1 * 1e6, 54.6 * 1e6, 219.3 * 1e6, 221.8 * 1e6, 219.4 * 1e6, 224.3 * 1e6};
  //Float_t ScaleFactor_55_NeqScalePerCent[nbins_cent] = {60.1 * 1e6, 54.6 * 1e6, 219.3 * 1e6, 221.8 * 1e6, 219.4 * 1e6, 221 * 1e6};
  Float_t ScaleFactor_65_NeqScalePerCent[nbins_cent] = {55.4 * 1e6, 50.4 * 1e6, 197 * 1e6, 201 * 1e6, 196 * 1e6, 189 * 1e6};
  Float_t ScaleFactor_80_NeqScalePerCent[nbins_cent] = {59 * 1e6, 51.9 * 1e6, 216 * 1e6, 221 * 1e6, 221 * 1e6, 247 * 1e6};
  //Float_t ScaleFactor_80_NeqScalePerCent[nbins_cent] = {59 * 1e6, 51.9 * 1e6, 216 * 1e6, 221 * 1e6, 221 * 1e6, 221 * 1e6};

  Float_t ScaleFactor_944Scale[nbins_cent] = {944 * 1e6 * 0.05,
					      944 * 1e6 * 0.05,
					      944 * 1e6 * 0.20,
					      944 * 1e6 * 0.20,
					      944 * 1e6 * 0.20,
					      944 * 1e6 * 0.20};

  // Float_t ScaleFactor_55_NeqScale[nbins_cent] = {54.7 * 1e6, 54.7 * 1e6, 218.8 * 1e6, 218.8 * 1e6, 218.8 * 1e6, 218.8 * 1e6}; 
  // Float_t ScaleFactor_65_NeqScale[nbins_cent] = {50.1 * 1e6, 50.1 * 1e6, 200.4 * 1e6, 200.4 * 1e6, 200.4 * 1e6, 200.4 * 1e6};
  // Float_t ScaleFactor_80_NeqScale[nbins_cent] = {53.83 * 1e6, 53.83 * 1e6, 215.3 * 1e6, 215.3 * 1e6, 215.3 * 1e6, 215.3 * 1e6};

  // Float_t ScaleFactor_55_NeqScale[nbins_cent] = {60.1 * 1e6, 54.7 * 1e6, 218.8 * 1e6, 218.8 * 1e6, 218.8 * 1e6, 218.8 * 1e6}; 
  // Float_t ScaleFactor_65_NeqScale[nbins_cent] = {55.4 * 1e6, 50.1 * 1e6, 200.4 * 1e6, 200.4 * 1e6, 200.4 * 1e6, 200.4 * 1e6};
  // Float_t ScaleFactor_80_NeqScale[nbins_cent] = {59.0 * 1e6, 53.83 * 1e6, 215.3 * 1e6, 215.3 * 1e6, 215.3 * 1e6, 215.3 * 1e6};

  // Float_t ScaleFactor_55_NeqScale[nbins_cent] = {54.5 * 1e6, 54.5 * 1e6, 218.8 * 1e6, 218.8 * 1e6, 218.8 * 1e6, 218.8 * 1e6}; 
  // Float_t ScaleFactor_65_NeqScale[nbins_cent] = {50.1 * 1e6, 50.1 * 1e6, 200.4 * 1e6, 200.4 * 1e6, 200.4 * 1e6, 200.4 * 1e6};
  // Float_t ScaleFactor_80_NeqScale[nbins_cent] = {53.8 * 1e6, 53.8 * 1e6, 215.3 * 1e6, 215.3 * 1e6, 215.3 * 1e6, 215.3 * 1e6};

  // Float_t ScaleFactor_80_NeqScale[nbins_cent] = {962.01 * 0.05 * 1e6,
  // 						 962.01 * 0.05 * 1e6,
  // 						 962.01 * 0.2 * 1e6,
  // 						 962.01 * 0.2 * 1e6,
  // 						 962.01 * 0.2 * 1e6,
  // 						 962.01 * 0.2 * 1e6};
  // Float_t ScaleFactor_65_NeqScale[nbins_cent] = {893.46 * 0.05 * 1e6,
  // 						 893.46 * 0.05 * 1e6,
  // 						 893.46 * 0.2 * 1e6,
  // 						 893.46 * 0.2 * 1e6,
  // 						 893.46 * 0.2 * 1e6,
  // 						 893.46 * 0.2 * 1e6};
  // Float_t ScaleFactor_55_NeqScale[nbins_cent] = {977.89 * 0.05 * 1e6,
  // 						 977.89 * 0.05 * 1e6,
  // 						 977.89 * 0.2 * 1e6,
  // 						 977.89 * 0.2 * 1e6,
  // 						 977.89 * 0.2 * 1e6,
  // 						 977.89 * 0.2 * 1e6};

  Float_t ScaleFactor_55_NeqScale[nbins_cent] = {53.6 * 1e6, 53.6 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6}; 
  Float_t ScaleFactor_65_NeqScale[nbins_cent] = {53.6 * 1e6, 53.6 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6}; 
  Float_t ScaleFactor_80_NeqScale[nbins_cent] = {53.6 * 1e6, 53.6 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6}; 
  
  
  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("akPu%dPFJetAnalyzer",radius);
  dir[3] = "akPu3CaloJetAnalyzer";
  dir[4] = "hiEvtAnalyzer";
  // dir[4] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "t",
    "HiTree"
    // , "jetObjTree"
  };

  for(int t = 0;t<N;t++){
    jetpbpb[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  evt_select = new TChain(Form("akPu%dJetAnalyzer/evtTree",radius));
  jet_select = new TChain(Form("akPu%dJetAnalyzer/jetTree",radius));

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;
    instr_Select>>filename_Select;

    if(printDebug)cout<<"HiForest filename = "<<filename_Forest.c_str()<<endl;

    jetpbpb[0]->Add(filename_Forest.c_str());
    jetpbpb[1]->Add(filename_Forest.c_str());
    jetpbpb[2]->Add(filename_Forest.c_str());
    jetpbpb[3]->Add(filename_Forest.c_str());
    jetpbpb[4]->Add(filename_Forest.c_str());
    jet_select->Add(filename_Select.c_str());
    evt_select->Add(filename_Select.c_str());
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[3]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[4]->GetEntries() << endl;
    if(printDebug)cout << "jet selection file " << jet_select->GetEntries() << endl;
    if(printDebug)cout << "event selection file" << evt_select->GetEntries() << endl;

  }
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[4]);
  jetpbpb[3]->AddFriend(jetpbpb[0]);
  jetpbpb[3]->AddFriend(jetpbpb[1]);
  jetpbpb[3]->AddFriend(jetpbpb[4]);
  
  jetpbpb[2]->AddFriend(evt_select);
  jetpbpb[3]->AddFriend(evt_select);

  // Forest files 
  int nref_F;
  float pt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float chMax_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];
  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  int jet55_F;
  int jet65_F;
  int jet80_F;
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet55_p_F;
  int jet65_p_F;
  int jet80_p_F;
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int hiNpix_F;
  int hiNpixelTracks_F;
  int hiBin_F;
  float hiHF_F;
  float hiZDC_F;
  float hiZDCplus_F;
  float hiZDCminus_F;
  int hiNtracks_F;
  int hiNtracksPtCut_F;
  int hiNtracksEtaCut_F;
  int hiNtracksEtaPtCut_F;
  int pcollisionEventSelection_F;
  int pHBHENoiseFilter_F;

  float calopt_F[1000];
  jetpbpb[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpbpb[4]->SetBranchAddress("evt",&evt_F);
  jetpbpb[4]->SetBranchAddress("run",&run_F);
  jetpbpb[4]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[4]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[4]->SetBranchAddress("hiHF", &hiHF_F);
  jetpbpb[4]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetpbpb[4]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jetpbpb[4]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jetpbpb[4]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jetpbpb[4]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jetpbpb[4]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);
  jetpbpb[4]->SetBranchAddress("hiZDC", &hiZDC_F);
  jetpbpb[4]->SetBranchAddress("hiZDCplus", &hiZDCplus_F);
  jetpbpb[4]->SetBranchAddress("hiZDCminus", &hiZDCminus_F);
  jetpbpb[4]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  jetpbpb[1]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpbpb[2]->SetBranchAddress("nref",&nref_F);
  jetpbpb[2]->SetBranchAddress("jtpt",pt_F);
  jetpbpb[2]->SetBranchAddress("jteta",eta_F);
  jetpbpb[2]->SetBranchAddress("jtphi",phi_F);
  jetpbpb[2]->SetBranchAddress("rawpt",rawpt_F);
  jetpbpb[2]->SetBranchAddress("jtpu",jtpu_F);
  jetpbpb[2]->SetBranchAddress("chargedMax",chMax_F);
  jetpbpb[2]->SetBranchAddress("chargedSum",chSum_F);
  jetpbpb[2]->SetBranchAddress("trackMax",trkMax_F);
  jetpbpb[2]->SetBranchAddress("trackSum",trkSum_F);
  jetpbpb[2]->SetBranchAddress("photonMax",phMax_F);
  jetpbpb[2]->SetBranchAddress("photonSum",phSum_F);
  jetpbpb[2]->SetBranchAddress("neutralMax",neMax_F);
  jetpbpb[2]->SetBranchAddress("neutralSum",neSum_F);
  jetpbpb[2]->SetBranchAddress("eSum",eSum_F);
  jetpbpb[2]->SetBranchAddress("eMax",eMax_F);
  jetpbpb[2]->SetBranchAddress("muSum",muSum_F);
  jetpbpb[2]->SetBranchAddress("muMax",muMax_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v1",&jet55_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v1",&jet65_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v1",&jet80_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  // event tree selection file:
  int hiBin_eS;
  int run_eS;
  int evt_eS;
  int lumi_eS;
  float vz_eS; 
  int isGoodEvent_eS;
  int nref_eS;
  int isMiMatch_eS[1000];
  int isClMatch_eS[1000];
  int isPFElecCut_eS[1000];
  int isTrackCut_eS[1000];
  int isMuCut_eS[1000];
  int index_eS[1000];
  float pt_eS[1000];
  float calopt_eS[1000];
  Double_t weight_eS;

  evt_select->SetBranchAddress("hiBin",&hiBin_eS);
  evt_select->SetBranchAddress("run_value",&run_eS);
  evt_select->SetBranchAddress("evt_value",&evt_eS);
  evt_select->SetBranchAddress("lumi_value",&lumi_eS);
  evt_select->SetBranchAddress("vz",&vz_eS);
  evt_select->SetBranchAddress("weight", &weight_eS);  
  evt_select->SetBranchAddress("isGoodEvt",&isGoodEvent_eS);
  evt_select->SetBranchAddress("nref",&nref_eS);
  evt_select->SetBranchAddress("index", index_eS);
  evt_select->SetBranchAddress("isMlMatch",isMiMatch_eS);
  evt_select->SetBranchAddress("isClMatch",isClMatch_eS);
  evt_select->SetBranchAddress("isPFElecCut",isPFElecCut_eS);
  evt_select->SetBranchAddress("isTrackCut",isTrackCut_eS);
  evt_select->SetBranchAddress("isMuCut",isMuCut_eS);
  evt_select->SetBranchAddress("pfpt",pt_eS);
  evt_select->SetBranchAddress("calopt",calopt_eS);

  // jet tree selection file:
  int hiBin_jS;
  int run_jS;
  int evt_jS;
  int lumi_jS;
  float vz_jS;
  int nref_jS;
  float pt_jS[1000];
  float eta_jS[1000];
  float eMax_jS[1000];
  float chSum_jS[1000];
  float phSum_jS[1000];
  float neSum_jS[1000];
  float muSum_jS[1000];
  float calopt_jS[1000];
  int   isCaloMatch_jS[1000];
  int   isMultiMatch_jS[1000];

  jet_select->SetBranchAddress("hiBin",&hiBin_jS);
  jet_select->SetBranchAddress("run_value",&run_jS);
  jet_select->SetBranchAddress("evt_value",&evt_jS);
  jet_select->SetBranchAddress("lumi_value",&lumi_jS);
  jet_select->SetBranchAddress("vz",&vz_jS);
  jet_select->SetBranchAddress("npf", &nref_jS);  
  jet_select->SetBranchAddress("pfpt", &pt_jS);  
  jet_select->SetBranchAddress("eMax", &eMax_jS);  
  jet_select->SetBranchAddress("pfeta", &eta_jS);  
  jet_select->SetBranchAddress("calopt", &calopt_jS);  
  jet_select->SetBranchAddress("isCaloMatch", &isCaloMatch_jS);  
  jet_select->SetBranchAddress("isMultiMatch", &isMultiMatch_jS);  
  jet_select->SetBranchAddress("chSum", &chSum_jS);  
  jet_select->SetBranchAddress("phSum", &phSum_jS);  
  jet_select->SetBranchAddress("neSum", &neSum_jS);  
  jet_select->SetBranchAddress("muSum", &muSum_jS);  

  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

TH1F *hpbpb_Jet80_eta_cent[nbins_cent],
	*hpbpb_Jet80_phi_cent[nbins_cent],
	*hpbpb_Jet80_eta_precut_cent[nbins_cent],
	*hpbpb_Jet80_phi_precut_cent[nbins_cent];
	
TH1F *hpbpb_Jet65_eta_cent[nbins_cent],
	*hpbpb_Jet65_phi_cent[nbins_cent],
	*hpbpb_Jet65_eta_precut_cent[nbins_cent],
	*hpbpb_Jet65_phi_precut_cent[nbins_cent];

TH1F *hpbpb_Jet55_eta_cent[nbins_cent],
	*hpbpb_Jet55_phi_cent[nbins_cent],
	*hpbpb_Jet55_eta_precut_cent[nbins_cent],
	*hpbpb_Jet55_phi_precut_cent[nbins_cent];
	
TH2F *hpbpb_Jet80_phieta_cent[nbins_cent],
	* hpbpb_Jet65_phieta_cent[nbins_cent],
	* hpbpb_Jet55_phieta_cent[nbins_cent],
	* hpbpb_Jet80_phieta_precut_cent[nbins_cent],
	* hpbpb_Jet65_phieta_precut_cent[nbins_cent],
	* hpbpb_Jet55_phieta_precut_cent[nbins_cent];

  // TH2F * hNpix_vs_hiBin = new TH2F("hNpix_vs_hiBin","",200, 0, 200, 120, 0, 60000);
  // TH2F * hNpixelTracks_vs_hiBin = new TH2F("hNpixelTracks_vs_hiBin","",200, 0, 200, 80, 0, 8000);
  // TH2F * hNTracks_vs_hiBin = new TH2F("hNTracks_vs_hiBin","",200, 0, 200, 100, 0, 5000);
  // TH2F * hNTracksPtCut_vs_hiBin = new TH2F("hNTracksPtCut_vs_hiBin","",200, 0, 200, 100, 0, 5000);
  // TH2F * hNTracksEtaCut_vs_hiBin = new TH2F("hNTracksEtaCut_vs_hiBin","",200, 0, 200, 100, 0, 5000);
  // TH2F * hNTracksEtaPtCut_vs_hiBin = new TH2F("hNTracksEtaPtCut_vs_hiBin","",200, 0, 200, 100, 0, 5000);

  TH2F * hNpix_vs_NJets_noCut = new TH2F("hNpix_vs_Njets_noCut","",50, 0, 100, 120, 0, 60000);
  TH2F * hNpix_vs_HF_noCut = new TH2F("hNpix_vs_HF_noCut","",200, 0, 8000, 120, 0, 60000);

  TH2F * hZDC_vs_HF_noCut = new TH2F("hZDC_vs_HF_noCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDCplus_vs_HF_noCut = new TH2F("hZDCplus_vs_HF_noCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDCminus_vs_HF_noCut = new TH2F("hZDCminus_vs_HF_noCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDC_vs_HF_cent0to0p5_noCut = new TH2F("hZDC_vs_HF_cent0to0p5_noCut","",200, 0, 8000, 800, 0, 480);
  // pile up cut equation
  // if (ZDC > 1200 - 4/15 * HF energy) continue;  
  TH2F * hZDCp_vs_HF_noCut = new TH2F("hZDCp_vs_HF_noCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDCp_vs_HF_cent0to0p5_noCut = new TH2F("hZDCp_vs_HF_cent0to0p5_noCut","",200, 0, 8000, 800, 0, 480);
  // pile up cut equation
  // select all of them:   
  TH2F * hZDCm_vs_HF_noCut = new TH2F("hZDCm_vs_HF_noCut","",200, 0, 8000, 300, 0, 120);
  TH2F * hZDCm_vs_HF_cent0to0p5_noCut = new TH2F("hZDCm_vs_HF_cent0to0p5_noCut","",200, 0, 8000, 300, 0, 120);
  // pile up cut equation
  // if (ZDCp > 250 - 1/20 * HF energy) continue;  
  TH2F * hNpixelTracks_vs_HF_noCut = new TH2F("hNpixelTracks_vs_HF_noCut","",200, 0, 8000, 80, 0, 8000);
  TH2F * hNTracks_vs_HF_noCut = new TH2F("hNTracks_vs_HF_noCut","",200, 0, 8000, 100, 0, 5000);
  TH1F * hJetComb_Npix_noCut = new TH1F("hJetComb_Npix_noCut","",120, 0, 60000);
  TH1F * hJetComb_vs_HF_noCut = new TH1F("hJetComb_vs_HF_noCut","",200, 0, 8000);
  TH1F * hJet80_vs_HF_noCut = new TH1F("hJet80_vs_HF_noCut","",200, 0, 8000);
  TH1F * hJet65_vs_HF_noCut = new TH1F("hJet65_vs_HF_noCut","",200, 0, 8000);
  TH1F * hJet55_vs_HF_noCut = new TH1F("hJet55_vs_HF_noCut","",200, 0, 8000);
  TH1F * hCentrality_fromHFBound_noCut = new TH1F("hCentrality_fromHFBound_noCut","",200, 0, 200);


  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************

  TH2F * hNpix_vs_NJets_SupernovaCut = new TH2F("hNpix_vs_Njets_SupernovaCut","",50, 0, 100, 120, 0, 60000);
  TH2F * hNpix_vs_HF_SupernovaCut = new TH2F("hNpix_vs_HF_SupernovaCut","",200, 0, 8000, 120, 0, 60000);

  TH2F * hZDC_vs_HF_SupernovaCut = new TH2F("hZDC_vs_HF_SupernovaCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDC_vs_HF_cent0to0p5_SupernovaCut = new TH2F("hZDC_vs_HF_cent0to0p5_SupernovaCut","",200, 0, 8000, 800, 0, 480);
  // pile up cut equation
  // if (ZDC > 1200 - 4/15 * HF energy) continue;  
  TH2F * hZDCp_vs_HF_SupernovaCut = new TH2F("hZDCp_vs_HF_SupernovaCut","",200, 0, 8000, 50, 0, 5);
  TH2F * hZDCp_vs_HF_cent0to0p5_SupernovaCut = new TH2F("hZDCp_vs_HF_cent0to0p5_SupernovaCut","",200, 0, 8000, 50, 0, 5);
  // pile up cut equation
  // select all of them:   
  TH2F * hZDCm_vs_HF_SupernovaCut = new TH2F("hZDCm_vs_HF_SupernovaCut","",200, 0, 8000, 300, 0, 120);
  TH2F * hZDCm_vs_HF_cent0to0p5_SupernovaCut = new TH2F("hZDCm_vs_HF_cent0to0p5_SupernovaCut","",200, 0, 8000, 300, 0, 120);
  // pile up cut equation
  // if (ZDCp > 250 - 1/20 * HF energy) continue;  
  TH2F * hNpixelTracks_vs_HF_SupernovaCut = new TH2F("hNpixelTracks_vs_HF_SupernovaCut","",200, 0, 8000, 80, 0, 8000);
  TH2F * hNTracks_vs_HF_SupernovaCut = new TH2F("hNTracks_vs_HF_SupernovaCut","",200, 0, 8000, 100, 0, 5000);
  TH1F * hJetComb_Npix_SupernovaCut = new TH1F("hJetComb_Npix_SupernovaCut","",120, 0, 60000);
  TH1F * hJetComb_vs_HF_SupernovaCut = new TH1F("hJetComb_vs_HF_SupernovaCut","",200, 0, 8000);
  TH1F * hJet80_vs_HF_SupernovaCut = new TH1F("hJet80_vs_HF_SupernovaCut","",200, 0, 8000);
  TH1F * hJet65_vs_HF_SupernovaCut = new TH1F("hJet65_vs_HF_SupernovaCut","",200, 0, 8000);
  TH1F * hJet55_vs_HF_SupernovaCut = new TH1F("hJet55_vs_HF_SupernovaCut","",200, 0, 8000);
  TH1F * hCentrality_fromHFBound_SupernovaCut = new TH1F("hCentrality_fromHFBound_SupernovaCut","",200, 0, 200);

  
  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************

  TH2F * hNpix_vs_NJets_PileUpCut = new TH2F("hNpix_vs_Njets_PileUpCut","",50, 0, 100, 120, 0, 60000);
  TH2F * hNpix_vs_HF_PileUpCut = new TH2F("hNpix_vs_HF_PileUpCut","",200, 0, 8000, 120, 0, 60000);
  
  TH2F * hZDC_vs_HF_PileUpCut = new TH2F("hZDC_vs_HF_PileUpCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDC_vs_HF_cent0to0p5_PileUpCut = new TH2F("hZDC_vs_HF_cent0to0p5_PileUpCut","",200, 0, 8000, 800, 0, 480);
  // pile up cut equation
  // if (ZDC > 1200 - 4/15 * HF energy) continue;  
  TH2F * hZDCp_vs_HF_PileUpCut = new TH2F("hZDCp_vs_HF_PileUpCut","",200, 0, 8000, 50, 0, 5);
  TH2F * hZDCp_vs_HF_cent0to0p5_PileUpCut = new TH2F("hZDCp_vs_HF_cent0to0p5_PileUpCut","",200, 0, 8000, 50, 0, 5);
  // pile up cut equation
  // select all of them:   
  TH2F * hZDCm_vs_HF_PileUpCut = new TH2F("hZDCm_vs_HF_PileUpCut","",200, 0, 8000, 300, 0, 120);
  TH2F * hZDCm_vs_HF_cent0to0p5_PileUpCut = new TH2F("hZDCm_vs_HF_cent0to0p5_PileUpCut","",200, 0, 8000, 300, 0, 120);
  // pile up cut equation
  // if (ZDCp > 250 - 1/20 * HF energy) continue;  
  TH2F * hNpixelTracks_vs_HF_PileUpCut = new TH2F("hNpixelTracks_vs_HF_PileUpCut","",200, 0, 8000, 80, 0, 8000);
  TH2F * hNTracks_vs_HF_PileUpCut = new TH2F("hNTracks_vs_HF_PileUpCut","",200, 0, 8000, 100, 0, 5000);
  TH1F * hJetComb_Npix_PileUpCut = new TH1F("hJetComb_Npix_PileUpCut","",120, 0, 60000);
  TH1F * hJetComb_vs_HF_PileUpCut = new TH1F("hJetComb_vs_HF_PileUpCut","",200, 0, 8000);
  TH1F * hJet80_vs_HF_PileUpCut = new TH1F("hJet80_vs_HF_PileUpCut","",200, 0, 8000);
  TH1F * hJet65_vs_HF_PileUpCut = new TH1F("hJet65_vs_HF_PileUpCut","",200, 0, 8000);
  TH1F * hJet55_vs_HF_PileUpCut = new TH1F("hJet55_vs_HF_PileUpCut","",200, 0, 8000);
  TH1F * hCentrality_fromHFBound_PileUpCut = new TH1F("hCentrality_fromHFBound_PileUpCut","",200, 0, 200);

  
  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************

  TH2F * hNpix_vs_NJets_ZDCLowCut = new TH2F("hNpix_vs_Njets_ZDCLowCut","",50, 0, 100, 120, 0, 60000);
  TH2F * hNpix_vs_HF_ZDCLowCut = new TH2F("hNpix_vs_HF_ZDCLowCut","",200, 0, 8000, 120, 0, 60000);
  
  TH2F * hZDC_vs_HF_ZDCLowCut = new TH2F("hZDC_vs_HF_ZDCLowCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDC_vs_HF_cent0to0p5_ZDCLowCut = new TH2F("hZDC_vs_HF_cent0to0p5_ZDCLowCut","",200, 0, 8000, 800, 0, 480);
  // pile up cut equation
  // if (ZDC > 1200 - 4/15 * HF energy) continue;  
  TH2F * hZDCp_vs_HF_ZDCLowCut = new TH2F("hZDCp_vs_HF_ZDCLowCut","",200, 0, 8000, 50, 0, 5);
  TH2F * hZDCp_vs_HF_cent0to0p5_ZDCLowCut = new TH2F("hZDCp_vs_HF_cent0to0p5_ZDCLowCut","",200, 0, 8000, 50, 0, 5);
  // pile up cut equation
  // select all of them:   
  TH2F * hZDCm_vs_HF_ZDCLowCut = new TH2F("hZDCm_vs_HF_ZDCLowCut","",200, 0, 8000, 300, 0, 120);
  TH2F * hZDCm_vs_HF_cent0to0p5_ZDCLowCut = new TH2F("hZDCm_vs_HF_cent0to0p5_ZDCLowCut","",200, 0, 8000, 300, 0, 120);
  // pile up cut equation
  // if (ZDCp > 250 - 1/20 * HF energy) continue;  
  TH2F * hNpixelTracks_vs_HF_ZDCLowCut = new TH2F("hNpixelTracks_vs_HF_ZDCLowCut","",200, 0, 8000, 80, 0, 8000);
  TH2F * hNTracks_vs_HF_ZDCLowCut = new TH2F("hNTracks_vs_HF_ZDCLowCut","",200, 0, 8000, 100, 0, 5000);
  TH1F * hJetComb_Npix_ZDCLowCut = new TH1F("hJetComb_Npix_ZDCLowCut","",120, 0, 60000);
  TH1F * hJetComb_vs_HF_ZDCLowCut = new TH1F("hJetComb_vs_HF_ZDCLowCut","",200, 0, 8000);
  TH1F * hJet80_vs_HF_ZDCLowCut = new TH1F("hJet80_vs_HF_ZDCLowCut","",200, 0, 8000);
  TH1F * hJet65_vs_HF_ZDCLowCut = new TH1F("hJet65_vs_HF_ZDCLowCut","",200, 0, 8000);
  TH1F * hJet55_vs_HF_ZDCLowCut = new TH1F("hJet55_vs_HF_ZDCLowCut","",200, 0, 8000);
  TH1F * hCentrality_fromHFBound_ZDCLowCut = new TH1F("hCentrality_fromHFBound_ZDCLowCut","",200, 0, 200);
  
  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************

  
  TH2F * hNpix_vs_NJets_AllCut = new TH2F("hNpix_vs_Njets_AllCut","",50, 0, 100, 120, 0, 60000);
  TH2F * hNpix_vs_HF_AllCut = new TH2F("hNpix_vs_HF_AllCut","",200, 0, 8000, 120, 0, 60000);
  
  TH2F * hZDC_vs_HF_AllCut = new TH2F("hZDC_vs_HF_AllCut","",200, 0, 8000, 800, 0, 480);
  TH2F * hZDC_vs_HF_cent0to0p5_AllCut = new TH2F("hZDC_vs_HF_cent0to0p5_AllCut","",200, 0, 8000, 800, 0, 480);
  // pile up cut equation
  // if (ZDC > 1200 - 4/15 * HF energy) continue;  
  TH2F * hZDCp_vs_HF_AllCut = new TH2F("hZDCp_vs_HF_AllCut","",200, 0, 8000, 50, 0, 5);
  TH2F * hZDCp_vs_HF_cent0to0p5_AllCut = new TH2F("hZDCp_vs_HF_cent0to0p5_AllCut","",200, 0, 8000, 50, 0, 5);
  // pile up cut equation
  // select all of them:   
  TH2F * hZDCm_vs_HF_AllCut = new TH2F("hZDCm_vs_HF_AllCut","",200, 0, 8000, 300, 0, 120);
  TH2F * hZDCm_vs_HF_cent0to0p5_AllCut = new TH2F("hZDCm_vs_HF_cent0to0p5_AllCut","",200, 0, 8000, 300, 0, 120);
  // pile up cut equation
  // if (ZDCp > 250 - 1/20 * HF energy) continue;  
  TH2F * hNpixelTracks_vs_HF_AllCut = new TH2F("hNpixelTracks_vs_HF_AllCut","",200, 0, 8000, 80, 0, 8000);
  TH2F * hNTracks_vs_HF_AllCut = new TH2F("hNTracks_vs_HF_AllCut","",200, 0, 8000, 100, 0, 5000);
  TH1F * hJetComb_Npix_AllCut = new TH1F("hJetComb_Npix_AllCut","",120, 0, 60000);
  TH1F * hJetComb_vs_HF_AllCut = new TH1F("hJetComb_vs_HF_AllCut","",200, 0, 8000);
  TH1F * hJet80_vs_HF_AllCut = new TH1F("hJet80_vs_HF_AllCut","",200, 0, 8000);
  TH1F * hJet65_vs_HF_AllCut = new TH1F("hJet65_vs_HF_AllCut","",200, 0, 8000);
  TH1F * hJet55_vs_HF_AllCut = new TH1F("hJet55_vs_HF_AllCut","",200, 0, 8000);
  TH1F * hCentrality_fromHFBound_AllCut = new TH1F("hCentrality_fromHFBound_AllCut","",200, 0, 200);

  
  // TH2F * hNTracksPtCut_vs_HF = new TH2F("hNTracksPtCut_vs_HF","",200, 0, 8000, 100, 0, 5000);
  // TH2F * hNTracksEtaCut_vs_HF = new TH2F("hNTracksEtaCut_vs_HF","",200, 0, 8000, 100, 0, 5000);
  // TH2F * hNTracksEtaPtCut_vs_HF = new TH2F("hNTracksEtaPtCut_vs_HF","",200, 0, 8000, 100, 0, 5000);


  TH1F *hpbpb_TrgObj80[nbins_cent];
  TH1F *hpbpb_TrgObj65[nbins_cent];
  TH1F *hpbpb_TrgObj55[nbins_cent];
  TH1F *hpbpb_TrgObjComb[nbins_cent];

  TH1F *hpbpb_944Scale_TrgObj80[nbins_cent];
  TH1F *hpbpb_944Scale_TrgObj65[nbins_cent];
  TH1F *hpbpb_944Scale_TrgObj55[nbins_cent];
  TH1F *hpbpb_944Scale_TrgObjComb[nbins_cent];
  
  TH1F *hpbpb_NeqScale_TrgObj80[nbins_cent];
  TH1F *hpbpb_NeqScale_TrgObj65[nbins_cent];
  TH1F *hpbpb_NeqScale_TrgObj55[nbins_cent];
  TH1F *hpbpb_NeqScale_TrgObjComb[nbins_cent];

  TH1F *hpbpb_NeqScalePerCent_TrgObj80[nbins_cent];
  TH1F *hpbpb_NeqScalePerCent_TrgObj65[nbins_cent];
  TH1F *hpbpb_NeqScalePerCent_TrgObj55[nbins_cent];
  TH1F *hpbpb_NeqScalePerCent_TrgObjComb[nbins_cent];

  TH1F *hpbpb_raw_TrgObj80[nbins_cent];
  TH1F *hpbpb_raw_TrgObj65[nbins_cent];
  TH1F *hpbpb_raw_TrgObj55[nbins_cent];
  TH1F *hpbpb_raw_TrgObjComb[nbins_cent];

  TH1F *hpbpb_anaBin_TrgObj80[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObj65[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObj55[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObjComb[nbins_cent];

  TH1F *hpbpb_JEC_TrgObj80[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj65[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj55[nbins_cent];
  TH1F *hpbpb_JEC_TrgObjComb[nbins_cent];

  TH1F *hpbpb_Smear_TrgObj80[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj65[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj55[nbins_cent];
  TH1F *hpbpb_Smear_TrgObjComb[nbins_cent];

  TH1F * hcentrality_Jet55[nbins_cent+2];
  TH1F * hcentrality_Jet55not65not80[nbins_cent+2];
  TH1F * hcentrality_Jet65[nbins_cent+2];
  TH1F * hcentrality_Jet65not80[nbins_cent+2];
  TH1F * hcentrality_Jet80[nbins_cent+2];
  TH1F * hcentrality_Jet55_Prescl[nbins_cent+2];
  TH1F * hcentrality_Jet55not65not80_Prescl[nbins_cent+2];
  TH1F * hcentrality_Jet55or65or80[nbins_cent+2];
  TH1F * hcentrality_Jet65_Prescl[nbins_cent+2];
  TH1F * hcentrality[nbins_cent+2];



  for(int i = 0;i<nbins_cent+2;++i){

    hcentrality[i] = new TH1F(Form("hcentrality_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55[i] = new TH1F(Form("hcentrality_Jet55_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55_Prescl[i] = new TH1F(Form("hcentrality_Jet55_Prescl_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55not65not80[i] = new TH1F(Form("hcentrality_Jet55not65not80_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55not65not80_Prescl[i] = new TH1F(Form("hcentrality_Jet55not65not80_Prescl_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55or65or80[i] = new TH1F(Form("hcentrality_Jet55or65or80_cent%d",i),"",200, 0, 200);
    hcentrality_Jet65[i] = new TH1F(Form("hcentrality_Jet65_cent%d",i),"",200, 0, 200);
    hcentrality_Jet65_Prescl[i] = new TH1F(Form("hcentrality_Jet65_Prescl_cent%d",i),"",200, 0, 200);
    hcentrality_Jet65not80[i] = new TH1F(Form("hcentrality_Jet65not80_cent%d",i),"",200, 0, 200);
    hcentrality_Jet80[i] = new TH1F(Form("hcentrality_Jet80_cent%d",i),"",200, 0, 200);

  }

  
  for(int i = 0;i<nbins_cent;++i){

    hpbpb_TrgObj80[i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObj65[i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObj55[i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObjComb[i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_raw_TrgObj80[i] = new TH1F(Form("hpbpb_raw_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_raw_TrgObj65[i] = new TH1F(Form("hpbpb_raw_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_raw_TrgObj55[i] = new TH1F(Form("hpbpb_raw_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_raw_TrgObjComb[i] = new TH1F(Form("hpbpb_raw_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_anaBin_TrgObj80[i] = new TH1F(Form("hpbpb_anaBin_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObj65[i] = new TH1F(Form("hpbpb_anaBin_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObj55[i] = new TH1F(Form("hpbpb_anaBin_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObjComb[i] = new TH1F(Form("hpbpb_anaBin_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    
    hpbpb_JEC_TrgObj80[i] = new TH1F(Form("hpbpb_JEC_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JEC_TrgObj65[i] = new TH1F(Form("hpbpb_JEC_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JEC_TrgObj55[i] = new TH1F(Form("hpbpb_JEC_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JEC_TrgObjComb[i] = new TH1F(Form("hpbpb_JEC_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_Smear_TrgObj80[i] = new TH1F(Form("hpbpb_Smear_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Smear_TrgObj65[i] = new TH1F(Form("hpbpb_Smear_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Smear_TrgObj55[i] = new TH1F(Form("hpbpb_Smear_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Smear_TrgObjComb[i] = new TH1F(Form("hpbpb_Smear_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

	hpbpb_Jet80_eta_cent[i] = new TH1F(Form("hpbpb_Jet80_etacent_R%d_cent%d",radius,i),Form("Eta distribution from Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5);
	hpbpb_Jet65_eta_cent[i] = new TH1F(Form("hpbpb_Jet65_etacent_R%d_cent%d",radius,i),Form("Eta distribution from Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5);
	hpbpb_Jet55_eta_cent[i] = new TH1F(Form("hpbpb_Jet55_etacent_R%d_cent%d",radius,i),Form("Eta distribution from Jet55 && !Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5);
  
	hpbpb_Jet80_phi_cent[i] = new TH1F(Form("hpbpb_Jet80_phicent_R%d_cent%d",radius,i),Form("Phi distribution from Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-pi,pi);
	hpbpb_Jet65_phi_cent[i] = new TH1F(Form("hpbpb_Jet65_phicent_R%d_cent%d",radius,i),Form("Phi distribution from Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-pi,pi);
	hpbpb_Jet55_phi_cent[i] = new TH1F(Form("hpbpb_Jet55_phicent_R%d_cent%d",radius,i),Form("Phi distribution from Jet55 && !Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-pi,pi);

	hpbpb_Jet80_eta_precut_cent[i] = new TH1F(Form("hpbpb_Jet80_etaprecut_cent_R%d_cent%d",radius,i),Form("Eta precut distribution from Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5);
	hpbpb_Jet65_eta_precut_cent[i] = new TH1F(Form("hpbpb_Jet65_etaprecut_cent_R%d_cent%d",radius,i),Form("Eta precut distribution from Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5);
	hpbpb_Jet55_eta_precut_cent[i] = new TH1F(Form("hpbpb_Jet55_etaprecut_cent_R%d_cent%d",radius,i),Form("Eta precut distribution from Jet55 && !Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5);
	hpbpb_Jet80_phi_precut_cent[i] = new TH1F(Form("hpbpb_Jet80_phiprecut_cent_R%d_cent%d",radius,i),Form("Phi precut distribution from Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-pi,pi);
	hpbpb_Jet65_phi_precut_cent[i] = new TH1F(Form("hpbpb_Jet65_phiprecut_cent_R%d_cent%d",radius,i),Form("Phi precut distribution from Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-pi,pi);
	hpbpb_Jet55_phi_precut_cent[i] = new TH1F(Form("hpbpb_Jet55_phiprecut_cent_R%d_cent%d",radius,i),Form("Phi precut distribution from Jet55 && !Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-pi,pi);
	
	hpbpb_Jet80_phieta_cent[i] = new TH2F(Form("hpbpb_Jet80_phieta_R%d_cent%d",radius,i),Form("Phi vs Eta distribution from Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5,30,-pi,pi);
	hpbpb_Jet65_phieta_cent[i] = new TH2F(Form("hpbpb_Jet65_phieta_R%d_cent%d",radius,i),Form("Phi vs Eta distribution from Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5,30,-pi,pi);
	hpbpb_Jet55_phieta_cent[i] = new TH2F(Form("hpbpb_Jet55_phieta_R%d_cent%d",radius,i),Form("Phi vs Eta distribution from Jet55 && !Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5,30,-pi,pi);
    hpbpb_Jet80_phieta_precut_cent[i] = new TH2F(Form("hpbpb_Jet80_phieta_precut_R%d_cent%d",radius,i),Form("Phi vs Eta precut distribution from Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5,30,-pi,pi);
	hpbpb_Jet65_phieta_precut_cent[i] = new TH2F(Form("hpbpb_Jet65_phieta_precut_R%d_cent%d",radius,i),Form("Phi vs Eta precut distribution from Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5,30,-pi,pi);
	hpbpb_Jet55_phieta_precut_cent[i] = new TH2F(Form("hpbpb_Jet55_phieta_precut_R%d_cent%d",radius,i),Form("Phi vs Eta precut distribution from Jet55 && !Jet65 && !Jet80 trigger R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),30,-2.5,2.5,30,-pi,pi);	
	
  }

  // get the HLT information per forest file.
  TH1F * hJet55_File = new TH1F("hJet55_File","",902, 0, 902);
  TH1F * hJet65_File = new TH1F("hJet65_File","",902, 0, 902);
  TH1F * hJet80_File = new TH1F("hJet80_File","",902, 0, 902);
  TH1F * hJet65not80_File = new TH1F("hJet65not80_File","",902, 0, 902);
  TH1F * hJet55not6580_File = new TH1F("hJet55not6580_File","",902, 0, 902);


  // ofstream myfile1;
  // myfile1.open (kFonametxt1.c_str());

  // ofstream myfile2;
  // myfile2.open (kFonametxt2.c_str());
  
  // now start the event loop for each file. 
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpbpb[0]->GetEntries();
  if(printDebug) nentries = 10;
  TRandom rnd; 

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jetpbpb[0]->GetEntry(nEvt);
    jetpbpb[1]->GetEntry(nEvt);
    jetpbpb[2]->GetEntry(nEvt);
    jetpbpb[4]->GetEntry(nEvt);
    jetpbpb[3]->GetEntry(nEvt);
    evt_select->GetEntry(nEvt);


    if(printDebug) cout<<"forest values = "<<hiBin_F<<", "<<evt_F<<", "<<run_F<<", "<<lumi_F<<", "<<vz_F<<endl;
    
    if(pHBHENoiseFilter_F == 0) continue; 
    if(pcollisionEventSelection_F == 0) continue; 
    if(fabs(vz_F)>15) continue;
    // if(!isGoodEvent_eS) continue; 


    

    int hfbin = findHFbin(hiHF_F)-1;
    if(hfbin == -1) continue;

    int cBin = findBin(hfbin);//tells us the centrality of the event.
    if(cBin == -1) continue;

    int jetCounter = 0;
    
    for(int g = 0;g<nref_F;g++){
      
      if(eta_F[g]>=-2 && eta_F[g]<2){ //to select inside 
	
	if(pt_F[g]>=50) jetCounter++;
	
      }//eta selection cut
      
    }// jet loop

    hZDCplus_vs_HF_noCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);
    hZDCminus_vs_HF_noCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

    if(run_F > 182323){

      //if((float)hiZDC_F/1350 > (1200 - (float)(4/15 * hiHF_F))) continue; 
      hZDC_vs_HF_noCut->Fill(hiHF_F, (float)hiZDC_F/1350);
      if(hfbin == 0) hZDC_vs_HF_cent0to0p5_noCut->Fill(hiHF_F, (float)hiZDC_F/1350);

    }

    if(run_F <= 182323){

      //if((float)hiZDCminus_F/1350 > (250 - (float)(1/20 * hiHF_F))) continue; 
    
      hZDCm_vs_HF_noCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);
      if(hfbin == 0) hZDCm_vs_HF_cent0to0p5_noCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

    }




    if(run_F > 182323){

      hZDCp_vs_HF_noCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);
      if(hfbin == 0) hZDCp_vs_HF_cent0to0p5_noCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

    }


    hCentrality_fromHFBound_noCut->Fill(hfbin);
	
    if(jet55_F == 1) hJet55_vs_HF_noCut->Fill(hiHF_F); 
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) {
      hJetComb_vs_HF_noCut->Fill(hiHF_F);
      hJetComb_Npix_noCut->Fill(hiNpix_F);
    }
    
    if(jet65_F == 1) hJet65_vs_HF_noCut->Fill(hiHF_F); 
    if(jet65_F == 1 && jet80_F == 0) {
      hJetComb_vs_HF_noCut->Fill(hiHF_F);
      hJetComb_Npix_noCut->Fill(hiNpix_F);
    }
    
    if(jet80_F == 1) {
      hJet80_vs_HF_noCut->Fill(hiHF_F); 
      hJetComb_vs_HF_noCut->Fill(hiHF_F);
      hJetComb_Npix_noCut->Fill(hiNpix_F);
    }


    hNpix_vs_HF_noCut->Fill(hiHF_F, hiNpix_F);
    hNpixelTracks_vs_HF_noCut->Fill(hiHF_F, hiNpixelTracks_F);
    hNTracks_vs_HF_noCut->Fill(hiHF_F, hiNtracks_F);
    hNpix_vs_NJets_noCut->Fill(jetCounter, hiNpix_F);

    // if(hiHF_F > 4500 && (float)hiZDC_F/1350 > 50)
    //   myfile2 << run_F << " " << lumi_F << " " << evt_F<<endl;

    
    // apply the correct supernova selection cut rejection here: 
    if(hiNpix_F < (38000 - 500*jetCounter)){
      
      if(run_F > 182323){

	//if((float)hiZDC_F/1350 > (1200 - (float)(4/15 * hiHF_F))) continue; 
	hZDC_vs_HF_SupernovaCut->Fill(hiHF_F, (float)hiZDC_F/1350);
	if(hfbin == 0) hZDC_vs_HF_cent0to0p5_SupernovaCut->Fill(hiHF_F, (float)hiZDC_F/1350);

      }

      if(run_F <= 182323){

	//if((float)hiZDCminus_F/1350 > (250 - (float)(1/20 * hiHF_F))) continue; 
    
	hZDCm_vs_HF_SupernovaCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);
	if(hfbin == 0) hZDCm_vs_HF_cent0to0p5_SupernovaCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

      }

      // if(hiZDCminus_F == 0){

      // 	hZDCp_vs_HF_SupernovaCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);
      // 	if(hfbin == 0) hZDCp_vs_HF_cent0to0p5_SupernovaCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);

      // }

      hCentrality_fromHFBound_SupernovaCut->Fill(hfbin);
    
	
      if(jet55_F == 1) hJet55_vs_HF_SupernovaCut->Fill(hiHF_F); 
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) {
	hJetComb_vs_HF_SupernovaCut->Fill(hiHF_F);
	hJetComb_Npix_SupernovaCut->Fill(hiNpix_F);
      }
    
      if(jet65_F == 1) hJet65_vs_HF_SupernovaCut->Fill(hiHF_F); 
      if(jet65_F == 1 && jet80_F == 0) {
	hJetComb_vs_HF_SupernovaCut->Fill(hiHF_F);
	hJetComb_Npix_SupernovaCut->Fill(hiNpix_F);
      }
    
      if(jet80_F == 1) {
	hJet80_vs_HF_SupernovaCut->Fill(hiHF_F); 
	hJetComb_vs_HF_SupernovaCut->Fill(hiHF_F);
	hJetComb_Npix_SupernovaCut->Fill(hiNpix_F);
      }


      hNpix_vs_HF_SupernovaCut->Fill(hiHF_F, hiNpix_F);
      hNpixelTracks_vs_HF_SupernovaCut->Fill(hiHF_F, hiNpixelTracks_F);
      hNTracks_vs_HF_SupernovaCut->Fill(hiHF_F, hiNtracks_F);
      hNpix_vs_NJets_SupernovaCut->Fill(jetCounter, hiNpix_F);

    }    


    // hNpix_vs_hiBin->Fill(hiBin_F, hiNpix_F);
    // hNpixelTracks_vs_hiBin->Fill(hiBin_F, hiNpixelTracks_F);
    // hNTracks_vs_hiBin->Fill(hiBin_F, hiNtracks_F);
    // hNTracksPtCut_vs_hiBin->Fill(hiBin_F, hiNtracksPtCut_F);
    // hNTracksEtaCut_vs_hiBin->Fill(hiBin_F, hiNtracksEtaCut_F);
    // hNTracksEtaPtCut_vs_hiBin->Fill(hiBin_F, hiNtracksEtaPtCut_F);
    // hNTracksPtCut_vs_HF->Fill(hiHF_F, hiNtracksPtCut_F);
    // hNTracksEtaCut_vs_HF->Fill(hiHF_F, hiNtracksEtaCut_F);
    // hNTracksEtaPtCut_vs_HF->Fill(hiHF_F, hiNtracksEtaPtCut_F);

    // if((hiZDCplus_F != 0) && ((float)hiZDC_F/500 >= (1200 - (float)4./15 * hiHF_F))) continue;
    // if((hiZDCplus_F == 0) && ((float)hiZDCminus_F/500 >= (250 - (float)1./20 * hiHF_F))) continue;    

    // if((run_F > 182323) && ((float)hiZDC_F/1350 >= (444.44 - (float)0.09876 * hiHF_F))) continue;
    // if((run_F <= 182323) && ((float)hiZDCminus_F/1350 >= (92.58 - (float)0.02057 * hiHF_F))) continue;        

    // if((run_F > 182323) && ((float)hiZDC_F/1350 >= (400 - (float)0.091 * hiHF_F))) continue;
    // if((run_F <= 182323) && ((float)hiZDCminus_F/1350 >= (65 - (float)0.0155 * hiHF_F))) continue;        

    if((float)hiZDCminus_F/1350 >= (90 - 0.0204 * hiHF_F)) continue;
    
    if(run_F > 182323){

      //if((float)hiZDC_F/1350 > (1200 - (float)(4/15 * hiHF_F))) continue; 
      hZDC_vs_HF_PileUpCut->Fill(hiHF_F, (float)hiZDC_F/1350);
      if(hfbin == 0) hZDC_vs_HF_cent0to0p5_PileUpCut->Fill(hiHF_F, (float)hiZDC_F/1350);

    }

    if(run_F <= 182323){

      //if((float)hiZDCminus_F/1350 > (250 - (float)(1/20 * hiHF_F))) continue; 
    
      hZDCm_vs_HF_PileUpCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);
      if(hfbin == 0) hZDCm_vs_HF_cent0to0p5_PileUpCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

    }

    // if(hiZDCminus_F == 0){

    //   hZDCp_vs_HF_PileUpCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);
    //   if(hfbin == 0) hZDCp_vs_HF_cent0to0p5_PileUpCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);

    // }
    


    hCentrality_fromHFBound_PileUpCut->Fill(hfbin);
    
	
    if(jet55_F == 1) hJet55_vs_HF_PileUpCut->Fill(hiHF_F); 
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) {
      hJetComb_vs_HF_PileUpCut->Fill(hiHF_F);
      hJetComb_Npix_PileUpCut->Fill(hiNpix_F);
    }
    
    if(jet65_F == 1) hJet65_vs_HF_PileUpCut->Fill(hiHF_F); 
    if(jet65_F == 1 && jet80_F == 0) {
      hJetComb_vs_HF_PileUpCut->Fill(hiHF_F);
      hJetComb_Npix_PileUpCut->Fill(hiNpix_F);
    }
    
    if(jet80_F == 1) {
      hJet80_vs_HF_PileUpCut->Fill(hiHF_F); 
      hJetComb_vs_HF_PileUpCut->Fill(hiHF_F);
      hJetComb_Npix_PileUpCut->Fill(hiNpix_F);
    }



    hNpix_vs_HF_PileUpCut->Fill(hiHF_F, hiNpix_F);
    hNpixelTracks_vs_HF_PileUpCut->Fill(hiHF_F, hiNpixelTracks_F);
    hNTracks_vs_HF_PileUpCut->Fill(hiHF_F, hiNtracks_F);
    hNpix_vs_NJets_PileUpCut->Fill(jetCounter, hiNpix_F);

    // if(hiHF_F >= 3200 && (float)hiZDC_F/1350 <= 2.5)
    //   myfile1 << run_F << " " << lumi_F << " " << evt_F<<endl;

    // if((hiZDCplus_F !=0 && hiZDCminus_F != 0) && hiHF_F >= 3200 && (float)hiZDC_F/1350 <= 1.5) continue;        
    // if(hiZDCplus_F == 0 && hiHF_F >= 3200 && (float)hiZDCminus_F/1350 <= 1.5) continue;
    
    if(run_F > 182323){

      //if(hiHF_F >= 3200 && hiZDC_F/1350 <= 5) continue;
      
      //if((float)hiZDC_F/1350 > (1200 - (float)(4/15 * hiHF_F))) continue; 
      hZDC_vs_HF_ZDCLowCut->Fill(hiHF_F, (float)hiZDC_F/1350);
      if(hfbin == 0) hZDC_vs_HF_cent0to0p5_ZDCLowCut->Fill(hiHF_F, (float)hiZDC_F/1350);

    }

    if(run_F <= 182323){

      //if((float)hiZDCminus_F/1350 > (250 - (float)(1/20 * hiHF_F))) continue; 
      //if(hiHF_F >= 3200 && hiZDCminus_F <= 5) continue;
    
      hZDCm_vs_HF_ZDCLowCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);
      if(hfbin == 0) hZDCm_vs_HF_cent0to0p5_ZDCLowCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

    }

    // if(hiZDCminus_F == 0){
      
    //   hZDCp_vs_HF_ZDCLowCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);
    //   if(hfbin == 0) hZDCp_vs_HF_cent0to0p5_ZDCLowCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);
    //   //continue;
    // }

    hCentrality_fromHFBound_ZDCLowCut->Fill(hfbin);
    
	
    if(jet55_F == 1) hJet55_vs_HF_ZDCLowCut->Fill(hiHF_F); 
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) {
      hJetComb_vs_HF_ZDCLowCut->Fill(hiHF_F);
      hJetComb_Npix_ZDCLowCut->Fill(hiNpix_F);
    }
    
    if(jet65_F == 1) hJet65_vs_HF_ZDCLowCut->Fill(hiHF_F); 
    if(jet65_F == 1 && jet80_F == 0) {
      hJetComb_vs_HF_ZDCLowCut->Fill(hiHF_F);
      hJetComb_Npix_ZDCLowCut->Fill(hiNpix_F);
    }
    
    if(jet80_F == 1) {
      hJet80_vs_HF_ZDCLowCut->Fill(hiHF_F); 
      hJetComb_vs_HF_ZDCLowCut->Fill(hiHF_F);
      hJetComb_Npix_ZDCLowCut->Fill(hiNpix_F);
    }



    hNpix_vs_HF_ZDCLowCut->Fill(hiHF_F, hiNpix_F);
    hNpixelTracks_vs_HF_ZDCLowCut->Fill(hiHF_F, hiNpixelTracks_F);
    hNTracks_vs_HF_ZDCLowCut->Fill(hiHF_F, hiNtracks_F);
    hNpix_vs_NJets_ZDCLowCut->Fill(jetCounter, hiNpix_F);


    
    // apply the correct supernova selection cut rejection here: 
    if(hiNpix_F >= (38000 - 500*jetCounter)){
      if(printDebug) cout<<"removed this supernova event"<<endl;
      continue;
    }



    if(run_F > 182323){

      //if((float)hiZDC_F/1350 > (1200 - (float)(4/15 * hiHF_F))) continue; 
      hZDC_vs_HF_AllCut->Fill(hiHF_F, (float)hiZDC_F/1350);
      if(hfbin == 0) hZDC_vs_HF_cent0to0p5_AllCut->Fill(hiHF_F, (float)hiZDC_F/1350);      

    }

    if(run_F <= 182323){

      //if((float)hiZDCminus_F/1350 > (250 - (float)(1/20 * hiHF_F))) continue; 
    
      hZDCm_vs_HF_AllCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);
      if(hfbin == 0) hZDCm_vs_HF_cent0to0p5_AllCut->Fill(hiHF_F, (float)hiZDCminus_F/1350);

    }

    // if(hiZDCminus_F == 0){

    //   hZDCp_vs_HF_AllCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);
    //   if(hfbin == 0) hZDCp_vs_HF_cent0to0p5_AllCut->Fill(hiHF_F, (float)hiZDCplus_F/1350);

    // }

    hCentrality_fromHFBound_AllCut->Fill(hfbin);
    
	
    if(jet55_F == 1) hJet55_vs_HF_AllCut->Fill(hiHF_F); 
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) {
      hJetComb_vs_HF_AllCut->Fill(hiHF_F);
      hJetComb_Npix_AllCut->Fill(hiNpix_F);
    }
    
    if(jet65_F == 1) hJet65_vs_HF_AllCut->Fill(hiHF_F); 
    if(jet65_F == 1 && jet80_F == 0) {
      hJetComb_vs_HF_AllCut->Fill(hiHF_F);
      hJetComb_Npix_AllCut->Fill(hiNpix_F);
    }
    
    if(jet80_F == 1) {
      hJet80_vs_HF_AllCut->Fill(hiHF_F); 
      hJetComb_vs_HF_AllCut->Fill(hiHF_F);
      hJetComb_Npix_AllCut->Fill(hiNpix_F);
    }

    hNpix_vs_HF_AllCut->Fill(hiHF_F, hiNpix_F);
    hNpixelTracks_vs_HF_AllCut->Fill(hiHF_F, hiNpixelTracks_F);
    hNTracks_vs_HF_AllCut->Fill(hiHF_F, hiNtracks_F);
    hNpix_vs_NJets_AllCut->Fill(jetCounter, hiNpix_F);

    hcentrality[cBin]->Fill(hfbin);
    hcentrality[nbins_cent+1]->Fill(hfbin);
    if(jet55_F == 1) hcentrality_Jet55[cBin]->Fill(hfbin);
    if(jet55_F == 1) hcentrality_Jet55[nbins_cent+1]->Fill(hfbin);
    if(jet55_F == 1) hcentrality_Jet55_Prescl[cBin]->Fill(hfbin, jet55_p_F);
    if(jet55_F == 1) hcentrality_Jet55_Prescl[nbins_cent+1]->Fill(hfbin, jet55_p_F);
    if(jet65_F == 1) hcentrality_Jet65[cBin]->Fill(hfbin);
    if(jet65_F == 1) hcentrality_Jet65[nbins_cent+1]->Fill(hfbin);
    if(jet65_F == 1) hcentrality_Jet65_Prescl[cBin]->Fill(hfbin, jet65_p_F);
    if(jet65_F == 1) hcentrality_Jet65_Prescl[nbins_cent+1]->Fill(hfbin, jet65_p_F);
    if(jet80_F == 1) hcentrality_Jet80[cBin]->Fill(hfbin);
    if(jet80_F == 1) hcentrality_Jet80[nbins_cent+1]->Fill(hfbin);
    if(jet55_F == 1 || jet65_F == 1 || jet80_F == 1) hcentrality_Jet55or65or80[cBin]->Fill(hfbin);
    if(jet55_F == 1 || jet65_F == 1 || jet80_F == 1) hcentrality_Jet55or65or80[nbins_cent+1]->Fill(hfbin);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80[cBin]->Fill(hfbin);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80[nbins_cent+1]->Fill(hfbin);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80_Prescl[cBin]->Fill(hfbin, jet55_p_F);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80_Prescl[nbins_cent+1]->Fill(hfbin, jet55_p_F);
    if(jet65_F == 1 && jet80_F == 0) hcentrality_Jet65not80[cBin]->Fill(hfbin);
    if(jet65_F == 1 && jet80_F == 0) hcentrality_Jet65not80[nbins_cent+1]->Fill(hfbin);

    if(jet55_F == 1) hJet55_File->Fill(startfile);
    if(jet65_F == 1) hJet65_File->Fill(startfile);
    if(jet80_F == 1) hJet80_File->Fill(startfile);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hJet55not6580_File->Fill(startfile);
    if(jet65_F == 1 && jet80_F == 0) hJet65not80_File->Fill(startfile);

    // jet_select->GetEntry(nGoodEvt);
    // ++nGoodEvt;

    if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
    //if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
    //if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;

    if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
    if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
    if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
    if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
    if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<", nref_jS = "<<nref_jS<<endl;

    if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;

    if(cBin == nbins_cent) continue;

    // //! Sort the jetTree jets according to pT
    // std::vector < Jet > vJet;
    // for(int jet2 = 0; jet2<nref_jS; ++jet2){
    //   //cout <<"\t \t jetTree *** "<< jet2 <<  ", pT " << pt_jS[jet2] <<  ", chSum : "<< chSum_jS[jet2] << endl;
    //   Jet ijet;
    //   ijet.id = jet2;
    //   ijet.pt = pt_jS[jet2];
    //   vJet.push_back(ijet);
    // }
    // std::sort (vJet.begin(), vJet.end(), compare_pt);
    // std::vector < Jet >::const_iterator itr;

    // int jet=0;
    // for(itr=vJet.begin(); itr!=vJet.end(); ++itr, ++jet){

    //   int jetLoc = (*itr).id;
    //   if(isMultiMatch_jS[jetLoc]) {
    // 	++itr;
    // 	jetLoc = (*itr).id;
    // 	if(itr == vJet.end())  break;
    //   }
    //   if(fabs(eta_jS[jetLoc]) > 2) continue;
    //   //if(isPFElecCut_eS[jet] != 1) continue;
    //   // if(isMiMatch_eS[jet]) continue;
    //   if(pt_jS[jetLoc] <15) continue;

    //   bool PFElecCut = false;

    //   Float_t Sumcand = chSum_jS[jetLoc] + phSum_jS[jetLoc] + neSum_jS[jetLoc] + muSum_jS[jetLoc];
    //   if(isCaloMatch_jS[jetLoc] == 1){
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.5 && calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.85 && eMax_jS[jetLoc]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_jS[jetLoc]/pt_jS[jetLoc] - (Float_t)9/7)) PFElecCut = true;
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.85) PFElecCut = true;
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.5 && eMax_jS[jetLoc]/Sumcand < 0.05) PFElecCut = true;
    //   }
    //   if(isCaloMatch_jS[jetLoc] == 0)
    // 	if(eMax_jS[jetLoc]/Sumcand < 0.05 ) PFElecCut = true;

    //   // if(!PFElecCut) continue;
      
    //   // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
    //   // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;

    Float_t wght = 1;       
    //   if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT = "<<calopt_F[index_eS[jet]]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;

	//precut loop
	
	for( int jet = 0; jet<nref_F; jet++ ){

	  
      if( fabs(eta_F[jet]) > 2.0 ) continue;
	  
      float recpt = pt_F[jet];
	  float eta = eta_F[jet];
	  float phi = phi_F[jet];
      
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0){
	//if(recpt > 140) continue;
	hpbpb_Jet55_eta_precut_cent[cBin]->Fill(eta, jet55_p_F* wght);
	hpbpb_Jet55_phi_precut_cent[cBin]->Fill(phi, jet55_p_F* wght);
	hpbpb_Jet55_phieta_precut_cent[cBin]->Fill(eta, phi, jet55_p_F* wght);
      }

      if(jet65_F == 1 && jet80_F == 0){
	//if(recpt > 140) continue;
	hpbpb_Jet65_eta_precut_cent[cBin]->Fill(eta, wght);
	hpbpb_Jet65_phi_precut_cent[cBin]->Fill(phi, wght);
	hpbpb_Jet65_phieta_precut_cent[cBin]->Fill(eta, phi, wght);
      }
      
      if(jet80_F == 1){
	hpbpb_Jet80_eta_precut_cent[cBin]->Fill(eta, wght);
	hpbpb_Jet80_phi_precut_cent[cBin]->Fill(phi, wght);
	hpbpb_Jet80_phieta_precut_cent[cBin]->Fill(eta, phi, wght);
      }
      
    }// jet loop
	
	
    for( int jet = 0; jet<nref_F; jet++ ){

	  
      if( fabs(eta_F[jet]) > 2.0 ) continue;
	  
      //if( isPFElecCut_eS[jet] != 1 ) continue;

      //if( chMax_F[jet] < 7 && trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;
      // if( trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;

      bool PFElecCut = false;

      Float_t Sumcand = chSum_F[jet] + phSum_F[jet] + neSum_F[jet] + muSum_F[jet];
      if(isClMatch_eS[jet] == 1){
    	if(calopt_eS[jet]/pt_F[jet] > 0.5 && calopt_eS[jet]/pt_F[jet] <= 0.85 && eMax_F[jet]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_eS[jet]/pt_F[jet] - (Float_t)9/7)) PFElecCut = true;
    	if(calopt_eS[jet]/pt_F[jet] > 0.85) PFElecCut = true;
    	if(calopt_eS[jet]/pt_F[jet] <= 0.5 && eMax_F[jet]/Sumcand < 0.05) PFElecCut = true;
      }
      if(isClMatch_eS[jet] == 0)
    	if(eMax_F[jet]/Sumcand < 0.05 ) PFElecCut = true;

      if(isPFElecCut_eS[jet] != PFElecCut && printDebug) 
	cout<<" pf e cut not same, run = "<<run_F<<" "<<run_eS<<", event = "<<evt_F<<" "<<evt_eS<<" , lumi = "<<lumi_F<<" "<<lumi_eS<<endl;
	  
      if(!PFElecCut) continue;

      // also need to cut on the high pT jets from the Jet 55 sample. unmatched PF jets with pfpt > 110 and calopt < 30 including -999, check on their high track content. at the moment dont worry about this.  
    
      float recpt = pt_F[jet];
	  float eta = eta_F[jet];
	  float phi = phi_F[jet];
      
	  if(recpt<=pfptmin)  continue; //pT Cut
	  
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0){
	//if(recpt > 140) continue;
	
	hpbpb_TrgObj55[cBin]->Fill(recpt, jet55_p_F* wght);
	hpbpb_raw_TrgObj55[cBin]->Fill(rawpt_F[jet], jet55_p_F* wght);
	hpbpb_anaBin_TrgObj55[cBin]->Fill(recpt, jet55_p_F* wght);
	hpbpb_JEC_TrgObj55[cBin]->Fill(recpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), jet55_p_F* wght);
	hpbpb_Smear_TrgObj55[cBin]->Fill(recpt + rnd.Gaus(0,1), jet55_p_F* wght);
	hpbpb_Jet55_eta_cent[cBin]->Fill(eta, jet55_p_F* wght);
	hpbpb_Jet55_phi_cent[cBin]->Fill(phi, jet55_p_F* wght);
	hpbpb_Jet55_phieta_cent[cBin]->Fill(eta, phi, jet55_p_F* wght);
      }

      if(jet65_F == 1 && jet80_F == 0){
	//if(recpt > 140) continue;
	hpbpb_TrgObj65[cBin]->Fill(recpt, wght);
	hpbpb_raw_TrgObj65[cBin]->Fill(rawpt_F[jet], wght);
	hpbpb_anaBin_TrgObj65[cBin]->Fill(recpt, wght);
	hpbpb_JEC_TrgObj65[cBin]->Fill(recpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj65[cBin]->Fill(recpt + rnd.Gaus(0,1), wght);
	hpbpb_Jet65_eta_cent[cBin]->Fill(eta, wght);
	hpbpb_Jet65_phi_cent[cBin]->Fill(phi, wght);
	hpbpb_Jet65_phieta_cent[cBin]->Fill(eta, phi, wght);
      }
      
      if(jet80_F == 1){
	hpbpb_TrgObj80[cBin]->Fill(recpt, wght);
	hpbpb_raw_TrgObj80[cBin]->Fill(rawpt_F[jet], wght);
	hpbpb_anaBin_TrgObj80[cBin]->Fill(recpt, wght);
	hpbpb_JEC_TrgObj80[cBin]->Fill(recpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj80[cBin]->Fill(recpt + rnd.Gaus(0,1), wght);
	hpbpb_Jet80_eta_cent[cBin]->Fill(eta, wght);
	hpbpb_Jet80_phi_cent[cBin]->Fill(phi, wght);
	hpbpb_Jet80_phieta_cent[cBin]->Fill(eta, phi, wght);
      }
      
    }// jet loop
    if(printDebug)cout<<endl;


  }// event loop

  
  for(int i = 0; i<nbins_cent; ++i){


    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj80[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj65[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj55[i]);
    
    divideBinWidth(hpbpb_TrgObjComb[i]);
    divideBinWidth(hpbpb_TrgObj80[i]);
    divideBinWidth(hpbpb_TrgObj65[i]);
    divideBinWidth(hpbpb_TrgObj55[i]);

    hpbpb_944Scale_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Clone(Form("hpbpb_HLT80_944Scale_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_944Scale_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Clone(Form("hpbpb_HLT65_944Scale_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_944Scale_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Clone(Form("hpbpb_HLT55_944Scale_R%d_20_eta_20_cent%d", radius, i));

    hpbpb_944Scale_TrgObj80[i]->Scale(1./ScaleFactor_944Scale[i]);
    hpbpb_944Scale_TrgObj65[i]->Scale(1./ScaleFactor_944Scale[i]);
    hpbpb_944Scale_TrgObj55[i]->Scale(1./ScaleFactor_944Scale[i]);
    hpbpb_944Scale_TrgObjComb[i] = (TH1F*)hpbpb_944Scale_TrgObj80[i]->Clone(Form("hpbpb_HLTComb_944Scale_R%d_20_eta_20_cent%d",radius, i));
    hpbpb_944Scale_TrgObjComb[i]->Add(hpbpb_944Scale_TrgObj65[i]);
    hpbpb_944Scale_TrgObjComb[i]->Add(hpbpb_944Scale_TrgObj55[i]);
    
    hpbpb_NeqScalePerCent_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Clone(Form("hpbpb_HLT80_NeqScalePerCent_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_NeqScalePerCent_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Clone(Form("hpbpb_HLT65_NeqScalePerCent_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_NeqScalePerCent_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Clone(Form("hpbpb_HLT55_NeqScalePerCent_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_NeqScalePerCent_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScalePerCent[i]);
    hpbpb_NeqScalePerCent_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScalePerCent[i]);
    hpbpb_NeqScalePerCent_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScalePerCent[i]);
    hpbpb_NeqScalePerCent_TrgObjComb[i] = (TH1F*)hpbpb_NeqScalePerCent_TrgObj80[i]->Clone(Form("hpbpb_HLTComb_NeqScalePerCent_R%d_20_eta_20_cent%d",radius, i));
    hpbpb_NeqScalePerCent_TrgObjComb[i]->Add(hpbpb_NeqScalePerCent_TrgObj65[i]);
    hpbpb_NeqScalePerCent_TrgObjComb[i]->Add(hpbpb_NeqScalePerCent_TrgObj55[i]);

    
    hpbpb_NeqScale_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Clone(Form("hpbpb_HLT80_NeqScale_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_NeqScale_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Clone(Form("hpbpb_HLT65_NeqScale_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_NeqScale_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Clone(Form("hpbpb_HLT55_NeqScale_R%d_20_eta_20_cent%d", radius, i));
    hpbpb_NeqScale_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_NeqScale_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_NeqScale_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    hpbpb_NeqScale_TrgObjComb[i] = (TH1F*)hpbpb_NeqScale_TrgObj80[i]->Clone(Form("hpbpb_HLTComb_NeqScale_R%d_20_eta_20_cent%d",radius, i));
    hpbpb_NeqScale_TrgObjComb[i]->Add(hpbpb_NeqScale_TrgObj65[i]);
    hpbpb_NeqScale_TrgObjComb[i]->Add(hpbpb_NeqScale_TrgObj55[i]);

    
    hpbpb_raw_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_raw_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_raw_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    
    hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj80[i]);
    hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj65[i]);
    hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj55[i]);

    divideBinWidth(hpbpb_raw_TrgObjComb[i]);
    divideBinWidth(hpbpb_raw_TrgObj80[i]);
    divideBinWidth(hpbpb_raw_TrgObj65[i]);
    divideBinWidth(hpbpb_raw_TrgObj55[i]);

    hpbpb_anaBin_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_anaBin_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_anaBin_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj80[i]);
    hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj65[i]);
    hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj55[i]);

    divideBinWidth(hpbpb_anaBin_TrgObjComb[i]);
    divideBinWidth(hpbpb_anaBin_TrgObj80[i]);
    divideBinWidth(hpbpb_anaBin_TrgObj65[i]);
    divideBinWidth(hpbpb_anaBin_TrgObj55[i]);

    hpbpb_Smear_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_Smear_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_Smear_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj80[i]);
    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj65[i]);
    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj55[i]);

    divideBinWidth(hpbpb_Smear_TrgObjComb[i]);

    hpbpb_JEC_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_JEC_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_JEC_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj80[i]);
    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj65[i]);
    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj55[i]);

    divideBinWidth(hpbpb_JEC_TrgObjComb[i]);
    
  }

  
  fout->Write();

  // myfile1.close();
  // myfile2.close();
  
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
