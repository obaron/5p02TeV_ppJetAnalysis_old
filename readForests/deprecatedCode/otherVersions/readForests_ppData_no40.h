// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests_ppData_no40.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppData_no40.exe


////////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TCut.h>
// utilities
#include <TMath.h>
#include <TRandom3.h>
#include <TStopwatch.h>
// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>


//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readForests_ppData_no40
const int defStartFile=50;
const int defEndFile=60; //inclusive boundary, runs 2 files by default
//const std::string defInFilelist = "filelists/5p02TeV_HighPtJet80_forests.txt";
const std::string defInFilelist = "filelists/5p02TeV_HighPtLowerJets_forests.txt";
const int defRadius=4;
const std::string defJetType="PF";
const std::string defOutputName = "readForests_ppData_no40_defOut.root";
const bool defDebugMode = true;
int readForests_ppData_no40(int startfile = defStartFile , int endfile = defEndFile ,
		       std::string inFilelist = defInFilelist , std::string outfile = defOutputName ,
		       int radius = defRadius , std::string jetType = defJetType , bool debugMode = defDebugMode );
const int readForestsArgCount=7+minArgs;

//// helper functions
float trigComb(bool *trg, int *pscl, float pt);
void divideBinWidth(TH1 *h);
float deltaphi(float phi1, float phi2);


//// CONSTANTS
// ---------------------------------------------------------------------------------------------------------------

// misc 
const char *etaWidth = (char*)"20_eta_20";

// binning arrays
const float ptbins[] = { 15., 30., 50., 80., 120., 170., 220., 300., 500. };
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one

const float JEC_ptbins[] = {
  17., 22., 27.,    //15-30
  33., 39., 47.,    //30-50
  55., 64., 74.,    //50-80
  84., 97., 114.,   //80-120
  133., 153.,      //120-170
  174., 196.,      //170-220
  220., 245., 272., //220-300
  300., 350., 400., //300-500
  550., 790., 1000. //500-inf
};
const int nbins_JEC_ptbins = sizeof(JEC_ptbins)/sizeof(float)-1;

const float etabins[] = {
  -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, 
  -3.839, -3.664, -3.489, -3.314, -3.139, 
  -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043, 
  -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, 
  -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, 
  -0.087, +0.000, +0.087, 
  +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, +0.879, +0.957, 
  +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, +1.653, +1.740, +1.830, +1.930, 
  +2.043, +2.172, +2.322, +2.500, +2.650, +2.853, +2.964, 
  +3.139, +3.314, +3.489, +3.664, +3.839, 
  +4.013, +4.191, +4.363, +4.538, +4.716, +4.889, +5.191
};
const int nbins_eta = sizeof(etabins)/sizeof(float)-1;

// string arrays

// variable names for QA Plots
const std::string var[] = {   
  "jtpt" ,  "rawpt",  //jets
  "jteta", "jtphi",   
  "trkMax", "trkSum", "trkHardSum", //jet constituents
  "chMax",  "chSum",  "chHardSum", 
  "phMax",  "phSum",  "phHardSum", 
  "neMax",  "neSum", 
  "eMax",   "eSum", 
  "muMax",  "muSum", 
  "Aj",     "xj" , "dphi", //dijet variables
  "leadJetPt", "subleadJetPt"
};
const int N_vars = sizeof(var)/sizeof(std::string);

//L1
const std::string L1BitStrings[]={//this array is a good idea
  //  "L1_SingleJet28_BptxAND",
  "L1_SingleJet40_BptxAND",
  "L1_SingleJet48_BptxAND",   
  "L1_SingleJet52_BptxAND"
};
const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);

const std::string L1PresclBranches[]={//this array is a good idea
  //  L1BitStrings[0]+"_Prescl",
  L1BitStrings[0]+"_Prescl",
  L1BitStrings[1]+"_Prescl",
  L1BitStrings[2]+"_Prescl"
};

// HLT                           
const std::string HLTJetType="Calo";//"PF";
const std::string HLTBitStrings[]={	
  //  "HLT_AK4"+HLTJetType+"Jet40_Eta5p1",
  "HLT_AK4"+HLTJetType+"Jet60_Eta5p1",			
  "HLT_AK4"+HLTJetType+"Jet80_Eta5p1",			
  "HLT_AK4"+HLTJetType+"Jet100_Eta5p1"    		
};						
//const int N_HLTBits=sizeof(HLTBitStrings)/sizeof(std::string);
const int N_HLTBits=N_L1Bits;

const std::string HLTPresclBranches[]={	
  //  HLTBitStrings[0]+"_v1_Prescl",
  HLTBitStrings[0]+"_v1_Prescl",
  HLTBitStrings[1]+"_v1_Prescl",
  HLTBitStrings[2]+"_v1_Prescl"
};			
			
const std::string HLTBranches[]={	
  //  HLTBitStrings[0]+"_v1",
  HLTBitStrings[0]+"_v1",
  HLTBitStrings[1]+"_v1",
  HLTBitStrings[2]+"_v1"
};						

// tree names+directories
const std::string treeNames[]={ 
  "GARBAGE ENTRY" , //use jet ana of choice later
  "hiEvtAnalyzer/HiTree" ,
  "skimanalysis/HltTree" ,
  "hltanalysis/HltTree" ,
  //  "hltobject/"+HLTBitStrings[0]+"_v" , 
  "hltobject/"+HLTBitStrings[0]+"_v" , 
  "hltobject/"+HLTBitStrings[1]+"_v" , 
  "hltobject/"+HLTBitStrings[2]+"_v"   
}; 
const int N_trees = sizeof(treeNames)/sizeof(std::string);


//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------


void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop 
    Float_t val = h->GetBinContent(i); 
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);     
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);    
    h->SetBinError(i,valErr);
  }//end nbinsX loop
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return;
}

float trigComb(bool *trgDec, int *treePrescl, double triggerPt){
  float weight_eS=0.;
  //  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. ) weight_eS = treePrescl[0];
  if(trgDec[0] && triggerPt>=60.  && triggerPt<80. ) weight_eS = treePrescl[0];
  if(trgDec[1] && triggerPt>=80.  && triggerPt<100.) weight_eS = treePrescl[1];
  if(trgDec[2] && triggerPt>=100. )                  weight_eS = treePrescl[2];
  return weight_eS;
}

float trigComb_sanityCheck(bool *trgDec, int *treePrescl, double triggerPt){
  std::cout<<"trigComb_sanityCheck Called..."<<std::endl<<std::endl;
  
  float weight_eS=0.;
  std::cout<<"before if statements..."<<std::endl;
  for(int k=0;k<N_L1Bits;k++){if(trgDec[k]) {
      std::cout<<"trgDec["<<k<<"]="<<trgDec[k]<<std::endl;
    std::cout<<"treePrescl["<<k<<"]="<<treePrescl[k]<<std::endl;
    std::cout<<"triggerPt="<<triggerPt<<std::endl;
    std::cout<<"weight_eS=="<<weight_eS<<std::endl<<std::endl;
    }}
  
  //  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. ) weight_eS = treePrescl[0];
  if(trgDec[0] && triggerPt>=60.  && triggerPt<80. ) weight_eS = treePrescl[0];
  if(trgDec[1] && triggerPt>=80.  && triggerPt<100.) weight_eS = treePrescl[1];
  if(trgDec[2] && triggerPt>=100. )                  weight_eS = treePrescl[2];

  std::cout<<"after if statements..."<<std::endl;
  for(int k=0;k<N_L1Bits;k++){if(trgDec[k]){
      std::cout<<"trgDec["<<k<<"]="<<trgDec[k]<<std::endl;
    std::cout<<"treePrescl["<<k<<"]="<<treePrescl[k]<<std::endl;
    std::cout<<"triggerPt="<<triggerPt<<std::endl;
    std::cout<<"weight_eS=="<<weight_eS<<std::endl<<std::endl;
    }}
  
  std::cout<<"trigComb_sanityCheck exiting..."<<std::endl<<std::endl;
  return weight_eS;
}

float deltaphi(float phi1, float phi2){
  float pi=TMath::Pi(); 
  float dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -= 2.*pi;
  return TMath::Abs(dphi);
}
