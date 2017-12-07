///////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
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


// I/O ------------------------
// env ------------------------

//CMSSW_BASE
const std::string CMSSW_BASE= 
  "/home/obaron/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string SCRATCH_BASE= 
  "/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";

const std::string input_dir=SCRATCH_BASE;//+"4.07.16_outputCondor/readForests_HIJetRECO_4-07-17/";
//const std::string readForests_condor_dir=CMSSW_BASE+"src/readForests/outputCondor";

const std::string ppData_inCondorDir=input_dir+"10.11.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_10-11-17_jetPlots_0.0eta2.0_rec56_HLTCalo/";
const std::string input_ppData_Filename=ppData_inCondorDir+"HighPtJetTrig_ak4PF-allFiles.root";
//const std::string input_ppData_Filename="HighPtJetTrig_ak4PF-allFiles.root";

const std::string ppMC_inCondorDir  =input_dir+"10.11.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-11-17_JERS_0.0eta2.0_rec56_gen32_HLTCalo/";
const std::string input_ppMC_Filename=ppMC_inCondorDir+"Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
//const std::string input_ppMC_Filename="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";

int main (int argc, char *argv[]){
  
  gStyle->SetOptFit(1111);
  
  std::cout<<" now opening Data File "<<std::endl<<input_ppData_Filename<<std::endl<<std::endl;
  TFile* finData = new TFile(input_ppData_Filename.c_str());
  
  TH1F *h_NEvents         = (TH1F*)finData->Get("NEvents_read");
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
  
  const Int_t NEvts=h_NEvents->GetEntries();
  const Int_t NEvtsPassvzCut=h_NEvents_vzCut->GetEntries();
  std::cout<<std::endl;
  std::cout<<"NEvents in dataset="<<NEvts<<std::endl;
  std::cout<<"NEvents passing vz/skim cut="<<NEvtsPassvzCut<<std::endl;
  std::cout<<std::endl;
  
  const double LumiEff_vz = (double) NEvtsPassvzCut/NEvts;//h_NEvents_vzCut->GetEntries()/h_NEvents->GetEntries(); 
  std::cout<<std::endl;
  std::cout<<"Lumi Efficiency post vz+skim cuts=NEvents_vzCut/NEvents ";
  std::cout<<"=LumiEff_vz ="<<LumiEff_vz<<std::endl;
  std::cout<<std::endl;

  const double intgrtdLumi=27.8*pow(10.,6.);// 25.8 pb^-1 to \microbarns^-1
  const double effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz; 
  std::cout<<std::endl;
  std::cout<<"effective Intgrted Lumi ="<<effIntgrtdLumi_vz<<std::endl;
  std::cout<<std::endl;

  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( "hWeightedVz" );
  TH1F *theRatio = (TH1F*)theDataEvtQAHist->Rebin(2,"theRatio");
  theRatio->Scale( 1/theRatio->GetBinWidth(1) );
  theRatio->Scale( 1/effIntgrtdLumi_vz );
  
  //theDataEvtQAHist->SetTitle (    h_Title.c_str() );
  //theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  //theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  //
  //theDataEvtQAHist->SetMarkerStyle(kDot);
  //theDataEvtQAHist->SetMarkerSize(1.1);
  //theDataEvtQAHist->SetMarkerColor( kBlack);
  //theDataEvtQAHist->SetLineColor( theRatioLineColor );
  //theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");
  TF1* fgaussData = NULL;
  TF1* fgaussMC = NULL;
  
  TCanvas* TCanMC = NULL;
  TCanvas* TCanWeight = NULL;
  
  std::cout<<" now opening MC File "<<std::endl<<input_ppMC_Filename<<std::endl<<std::endl;
  TFile* finMC = new TFile(input_ppMC_Filename.c_str());  
  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( "hpthatWeightedVz" );
  theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(1) );
  theMCEvtQAHist->Scale( theRatio->Integral()/theMCEvtQAHist->Integral() );
  
  TCanMC = new TCanvas("TCanMC","cMC",600,600);   
  TCanWeight = new TCanvas("TCanWeight","Weight: Data/MC",600,600);   
    
  //TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("theDataHistClone"); I'm going to replace this with the rebin which should make its own clone
  double norm = theRatio->GetMaximumStored();
  
  TH1F* binWeight = new TH1F("binWeight","Bin-based Weight",500,-25,25);
  TH1F* fnWeight = new TH1F("fnWeight","Fn-based Weight",500,-25,25);
  
  //make fit functions
  fgaussData = new TF1("fgaussData","gaus", -24, 24);
  fgaussData->SetParameters(norm, 0.9999, 0.15); 

  fgaussMC = new TF1("fgaussMC","gaus", -24, 24);
  fgaussMC->SetParameters(norm, 0.9999, 0.15); 
  
  int fitstatusdata = 0; 
  fitstatusdata = theRatio->Fit(fgaussData, "R");
  std::cout<< "Data Fit Status: "<< fitstatusdata<< ", Fit Error: "<< fgaussData->GetParError(1)<< std::endl;
  
  int fitstatusMC = 0; 
  fitstatusMC = theMCEvtQAHist->Fit(fgaussMC, "R");
  std::cout<< "MC Fit Status: "<< fitstatusMC<< ", Fit Error: "<< fgaussMC->GetParError(1)<< std::endl;
  
  //Compare fit to histogram
  TCanMC->cd();
  theMCEvtQAHist->Draw();
  fgaussMC->Draw("same");
  TCanMC->Print("MCgaussfit.png","png");
  
  theRatio->Divide(theMCEvtQAHist);
  //theRatio->Draw();  
  //theRatio->SetLineColor( altRatioLineColor1 );
  //theEvtQALeg->AddEntry(theRatio,"MC not vz-weighted","lp");
  
  Float_t theVzBinWidth=theRatio->TH1::GetBinWidth(1);
  Float_t xLow = theRatio->TH1::GetBinLowEdge(1), xHigh=theRatio->TH1::GetBinLowEdge(500);
  Float_t NvzWeightBins_F=(xHigh-xLow)/(theVzBinWidth);
  Int_t NvzWeightBins = theRatio->TH1::GetNbinsX();//should be 60
  
  std::cout<<"float-division says NBins="<<NvzWeightBins_F<<" (exactly)"<<std::endl;
  std::cout<<"int-typecast says NBins="<<NvzWeightBins<<std::endl;
  
  std::cout<<"now grabbing vzWeights for "<<NvzWeightBins<<" bins for ( "<<xLow<<"< vz <"<<xHigh<<" )"<<std::endl;
  std::cout<<std::endl;
  
  
  for (int i=1;i<=NvzWeightBins;++i){//binsX loop

	Float_t hist_xLow = theRatio->TH1::GetBinLowEdge(i);
	std::cout<<"Low Bin Edge = "<<hist_xLow<<std::endl;
    Float_t vzWeight = theRatio->TH1::GetBinContent(i);    //TH1 bin counting starts at i=1?! why?!
	binWeight->SetBinContent(i,vzWeight);
	//function fit ratio weights //No no no - this needs to be the x value, not the bin content! //so do I want center, low edge, or high edge? Or something else?
	Double_t gaussMC = fgaussMC->Eval(theMCEvtQAHist->GetBinLowEdge(i)+(theMCEvtQAHist->GetBinWidth(i))/2);
	Double_t gaussData = fgaussData->Eval(theMCEvtQAHist->GetBinLowEdge(i)+(theMCEvtQAHist->GetBinWidth(i))/2);
	Double_t gaussFit = (gaussData/gaussMC);
	//Double_t test = 0.4;
	fnWeight->SetBinContent(i,gaussFit);
		
    if(theRatio->GetBinContent(i)<=0.)
      std::cout<<"warning! bin content in data hist zero (numerator)"<<std::endl;
	  std::cout<<"bin content = "<<theRatio->GetBinContent(i)<<std::endl;
    if(i%5==0) {
      
      std::cout<<"i=="<<i<<", vzWeight="<<vzWeight <<" , vzLow="<<xLow<<std::endl;
	  std::cout<<"MC function LowEdge = "<<theMCEvtQAHist->GetBinLowEdge(i)<<std::endl;
	  std::cout<<"Function fit weight="<<gaussFit<<std::endl;
	}
    else{
      std::cout<<"i=="<<i<<", vzWeight="<<vzWeight<< std::endl;
	  std::cout<<"Function fit weight="<<gaussFit<<std::endl;
	}
    std::cout<<std::endl;
    //Float_t leftSideOfBin=xLow+(i)*theVzBinWidth;
    //Float_t rightSideOfBin=xLow+(i+1)*theVzBinWidth;
    //std::cout<<"i=="<<i<<", "<<leftSideOfBin<<"<vz<="<<rightSideOfBin<<", vzWeight="<<vzWeight<< std::endl;
    xLow=xLow+theRatio->TH1::GetBinWidth(1);
    
    
  }
  //vz from function loop - testing. I think this needs its own loop.
  	  

  TH1F* fweightRatio = (TH1F*)binWeight->Clone("fweightRatio");
  fweightRatio->Divide(fnWeight);
  TCanWeight->cd();
  fweightRatio->Draw();
  //fnWeight->Draw();
  TCanWeight->Print("WeightRatio.png","png");
/*
//Might not be a good approach:

  //make ratio of fit function
  fgaussratio = new TF1("fgaussratio",fgaussData/fgaussMC,-2,2);
*/
  
  
  return 0 ;
}
