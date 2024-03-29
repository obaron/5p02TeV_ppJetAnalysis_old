#include "printPlots.h"


// code/job switches ------------------------

//options
const bool debugMode=true, doEventCounts=true, doJetIDPlots=true;

//draw switches
const bool drawEvtQAPlots=true;
const bool drawJetQAPlots=true;
const bool drawJetConstituentPlots=drawJetQAPlots&&true, drawDijetPlots=drawJetQAPlots&&false;
const bool drawJetRapBinsPlot=false;//, drawGENJetRapBinsPlot=true;  


//// hist painting ------------------------
//
//// line colors
//const int theDataOverlayLineColor=1, theMCOverlayLineColor=1;//, altOverlayLineColor=3; 
//const int theRatioLineColor=1;//,altRatioLineColor1=8, altRatioLineColor2=7;
//
//// marker colors
//const int theDataOverlayMarkerColor=2, theMCOverlayMarkerColor=4;//,theRatioMarkerColor=9;












// the macro ------------------------
int printPlots_jetPlots(const std::string input_ppData_condorDir , const std::string input_ppMC_condorDir , 
			const std::string output_PDFname_base){
  

  globalHistStyle();
  
  ////format of the filename is always HighPt{filelist}_ak{3,4}{PF,Calo}-allFiles.root
  //std::size_t jobTypePos=input_ppMC_condorDir.find("__")+2;  
  //std::string jobType;
  //if(jobTypePos!=std::string::npos) jobType=input_ppMC_condorDir.substr( jobTypePos);
  //else                              jobType="noJobType";
  //if(debugMode)std::cout<<"jobType string is = "<<jobType<<std::endl;
  
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=input_ppData_condorDir.find("_ak")+3;  
  const std::string radiusInt= input_ppData_condorDir.substr( radPos,1);
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=input_ppData_condorDir.find("Jets");;
  //if(debugMode)std::cout<<"radPos is = "<<radPos<<std::endl;
  //if(debugMode)std::cout<<"jetsPos is = "<<jetsPos<<std::endl;
  //if(debugMode)std::cout<<"jetTypePos is ="<<jetTypePos<<std::endl;
  
  
  
  //make final jetType strings
  const std::string jetType=input_ppData_condorDir.substr( jetTypePos,(jetsPos-jetTypePos) );
  const std::string fullJetType="ak"+radiusInt+jetType;
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  

  //put together input file strings
  //const std::string input_ppData_Filename="ppData_HighPtJet80_" +fullJetType+ "Jets-allfiles.root";
  const std::string input_ppData_Filename="ppData_HighPtJet80andLower_" +fullJetType+ "Jets-allfiles.root";
  const std::string input_ppMC_Filename  ="ppMC_Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "Jets-allfiles.root";
  const std::string ppData_fullFilename=inputDir+input_ppData_condorDir+input_ppData_Filename;
  const std::string ppMC_fullFilename  =inputDir+input_ppMC_condorDir+input_ppMC_Filename;
  
  




  // OPEN INPUT SECTION
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;  

  
  
  TFile *finData=NULL;  
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppData_condorDir<<std::endl<<std::endl;
  finData = new TFile(ppData_fullFilename.c_str(), "READ");      
  if(!finData){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_Filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);    
  }  
  
  TFile *finMC=NULL;
  std::cout<<" now opening ppMC: "<<std::endl<<input_ppMC_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppMC_condorDir<<std::endl<<std::endl;  
  finMC = new TFile(ppMC_fullFilename.c_str(), "READ");      
  if(!finMC){
    std::cout << " MC file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppMC_condorDir =" << input_ppMC_condorDir << std::endl;
    std::cout << "input_ppMC_filename  =" << input_ppMC_Filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);    
  }
  
  
  
  
 



  // GET OUTPUT PDF FILE READY
  std::string thePDFFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetPlots.pdf";
  std::string open_thePDFFileName=thePDFFileName+"[";    std::string close_thePDFFileName=thePDFFileName+"]";  
  
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1000, 800);
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  //temp_canvOpen->UseCurrentStyle();
  temp_canvOpen->Close();  
  

  // evtcounts/effective integrated luminosity ----------------------  
  long double theLumi;
  if(doEventCounts){
    printDataEventCountReport((TFile*) finData);
    //getEventCounts( (TFile*)finData, true );
    //getEventCounts( (TFile*)finMC,   false);
    theLumi=computeEffLumi( (TFile*) finData);   }
  else {
    std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
    theLumi=intgrtdLumi;}
  


  // evt plots ----------------------
  if(drawEvtQAPlots){
    
    std::cout<<std::endl<<" drawing evt QA Plots now! "<<std::endl<<std::endl;
    
    printMCEvtQAHist( (TFile*) finMC   , "hpthat" ,
		      (std::string) thePDFFileName );
    
    printMCEvtQAHist( (TFile*) finMC   , "hWeightedpthat" ,
		      (std::string) thePDFFileName );
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVz" , 
		       (TFile*) finMC   , "hWeightedVz" ,
		       (int) 10, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVy" , 
		       (TFile*) finMC   , "hWeightedVy" ,
		       (int) 5, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVx" , 
		       (TFile*) finMC   , "hWeightedVx" ,
		       (int) 5, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVr" , 
		       (TFile*) finMC   , "hWeightedVr" ,
		       (int) 5, (std::string) thePDFFileName  , (long double) theLumi  ) ;

        
    //printEvtNrefQAHist( (TFile*) finData , "hWNref" , 
    //			(TFile*) finMC   , "hWNref" ,
    //			(std::string) thePDFFileName  , (std::string) fullJetType, 
    //			(long double) theLumi  );
    //
    //printEvtNrefQAHist( (TFile*) finData , "hWjetsPEvt" , 
    //			(TFile*) finMC   , "hWjetsPEvt" ,
    //			(std::string) thePDFFileName  , (std::string) fullJetType, 
    //			(long double) theLumi  );
    //
    //printEvtNrefQAHist( (TFile*) finData , "hLeadJetPt" , 
    //			(TFile*) finMC   , "hLeadJetPt" ,
    //			(std::string) thePDFFileName  , (std::string) fullJetType, 
    //			(long double) theLumi  );

    
    
  }  
  else std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl; 
  
  
  
  // jet plots ----------------------  
  if(drawJetQAPlots){
    std::cout<<" drawing jet QA Plots..."<<std::endl;    
    
    std::string jetIDInt;
    if(doJetIDPlots)jetIDInt="1";
    else jetIDInt="0";  
    
    //TH1s
    for(int j=0;j<N_vars;j++){ 
      std::cout<<std::endl;
      if(debugMode)std::cout<<std::endl<<" var ="<<var[j]<<", j="<<j<<std::endl;
      
      //constituent/dijet skips
      bool skipConstitPlot= (!drawJetConstituentPlots && j>=jetConstits_varStart && j<dijet_varStart);
      bool skipDijetPlot= (!drawDijetPlots && j>=dijet_varStart);
      
      //skip plot
      bool skipPlot=skipConstitPlot||skipDijetPlot;
      if(skipPlot) {
	if(debugMode)std::cout<<"skipping jet plot for "<<var[j]<<std::endl;
	continue;}
      
      //printPlot      
      std::string inHistName="hJetQA_"+jetIDInt+"wJetID_"+var[j];          
      printJetQAHist( (TFile*) finData , (TFile*) finMC   ,  (int) j, (bool)doJetIDPlots, 
		      (std::string) inHistName , (std::string) thePDFFileName  , (std::string) fullJetType, 
		      (long double) theLumi  );      
    }
    
  }
  
  
  
  
  
  // jet plots ----------------------  
  if(drawJetRapBinsPlot){
    
    std::cout<<" drawing dual-diff xsec plot(s)..."<<std::endl;
    
    //print dual diff xsec w/ data and MC (RECO)
    printJetSpectraRapHists( (TFile*) finData , (TFile*) finMC   ,  (bool)doJetIDPlots, 
			     (std::string) thePDFFileName  , (std::string) fullJetType, 
			     (long double) theLumi  );   
    
    //print dual diff xsec w/ data only (RECO)
    printJetSpectraRapHists( (TFile*) finData , NULL   ,  (bool)doJetIDPlots, 
    			     (std::string) thePDFFileName  , (std::string) fullJetType, 
			     (long double) theLumi  );   
    
    //print dual diff xsec w/ MC only (RECO v GEN)
    printMCJetSpectraRapHists( (TFile*) finMC   ,  (bool)doJetIDPlots, 
			       (std::string) thePDFFileName  , (std::string) fullJetType);

  }
  





  
  
  
  
  
  


  
  
  // for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1200, 600);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();  
  
  return 0;
  
}// end printplots

				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=4 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetPlots.exe <ppDataCondorDir> <ppMCCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  //if(argc==1)rStatus=printPlots_jetPlots();
  if(argc==4) rStatus=printPlots_jetPlots( (const std::string) argv[1], (const std::string) argv[2],
					(const std::string) argv[3]);
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
   
}



///////////////////////
//// PROBLEM SECTION, do NOT delete comments or other commented code
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//
//if(debugMode)std::cout<<std::endl<<"N_trigs="<<N_trigs<<" and j="<<j<<std::endl;
//if(debugMode)std::cout<<"j<(N_trigs-1)="<<(bool)(j<(N_trigs-1))<<std::endl;
//if(debugMode)std::cout<<"j<N_trigs="<<(bool)(j<(N_trigs))<<std::endl;
//if(debugMode)std::cout<<"j<(N_trigs+1)="<<(bool)(j<(N_trigs+1))<<std::endl;	
//
////assert(j!=N_trigs);// NEVER FIRES
////assert(j<N_trigs); // NEVER FIRES
//
////open the hists + do scaling 
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//
//if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
//  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
//  continue;}  
//else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
//
////if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
////  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
////  continue;}	
////else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
//  
////open the hists + do scaling 
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//// PROBLEM SECTION
///////////////////////








//const bool drawJetTrigQAPlots=true, drawJetRapBinsPlot=false;
//const bool drawDataMCOverlaysInput=false; 		 //this should always be false, but needs to be cleaned up later.
//const bool drawDataMCOverlays   = drawDataMCOverlaysInput;
//const bool drawDataMCRatios     = true;   //!drawDataMCOverlays;
//const bool drawTrigDataOverlays = false;  //drawDataMCOverlays && drawJetTrigQAPlots;
//const bool drawTrigDataRatios   = false; //drawDataMCRatios   && drawJetTrigQAPlots;
//const int theTrigOverlayLineColor[]  ={  1,  1,  1,  1,  1,  1 };
//const int theTrigCombMarkerColor=1, altTrigCombMarkerColor=12;
//const int theTrigOverlayMarkerColor[]={  2,  3,  6,  7,  1,  4 };
//const int theRapOverlayMarkerColor[] ={  2,  3,  6,  7,  1,  4 , 8, 9};
//const int theTrigOverlayMarker[]     ={ 20, 20, 20, 20, 20, 32 };
//const std::string crossSectionAxTitle="#sigma (#mub)";
//const std::string AUAxTitle="A.U.";
//const std::string ratioTitle="MC/Data";
