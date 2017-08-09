#!/bin/bash

#echo ""
#echo ""
#echo "compiling printPlots_MCJEC"
#rootcompile printPlots_MCJEC.C
#echo ""
#echo "Running..."
#echo ""
#echo ""
#./printPlots_MCJEC.exe  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-15-16_MCJEC_jtID0v2/ jetID0v2
#./printPlots_MCJEC.exe  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-15-16_MCJEC_jtIDTv2/ jetIDTv2

echo ""
echo ""
echo "compiling printPlots_jetPlots"
echo ""
echo ""
rootcompile printPlots_jetPlots.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""
./printPlots_jetPlots.exe  ppData_HighPtJetTrig_ak4PFJets_12-15-16_jetPlots_jtIDv2/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-15-16_jetPlots_jtIDv2/ jetIDv2_eta2 0

./printPlots_jetPlots.exe  ppData_HighPtJetTrig_ak4PFJets_12-15-16_jetPlots_jtIDv2/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-15-16_jetPlots_jtIDv2/ jetIDv2_eta2 1







#echo ""
#echo ""
#echo "compiling eta check, genpt 30 to 50..."
#rootcompile printPlots_MCJEC_genpt30to50.C
#echo ""
#echo "done compiling. Running..."
#echo ""
#echo ""
#./printPlots_MCJEC_genpt30to50.exe  12.09.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-09-16_MCJEC_jtID0 jtID0 
#
#
#echo ""
#echo ""
#echo "compiling eta check, genpt 150 to 200..."
#rootcompile printPlots_MCJEC_genpt150to200.C
#echo ""
#echo "done compiling. Running..."
#echo ""
#echo ""
#./printPlots_MCJEC_genpt150to200.exe  12.09.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-09-16_MCJEC_jtID0 jtID0 
#
#
#
#echo ""
#echo ""
#echo "compiling printPlots_MCJEC_MCEff"
#rootcompile printPlots_MCJEC_MCEff.C
#echo ""
#echo "done compiling. Running..."
#echo ""
#echo ""
#./printPlots_MCJEC_MCEff.exe  12.09.16_outputCondor!/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-09-16_MCJEC_jtID0 jtID0

 