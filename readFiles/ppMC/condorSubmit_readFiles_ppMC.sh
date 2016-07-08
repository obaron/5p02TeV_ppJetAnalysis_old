#!/bin/bash

echo ""
# error conditions 
if [[ $# -ne 5 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_readFiles_ppMC.sh <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <debug>"
    echo "set <NJobs>=-1 to run the entire filelist with jobs of size <NFilesPerJob>"
    return 1
#elif [[ ! $3 =~ ^-?[0-9]+$ ]] # check integer input
#then
#    echo "non integer <startFilePos>, exit"
#    return 1
#elif [[ $3 -lt 0 ]] # check for valid startFilePos
#then
#    echo "bad third argument"
#    echo "0 <= <startFilePos> <= nFiles-1"
#    return 1
fi

# input arguments to submit script
NJobs=$1
NFilesPerJob=$2
startFilePos=$3
filelistIn=$4  #echo "filelistIn is ${filelistIn}" #debug
debug=$5

# one condor job submit per NFilesPerJob until we submit NJobs
nFiles=`wc -l < $filelistIn`
if [[ $NJobs -eq -1 ]]
then
    echo "submitting jobs for all files in list"
    NJobs=$(( $nFiles / $NFilesPerJob ))
    if [[ $(( $nFiles % $NFilesPerJob ))  -gt 0 ]]
    then
    NJobs=$(( $NJobs + 1 ))
    fi
    startFilePos=0
elif [[ $startFilePos -ge $nFiles ]]
then
    echo "bad <startFilePos>, exit"
    echo "0 <= <startFilePos> < nFiles-1"
    return 1
elif [[ ! $startFilePos =~ ^-?[0-9]+$ ]]
then
    echo "non integer <startFilePos>, exit"
    return 1
elif [[ $startFilePos -lt 0 ]]
then
    echo "bad third argument"
    echo "0 <= <startFilePos> <= nFiles-1"
    return 1
fi

# some debug info, just in case
NFilesRequested=$(( $NJobs * $NFilesPerJob ))
echo "require ${NFilesRequested} files for ${NJobs} jobs"
echo "# of files in list: ${nFiles}"
echo "starting at file position ${startFilePos}..."
echo ""

# additional inputs to the run script and .exe, these don't change too much
radius=4
jetType="PF"
# old and shoul consider removing
destination="/mnt/hadoop/cms/store/user/ilaflott/5p02TeV_ppJetAnalysis/readFiles/ppMC" 

# create output folder/logfileNames with name based on filelist
filelist=${filelistIn##*/} #
echo "filelist is ${filelist}"
filelistTitle=${filelist%_*} #
echo "filelistTitle is ${filelistTitle}"
energy=${filelistTitle%%_*} #
echo "energy is ${energy}"
trig=${filelistTitle#*_} #
echo "trig is ${trig}"
dirName="readFiles_ppMC_${energy}_${trig}_$(date +"%Y-%m-%d__%H_%M")" #
echo "dirName is ${dirName}"
outName="${trig}_ak${radius}${jetType}" #
echo "outName is ${outName}"
logFileDir="${PWD}/outputCondor/${dirName}" #
echo "logFileDir is ${logFileDir}"
if [ -d "${logFileDir}" ]; then
    rm -rf "${logFileDir}"
fi
mkdir $logFileDir
echo "log files in outputCondor/${dirName}"

# cmsenv for condor
echo "cmsenv'ing..."
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src
cmsenv
cd -

# compile code executable, same as rootcompile in my .bashrc
echo "compiling..."
g++ readFiles_ppMC.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readFiles_ppMC.exe || return 1
cp readFiles_ppMC.* "${logFileDir}"
cp condorRun_readFiles_ppMC.sh "${logFileDir}"
cp ${filelistIn} "${logFileDir}"
cd ${logFileDir}

### CREATE NAMES AND FILES, THEN SUBMIT ###
NthJob=0
filelistLength=$nFiles
lastFilePos=$(( $filelistLength-1 ))
startfile=0
endfile=0
while [ $NthJob -lt $NJobs ]
do 
    echo ""
    echo "SPLITTING FILES FOR JOB #${NthJob}"

    # start/end file 
    if [[ $NthJob -le 0 ]]
    then
	startfile=$startFilePos
	endfile=$(( $startfile + $NFilesPerJob ))
	endfile=$(( $endfile - 1 ))
    else 
	startfile=$(( $endfile + 1 ))
	endfile=$(( $startfile + $NFilesPerJob ))
	endfile=$(( $endfile - 1 ))
    fi
    #echo "startfile is ${startfile}" 
    # check; end of filelist
    if [[ $endfile -ge $lastFilePos ]] 
    then 
	echo "end of filelist!"
	let endfile=$lastFilePos
	let NthJob=$(( $NJobs - 1 )) 
    fi
    #echo "endfile is ${endfile}"     

    
    # for next job
    NthJob=$(($NthJob + 1))

    # define output names for job submission
    fileRange="${startfile}to${endfile}"
    Error="${outName}-${fileRange}.err"
    Output="${outName}-${fileRange}.out"
    Log="${outName}-${fileRange}.log"
    outfile="${outName}-${fileRange}.root"
    
    # create the condor submit file
    cat > ${logFileDir}/subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
Executable     = condorRun_readFiles_ppMC.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $startfile $endfile $filelist $outfile $radius $jetType $debug
Input          = /dev/null
Error          = ${logFileDir}/$Error
Output         = ${logFileDir}/$Output
Log            = ${logFileDir}/$Log
# get the environment (path, etc.)
GetEnv         = True
# prefer to run on fast, 64 bit computers
Rank           = kflops
Requirements   = Arch == "X86_64"
should_transfer_files   = YES
transfer_input_files = ${filelist},readFiles_ppMC.exe
when_to_transfer_output = ON_EXIT
Queue
EOF

    # submit the job defined in the above submit file
    echo "running readFiles_ppMC on files #${startfile} to #${endfile}"
    condor_submit ${logFileDir}/subfile    
    sleep 1s #my way of being nicer to condor, not sure it really matters but i'm paranoid
done

cd -
echo "done."
return