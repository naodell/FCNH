#!/bin/sh

### Required parameters #####

DATANAME=$1
COUNT=$2

### Specify addtional arguments here ####
SUFFIX=$3
SELECTION=$4
PERIOD=$5

echo ${_CONDOR_SCRATCH_DIR}
### Transfer files, prepare directory ###
source /etc/bashrc prod
export OSG_APP=/software/tier3/osg
export SCRAM_ARCH=slc5_amd64_gcc462
source /software/tier3/osg/cmsset_default.sh

scram p CMSSW CMSSW_5_3_10
cd ./CMSSW_5_3_10/src
cmsenv 

cp ../../source.tar.gz .
tar -xzf source.tar.gz
cd Analysis_CMS/fcncAnalysis
cp ../../../../input_${DATANAME}_${COUNT}.txt input.txt
rm histos/*root

root -l -b -q 'run.C(1e9, "'$SUFFIX' '$SELECTION' '$PERIOD'")'

### Copy output and cleanup ###
rename .root _${DATANAME}_$COUNT.root histos/fcncHistograms*

cp histos/fcncHistograms_*.root ${_CONDOR_SCRATCH_DIR}
