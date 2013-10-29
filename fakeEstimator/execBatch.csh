#!/bin/csh

### Required parameters #####

set dataName  = $1
set count     = $2

### Specify addtional arguments here ####
set suffix    = $3
set selection = $4
set period    = $5

### Transfer files, prepare directory ###
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_5_3_2
cd ./CMSSW_5_3_2/src
cmsenv 

#cd ${_CONDOR_SCRATCH_DIR}

cp ../../source.tar.gz .
tar -xzf source.tar.gz
cd Analysis_CMS/fakeEstimator
cp ../../../../input_${dataName}_${count}.txt input.txt

root -l -b -q 'run.C(1e9, "'$suffix' '$selection' '$period'")'

### Copy output and cleanup ###
cp fakeHistograms.root ${_CONDOR_SCRATCH_DIR}/fakeHistograms_${dataName}_${count}.root
