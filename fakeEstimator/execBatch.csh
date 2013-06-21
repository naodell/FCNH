#!/bin/csh

### Required parameters #####

set outDir    = $1
set count     = $2
set dataName  = $3

### Specify addtional arguments here ####
set suffix    = $4
set selection = $5
set period    = $6

### Transfer files, prepare directory ###
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_5_3_2
cd CMSSW_5_3_2/src
cmsenv 

cp ../../source.tar.gz .
tar -xzf source.tar.gz
cd analysis/fakeEstimator
cp ../../../../input_${dataName}_${count}.txt input.txt
rm histos/*root

sed -i "s/SUFFIX/${suffix}/g" fakeAnalyzer.C
root -l -b -q run.C

### Copy output and cleanup ###
rename .root _${dataName}_$count.root histos/fakeHistograms*
cp histos/fakeHistograms_*.root $outDir/.
