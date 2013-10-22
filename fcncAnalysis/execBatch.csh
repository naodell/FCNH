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
cd Analysis_CMS/fcncAnalysis
cp ../../../../input_${dataName}_${count}.txt input.txt
rm histos/*root

root -l -b -q 'run.C(1e9, "'$suffix' '$selection' '$period'")'

### Copy output and cleanup ###
rename .root _${dataName}_$count.root histos/fcncHistograms*

#cp histos/fcncHistograms_*.root $outDir/.
cp histos/fcncHistograms_*.root ${_CONDOR_SCRATCH_DIR}
