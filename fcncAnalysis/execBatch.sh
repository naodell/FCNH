#!/bin/sh

### Required parameters #####

dataName  = $1
count 	  = $2

### Specify addtional arguments here ####
suffix    = $3
selection = $4
period    = $5

### Transfer files, prepare directory ###
source /etc/bashrc prod
export OSG_APP=/software/tier3/osg
export SCRAM_ARCH=slc5_amd64_gcc462
source /software/tier3/osg/cmsset_default.sh

scram p CMSSW CMSSW_5_3_2
cd ./CMSSW_5_3_2/src
cmsenv 

#cd ${_CONDOR_SCRATCH_DIR}

cp ../../source.tar.gz .
tar -xzf source.tar.gz
cd Analysis_CMS/fcncAnalysis
cp ../../../../input_${dataName}_${count}.txt input.txt
rm histos/*root

./fcncLocal.csh $suffix $selection $period

### Copy output and cleanup ###
rename .root _${dataName}_$count.root histos/fcncHistograms*

#cp histos/fcncHistograms_*.root $outDir/.
cp histos/fcncHistograms_*.root ${_CONDOR_SCRATCH_DIR}
