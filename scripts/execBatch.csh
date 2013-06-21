#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_4_1
cd CMSSW_4_4_1/src
cmsenv 
#### Leave this blank #######

#############################


set srcDir    = $1
set outDir    = $2
set count     = $3
set dataName  = $4

### Specify addtional arguments here ####
set suffix    = $5
set trigger   = $6
set selection = $7
set period    = $8

### Transfer files, prepare directory ###
mkdir data
mkdir plugins
mkdir printouts

cp $srcDir/higgsAnalyzer_Template.C ./higgsAnalyzer.C
cp $srcDir/higgsAnalyzer.h .
cp $srcDir/../src/TC*.cc .
cp $srcDir/../src/TC*.h .
cp $srcDir/data/*root data/.
cp $srcDir/plugins/*.* plugins/.
cp $srcDir/mvaWeights .

### Do analysis stuff ###
sed -i "s/SUFFIX/$suffix/g" higgsAnalyzer.C
sed -i "s/TRIGGER/$trigger/g" higgsAnalyzer.C
sed -i "s/SELECTION/$selection/g" higgsAnalyzer.C
sed -i "s/PERIOD/$period/g" higgsAnalyzer.C

cat input.txt

cat > run.C << +EOF

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void run() {

    gROOT->LoadMacro("TCJet.cc+");
    gROOT->LoadMacro("TCMET.cc+");
    gROOT->LoadMacro("TCGenJet.cc+");
    gROOT->LoadMacro("TCGenParticle.cc+");
    gROOT->LoadMacro("TCElectron.cc+");
    gROOT->LoadMacro("TCMuon.cc+");
    gROOT->LoadMacro("TCTau.cc+");
    gROOT->LoadMacro("TCPhoton.cc+");
    gROOT->LoadMacro("TCPrimaryVtx.cc+");
    gROOT->LoadMacro("TCTriggerObject.cc+");

    //analysis plugins (selectors, utiltities, etc.)
    gROOT->LoadMacro("plugins/WeightUtils.cc+");
    gROOT->LoadMacro("plugins/TriggerSelector.cc+");

    TChain* fChain = new TChain("ntupleProducer/eventTree");

    ifstream sourceFiles("input.txt");
    char line[1000];
    while (sourceFiles >> line) {
        fChain->Add(line);      
    }
    sourceFiles.close();

    TStopwatch timer;
    timer.Start();

    fChain->Process("higgsAnalyzer.C+");
}

+EOF

root -l -b -q run.C

### Copy output and cleanup ###
cp higgsHistograms.root $outDir/higgsHistograms_${dataName}_$count.root
rm -r ./*
