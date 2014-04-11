#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void run(Long64_t nEntries = 1e5, string args = "TEST muon 2012") {

    //container classes
    gROOT->LoadMacro("../src/TCPhysObject.cc+");
    gROOT->LoadMacro("../src/TCTrack.cc+");
    gROOT->LoadMacro("../src/TCEGamma.cc+");
    gROOT->LoadMacro("../src/TCJet.cc+");
    gROOT->LoadMacro("../src/TCMET.cc+");
    gROOT->LoadMacro("../src/TCElectron.cc+");
    gROOT->LoadMacro("../src/TCMuon.cc+");
    gROOT->LoadMacro("../src/TCTau.cc+");
    gROOT->LoadMacro("../src/TCPhoton.cc+");
    gROOT->LoadMacro("../src/TCGenJet.cc+");
    gROOT->LoadMacro("../src/TCGenParticle.cc+");
    gROOT->LoadMacro("../src/TCPrimaryVtx.cc+");
    gROOT->LoadMacro("../src/TCTriggerObject.cc+");

    //analysis plugins (selectors, utiltities, etc.)
    gROOT->LoadMacro("../plugins/HistManager.cc+");
    gROOT->LoadMacro("../plugins/EGammaMvaEleEstimator.cc+");
    gROOT->LoadMacro("../plugins/rochcor2012jan22.C+");
    gROOT->LoadMacro("../plugins/WeightUtils.cc+");
    gROOT->LoadMacro("../plugins/TriggerSelector.cc+");
    gROOT->LoadMacro("../plugins/Selector.cc+");

    TChain* fChain = new TChain("ntupleProducer/eventTree");

    ifstream sourceFiles("input.txt");
    char line[2048];
    int  count = 0;

    while (sourceFiles >> line) {
        fChain->Add(line);      
        ++count;
    }
    cout << count << " files added!"<<endl;
    sourceFiles.close();

    TStopwatch timer;
    timer.Start();

    fChain->Process("fakeAnalyzer.C+", args.c_str(), nEntries, 0);

    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() << endl;
    cout << "RealTime : " << timer.RealTime() << endl;
    cout << "\n";
}

