#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void run(Long64_t nEntries = 1e6, string suffix = "DATA_MUON", string selection = "3l") 
{

    //container classes
    gROOT->LoadMacro("../src/TCPhysObject.cc+");

    TChain* fChain = new TChain(("tree" + selection + "_" + suffix).c_str());

    fChain->Add("~/eos/histos/fcnh_cut1_2012_20131217_205058.root");      

    TStopwatch timer;
    timer.Start();

    fChain->Process("analyzer.C+", (selection + "_" + suffix).c_str(), nEntries, 0);

    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() << endl;
    cout << "RealTime : " << timer.RealTime() << endl;
    cout << "\n";
}

