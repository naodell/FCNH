#ifndef fakeAnalyzer_h
#define fakeAnalyzer_h

// c++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <sstream> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <bitset>

// ROOT libraries
#include <TROOT.h>
#include <TChain.h>
#include <TSelector.h>
#include <TH2.h>
#include <TStyle.h>
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TProfile.h"
#include "TMatrix.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TCanvas.h"

// object definitions 
#include "../interface/TCPhysObject.h"
#include "../interface/TCJet.h"
#include "../interface/TCMET.h"
#include "../interface/TCElectron.h"
#include "../interface/TCMuon.h"
#include "../interface/TCTau.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCGenJet.h"
#include "../interface/TCGenParticle.h"
#include "../interface/TCPrimaryVtx.h"
#include "../interface/TCTriggerObject.h"

// plugins
#include "../plugins/WeightUtils.h"
#include "../plugins/Selector.h"
#include "../plugins/TriggerSelector.h"
#include "../plugins/HistManager.h"

using namespace std;

typedef vector<TCPhysObject> vObj;

class fakeAnalyzer : public TSelector {

    private:

        TFile* histoFile;
        TTree* thisTree;

        // Utilities and selectors
        WeightUtils     *weighter;
        HistManager     *histManager;
        Selector        *selector;
        TriggerSelector *triggerSelector;

        // Useful global variables
        string  suffix;
        string  selection;
        string  period;

        bool    zTagged;
        bool    ossfTagged;

        TVector3*       selectedVtx;
        TLorentzVector  ZP4;
        TCPhysObject    tag, lep1, lep2, lep3;
        //TCJet           probeJet;
        //vector<TCJet>   jets;

        Float_t dileptonMass;

    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types                                                                                                               
        TClonesArray    *recoJets;
        //TClonesArray    *recoJPT;
        TClonesArray    *recoElectrons;
        TClonesArray    *recoMuons;
        TClonesArray    *recoTaus;
        TClonesArray    *recoPhotons;
        TCMET           *recoMET;
        TClonesArray    *genJets;
        TClonesArray    *genParticles;
        TClonesArray    *primaryVtx;
        TVector3        *beamSpot;

        Int_t           nPUVertices;
        Float_t         nPUVerticesTrue;
        Bool_t          isRealData;
        UInt_t          runNumber;
        ULong64_t       eventNumber;
        UInt_t          lumiSection;
        UInt_t          bunchCross;
        Float_t         ptHat;
        Float_t         qScale;
        Float_t         rhoFactor;
        Float_t         rho25Factor;
        Float_t         rhoMuFactor;
        ULong64_t       triggerStatus;
        UInt_t          hltPrescale[64];

        Bool_t          NoiseFilters_isScraping;
        Bool_t          NoiseFilters_isNoiseHcalHBHE;
        Bool_t          NoiseFilters_isNoiseHcalLaser;
        Bool_t          NoiseFilters_isNoiseEcalTP;
        Bool_t          NoiseFilters_isNoiseEcalBE;
        Bool_t          NoiseFilters_isCSCTightHalo;
        Bool_t          NoiseFilters_isCSCLooseHalo;

        // List of branches
        TBranch        *b_recoJets;   //!
        TBranch        *b_recoElectrons;   //!
        TBranch        *b_recoMuons;   //!
        TBranch        *b_recoTaus;   //!
        TBranch        *b_recoPhotons;   //!
        TBranch        *b_recoMET;   //!
        TBranch        *b_genJets;   //!
        TBranch        *b_genParticles;   //!
        TBranch        *b_primaryVtx;   //!
        TBranch        *b_beamSpot;   //!
        TBranch        *b_nPUVertices;   //!
        TBranch        *b_nPUVerticesTrue;   //! 
        TBranch        *b_isRealData;   //!
        TBranch        *b_runNumber;   //!
        TBranch        *b_eventNumber;   //!  
        TBranch        *b_lumiSection;   //! 
        TBranch        *b_bunchCross;   //!
        TBranch        *b_ptHat;   //!   
        TBranch        *b_qScale;   //!
        TBranch        *b_rhoFactor;   //! 
        TBranch        *b_rho25Factor;   //!
        TBranch        *b_rhoMuFactor;   //!
        TBranch        *b_triggerStatus;   //!
        TBranch        *b_hltPrescale;   //!
        TBranch        *b_NoiseFilters;   //!

        //For counting events
        unsigned          eventCount[16];

        fakeAnalyzer(TTree * /*tree*/ =0) { }
        virtual ~fakeAnalyzer() { }
        virtual int     Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        //virtual void    SlaveBegin(TTree *tree) { TString option = GetOption();};
        virtual void    Init(TTree *tree);
        virtual bool    Notify();
        virtual bool    Process(Long64_t entry);
        virtual int     GetEntry(Long64_t entry, int getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        //virtual void    SlaveTerminate() {};
        virtual void    Terminate();

        virtual void    FillDenominatorHists(TCPhysObject&, vector<TCJet>&);
        virtual void    FillNumeratorHists(TCPhysObject&, vector<TCJet>&);
        virtual void    FillClosureHists(TCPhysObject&, vector<TCJet>&);
        virtual void    FillJetHists(TCPhysObject&, vector<TCJet>&, string);
        virtual bool    GenProbeMatcher(TCPhysObject&, vector<TCGenParticle>&);
        virtual void    DoZTag(vObj&, float);
        virtual bool    CheckQCD2lCR(vector<TCJet>&, TCPhysObject&);
        virtual bool    CheckZPlusJetCR(TCPhysObject&);
        virtual float   CalculateTransMass(TCPhysObject&, TCMET&);
        virtual void    SetEffWeights(vObj&, vector<TCJet>&, string);

        virtual string  str(int i) {return static_cast<ostringstream*>( &(ostringstream() << i) )->str();}


        ClassDef(fakeAnalyzer,0);
};

#endif

#ifdef fakeAnalyzer_cxx
void fakeAnalyzer::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    recoJets = 0;
    recoMET = 0;
    genJets = 0;
    genParticles = 0;
    primaryVtx = 0;
    recoMuons = 0;
    recoElectrons = 0;
    recoPhotons = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;

    thisTree    = tree; 
    fChain      = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
    fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
    fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
    fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
    fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
    fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
    fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
    fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);

    //fChain->SetBranchAddress("beamSpot", &beamSpot, &b_beamSpot);
    fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
    fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
    fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
    fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);
    fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
    fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
    fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
    fChain->SetBranchAddress("rho25Factor", &rho25Factor, &b_rho25Factor);
    fChain->SetBranchAddress("rhoMuFactor", &rhoMuFactor, &b_rhoMuFactor);
    fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
    fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
    fChain->SetBranchAddress("NoiseFilters", &NoiseFilters_isScraping, &b_NoiseFilters);

}

bool fakeAnalyzer::Notify()
{

    UInt_t initEvents;
    TFile   *inFile     = thisTree->GetCurrentFile();
    TTree   *jobTree    = (TTree*)inFile->Get("ntupleProducer/jobTree");
    TBranch *nEvents    =  jobTree->GetBranch("nEvents");

    cout << "Analyzing file " << inFile->GetName() << endl;

    nEvents->SetAddress(&initEvents);
    nEvents->GetEntry(0);

    eventCount[0] += initEvents;

    return kTRUE;
    return kTRUE;
}

#endif 
