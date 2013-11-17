#ifndef fcncAnalyzer_h
#define fcncAnalyzer_h

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
#include "TFile.h"
#include "TProfile.h"
#include "TMatrix.h"
#include "TMatrixT.h"
#include "TRandom3.h"

// MVA tools
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

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
#include "../plugins/Selector.h"
#include "../plugins/TriggerSelector.h"
#include "../plugins/WeightUtils.h"
#include "../plugins/HistManager.h"

using namespace std;

const unsigned short N_CATEGORIES = 16;
const string categoryNames[] = {
    //inclusive
    "inclusive",
    "ss_inclusive",
    "os_inclusive",
    "3l_inclusive",

    // flavor
    "ss_ee",
    "ss_emu",
    "ss_mumu",
    "os_ee",
    "os_emu",
    "os_mumu",
    "3l_eee",
    "3l_eemu",
    "3l_emumu",
    "3l_mumumu",

    // WH 
    "3l_SSSF",
    "3l_OSSF",

    // geometric
    //"ss_endcap",
    //"ss_mixed",
    //"ss_barrel",
    //"os_endcap",
    //"os_mixed",
    //"os_barrel",
    //"3l_endcap",
    //"3l_2end1bar",
    //"3l_1end2bar",
    //"3l_barrel",
};

const unsigned short N_CUTS = 8;
const string cutNames[] = {"preselection", "Z veto",  "MET selection", "jet selection", "BDT", "WZ_CR", "ttbar_CR", "ttZ_CR"};

typedef vector<TCPhysObject> vObj;

class fcncAnalyzer : public TSelector {

    private:

        TFile* histoFile[N_CUTS];
        TTree* thisTree;
        ofstream fout[8];

        // Job setup
        string  suffix;
        string  selection;
        string  period;

        // Random number generator
        TRandom3* rnGenerator;

        // Utilities and selectors
        WeightUtils     *weighter;
        HistManager     *histManager;
        Selector        *selector;
        TriggerSelector *triggerSelector;

        //Event variables
        bool zTagged;
        bool ossfTagged;
        string subdir;

        TLorentzVector dileptonP4;
        TLorentzVector lep1P4, lep2P4, lep3P4; // If event is zTagged, lep1 and lep2 are associated to the z

        Float_t     MET;
        Float_t     metPhi;
        Float_t     METLD;
        Float_t     HT;
        Float_t     HTs;
        Float_t     MHT;
        Float_t     MT;

        Int_t       lepMult;
        Int_t       bJetMult;
        Int_t       jetMult;

        bitset<18>  evtCategory;
        UInt_t      flavorCat;
        UInt_t      chargeCat;

        Float_t     evtWeight;

        // trees for lepton mva
        TTree*  muTree;
        TTree*  eleTree;

        Float_t sip3d;
        Float_t chPFIso, neuPFIso;
        Float_t drLepJet, ptRatioLepJet, btagLepJet;
        Float_t dxy, dz; // muon specific
        Float_t mva; // electron specific
        Int_t   missHits; 


        // Simple ntuples for MVA
        TTree*      mvaTree;

        Float_t     dileptonMassOS;
        Float_t     trileptonMass;
        Float_t     dileptonDROS;

        Float_t     lep1Pt, lep2Pt, lep3Pt;
        Float_t     lep1Eta, lep2Eta, lep3Eta;
        Float_t     lep1Phi, lep2Phi, lep3Phi;
        Float_t     bJetPt, bJetEta, bJetPhi;

        // MVA reader and modified input variables
        TMVA::Reader* mvaReader;
        Float_t     f_flavorCat;
        Float_t     f_chargeCat;
        Float_t     f_bJetMult;
        Float_t     f_jetMult;

    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain

        // Declaration of leaf types                                                                                                               
        TClonesArray    *recoJets;
        //TClonesArray    *recoJPT;
        TClonesArray    *recoElectrons;
        TClonesArray    *recoMuons;
        TClonesArray    *recoPhotons;
        TCMET           *recoMET;
        TClonesArray    *triggerObjects;
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
        TBranch        *b_recoPhotons;   //!
        TBranch        *b_recoMET;   //!
        TBranch        *b_triggerObjects;   //!
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
        //TBranch        *b_evtWeight;   //! 
        TBranch        *b_rhoFactor;   //! 
        TBranch        *b_rho25Factor;   //!
        TBranch        *b_rhoMuFactor;   //!
        TBranch        *b_triggerStatus;   //!
        TBranch        *b_hltPrescale;   //!
        TBranch        *b_NoiseFilters;   //!

        //For counting events
        int          eventCount[16];
        //For counting weighted events
        float        eventCountWeighted[16];

        fcncAnalyzer(TTree * /*tree*/ =0) { }
        virtual ~fcncAnalyzer() { }
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

        virtual bool    AnalysisSelection(vObj, vector<TCJet>, vector<TCJet>, vector<TCJet>, TVector3, string);
        virtual void    GetFakeBG(vObj, vObj, vector<TCJet>, vector<TCJet>, vector<TCJet>, TVector3);

        // Plot methods
        virtual void    MakePlots(vObj, vector<TCJet>, vector<TCJet>, TCMET, TVector3, unsigned);
        virtual void    MakeQMisIDPlots(vObj);
        virtual void    Make4lPlots(vObj, TCMET); //, vector<TCJet>, vector<TCJet>);
        virtual void    LeptonPlots(vObj, vector<TCJet>, vector<TCJet>, TVector3);
        virtual void    JetPlots(vector<TCJet>, vector<TCJet>);
        virtual void    MetPlots(TCMET, vObj);
        virtual void    DileptonPlots2D(vObj);
        virtual void    GenPlots(vector<TCGenParticle>, vObj);
        virtual void    MiscPlots();
        virtual void    FillYieldHists(string, float, unsigned);

        // Set/Get methods
        virtual void    SetEventCategory(vObj);
        virtual void    SetVarsMVA(vObj, vector<TCJet>, vector<TCJet>);
        virtual void    FillLepMVA(vector<TCMuon>, vector<TCElectron>, vector<TCJet>, TVector3);
        virtual void    SetEventVariables(vObj, vector<TCJet>, vector<TCJet>, TCMET);
        virtual void    SetYields(unsigned);
        virtual int     GetHistCategory(unsigned);

        // helper functions
        virtual string  str(int i) {return static_cast<ostringstream*>( &(ostringstream() << i) )->str();}
        virtual TLorentzVector  CalculateNuP4(TLorentzVector, TCMET);
        virtual float           CalculateTransMass(TCPhysObject, TCMET);
        virtual bool            CosmicMuonFilter(TCPhysObject, TCPhysObject);
        virtual float           CalculateFourLeptonMass(vObj);
        //virtual float           CalculateChi2Mass(vObj, vObj, vObj, TCMET);

        ClassDef(fcncAnalyzer,0);
};

#endif

#ifdef fcncAnalyzer_cxx
void fcncAnalyzer::Init(TTree *tree)
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
    triggerObjects = 0;

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
    fChain->SetBranchAddress("triggerObjects", &triggerObjects, &b_triggerObjects);
    fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
    fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
    fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);

    //fChain->SetBranchAddress("beamSpot", &beamSpot, &b_beamSpot);
    fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
    fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
    fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
    fChain->SetBranchAddress("rho25Factor", &rho25Factor, &b_rho25Factor);
    fChain->SetBranchAddress("rhoMuFactor", &rhoMuFactor, &b_rhoMuFactor);

    fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
    fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);

    fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
    fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
    fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
    fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
    fChain->SetBranchAddress("NoiseFilters", &NoiseFilters_isScraping, &b_NoiseFilters);
}

bool fcncAnalyzer::Notify()
{
    UInt_t initEvents;
    TFile   *inFile  = thisTree->GetCurrentFile();
    TTree   *jobTree = (TTree*)inFile->Get("ntupleProducer/jobTree");
    TBranch *nEvents = jobTree->GetBranch("nEvents");

    cout << "Analyzing file " << inFile->GetName() << endl;

    nEvents->SetAddress(&initEvents);
    nEvents->GetEntry(0);

    evtWeight = initEvents;
    evtCategory.reset();

    SetYields(0);

    return kTRUE;
}

#endif 
