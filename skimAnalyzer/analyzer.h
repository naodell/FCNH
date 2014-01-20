//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 15 18:03:53 2014 by ROOT version 5.32/00
// from TTree tree3l_DATA_ELECTRON/Tree for cut MVA
// found on file: ../fcncAnalysis/combined_histos/fcnh_cut1_2012_20131217_205058.root
//////////////////////////////////////////////////////////

#ifndef analyzer_h
#define analyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1D.h>
#include "../interface/TCPhysObject.h"

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class analyzer : public TSelector {

    TFile* histoFile;

    TH1D*   h1_Met;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Float_t         evtWeight;
   Int_t           flavorCat;
   Int_t           chargeCat;
   Float_t         met;
   Float_t         metPhi;
   Float_t         HT;
   Float_t         MT;
   Int_t           jetMult;
   Int_t           bJetMult;
   Float_t         trileptonMass;
   Float_t         dileptonMassOS;
   Float_t         dileptonDROS;
   TClonesArray    *jets;
   TClonesArray    *leptons;

   // List of branches
   TBranch        *b_evtWeight;   //!
   TBranch        *b_flavorCat;   //!
   TBranch        *b_chargeCat;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MT;   //!
   TBranch        *b_jetMult;   //!
   TBranch        *b_bJetMult;   //!
   TBranch        *b_trileptonMass;   //!
   TBranch        *b_dileptonMassOS;   //!
   TBranch        *b_dileptonDROS;   //!
   TBranch        *b_jets;   //!
   TBranch        *b_leptons;   //!

   analyzer(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~analyzer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree) { TString option = GetOption();}
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate() {};
   virtual void    Terminate();

   ClassDef(analyzer,0);
};

#endif

#ifdef analyzer_cxx
void analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jets = 0;
   leptons = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("flavorCat", &flavorCat, &b_flavorCat);
   fChain->SetBranchAddress("chargeCat", &chargeCat, &b_chargeCat);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MT", &MT, &b_MT);
   fChain->SetBranchAddress("jetMult", &jetMult, &b_jetMult);
   fChain->SetBranchAddress("bJetMult", &bJetMult, &b_bJetMult);
   fChain->SetBranchAddress("trileptonMass", &trileptonMass, &b_trileptonMass);
   fChain->SetBranchAddress("dileptonMassOS", &dileptonMassOS, &b_dileptonMassOS);
   fChain->SetBranchAddress("dileptonDROS", &dileptonDROS, &b_dileptonDROS);
   fChain->SetBranchAddress("jets", &jets, &b_jets);
   fChain->SetBranchAddress("leptons", &leptons, &b_leptons);
}

Bool_t analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef analyzer_cxx
