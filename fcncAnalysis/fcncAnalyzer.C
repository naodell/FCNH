#define fcncAnalyzer_cxx

#include "fcncAnalyzer.h"

using namespace std;

/////////////////////////////
//Specify parameters here. //
/////////////////////////////


const bool      doPrintout  = false;
const bool      doGenPrint  = false;

// MVA switches
const bool      doMVACut    = true;
const bool      doMVATree   = false;
const bool      doLepTree   = false;

// Data-driven BG estimation switches
bool doQFlips = true;
bool doFakes  = true;


/////////////////
//Analysis cuts//
/////////////////


const float   jetPtCut[]        = {30., 15.};
const float   muPtCut[]         = {10., 3.};
const float   elePtCut[]        = {10., 7.};
const float   phoPtCut[]        = {10., 10.};
const float   leptonPtCut[]     = {20., 10.};
const float   metCut[]          = {30., 0.};
const float   htCut[]           = {13., 14.};
const float   bJetVeto          = 1e9;

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
bool BTagSortCondition(TCJet j1, TCJet j2) {return (j1.BDiscriminatorMap("CSV") > j2.BDiscriminatorMap("CSV"));} 

void fcncAnalyzer::Begin(TTree* tree) 
{
    // Job config
    TObjArray *args = (TObjArray*)fOption.Tokenize(" ");

    suffix      = (string)((TObjString*)args->At(0))->GetString();
    selection   = (string)((TObjString*)args->At(1))->GetString();
    period      = (string)((TObjString*)args->At(2))->GetString();

    cout << suffix << " " << selection << " " << period << endl;

    // Get trigger names from jobTree
    vector<string>* triggerNames = 0;
    TFile *inFile   = tree->GetCurrentFile();
    TTree *jobTree  = (TTree*)inFile->Get("ntupleProducer/jobTree");

    jobTree->SetBranchAddress("triggerNames", &triggerNames);
    jobTree->GetEntry();

    // Initialize utilities and selectors here //
    selector        = new Selector(muPtCut, elePtCut, jetPtCut, phoPtCut);
    weighter        = new WeightUtils(suffix, period, selection, isRealData);
    triggerSelector = new TriggerSelector(selection, period, *triggerNames, true);

    // Initialize histograms //
    TString option = GetOption();
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    histManager  = new HistManager();
    subdir = suffix;

    for (unsigned iCut = 0; iCut < N_CUTS; ++iCut) {
        string index = str(iCut+1);

        histoFile[iCut] = new TFile(("histos/fcncHistograms_cut" + index + ".root").c_str(), "RECREATE");
        histManager->AddFile(histoFile[iCut]);
        histManager->SetFileNumber(iCut);

        histoFile[iCut]->mkdir("TESTS", "TESTS");
        histoFile[iCut]->GetDirectory("TESTS")->mkdir(subdir.c_str());

        for (unsigned i = 0; i < N_CATEGORIES; ++i) { 

            histoFile[iCut]->mkdir(categoryNames[i].c_str(), categoryNames[i].c_str());
            histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(suffix.c_str(), suffix.c_str());

            if (doQFlips && (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "TEST")) 
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir("QFlips", "QFlips");
            else
                doQFlips = false;

            if (doFakes && (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "DATA_MUON" || suffix == "TEST")) 
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir("Fakes", "Fakes");
            else
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(("Fakes_"+suffix).c_str(), ("Fakes_"+suffix).c_str());
        }
    }

    // Initialize lepton MVA
    if (doLepTree) {
        histoFile[0]->cd();
        muTree = new TTree(("muTree_" + suffix).c_str(), "Tree for lepton MVA");
        muTree->Branch("sip3d", &sip3d, "sip3d/F");
        muTree->Branch("chPFIso", &chPFIso, "chPFIso/F");
        muTree->Branch("neuPFIso", &neuPFIso, "neuPFIso/F");
        muTree->Branch("drLepJet", &drLepJet, "drLepJet/F");
        muTree->Branch("ptRatioLepJet", &ptRatioLepJet, "ptRatioLepJet/F");
        muTree->Branch("btagLepJet", &btagLepJet, "btagLepJet/F");
        muTree->Branch("dxy", &dxy, "dxy/F");
        muTree->Branch("dz", &dz, "dz/F");


        eleTree = new TTree(("eleTree_" + suffix).c_str(), "Tree for lepton MVA");
        eleTree->Branch("sip3d", &sip3d, "sip3d/F");
        eleTree->Branch("chPFIso", &chPFIso, "chPFIso/F");
        eleTree->Branch("neuPFIso", &neuPFIso, "neuPFIso/F");
        eleTree->Branch("drLepJet", &drLepJet, "drLepJet/F");
        eleTree->Branch("ptRatioLepJet", &ptRatioLepJet, "ptRatioLepJet/F");
        eleTree->Branch("btagLepJet", &btagLepJet, "btagLepJet/F");
        eleTree->Branch("mva", &mva, "mva/F");
        eleTree->Branch("missHits", &missHits, "missHits/I");
    }

    // Initialize pass tree for MVA input //
    selJets     = new TClonesArray("TCPhysObject");
    selLeptons  = new TClonesArray("TCPhysObject");
    if (doMVATree) {
        histoFile[0]->cd();
        // make base tree for 3l and ss selections
        tree3l = new TTree(("tree3l_" + suffix).c_str(), "Tree for cut MVA");
        tree3l->Branch("evtWeight", &evtWeight, "evtWeight/F");
        tree3l->Branch("flavorCat", &flavorCat, "flavorCat/I");
        tree3l->Branch("chargeCat", &chargeCat, "chargeCat/I");
        tree3l->Branch("met", &MET, "met/F");
        tree3l->Branch("metPhi", &metPhi, "metPhi/F");
        tree3l->Branch("HT", &HT, "HT/F");
        tree3l->Branch("MT", &MT, "MT/F");
        tree3l->Branch("jetMult", &jetMult, "jetMult/I");
        tree3l->Branch("bJetMult", &bJetMult, "bJetMult/I");
        tree3l->Branch("trileptonMass", &trileptonMass, "trileptonMass/F");
        tree3l->Branch("dileptonMassOS", &dileptonMassOS, "dileptonMassOS/F");
        tree3l->Branch("dileptonDROS", &dileptonDROS, "dileptonDROS/F");
        tree3l->Branch("jets",&selJets, 6400, 0);
        tree3l->Branch("leptons",&selLeptons, 6400, 0);

        treeSS = new TTree(("treeSS_" + suffix).c_str(), "Tree for cut MVA");
        treeSS->Branch("evtWeight", &evtWeight, "evtWeight/F");
        treeSS->Branch("flavorCat", &flavorCat, "flavorCat/I");
        treeSS->Branch("chargeCat", &chargeCat, "chargeCat/I");
        treeSS->Branch("met", &MET, "met/F");
        treeSS->Branch("metPhi", &metPhi, "metPhi/F");
        treeSS->Branch("HT", &HT, "HT/F");
        treeSS->Branch("MT", &MT, "MT/F");
        treeSS->Branch("jetMult", &jetMult, "jetMult/I");
        treeSS->Branch("bJetMult", &bJetMult, "bJetMult/I");
        treeSS->Branch("dileptonMass", &dileptonMassOS, "dileptonMass/F");
        treeSS->Branch("dileptonDR", &dileptonDROS, "dileptonDR/F");
        treeSS->Branch("jets",&selJets, 6400, 0);
        treeSS->Branch("leptons",&selLeptons, 6400, 0);

        if (doQFlips && (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "TEST")) {
            treeQFlips = new TTree("treeSS_QFlips", "Tree for same-sign cut MVA");
            treeQFlips->Branch("evtWeight", &evtWeight, "evtWeight/F");
            treeQFlips->Branch("flavorCat", &flavorCat, "flavorCat/I");
            treeQFlips->Branch("chargeCat", &chargeCat, "chargeCat/I");
            treeQFlips->Branch("met", &MET, "met/F");
            treeQFlips->Branch("metPhi", &metPhi, "metPhi/F");
            treeQFlips->Branch("HT", &HT, "HT/F");
            treeQFlips->Branch("MT", &MT, "MT/F");
            treeQFlips->Branch("jetMult", &jetMult, "jetMult/I");
            treeQFlips->Branch("bJetMult", &bJetMult, "bJetMult/I");
            treeQFlips->Branch("dileptonMass", &dileptonMassOS, "dileptonMass/F");
            treeQFlips->Branch("dileptonDR", &dileptonDROS, "dileptonDR/F");
        }

        if (doFakes && (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "DATA_MUON" || suffix == "TEST")) {
            treeFakes3l = new TTree("tree3l_Fakes", "Tree for 3l cut MVA");
            treeFakes3l->Branch("evtWeight", &evtWeight, "evtWeight/F");
            treeFakes3l->Branch("flavorCat", &flavorCat, "flavorCat/I");
            treeFakes3l->Branch("chargeCat", &chargeCat, "chargeCat/I");
            treeFakes3l->Branch("met", &MET, "met/F");
            treeFakes3l->Branch("metPhi", &metPhi, "metPhi/F");
            treeFakes3l->Branch("HT", &HT, "HT/F");
            treeFakes3l->Branch("MT", &MT, "MT/F");
            treeFakes3l->Branch("jetMult", &jetMult, "jetMult/I");
            treeFakes3l->Branch("bJetMult", &bJetMult, "bJetMult/I");
            treeFakes3l->Branch("trileptonMass", &trileptonMass, "trileptonMass/F");
            treeFakes3l->Branch("dileptonMassOS", &dileptonMassOS, "dileptonMassOS/F");
            treeFakes3l->Branch("dileptonDROS", &dileptonDROS, "dileptonDROS/F");

            treeFakesSS = new TTree("treeSS_Fakes", "Tree for same-sign cut MVA");
            treeFakesSS->Branch("evtWeight", &evtWeight, "evtWeight/F");
            treeFakesSS->Branch("flavorCat", &flavorCat, "flavorCat/I");
            treeFakesSS->Branch("chargeCat", &chargeCat, "chargeCat/I");
            treeFakesSS->Branch("met", &MET, "met/F");
            treeFakesSS->Branch("metPhi", &metPhi, "metPhi/F");
            treeFakesSS->Branch("HT", &HT, "HT/F");
            treeFakesSS->Branch("MT", &MT, "MT/F");
            treeFakesSS->Branch("jetMult", &jetMult, "jetMult/I");
            treeFakesSS->Branch("bJetMult", &bJetMult, "bJetMult/I");
            treeFakesSS->Branch("dileptonMass", &dileptonMassOS, "dileptonMass/F");
            treeFakesSS->Branch("dileptonDR", &dileptonDROS, "dileptonDR/F");
        }
    }

    if (doMVACut) {
<<<<<<< HEAD
        //string mva3lCats[4] = {"eee", "eemu", "emumu", "mumumu"};
        string mva3lCats[4] = {"inclusive"};
<<<<<<< HEAD
        for (unsigned i = 0; i < 1; ++i) {
=======
        string mva3lCats[4] = {"eee", "eemu", "emumu", "mumumu"};
        for (unsigned i = 0; i < 4; ++i) {
>>>>>>> parent of 72355eb... Minor changes
=======
        for (unsigned i = 0; i < 4; ++i) {
>>>>>>> parent of 17d9010... Adding new BDT weight files
            TMVA::Reader* mvaReader = new TMVA::Reader("!Color:!Silent");

            mvaReader->AddVariable("met", &MET);
            mvaReader->AddVariable("HT", &HT);
            mvaReader->AddVariable("MT", &MT);
            mvaReader->AddVariable("bJetMult", &f_bJetMult);
            mvaReader->AddVariable("dileptonMassOS", &dileptonMassOS);
            //mvaReader->AddVariable("dileptonDROS", &dileptonDROS);
            //mvaReader->AddVariable("trileptonMass", &trileptonMass);

            mvaReader->AddSpectator("flavorCat", &f_flavorCat);
            mvaReader->AddSpectator("chargeCat", &f_chargeCat);
            mvaReader->AddSpectator("jetMult", &f_jetMult);
            mvaReader->AddSpectator("evtWeight", &evtWeight);

            mvaReader->BookMVA("test", ("../data/weights/20131210_001709/TMVAClassification_3l_" + mva3lCats[i] + "_BDT.weights.xml").c_str());
            mva3lReader.push_back(mvaReader);
        }

<<<<<<< HEAD
        string mvaSSCats[3] = {"inclusive"};
<<<<<<< HEAD
        for (unsigned i = 0; i < 1; ++i) {
=======
        string mvaSSCats[3] = {"ee", "emu", "mumu"};
        for (unsigned i = 0; i < 3; ++i) {
>>>>>>> parent of 72355eb... Minor changes
=======
        for (unsigned i = 0; i < 3; ++i) {
>>>>>>> parent of 17d9010... Adding new BDT weight files
            TMVA::Reader* mvaReader = new TMVA::Reader("!Color:!Silent");

            mvaReader->AddVariable("met", &MET);
            mvaReader->AddVariable("HT", &HT);
            mvaReader->AddVariable("MT", &MT);
            mvaReader->AddVariable("bJetMult", &f_bJetMult);
            mvaReader->AddVariable("dileptonMass", &dileptonMassOS);
            //mvaReader->AddVariable("dileptonDR", &dileptonDROS);

            mvaReader->AddSpectator("flavorCat", &f_flavorCat);
            mvaReader->AddSpectator("chargeCat", &f_chargeCat);
            mvaReader->AddSpectator("jetMult", &f_jetMult);
            mvaReader->AddSpectator("evtWeight", &evtWeight);

            mvaReader->BookMVA("test", ("../data/weights/20131210_001709/TMVAClassification_SS_" + mvaSSCats[i] + "_BDT.weights.xml").c_str());
            mvaSSReader.push_back(mvaReader);
        }
    }

    // initialize some global variables
    for (unsigned i = 0; i < 16; ++i) {
        eventCount[i] = 0;
        eventCountWeighted[i] = 0;
    }

    evtWeight   = 1.;
}

bool fcncAnalyzer::Process(Long64_t entry)
{  
    GetEntry(entry);

    // Reset event information
    selector->PurgeObjects();
    evtCategory.reset();
    evtWeight = 1.;
    subdir = suffix;

    if (eventCount[1] == 0) {
        weighter->SetDataBit(isRealData);
        triggerSelector->SetDataBit(isRealData);
        selector->SetDataBit(isRealData);
    }

    if (eventCount[1] == 0)
        cout << "Starting analysis..." << endl;
    else if (eventCount[1] % (int)1e4 == 0) 
        cout << eventCount[5] << " events passed of " << eventCount[1] << " checked!" << endl;

    SetYields(1);

    //////////////////
    //Trigger status//
    //////////////////


    for(int i = 0; i < 64; ++i) {
        if (triggerStatus & ULong64_t(0x1) << i) {
            histManager->Fill1DHist(i+1, "h1_TriggerStatus", "Triggers", 64, 0.5, 64.5);
        }
    } 

    bool triggerPass = false;
    triggerPass = triggerSelector->SelectTrigger(triggerStatus, hltPrescale);

    // Double electron workaround.  Gets rid of hopelessly prescaled events fo July 20-26, 2011
    //if (selection == "electron" && (runNumber > 171200 && runNumber < 171600)) return kTRUE;

    if (!triggerPass) 
        return kTRUE;
    else
        SetYields(2);

    vstring passNames = triggerSelector->GetPassNames();

    //for (unsigned i = 0; i < passNames.size(); ++i) cout << passNames[i] << endl;
    //cout << endl;

    if (passNames.size() == 0) passNames.push_back("NULL");


    ////////////////////////////
    //Check the event vertices//
    ////////////////////////////


    if (!isRealData) {
        histManager->SetFileNumber(0);
        histManager->SetDirectory("inclusive/" + subdir);
        histManager->SetWeight(1);

        histManager->Fill1DHist(nPUVertices,
                "h1_SimVertexMult", "Multiplicity of simulated vertices", 100, 0.5, 100.5);
        histManager->Fill1DHist(nPUVerticesTrue,
                "h1_SimVertexMultTrue", "True simulated PU", 500, 0., 100.);
    }

    selector->PVSelector(primaryVtx);

    if (selector->GetSelectedPVs().size() < 1) 
        return kTRUE;
    else
        SetYields(3);

    TVector3 selectedVtx = *selector->GetSelectedPVs()[0];


    //////////////////
    // Data quality //
    //////////////////


    if (
            isRealData 
            && ( NoiseFilters_isCSCTightHalo || NoiseFilters_isNoiseHcalHBHE || NoiseFilters_isScraping )
       ) return kTRUE;
    else
        SetYields(4);

    /////////////////////////////
    // Get gen level particles //
    /////////////////////////////

    vector<TCGenParticle> higgs, dubyas, Zeds, gElectrons, gMuons, gTaus, gLeptons;
    vector<TCGenJet> gJets;

    if (!isRealData) {
        selector->GenParticleSelector(genParticles, 25, 3, "Higgs");
        selector->GenParticleSelector(genParticles, 24, 3, "Dubya");
        selector->GenParticleSelector(genParticles, 23, 3, "Zeds");
        selector->GenParticleSelector(genParticles, 11, 3, "electrons");
        selector->GenParticleSelector(genParticles, 13, 3, "muons");
        selector->GenParticleSelector(genParticles, 15, 3, "taus");

        higgs       = selector->GetSelectedGenParticles("Higgs");
        dubyas      = selector->GetSelectedGenParticles("Dubya");
        Zeds        = selector->GetSelectedGenParticles("Zeds");
        gElectrons  = selector->GetSelectedGenParticles("electrons");
        gMuons      = selector->GetSelectedGenParticles("muons");
        gTaus       = selector->GetSelectedGenParticles("taus");

        gLeptons.insert(gLeptons.end(), gElectrons.begin(), gElectrons.end());
        gLeptons.insert(gLeptons.end(), gMuons.begin(), gMuons.end());
        gLeptons.insert(gLeptons.end(), gTaus.begin(), gTaus.end());

        selector->GenJetSelector(genJets);
        gJets = selector->GetSelectedGenJets();
    }

    if (doGenPrint) {
        // Higgs
        for (unsigned i = 0; i < higgs.size(); ++i)
            if (higgs[i].GetStatus() == 3) 
                cout << higgs[i].GetStatus() << ", " << higgs[i].M() << ", " << higgs[i].Mother() << ", higgs" << endl;

        // Vector bosons
        for (unsigned i = 0; i < dubyas.size(); ++i)
            if (dubyas[i].GetStatus() == 3) 
                cout << "\t status = " << dubyas[i].GetStatus() <<  ", mass = " << dubyas[i].M() << ", pt = " << dubyas[i].Pt() << ", mother = " << dubyas[i].Mother() << ", dubyas" << endl;
        for (unsigned i = 0; i < Zeds.size(); ++i)
            if (Zeds[i].GetStatus() == 3) 
                cout << "\t status = " << Zeds[i].GetStatus() << ", mass = " << Zeds[i].M() << ", pt = " << Zeds[i].Pt() << ", mother = " << Zeds[i].Mother() << ", Zeds" << endl;

        // leptons
        if ((gElectrons.size() + gMuons.size() + gTaus.size()) > 0) {
            for (unsigned i = 0; i < gElectrons.size(); ++i)
                cout << "\t\t" << gElectrons[i].GetStatus() << ", " << gElectrons[i].Pt() << ", " << gElectrons[i].Eta() << ", " << gElectrons[i].Mother() << ", " << gElectrons[i].Grandmother() << ", electrons" << endl;
            for (unsigned i = 0; i < gMuons.size(); ++i)
                cout << "\t\t" << gMuons[i].GetStatus() << ", " << gMuons[i].Pt() << ", " << gMuons[i].Eta() << ", " << gMuons[i].Mother() << ", " << gMuons[i].Grandmother() << ", muons" << endl;
            for (unsigned i = 0; i < gTaus.size(); ++i)
                cout << "\t\t" << gTaus[i].GetStatus() << ", " << gTaus[i].Pt() << ", " << gTaus[i].Eta() << ", " << gTaus[i].Mother() << ", " << gTaus[i].Grandmother() << ", taus" << endl;
        }

        cout << "\n" << endl;
    }


    //////////////////////
    // object selection //
    //////////////////////

    // Set input variables
    selector->SetRho(rho25Factor);

    // Run selectors
    // N.B. The order that the selectors are run matters.  This
    // is because overlap between leptons and jets is checked.
    // If there is overlap between an electron and muon, it is 
    // assumed that the muon is faking the electron.  Jets that 
    // ovelap with leptons are assumed to be leptons, but a collection
    // of the overlap jets are saved to discriminate against b-jets
    // that fake muons.

    selector->MuonSelector(recoMuons);
    selector->ElectronSelector(recoElectrons);
    selector->JetSelector(recoJets);

    // Get analysis leptons
    vector<TCMuon>      muons       = selector->GetSelectedMuons("tight");
    vector<TCElectron>  electrons   = selector->GetSelectedElectrons("tight");

    vObj leptons;
    leptons.insert(leptons.end(), muons.begin(), muons.end());
    leptons.insert(leptons.end(), electrons.begin(), electrons.end());

    // Get loose leptons
    vObj extraLeptons;
    vector<TCElectron>  looseElectrons  = selector->GetSelectedElectrons("loose");
    vector<TCMuon>      looseMuons      = selector->GetSelectedMuons("loose");
    extraLeptons.insert(extraLeptons.end(), looseElectrons.begin(), looseElectrons.end());
    extraLeptons.insert(extraLeptons.end(), looseMuons.begin(), looseMuons.end());

    // Get overlap electrons
    vector<TCElectron> olElectrons = selector->GetSelectedElectrons("tight_overlap");

    // Get fakeable leptons
    vObj fakeables;
    vector<TCElectron>  fakeableElectrons  = selector->GetSelectedElectrons("fakeable");
    vector<TCMuon>      fakeableMuons      = selector->GetSelectedMuons("fakeable");
    fakeables.insert(fakeables.end(), fakeableElectrons.begin(), fakeableElectrons.end());
    fakeables.insert(fakeables.end(), fakeableMuons.begin(), fakeableMuons.end());

    // Get jets
    vector<TCJet> allJets;
    vector<TCJet> jets      = selector->GetSelectedJets("tight");
    vector<TCJet> bJetsM    = selector->GetSelectedJets("bJetsMedium");
    vector<TCJet> bJetsL    = selector->GetSelectedJets("bJetsLoose");
    vector<TCJet> fwdJets   = selector->GetSelectedJets("forward");
    vector<TCJet> muJets    = selector->GetSelectedJets("muJets");
    vector<TCJet> eleJets   = selector->GetSelectedJets("eleJets");

    jets.insert(jets.end(), fwdJets.begin(), fwdJets.end());
    //jets.insert(jets.end(), bJetsL.begin(), bJetsL.end());

    allJets.insert(allJets.end(), jets.begin(), jets.end());
    allJets.insert(allJets.end(), bJetsM.begin(), bJetsM.end());
    allJets.insert(allJets.end(), muJets.begin(), muJets.end());
    allJets.insert(allJets.end(), eleJets.begin(), eleJets.end());

    // Order collections by pt
    sort(extraLeptons.begin(), extraLeptons.end(), P4SortCondition);
    sort(jets.begin(), jets.end(), P4SortCondition);
    sort(allJets.begin(), allJets.end(), P4SortCondition);
    sort(bJetsM.begin(), bJetsM.end(), BTagSortCondition);
    sort(bJetsL.begin(), bJetsL.end(), BTagSortCondition);
    sort(leptons.begin(), leptons.end(), P4SortCondition);

    // Fill lepton mva tree
    if (doLepTree)
        FillLepMVA(selector->GetSelectedMuons("premva"), selector->GetSelectedElectrons("premva"), allJets, selectedVtx);

    histManager->SetFileNumber(0);
    histManager->SetDirectory("inclusive/" + subdir);
    histManager->SetWeight(1);

    //!!!!!!!!!!!!!!!!!!!!!!!!//
    //                        //
    //  Overlap Jets/Leptons  //
    //                        //
    //!!!!!!!!!!!!!!!!!!!!!!!!//

    histManager->Fill1DHist(muJets.size() + eleJets.size(), 
            "h1_OverlapJetMult", "(e/#mu)-jet multiplicity;N_{jets};Entries / bin", 5, -0.5, 4.5);

    if (muons.size() >= 1)
        histManager->Fill1DHist(olElectrons.size(), "h1_OverlapEleMu", ";(e/#mu) multiplicity;Entries", 4, -0.5, 3.5);

    if (olElectrons.size() > 0) return kTRUE;


    //!!!!!!!!!!!!!!!!!!!!!!!!//
    //                        //
    //  Event categorization  //
    //  and weighting....     //
    //                        //
    //!!!!!!!!!!!!!!!!!!!!!!!!//


    histManager->Fill1DHist(leptons.size(), "h1_LeptonMult", "lepton multiplicity;N_{leptons};Events / bin", 6, -0.5, 5.5);

    if (!isRealData) {
        GenPlots(gLeptons, leptons);
    }


    if (leptons.size() == 1) {

        //!!! Single leptons just for fakes !!!//
        if (leptons[0].Pt() < leptonPtCut[0]) 
            return kTRUE;

    } else if (leptons.size() == 2) {

        // Electron charge misid control region needs a lower lead lepton pt cut //
        if (
                leptons[0].Pt() > 10. && leptons[1].Pt() > 10.  // pt cut
                && leptons[0].Type() == "electron" && leptons[1].Type() == "electron" // get electrons
                && (fabs((leptons[0] + leptons[1]).M() - 91.2) < 15) // Z mass window
           ) 
            MakeQMisIDPlots(leptons);

        //!!! Dilepton selection !!!//
        if ( leptons[0].Pt() < leptonPtCut[0] || leptons[1].Pt() < leptonPtCut[1]) 
            return kTRUE;

    } else if (leptons.size() == 3) {

        //!!! Trilepton selection !!!//
        if ( 
                leptons[0].Pt() < leptonPtCut[0] 
                || leptons[1].Pt() < leptonPtCut[1] 
                || leptons[2].Pt() < leptonPtCut[1]
                || fabs(leptons[0].Charge() + leptons[1].Charge() + leptons[2].Charge()) != 1
           ) 
            return kTRUE;

    } else if (leptons.size() == 4) {

        //!!! Tetralepton selection !!!//
        if (leptons[0].Pt() < leptonPtCut[0] || leptons[1].Pt() < leptonPtCut[1]) 
            return kTRUE;

    } else if (leptons.size() > 4)
        return kTRUE;

    if (leptons.size() > 1) { // Only do signal extraction if there are at least two leptons

        //!! low mass resonance rejection !!//
        bool lowMassOS  = false;
        bool isCosmics  = false;

        for (unsigned i = 1; i < leptons.size(); ++i) {
            for (unsigned j = 0; j < i; ++j) {
                if ((leptons[i] + leptons[j]).M() < 12) 
                    lowMassOS = true;

                if (isRealData && leptons[i].Type() == "muon" && leptons[j].Type() == "muon")
                    if (CosmicMuonFilter(leptons[i], leptons[j]))
                        isCosmics = true;
            }
        }

        if (lowMassOS || isCosmics) return kTRUE;


        //!!!!!!!!!!!!!!!!!!!!!!!!!!//
        // End of preselection cuts //
        //!!!!!!!!!!!!!!!!!!!!!!!!!!//


        weighter->SetObjects(leptons, jets, nPUVerticesTrue, passNames[0]);
        evtWeight *= weighter->GetTotalWeight();
        histManager->SetWeight(evtWeight);

        SetEventCategory(leptons);
        SetEventVariables(leptons, jets, bJetsM, *recoMET); 

        // Pre-selection mc-truth plots
        if (!isRealData) {
            unsigned cat = 1 + ((evtCategory.to_ulong() >> 2) & 0x3);
            histManager->SetDirectory(categoryNames[cat] + "/" + subdir);
            GenPlots(gLeptons, leptons);

            cat = GetHistCategory(2) - 10;
            histManager->SetDirectory(categoryNames[cat] + "/" + subdir);
            GenPlots(gLeptons, leptons);
        }

        AnalysisSelection(leptons, jets, bJetsM, bJetsL, selectedVtx, suffix);

        if (doQFlips) {

            // Electron charge flip background estimation: applies weights based on
            // charge flip probability in Z->ee events.
            if (
                    leptons.size() == 2 
                    && leptons[0].Charge() != leptons[1].Charge()
                    && (leptons[0].Type() == "electron" || leptons[1].Type() == "electron")
               ) { 

                Float_t qFlipWeight = weighter->GetQFlipWeight();

                evtWeight *= qFlipWeight;
                vObj flipLeptons = leptons;

                if (runNumber%2 == 0)
                    flipLeptons[0].SetCharge(flipLeptons[1].Charge());
                else
                    flipLeptons[1].SetCharge(flipLeptons[0].Charge());

                SetEventCategory(flipLeptons);
                SetEventVariables(flipLeptons, jets, bJetsM, *recoMET); 
                AnalysisSelection(flipLeptons, jets, bJetsM, bJetsL, selectedVtx, "QFlips");

                evtWeight /= qFlipWeight; // Remove charge flip weight
            }
        }
    } //Allow for the case of one or no leptons for fakes

    if (doFakes) {

        // Do clean-up that's only done for analysis leptons 
        // remove fakeable from jet collection and make sure fakeable
        // object does not overlap with passing object
        bool leptonMatched      = false;
        bool lowMassResonance   = false;
        vObj matchedFakeables, unmatchedFakeables;

        for (unsigned i = 0; i < fakeables.size(); ++i) {

            for (unsigned j = 0; j < leptons.size(); ++j) {
                if (fakeables[i].DeltaR(leptons[j]) < 0.1) {
                    fakeables.erase(fakeables.begin() + i);
                    --i;
                    leptonMatched = true;
                    break;
                }

                if ((fakeables[i] + leptons[j]).M() < 10.) {
                    lowMassResonance = true;
                    break;
                }
            }

            if (leptonMatched || lowMassResonance) continue;

            bool fakeMatched        = false;
            for (unsigned j = 0; j < fakeables.size(); ++j) {
                if (i == j) continue;

                if (fakeables[i].DeltaR(fakeables[j]) < 0.1) {
                    fakeMatched = true;
                    matchedFakeables.push_back(fakeables[i]);

                } else if ((fakeables[i] + fakeables[j]).M() < 10.) {
                    lowMassResonance = true;
                    break;
                }
            }

            if (!fakeMatched) 
                unmatchedFakeables.push_back(fakeables[i]);
        }

        vector<TCJet> fJets     = selector->GetSelectedJets("tight_NoFakes");
        vector<TCJet> fBJetsM   = selector->GetSelectedJets("bJetsMedium_NoFakes");
        vector<TCJet> fBJetsL   = selector->GetSelectedJets("bJetsLoose_NoFakes");

        //fJets.insert(fJets.end(), fBJetsL.begin(), fBJetsL.end());

        sort(fJets.begin(), fJets.end(), P4SortCondition);
        sort(fBJetsM.begin(), fBJetsM.end(), BTagSortCondition);
        sort(fBJetsL.begin(), fBJetsL.end(), BTagSortCondition);

        if (!lowMassResonance) {

            histManager->SetFileNumber(0);
            histManager->SetDirectory("inclusive/" + suffix);
            histManager->Fill1DHist(matchedFakeables.size()/2., "h1_fakeableOverlapMult", ";e #mu overlap pairs;Entries", 3, -0.5, 2.5);

            if (matchedFakeables.size() == 0) {
                GetFakeBG(leptons, fakeables, fJets, fBJetsM, fBJetsL, selectedVtx);
            } /*else if (matchedFakeables.size() == 2) {
                evtWeight = evtWeight/(matchedFakeables.size());
                for (unsigned i = 0; i < matchedFakeables.size(); ++i) {
            //cout << matchedFakeables[i].Pt() << ", " << matchedFakeables[i].Eta() << endl;
            vObj tmpFakeables = unmatchedFakeables;
            tmpFakeables.push_back(matchedFakeables[i]);
            GetFakeBG(leptons, tmpFakeables, fJets, fBJetsM, fBJetsL, selectedVtx);
            }
            //cout << endl;
            evtWeight *= matchedFakeables.size();
            }*/
        }
    }

    return kTRUE;
}

void fcncAnalyzer::Terminate()
{
    cout<<"\nRunning over "<<suffix<<" dataset with "<<selection<<" selection for "<<period<<" data."<<"\n"<<endl;
    cout<<"| CUT DESCRIPTION                    |\t" << "\t|"<<endl;
    cout<<"| Initial number of events:          |\t" << eventCountWeighted[0] << "\t|"<<endl;
    cout<<"| Number of events ntuplized:        |\t" << eventCount[1]  << "\t|\t" << eventCountWeighted[1] << "\t|"<<endl;
    cout<<"| Pass HLT selection:                |\t" << eventCount[2]  << "\t|\t" << eventCountWeighted[2] << "\t|"<<endl;
    cout<<"| Good PV:                           |\t" << eventCount[3]  << "\t|\t" << eventCountWeighted[3] << "\t|"<<endl;
    cout<<"| Data quality bits:                 |\t" << eventCount[4]  << "\t|\t" << eventCountWeighted[4] << "\t|"<<endl;
    cout<<"| Lepton selection:                  |\t" << eventCount[5]  << "\t|\t" << eventCountWeighted[5] << "\t|"<<endl;

    // FCNH selection //
    cout<<"| Z veto:                            |\t" << eventCount[6]  << "\t|\t" << eventCountWeighted[6] << "\t|"<<endl;
    cout<<"| At least two jets:                 |\t" << eventCount[7]  << "\t|\t" << eventCountWeighted[7] << "\t|"<<endl;
    cout<<"| MET cut:                           |\t" << eventCount[8]  << "\t|\t" << eventCountWeighted[8] << "\t|"<<endl;
    cout<<"| HT cut:                            |\t" << eventCount[9]  << "\t|\t" << eventCountWeighted[9] << "\t|"<<endl;
    cout<<"| BDT:                               |\t" << eventCount[15]  << "\t|\t" << eventCountWeighted[15] << "\t|"<<endl;


    // Control regions //
    cout<<"\nControl region event yields."<<"\n"<<endl;
    cout<<"| WZ:                                |\t" << eventCount[10]  << "\t|\t" << eventCountWeighted[10] << "\t|"<<endl;
    cout<<"| ttbar:                             |\t" << eventCount[11]  << "\t|\t" << eventCountWeighted[11] << "\t|"<<endl;
    cout<<"| ttZ:                               |\t" << eventCount[12]  << "\t|\t" << eventCountWeighted[12] << "\t|"<<endl;
    cout<<"| low eta:                           |\t" << eventCount[13]  << "\t|\t" << eventCountWeighted[13] << "\t|"<<endl;
    cout<<"| ZZ:                                |\t" << eventCount[14]  << "\t|\t" << eventCountWeighted[14] << "\t|"<<endl;


    //for (int i = 0; i < 8; ++i) fout[i].close();

    // Set alphanumeric bins for charge and flavor histograms
    string chLabels[12] = {"--", "-+", "+-", "++", "---", "--+", "-+-", "-++", "+--", "+-+", "++-", "+++"};
    string flLabels[12] = {"ee", "e#mu", "#mue", "#mu#mu", "eee", "ee#mu", "e#mue", "e#mu#mu", "#muee", "#mue#mu", "#mu#mue", "#mu#mu#mu"};

    for (unsigned i = 0; i < N_CATEGORIES; ++i) {
        histManager->SetDirectory(categoryNames[i] + "/" + suffix);
        histManager->SetBinNames1D("h1_LeptonCharge", chLabels, 12);
        histManager->SetBinNames1D("h1_LeptonFlavor", flLabels, 12);
        histManager->SetBinNames2D("h2_LepChargeVsFlavor", chLabels, 12, flLabels, 12);
    }

    // Clean-up and write to output files
    selector->Delete();
    triggerSelector->Delete();
    weighter->Delete();
    histManager->Delete();

    // Close ntuple file
    for (unsigned i = 0; i < N_CUTS; ++i) {
        histoFile[i]->Write();
        histoFile[i]->Close();  
    }
}

bool fcncAnalyzer::AnalysisSelection(vObj leptons, vector<TCJet> jets, vector<TCJet> bJetsM, vector<TCJet> bJetsL, TVector3 PV, string histDir)
{
    subdir = histDir;

    // ZZ control region //
    if (leptons.size() == 4) {
        if ( bJetsM.size() == 0) {
            Make4lPlots(leptons, *recoMET);
            SetYields(14);
        }
        return kTRUE; 
    }

    //!!!!!!!!!!!!!!!!!//
    // Control regions //
    //!!!!!!!!!!!!!!!!!//


    // WZ control region //
    if (
            zTagged 
            && leptons.size() == 3 
            && bJetsM.size() == 0
            && jets.size() > 1
            && METLD > 0.3
       ) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 6);
        SetYields(10);
    }

    // ttbar control region //
    if (
            leptons.size() == 2 
            && (bJetsM.size() == 1 || bJetsL.size() == 2)
            && leptons[0].Type() != leptons[1].Type()
            && leptons[0].Charge() != leptons[1].Charge()
            && MET > 30
       ) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 7);
        SetYields(11);
    }

    // ttZ control region //
    if (
            leptons.size() >= 3 
            && zTagged
            && (bJetsM.size() >= 1 || bJetsL.size() >= 2)
            && METLD > 0.2 
       ) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 8);
        SetYields(12);
    }

    // Z+fake control region //
    if (
            zTagged 
            && leptons.size() == 3 
            && METLD < 0.3
       ) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 9);
        SetYields(13);
    }


    //!!!!!!!!!!!!!//
    // Do MVA Trees//
    //!!!!!!!!!!!!!//


    if (leptons.size() == 3) {
        // Fill MVA trees //
        SetVarsMVA(leptons, bJetsM, jets);
        if (doMVATree) {
            if (histDir.substr(0, 5) == "Fakes" && isRealData && doFakes) {
                //cout << jetMult << endl;
                treeFakes3l->Fill();
            } else {
                tree3l->Fill();
            }
        }
    } else if (leptons.size() == 2 && leptons[0].Charge() == leptons[1].Charge()) {
        // Fill MVA trees //
        SetVarsMVA(leptons, bJetsM, jets);
        if (doMVATree) {
            if (histDir.substr(0, 5) == "Fakes" && isRealData && doFakes)
                treeFakesSS->Fill();
            else if (histDir == "QFlips" && isRealData && doQFlips)
                treeQFlips->Fill();
            else
                treeSS->Fill();
        }
    }


    //!!!!!!!!!!!!!!!!!!!!!!//
    //                      //
    //  Analysis selection  //
    //  Cut n' Count!!!     //
    //                      //
    //!!!!!!!!!!!!!!!!!!!!!!//


    // Preselection Plots
    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 0);
    SetYields(5);


    //!! Z-veto !!//
    if (
            leptons.size() == 2 
            && leptons[0].Charge() == leptons[1].Charge()
            && leptons[0].Type() == leptons[1].Type()
            && fabs((leptons[0] + leptons[1]).M() - 91.2) < 15. 
       ) 
        return true;
    else if (leptons.size() == 3 && (zTagged || (dileptonMassOS > 40 && fabs(trileptonMass - 91.2) < 7.5))) 
        return true;

    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 1);
    SetYields(6);


    //!! Require at least one b-jet !!//
    if (bJetsM.size() + jets.size() <= 1) return true;
    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 2);
    SetYields(7);

    if (jetMult > 1 && HT < 60)
        cout << HT << endl;


    //!! Do mva selection !!//
    if (doMVACut) {
        float mvaValue = -99.;
        float mvaCut = -99.;
        if (leptons.size() == 3) {
            if (flavorCat == 5) {
                mvaValue = mva3lReader[0]->EvaluateMVA("test");
                mvaCut   = -0.1578;
            } else if (flavorCat == 6 || flavorCat == 7 || flavorCat == 9) {
                mvaValue = mva3lReader[1]->EvaluateMVA("test");
                mvaCut   = -0.0289;
            } else if (flavorCat == 8 || flavorCat == 10 || flavorCat == 11) {
                mvaValue = mva3lReader[2]->EvaluateMVA("test");
                mvaCut   = -0.0854;
            } else if (flavorCat == 12) {
                mvaValue = mva3lReader[3]->EvaluateMVA("test");
                mvaCut   = -0.1532;
            }

            histManager->SetFileNumber(1);
            histManager->SetDirectory("3l_inclusive/" + subdir);
            histManager->Fill1DHist(mvaValue, "h1_BDT", "BDT value;Entries / bin;BDT", 36, -1., 0.2);

        } else if (leptons.size() == 2 && leptons[0].Charge() == leptons[1].Charge()) {
            if (flavorCat == 1) {
                mvaValue = mvaSSReader[0]->EvaluateMVA("test");
                mvaCut   = 0.1390;
            } else if (flavorCat == 2 || flavorCat == 3) {
                mvaValue = mvaSSReader[1]->EvaluateMVA("test");
                mvaCut   = -0.2249;
            } else if (flavorCat == 4) {
                mvaValue = mvaSSReader[2]->EvaluateMVA("test");
                mvaCut   = -0.499;
            }

            histManager->SetFileNumber(1);
            histManager->SetDirectory("ss_inclusive/" + subdir);
            histManager->Fill1DHist(mvaValue, "h1_BDT", "BDT value;Entries / bin;BDT", 36, -1., 0.2);
        }

        if (mvaValue > mvaCut) {
            MakePlots(leptons, jets, bJetsM, *recoMET, PV, 5);
            SetYields(15);
        }
    }

    //!! MET cut !!//
    if (leptons.size() == 2){
        if (leptons[0].Charge() == leptons[1].Charge()) 
            if (recoMET->Mod() < metCut[0])
                return true;
    } else if (leptons.size() == 3) {
        if (recoMET->Mod() < metCut[1]) 
            return true;
    }
    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 3);
    SetYields(8);

    //!! HT cut !!//
    if (leptons.size() == 2){
        if (leptons[0].Charge() == leptons[1].Charge()) 
            if (sqrt(HT) < htCut[0])
                return true;
    } else if (leptons.size() == 3) {
        if (sqrt(HT) < htCut[1]) 
            return true;
    }
    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 4);
    SetYields(9);

    return true;
}

void fcncAnalyzer::GetFakeBG(vObj leptons, vObj fakeables, vector<TCJet> jets, vector<TCJet> bJetsM, vector<TCJet> bJetsL, TVector3 PV)
{
    // Do application of fake rates here.  This is done for the case of
    // ppf, pff, pf and ff events.  
    if (
<<<<<<< HEAD
            (fakeables.size() == 1 && (leptons.size() == 2 || leptons.size() == 1))
            || (fakeables.size() == 2 && (leptons.size() == 1 || leptons.size() == 0)) 
       ) {

        string fakeCat = GetFakeCategory(fakeables);
        cout << fakeCat << endl;

        Float_t fakeWeight1 = 1.;
        Float_t fakeWeight2 = 1.;
        if (fakeables.size() == 2) {
            fakeWeight1 = weighter->GetFakeWeight(fakeables[0], "QCD2l");
            fakeWeight2 = weighter->GetFakeWeight(fakeables[1], "QCD2l");
            //cout << fakeables[0].Type() << ", " << fakeables[1].Type() << ", " << fakeWeight1 << fakeWeight2 << endl;
        } else if (fakeables.size() == 1) 
            fakeWeight1 = weighter->GetFakeWeight(fakeables[0], "QCD2l");
            //cout << fakeables.size() << ", " << fakeables[0].Type() << ", " << fakeWeight1 << endl;

        Float_t fakeWeight = fakeWeight1*fakeWeight2;
=======
            (fakeables.size() == 1 && ((leptons.size() == 2 && leptons[0].Charge() != leptons[1].Charge()) || leptons.size() == 1))
            || (fakeables.size() == 2 && leptons.size() < 2) 
            //|| (fakeables.size() == 3 && leptons.size() == 0)
       ) {

        Float_t fakeWeight = weighter->GetFakeWeight(fakeables, "QCD2l");
        //cout << fakeWeight << endl;
>>>>>>> parent of 17d9010... Adding new BDT weight files

        if (fakeWeight <= 0)
            return;
        else
            evtWeight *= fakeWeight;

<<<<<<< HEAD
        // Fake rate fudge 
        if (fakeables.size() == 1 && leptons.size() == 2) 
            if (fakeables[0].Type() == "electron" && leptons[0].Type() == leptons[1].Type()) 
                fakeWeight *= 0.4;


=======
>>>>>>> parent of 17d9010... Adding new BDT weight files
        vObj leptonsPlusFakes = leptons;
        leptonsPlusFakes.insert(leptonsPlusFakes.end(), fakeables.begin(), fakeables.end());
        sort(leptonsPlusFakes.begin(), leptonsPlusFakes.end(), P4SortCondition);

        SetEventCategory(leptonsPlusFakes);
        SetEventVariables(leptonsPlusFakes, jets, bJetsM, *recoMET); 

<<<<<<< HEAD
        for (unsigned i = 0; i < leptonsPlusFakes.size(); ++i) {
            cout << leptonsPlusFakes[i].IsFake() << ", " << leptonsPlusFakes[i].Type() << "\t";
        }
        cout << endl;
=======
        if (fakeables.size() >= 1) {
            unsigned flCategory = GetHistCategory(2) - 10;
            histManager->SetFileNumber(0);
            histManager->SetDirectory(categoryNames[flCategory] + "/" + suffix);
            histManager->Fill1DHist(recoMET->DeltaPhi(fakeables[0].P2()),
                    "h1_MetFakeableDeltaPhi", "#Delta#phi(fakeable, MET);#Delta#phi(fakeable, MET);Entries / bin", 36, 0., TMath::Pi());

            vector<TCJet> muFakeJets    = selector->GetSelectedJets("muFakes");
            vector<TCJet> eleFakeJets   = selector->GetSelectedJets("eleFakes");

            for (unsigned i = 0; i < muFakeJets.size(); ++i) {
                histManager->Fill1DHist(muFakeJets[i].BDiscriminatorMap("CSV"),
                        "h1_MatchedMuJetBDiscr", "matched #mu-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
            }
            for (unsigned i = 0; i < eleFakeJets.size(); ++i) {
                histManager->Fill1DHist(eleFakeJets[i].BDiscriminatorMap("CSV"),
                        "h1_MatchedEleJetBDiscr", "matched e-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
            }

        }

        //for (unsigned i = 0; i < leptonsPlusFakes.size(); ++i) {
        //    for (unsigned j = 0; j < jets.size(); ++j) {
        //        if (leptonsPlusFakes[i].DeltaR(jets[j]) < 0.5) 
        //            cout << leptonsPlusFakes[i].DeltaR(jets[j]) << endl;
        //    }
        //}
>>>>>>> parent of 17d9010... Adding new BDT weight files

        // Enforce same-sign dilepton/trilepton selection with fake leptons
        if (leptonsPlusFakes.size() == 2) { 
            if (
                    leptonsPlusFakes[0].Charge() == leptonsPlusFakes[1].Charge()
                    && leptonsPlusFakes[0].Pt() > leptonPtCut[0] && leptonsPlusFakes[1].Pt() > leptonPtCut[1]
               ) {
                if (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "DATA_MUON" || suffix == "TEST") 
                    AnalysisSelection(leptonsPlusFakes, jets, bJetsM, bJetsL, PV, "Fakes");
                else                                                  
                    AnalysisSelection(leptonsPlusFakes, jets, bJetsM, bJetsL, PV, "Fakes_"+suffix);
            }
        } else if ( leptonsPlusFakes.size() == 3) {
            if (
                    leptonsPlusFakes[0].Pt() > leptonPtCut[0] 
                    && leptonsPlusFakes[1].Pt() > leptonPtCut[1]
                    && leptonsPlusFakes[2].Pt() > leptonPtCut[1]
                    && fabs(leptonsPlusFakes[0].Charge() + leptonsPlusFakes[1].Charge() + leptonsPlusFakes[2].Charge()) == 1
               ) {
                if (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "DATA_MUON" || suffix == "TEST") 
                    AnalysisSelection(leptonsPlusFakes, jets, bJetsM, bJetsL, PV, "Fakes");
                else                                                  
                    AnalysisSelection(leptonsPlusFakes, jets, bJetsM, bJetsL, PV, "Fakes_"+suffix);
            }
        }
        //cout << fakeWeight << "\t" << leptonsPlusFakes.size() << endl;
        evtWeight /= fakeWeight; // Remove fake weight from event weight

    } else
        return;
}

void fcncAnalyzer::MakePlots(vObj leptons, vector<TCJet> jets, vector<TCJet> bJets, TCMET met, TVector3 PV, unsigned cutLevel)
{

    // Fill histograms that correspond to event category
    // and histograms that are inclusive for 2 and 3 lepton
    // selections separately

    histManager->SetFileNumber(cutLevel);

    //cout << "\n" << endl;

    for (int i = 0; i < 4; ++i) { 
        unsigned histCategory = 0;

        switch (i) {
            case 1:
                // inclusive lepton categories
                histCategory = 1 + ((evtCategory.to_ulong() >> 2) & 0x3);
                break;
            case 2:
                // flavor categories
                histCategory = GetHistCategory(2) - 10;
                break;
            case 3:
                // WH categories
                histCategory = GetHistCategory(3);
                break;
            case 4:
                // eta categories
                histCategory = GetHistCategory(2);
                break;
            default:
                histCategory = 0;
        }

        if (i != 0 && histCategory == 0) continue;

        histManager->SetDirectory(categoryNames[histCategory] + "/" + subdir);

        LeptonPlots(leptons, jets, bJets, PV);
        MetPlots(met, leptons);
        JetPlots(jets, bJets);
        DileptonPlots2D(leptons);
        MiscPlots();

        histManager->SetWeight(1);
        histManager->Fill1DHist(primaryVtx->GetSize(),
                "h1_PvMultUnweighted", "Multiplicity of PVs", 51, -0.5, 50.);
        histManager->SetWeight(evtWeight);
    }
}

void fcncAnalyzer::LeptonPlots(vObj leptons, vector<TCJet> jets, vector<TCJet> bJets, TVector3 PV)
{

    unsigned centralCount = 0;

    for (unsigned i = 0; i < leptons.size(); ++i) {
        string index = str(i+1);

        histManager->Fill1DHist(leptons[i].Pt(),
                "h1_Lepton" + index + "Pt", "p_{T} lepton " + index + " ;p_{T}^{l" + index + "} (GeV);Entries / 5 GeV", 70, 0., 350.);
        histManager->Fill1DHist(leptons[i].Eta(),
                "h1_Lepton" + index + "Eta", "#eta lepton " + index + ";#eta^{l" + index + "};Entries / bin", 50, -2.5, 2.5);
        histManager->Fill1DHist(leptons[i].Phi(),
                "h1_Lepton" + index + "Phi", "#phi lepton " + index + ";#phi^{l" + index + "};Entries / bin", 36, -TMath::Pi(), TMath::Pi());


        histManager->Fill1DHist(leptons[i].Dxy(&PV), 
                "h1_Lepton" + index + "dxy", "d_{xy} leptons " + index + ";d_{xy} (cm);Entries / bin", 100., -.02, 0.02);
        histManager->Fill1DHist(leptons[i].Dz(&PV), 
                "h1_Lepton" + index + "dz", "d_{z} leptons " + index + ";d_{z} (cm);Entries / bin", 100., -0.15, 0.15);

        if (leptons[i].Type() == "electron") {
            histManager->Fill1DHist(leptons[i].Pt(),
                    "h1_ElectronPt", "p_{T} Electron;p_{T,e} (GeV);Entries / 5 GeV", 70, 0., 350.);
            histManager->Fill1DHist(leptons[i].Eta(),
                    "h1_ElectronEta", "#eta Electron;#eta_{e};Entries / bin", 50, -2.5, 2.5);
            histManager->Fill1DHist(leptons[i].Phi(),
                    "h1_ElectronPhi", "#phi Electron;#phi_{e};Entries / bin", 36, -TMath::Pi(), TMath::Pi());

            histManager->Fill1DHist(leptons[i].Dxy(&PV), 
                    "h1_ElectronDxy", "electron d_{xy};d_{xy} (cm);Entries / bin", 100., -.02, 0.02);
            histManager->Fill1DHist(leptons[i].Dz(&PV), 
                    "h1_ElectronDz", "electron d_{z};d_{z} (cm);Entries / bin", 100., -0.15, 0.15);

        } else if (leptons[i].Type() == "muon") {
            histManager->Fill1DHist(leptons[i].Pt(),
                    "h1_MuonPt", "p_{T} muon;p_{T,#mu} (GeV);Entries / 5 GeV", 70, 0., 350.);
            histManager->Fill1DHist(leptons[i].Eta(),
                    "h1_MuonEta", "#eta muon;#eta_{#mu};Entries / bin", 50, -2.5, 2.5);
            histManager->Fill1DHist(leptons[i].Phi(),
                    "h1_MuonPhi", "#phi muon;#phi_{#mu};Entries / bin", 36, -TMath::Pi(), TMath::Pi());

            histManager->Fill1DHist(leptons[i].Dxy(&PV), 
                    "h1_MuonDxy", "muon d_{xy};d_{xy};Entries / bin", 100., -.02, 0.02);
            histManager->Fill1DHist(leptons[i].Dz(&PV), 
                    "h1_MuonDz", "muon d_{z};d_{z};Entries / bin", 100., -0.15, 0.15);
        }

        if (fabs(leptons[i].Eta()) < 1.) 
            ++centralCount;

        if (bJets.size() > 0) {
            histManager->Fill1DHist(fabs(leptons[i].DeltaPhi(bJets[0])),
                    "h1_Lepton" + index + "BJetDeltaPhi", "#Delta #phi;#Delta #phi_{l" + index + ",b};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(leptons[i].Eta() - bJets[0].Eta()),
                    "h1_Lepton" + index + "BJetDeltaEta", "#Delta #eta; #eta_{l" + index + ",b};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(leptons[i].DeltaR(bJets[0])),
                    "h1_Lepton" + index + "BJetDeltaR", "#Delta R;#Delta R_{l" + index + ",b};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(leptons[i].Pt() - bJets[0].Pt())/(leptons[i].Pt() + bJets[0].Pt()),
                    "h1_Lepton" + index + "BJetDeltaPt", "#Delta p_{T,l" + index + "b}/#Sigma p_{T,l" + index + "b};#Delta p_{T,l" + index + "b}/#Sigma p_{T,l" + index + "b};Entries / bin", 100, 0., 1.);
        }

        if (jets.size() > 0) {
            histManager->Fill1DHist(fabs(leptons[i].DeltaPhi(jets[0])),
                    "h1_Lepton" + index + "JetDeltaPhi", "#Delta #phi;#Delta #phi_{l" + index + ",b};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(leptons[i].Eta() - jets[0].Eta()),
                    "h1_Lepton" + index + "JetDeltaEta", "#Delta #eta; #eta_{l" + index + ",b};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(leptons[i].DeltaR(jets[0])),
                    "h1_Lepton" + index + "JetDeltaR", "#Delta R;#Delta R_{l" + index + ",b};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(leptons[i].Pt() - jets[0].Pt())/(leptons[i].Pt() + jets[0].Pt()),
                    "h1_Lepton" + index + "JetDeltaPt", "#Delta p_{T,l" + index + "b}/#Sigma p_{T,l" + index + "b};#Delta p_{T,l" + index + "b}/#Sigma p_{T,l" + index + "b};Entries / bin", 100, 0., 1.);
        }

        for (unsigned j = 0; j < i; ++j) {
            string jndex = str(j+1);

            histManager->Fill1DHist((leptons[i] + leptons[j]).M(),
                    "h1_DileptonMass" + index + jndex, "dilepton M_{" + index + jndex + "};M_{" + index + jndex + "} (GeV/c^{2});Entries / 5 GeV", 80, 0., 400.);
            histManager->Fill1DHist((leptons[i] + leptons[j]).Mt(),
                    "h1_DileptonTransMass" + index + jndex, "dilepton M_{T," + index + jndex + "};M_{T," + index + jndex + "} (GeV/c^{2});Entries / 5 GeV", 100, 0., 500.);
            histManager->Fill1DHist((leptons[i] + leptons[j]).Pt(),
                    "h1_DileptonQt" + index + jndex, "dilepton q_{T," + index + jndex + "};q_{T}^{" + index + jndex + "} (GeV);Entries / 5 GeV", 100, 0., 500.);

            histManager->Fill1DHist(fabs(leptons[i].DeltaPhi(leptons[j])),
                    "h1_DileptonDeltaPhi" + index + jndex, "dilepton #Delta #phi_{" + index + jndex + "};#Delta #phi_{" + index + jndex + "};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(leptons[i].Eta() - leptons[j].Eta()),
                    "h1_DileptonDeltaEta" + index + jndex, "dilepton #Delta #eta_{" + index + jndex + "};#Delta #eta_{" + index + jndex + "};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(leptons[i].DeltaR(leptons[j])),
                    "h1_DileptonDeltaR" + index + jndex, "dilepton #Delta R_{" + index + jndex + "};#Delta R_{" + index + jndex + "};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(leptons[i].Pt() - leptons[j].Pt())/(leptons[i].Pt() + leptons[j].Pt()),
                    "h1_DileptonDeltaPt" + index + jndex, "dilepton #Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};#Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};Entries / bin", 50, 0., 1.);
            histManager->Fill1DHist(fabs((leptons[i] + leptons[j]).Pt())/(leptons[i].Pt() + leptons[j].Pt()),
                    "h1_DileptonBalance" + index + jndex, "dilepton #Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};#Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};Entries / bin", 50, 0., 1.);

        }
    }

    // OSSF dilepton pair variables
    if (ossfTagged) {
        histManager->Fill1DHist(dileptonP4.M(),
                "h1_DileptonOSMass", "OS dilepton M;M_{OS};Entries / 4 GeV", 100, 0., 400.);
        histManager->Fill1DHist(dileptonP4.Mt(),
                "h1_DileptonOSTransMass", "OS dilepton MT;MT_{OS};Entries / 5 GeV", 100, 0., 500.);
        histManager->Fill1DHist(dileptonP4.Pt(),
                "h1_DileptonOSQt", "dilepton q_{T,OS};q_{T}^{OS};Entries / 5 GeV", 100, 0., 500.);
        histManager->Fill1DHist(fabs(lep1P4.DeltaPhi(lep2P4)),
                "h1_DileptonOSDeltaPhi", "dilepton #Delta #phi_{OS};#Delta #phi_{OS};Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(lep1P4.Eta() - lep2P4.Eta()),
                "h1_DileptonOSDeltaEta", "dilepton #Delta #eta_{OS};#Delta #eta_{OS};Entries / bin", 60, 0., 6.);
        histManager->Fill1DHist(fabs(lep2P4.DeltaR(lep1P4)),
                "h1_DileptonOSDeltaR", "dilepton #Delta R_{OS};#Delta R_{OS};Entries / bin", 70, 0., 7.);
        histManager->Fill1DHist(fabs(lep1P4.Pt() - lep2P4.Pt())/(lep1P4.Pt() + lep2P4.Pt()),
                "h1_DileptonOSDeltaPt", "dilepton #Delta p_{T, OS}/#Sigma p_{T, OS};#Delta p_{T, OS}/#Sigma p_{T, OS};Entries / bin", 50, 0., 1.);
        histManager->Fill1DHist(dileptonP4.Pt()/(lep1P4.Pt() + lep2P4.Pt()),
                "h1_DileptonOSBalance", "dilepton #Delta p_{T, OS}/#Sigma p_{T, OS};#Delta p_{T, OS}/#Sigma p_{T, OS};Entries / bin", 50, 0., 1.);

        if (bJets.size() > 0) {
            histManager->Fill1DHist(fabs(dileptonP4.DeltaPhi(bJets[0])),
                    "h1_DileptonBJetDeltaPhi", "#Delta #phi;#Delta #phi_{ll,b};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(dileptonP4.Eta() - bJets[0].Eta()),
                    "h1_DileptonBJetDeltaEta", "#Delta #eta; #eta_{ll,b};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(dileptonP4.DeltaR(bJets[0])),
                    "h1_DileptonBJetDeltaR", "#Delta R;#Delta R_{ll,b};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(dileptonP4.Pt() - bJets[0].Pt())/(dileptonP4.Pt() + bJets[0].Pt()),
                    "h1_DileptonBJetDeltaPt", "#Delta p_{T,llb}/#Sigma p_{T,llb};#Delta p_{T,llb}/#Sigma p_{T,llb};Entries / bin", 50, 0., 1.);
        }

        if (jets.size() > 0) {
            histManager->Fill1DHist(fabs(dileptonP4.DeltaPhi(jets[0])),
                    "h1_DileptonJetDeltaPhi", "#Delta #phi;#Delta #phi_{ll,b};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(dileptonP4.Eta() - jets[0].Eta()),
                    "h1_DileptonJetDeltaEta", "#Delta #eta; #eta_{ll,b};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(dileptonP4.DeltaR(jets[0])),
                    "h1_DileptonJetDeltaR", "#Delta R;#Delta R_{ll,b};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(dileptonP4.Pt() - jets[0].Pt())/(dileptonP4.Pt() + jets[0].Pt()),
                    "h1_DileptonJetDeltaPt", "#Delta p_{T,llb}/#Sigma p_{T,llb};#Delta p_{T,llb}/#Sigma p_{T,llb};Entries / bin", 100, 0., 1.);
        }
    }

    if (leptons.size() == 3) {

        TLorentzVector trileptonP4 = leptons[0] + leptons[1] + leptons[2];
        histManager->Fill1DHist(trileptonP4.M(),
                "h1_TrileptonMass", "M_{lll};M_{lll};Entries / 5 GeV", 60, 0., 300.);
        histManager->Fill1DHist(trileptonP4.Pt(),
                "h1_TrileptonPt", "p_{T,3l};p_{T,3l};Entries / 5 GeV", 40, 0., 400.);

        histManager->Fill1DHist(dileptonP4.DeltaR(lep3P4), 
                "h1_DileptonLepDeltaR", "#Delta R(OS,l3);#Delta R(ll,l);Entries / bin", 70, 0., 7.);
        histManager->Fill1DHist(fabs(dileptonP4.DeltaPhi(lep3P4)), 
                "h1_DileptonLepDeltaPhi", "#Delta #phi(OS,l3);#Delta #phi(ll,l);Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(dileptonP4.Eta() - lep3P4.Eta()), 
                "h1_DileptonLepDeltaEta", "#Delta #eta(OS,l3);#Delta #eta(ll,l);Entries / bin", 60, 0., 6.);
        histManager->Fill1DHist(MT,
                "h1_Lep3MetMT", ";MT_{l3,MET};Entries / 5 GeV", 60, 0., 300.);
    }

    if (leptons.size() == 2) {
        if (leptons[0].Charge() == leptons[1].Charge()){
            histManager->Fill1DHist(leptons[0].Charge()*int(jets.size()),
                    "h1_JetMultCharge", "q_{ll}*N_{jets} for ss dileptons;q_{ll}*N_{jets};Entries / bin", 21, -10.5, 10.5);
        }
    }

    for (unsigned i = 0; i < bJets.size(); ++i) {
        if (fabs(bJets[i].Eta()) < 1.)
            ++centralCount;
    }
    for (unsigned i = 0; i < jets.size(); ++i) {
        if (fabs(jets[i].Eta()) < 1.)
            ++centralCount;
    }

    histManager->Fill1DHist(centralCount,
            "h1_Centrality", "Event Centrality;centrality;Entries / bin", 10, -0.5, 9.5);
    histManager->Fill1DHist(HT,
            "h1_HT", "HT;HT;Entries / 10 GeV", 75, 0., 1500.);
    histManager->Fill1DHist(sqrt(HT),
            "h1_HT", "HT;HT;Entries / bin", 80, 0., 40.);
    histManager->Fill1DHist(HTs + MET,
            "h1_HTs", "HT_{s};HT_{s};Entries / 10 GeV", 100, 0., 2000.);
    histManager->Fill1DHist(MHT/HTs,
            "h1_MHTOverHTs", "MHT/HT_{s};MHT/HT_{s};Entries / bin", 50, 0., 1.); 
    histManager->Fill1DHist(MHT,
            "h1_MHT", "MHT;MHT;Entries / bin", 35, 0., 350.); 
    histManager->Fill1DHist(MHT - MET,
            "h1_MHT-MET", "MHT - MET;MHT - MET;Entries / 6 GeV", 50, -150., 150.); 
    histManager->Fill1DHist(METLD,
            "h1_METLD", "METLD;METLD;Entries / bin", 100, 0., 1.); 

    histManager->Fill2DHist(HT, MET,
            "h2_metVsHt", "MET vs HT;HT;MET", 50, 0., 1000., 35, 0., 350.); 
    histManager->Fill2DHist(HTs, MET,
            "h2_metVsHts", "MET vs HTs;HTs;MET", 50, 0., 1000., 35, 0., 350.); 
    histManager->Fill2DHist(sqrt(HT), MET,
            "h2_metVsSqrtHt", "MET vs #sqrt{HT};#sqrt{HT};MET", 50, 0., 40., 35, 0., 350.); 
}


void fcncAnalyzer::MetPlots(TCMET met, vObj leptons)
{
    histManager->Fill1DHist(met.Mod(), 
            "h1_Met", "MET;MET;Entries / 10 GeV", 35, 0., 350.);
    histManager->Fill1DHist(met.Phi() - TMath::Pi(),
            "h1_MetPhi", "#phi MET;#phi;Entries / 0.087 rad", 36, -TMath::Pi(), TMath::Pi());
    histManager->Fill1DHist(met.SumEt(),
            "h1_MetSumEt", "#Sigma E_{T} of MET;#Sigma E_{T};Entries / 20 GeV", 75, 0., 2600.);

    float dPhiMin = TMath::Pi();
    unsigned iLep = 0;

    for (unsigned i = 0; i < leptons.size(); ++i) {
        string index = str(i+1);

        if (fabs(met.DeltaPhi(leptons[i].P2())) < dPhiMin) {
            dPhiMin = fabs(met.DeltaPhi(leptons[i].P2()));
            iLep = i;
        }

        histManager->Fill1DHist(fabs(met.DeltaPhi(leptons[i].P2())),
                "h1_MetLepton" + index + "DeltaPhi", "#Delta#phi(l" + index + ", MET);#Delta#phi(l" + index + ", MET);Entries / bin", 36, 0., TMath::Pi());

    }

    histManager->Fill1DHist(fabs(met.DeltaPhi(leptons[iLep].P2())),
            "h1_MetLepDeltaPhiMin", "min(#Delta#phi(l, MET));min(#Delta#phi(l, MET));Entries / bin", 36, 0., TMath::Pi());
    histManager->Fill1DHist(iLep + 1,
            "h1_nearLepIndex", "Index of lepton nearest to MET in phi;iLep;Entries / bin", 4, 0.5, 4.5);


    if (dPhiMin < TMath::Pi()/2.)
        histManager->Fill1DHist(met.Mod()*sin(dPhiMin), 
                "h1_ProjectedMet", "Projected MET;projMET (GeV);Entries / 10 GeV", 35, 0., 350.);
    else 
        histManager->Fill1DHist(met.Mod(), 
                "h1_ProjectedMet", "Projected MET;projMET (GeV);Entries / 10 GeV", 35, 0., 350.);
}


void fcncAnalyzer::JetPlots(vector<TCJet> jets, vector<TCJet> bJets)
{
    float ptBins[] = {15, 30, 40, 50, 60, 70, 90, 120, 160, 200, 1000};

    histManager->Fill1DHist(jets.size(),
            "h1_JetMult", "Multiplicity of jets;N_{jets};Entries / bin", 11, -0.5, 10.5);
    histManager->Fill1DHist(bJets.size(),
            "h1_BJetMult", "Multiplicity of b-jets;N_{b-jets};Entries / bin", 6, -0.5, 5.5);
    histManager->Fill1DHist(jets.size() + bJets.size(),
            "h1_AllJetMult", "Multiplicity of all jets;N_{jets};Entries / bin", 11, -0.5, 10.5);

    for (unsigned i = 0; i < jets.size(); ++i) {
        string index = str(i+1);

        histManager->Fill1DHist(jets[i].BDiscriminatorMap("CSV"),
                "h1_Jet" + index + "BDiscr", "jet " + index + " b-discriminator;CSV (jet " + index + ");Entries / 0.1", 46, -4., 0.6);
        histManager->Fill1DHist(jets[i].Pt(),
                "h1_Jet" + index + "Pt", "p_{T} of jet" + index + ";p_{T}^{j" + index + "};Entries / 10 GeV", 59, 10., 600.);
        histManager->Fill1DHist(jets[i].Eta(),
                "h1_Jet" + index + "Eta", "#eta of jet" + index + ";#eta^{j" + index + "};Entries / bin", 100, -5., 5.);
        histManager->Fill1DHist(jets[i].Phi(),
                "h1_Jet" + index + "Phi", "#phi of jet" + index + ";#phi^{j" + index + "};Entries / bin", 36, -TMath::Pi(), TMath::Pi());

        histManager->SetWeight(1);
        if (abs(jets[i].JetFlavor()) == 5) { // b-jets that are not correctly tagged
            histManager->Fill1DHistUnevenBins(jets[i].Pt(),
                    "h1_BTruthDenomPt", "b flavor jet p_{T};p_{T}", 10, ptBins);
        } else if (abs(jets[i].JetFlavor()) == 4) {
            histManager->Fill1DHistUnevenBins(jets[i].Pt(),
                    "h1_CTruthNumerPt", "c flavor jet p_{T};p_{T}", 10, ptBins);
            histManager->Fill1DHistUnevenBins(jets[i].Pt(),
                    "h1_CTruthDenomPt", "c flavor jet p_{T};p_{T}", 10, ptBins);
        } else if (abs(jets[i].JetFlavor()) != 0) {
            histManager->Fill1DHistUnevenBins(jets[i].Pt(),
                    "h1_MistagDenomPt", "mistag jet p_{T};p_{T}", 10, ptBins);
        }
        histManager->SetWeight(evtWeight);
    }

    for (unsigned i = 0; i < bJets.size(); ++i) {
        string index = str(i+1);

        histManager->Fill1DHist(bJets[i].BDiscriminatorMap("CSV"),
                "h1_BJet" + index + "BDiscr", "jet " + index + " b-discriminator;CSV (jet " + index + ");Entries / 0.1", 70, 0.5, 1.2);
        histManager->Fill1DHist(bJets[i].Pt(),
                "h1_BJet" + index + "Pt", "p_{T} of b-jet " + index + ";p_{T}^{b" + index + "};Entries / 10 GeV", 59, 10., 600.);
        histManager->Fill1DHist(bJets[i].Eta(),
                "h1_BJet" + index + "Eta", "#eta of b-jet " + index + ";#eta^{b" + index + "};Entries / bin", 50, -2.5, 2.5);
        histManager->Fill1DHist(bJets[i].Phi(),
                "h1_BJet" + index + "Phi", "#phi of b-jet " + index + ";#phi^{b" + index + "};Entries / bin", 36, -TMath::Pi(), TMath::Pi());


        histManager->SetWeight(1);
        if (abs(bJets[i].JetFlavor()) == 5) { // Correctly tagged b-jets
            histManager->Fill1DHistUnevenBins(bJets[i].Pt(),
                    "h1_BTruthNumerPt", "b flavor jet p_{T};p_{T}", 10, ptBins);
            histManager->Fill1DHistUnevenBins(bJets[i].Pt(),
                    "h1_BTruthDenomPt", "b flavor jet p_{T};p_{T}", 10, ptBins);
        } else if (abs(bJets[i].JetFlavor()) == 4) {
            histManager->Fill1DHistUnevenBins(bJets[i].Pt(),
                    "h1_CTruthDenomPt", "c flavor jet p_{T};p_{T}", 10, ptBins);
        } else if (abs(bJets[i].JetFlavor()) != 0) { // misidentified b-jets (b-tagged light jets)
            histManager->Fill1DHistUnevenBins(bJets[i].Pt(),
                    "h1_MistagDenomPt", "b flavor jet p_{T};p_{T}", 10, ptBins);
            histManager->Fill1DHistUnevenBins(bJets[i].Pt(),
                    "h1_MistagNumerPt", "b flavor jet p_{T};p_{T}", 10, ptBins);
        }
        histManager->SetWeight(evtWeight);
    }

    // Angular correlations between leading jets //
    if (bJets.size() > 0 && jets.size() > 0) {
        histManager->Fill1DHist(bJets[0].DeltaPhi(jets[0]), 
                "h1_JetBJetDeltaPhi", "#Delta#phi(j_{1}, b_{1});#Delta#phi(j_{1}, b_{1});Entries / bin", 36, 0., TMath::Pi()); 
        histManager->Fill1DHist(fabs(bJets[0].Eta() - jets[0].Eta()), 
                "h1_JetBJetDeltaEta", "#Delta#eta(j_{1}, b_{1});#Delta#eta(j_{1}, b_{1});Entries / bin", 50, 0., 9.); 
        histManager->Fill1DHist(bJets[0].DeltaR(jets[0]), 
                "h1_JetBJetDeltaR", "#Delta R(j_{1}, b_{1});#Delta R(j_{1}, b_{1});Entries / bin", 70, 0., 7.);
    }
}


void fcncAnalyzer::DileptonPlots2D(vObj leptons)
{

    for (unsigned i = 1; i < leptons.size(); ++i) {
        string index = str(i+1);

        for (unsigned j = 0; j < i; ++j) {
            string jndex = str(j+1);

            if (
                    leptons.size() == 3 
                    && leptons[i].Charge()  != leptons[j].Charge()
                    && leptons[i].Type()    == leptons[j].Type()
               ) {

                unsigned k = 3-i-j;
                string kndex = str(k+1);

                // Dalitz Plots //
                histManager->Fill2DHist((leptons[i] + leptons[j]).M(), (leptons[j] + leptons[k]).M(),
                        "h2_DileptonM" + jndex + kndex + "VsM" + index + jndex, "M_{" + jndex + kndex + "} vs M_{" + index + jndex + "};M_{" + jndex + kndex + "};M_{" + index + jndex + "}", 30, 0., 300., 30, 0., 300.);
                histManager->Fill2DHist(pow((leptons[i] + leptons[j]).M(), 2), pow((leptons[j] + leptons[k]).M(), 2),
                        "h2_DalitzM" + jndex + kndex + "VsM" + index + jndex, "M^{2}_{" + jndex + kndex + "} vs M^{2}_{" + index + jndex + "};M^{2}_{" + jndex + kndex + "};M^{2}_{" + index + jndex + "}", 30, 0., 90000., 30, 0., 90000.);


                histManager->Fill2DHist(leptons[i].DeltaR(leptons[j]), (leptons[i] + leptons[j]).Pt(),
                        "h2_DileptonQtVsDeltaROS", "q_{T} vs #Delta R;#Delta R;q_{T}", 70, 0., 7., 50, 0., 500.);
                histManager->Fill2DHist(leptons[i].DeltaR(leptons[j]), (leptons[i] + leptons[j]).M(),
                        "h2_DileptonMVsDeltaROS", "M_{OS} vs #Delta R;#Delta R;M_{OS}", 70, 0., 7., 40, 0., 400.);
                histManager->Fill2DHist((leptons[i] + leptons[j]).M(), (leptons[0] + leptons[1] + leptons[2]).M(),
                        "h2_TrileptonMVsDileptonMOS", "M_{lll} vs M_{OS};M_{OS};M_{lll}", 60, 0., 300., 60, 0., 300.);
            }
        } 
    }
}

void fcncAnalyzer::MiscPlots()
{
    histManager->Fill1DHist(chargeCat, 
            "h1_LeptonCharge", "lepton charge categories;charge category;Entries / bin", 12, 0.5, 12.5);
    histManager->Fill1DHist(flavorCat, 
            "h1_LeptonFlavor", "lepton flavor categories;flavor category;Entries / bin", 12, 0.5, 12.5);
    histManager->Fill2DHist(chargeCat, flavorCat, 
            "h2_LepChargeVsFlavor", "lepton flavor vs. charge;charge category;flavor category", 12, 0.5, 12.5, 12, 0.5, 12.5);


    histManager->Fill1DHist(evtWeight,
            "h1_EventWeight", "event weight", 100, 0., 3.);
    histManager->Fill1DHist(primaryVtx->GetSize(),
            "h1_PvMult", "Multiplicity of PVs", 51, -0.5, 50.);
    //histManager->Fill1DHist(primaryVtx[0].Z(),
    //        "h1_PvZ", "z_{PV};z_{PV};Entries / bin" 50, 0.5, 50.5);

    // Histograms for systematic errors
    histManager->Fill1DHist(weighter->GetFakeUncertainty(),
            "h1_FakeWeightUncertainty", "fake rate error;#sigma_{fake};Entries / bin", 40, 0., 0.1);
}

void fcncAnalyzer::MakeQMisIDPlots(vObj electrons)
{
    histManager->SetFileNumber(0);
    histManager->SetDirectory("inclusive/" + subdir);

    float ptBins[] = {0., 35., 75., 150.};
    float etaBins[] = {0., 0.8, 1.479, 2.5};

    unsigned iEta1, iPt1, iEta2, iPt2;

    // Set iEta bins for leading and trailing electrons
    if (fabs(electrons[0].Eta()) < 0.8)
        iEta1 = 0;
    else if (fabs(electrons[0].Eta()) > 0.8 && fabs(electrons[0].Eta()) < 1.479)
        iEta1 = 1;
    else if (fabs(electrons[0].Eta()) > 1.479)
        iEta1 = 2;

    if (fabs(electrons[1].Eta()) < 0.8)
        iEta2 = 0;
    else if (fabs(electrons[1].Eta()) > 0.8 && fabs(electrons[1].Eta()) < 1.479)
        iEta2 = 1;
    else if (fabs(electrons[1].Eta()) > 1.479)
        iEta2 = 2;


    // Set iPt bins for leading and trailing electrons
    if (electrons[0].Pt() < 35.)
        iPt1 = 1;
    else if (electrons[0].Pt() > 35 && electrons[0].Pt() < 75)
        iPt1 = 2;
    else if (electrons[0].Pt() > 75)
        iPt1 = 3;

    if (electrons[1].Pt() < 35.)
        iPt2 = 1;
    else if (electrons[1].Pt() > 35 && electrons[1].Pt() < 75)
        iPt2 = 2;
    else if (electrons[1].Pt() > 75)
        iPt2 = 3;

    //cout << "===========================" << endl;
    //cout << iPt1 << ", " << iEta1 << "\t\t" << electrons[0].Pt() << ", " << electrons[0].Eta() << "\t\t" << 3*iEta1 + iPt1 << endl;
    //cout << iPt2 << ", " << iEta2 << "\t\t" << electrons[1].Pt() << ", " << electrons[1].Eta() << "\t\t" << 3*iEta2 + iPt2 << endl;
    //cout << "===========================" << endl;

    if (electrons[0].Charge() == electrons[1].Charge()) {
        histManager->Fill2DHistUnevenBins(electrons[0].Pt(), electrons[0].Eta(),
                "h2_LeadElecQMisIDNumer", "lead e charge misID (numerator);p_{T};#eta", 3, ptBins, 3, etaBins); 
        histManager->Fill2DHistUnevenBins(electrons[1].Pt(), electrons[1].Eta(),
                "h2_TrailingElecQMisIDNumer", "trailing e charge misID (numerator);p_{T};#eta", 3, ptBins, 3, etaBins); 

        histManager->Fill2DHist(3*iEta1 + iPt1, 3*iEta2 + iPt2,
                "h2_DileptonQMisIDNumer", "e charge misID (numerator);e_{leading};e_{trailing}", 9, 0.5, 9.5, 9, 0.5, 9.5);
    }

    if (electrons[0].Charge() != electrons[1].Charge()) {
        histManager->Fill2DHistUnevenBins(electrons[0].Pt(), electrons[0].Eta(),
                "h2_LeadElecQMisIDDenom", "lead e charge misID (denominator);p_{T};#eta", 3, ptBins, 2, etaBins); 
        histManager->Fill2DHistUnevenBins(electrons[1].Pt(), electrons[1].Eta(),
                "h2_TrailingElecQMisIDDenom", "trailing e charge misID (denominator);p_{T};#eta", 3, ptBins, 2, etaBins); 

        histManager->Fill2DHist(3*iEta1 + iPt1, 3*iEta2 + iPt2,
                "h2_DileptonQMisIDDenom", "e charge misID (denominator);e_{leading};e_{trailing}", 9, 0.5, 9.5, 9, 0.5, 9.5);
    }
}

void fcncAnalyzer::Make4lPlots(vObj leptons, TCMET met)//, vector<TCJet> jets, vector<TCJet> bJets) 
{
    histManager->SetFileNumber(0);
    histManager->SetDirectory("inclusive/" + subdir);

    TLorentzVector tetraLeptonP4 = leptons[0] + leptons[1] + leptons[2] + leptons[3];
    float pt4l = leptons[0].Pt() + leptons[1].Pt() + leptons[2].Pt() + leptons[3].Pt();

    histManager->Fill1DHist(tetraLeptonP4.M(), 
            "h1_4lMass", "M_{4l};M_{4l};Entries / 5 GeV", 40, 50., 250.);
    histManager->Fill1DHist(tetraLeptonP4.Pt(), 
            "h1_4lPt", "p_{T,4l};p_{T,4l};Entries / 5 GeV", 20, 0., 500.);
    histManager->Fill1DHist(pt4l, 
            "h1_4lSumPt", "#Sigma p_{T,4l};#Sigma p_{T,4l};Entries / 5 GeV", 20, 0., 500.);

    histManager->Fill1DHist(met.Mod(), 
            "h1_4lMet", "MET (4l);MET;Entries / 10 GeV", 35, 0., 350.);

    //if (fabs(tetraLeptonP4.M() - 91.2) < 15) {
    //    histManager->Fill1DHist(leptons[3], 
    //            "h1_Z4lLep4Pt", "MET (4l);MET;Entries / 10 GeV", 35, 0., 350.);
    //}

}
void fcncAnalyzer::GenPlots(vector<TCGenParticle> gen, vObj leptons)
{

    float    effBins[]  = {10., 20., 30., 45., 65., 90., 140., 200., 400.};

    if (gen.size() == 3) {
        histManager->Fill1DHist((gen[0] + gen[1] + gen[2]).M(),
                "h1_GenTrileptonMass", "gen trilepton mass;M_{lll};Entries / 5 GeV", 60, 0., 300.);
    }

    unsigned higgsLepCount  = 0;
    //unsigned topLepCount    = 0;

    for (unsigned i = 0; i < gen.size(); ++i) {

        if (gen.size() == 3) {
            if (gen[i].Grandmother() == 25) {
                string index = str(higgsLepCount + 1);
                histManager->Fill1DHist(gen[i].Pt(),
                        "h1_GenHiggsLeptonPt" + index, "Higgs lepton " + index + " p_{T};p_{T," + index + "};Entries / 5 GeV", 40, 0., 200.);
                histManager->Fill1DHist(gen[i].Eta(),
                        "h1_GenHiggsLeptonEta" + index, "Higgs lepton " + index + " #eta;#eta_{" + index + "};Entries / bin", 50, -5., 5.);

                ++higgsLepCount;

            }

            if (fabs(gen[i].Grandmother()) == 6) {
                histManager->Fill1DHist(gen[i].Pt(),
                        "h1_TopLeptonPt", "top lepton p_{T};p_{T};Entries / 5 GeV", 40, 0., 200.);
                histManager->Fill1DHist(gen[i].Eta(),
                        "h1_TopLeptonEta", "top lepton #eta;#eta;Entries / bin", 50, -5., 5.);
            }
        }



        for (unsigned j = 0; j < leptons.size(); ++j) {
            bool eMatched   = false;
            bool muMatched  = false;

            // Match electrons
            if (fabs(gen[i].GetPDGId()) == 11 
                    && leptons[j].Type() == "electron"
                    && leptons[j].DeltaR(gen[i]) < 0.01) {

                short qxq = gen[i].Charge()*leptons[j].Charge(); 
                eMatched = true;

                histManager->Fill1DHist(qxq,
                        "h1_GenEleChargeMisId", ";q_{gen} x q_{reco};Electrons / bin", 3, -1.5, 1.5);
                histManager->Fill1DHist((gen[i].Pt() - leptons[j].Pt())/gen[i].Pt(),
                        "h1_GenEleBalance", ";(p_{T,gen} - p_{T,reco})/p{T,gen};Electrons / bin", 50, -2.5, 2.5);
                histManager->Fill1DHist(gen[i].DeltaR(leptons[j]),
                        "h1_GenEleDeltaR", ";#Delta_{R}(e_{gen},e_{reco});Electrons / bin", 30, 0, 0.3);

                if (qxq == -1) {
                    histManager->Fill1DHist(gen[i].Pt(),
                            "h1_GenEleMisIdPt", "p_{T} of charge mis-id e;p_{T};Electrons / 10 GeV", 29, 10., 300.);
                    histManager->Fill1DHist(gen[i].Eta(),
                            "h1_GenEleMisIdEta", "#eta of charge mis-id e;#eta;Electrons / bin", 25, -2.5, 2.5);
                }
            }


            if (fabs(gen[i].GetPDGId()) == 11 && eMatched)
                histManager->FillProfile(gen[i].Pt(), 1, "p1_GenEleEffPt", ";p_{T,gen};Efficiency", 8, effBins);

            // Match muons
            if (fabs(gen[i].GetPDGId()) == 13 
                    && leptons[j].Type() == "muon"
                    && leptons[j].DeltaR(gen[i]) < 0.01) {

                short qxq = gen[i].Charge()*leptons[j].Charge(); 
                muMatched = true;

                histManager->Fill1DHist(qxq,
                        "h1_GenMuChargeMisId", ";q_{gen} x q_{reco};Muons / bin", 3, -1.5, 1.5);
                histManager->Fill1DHist((gen[i].Pt() - leptons[j].Pt())/gen[i].Pt(),
                        "h1_GenMuBalance", ";(p_{T,gen} - p_{T,reco})/p{T,gen};Muons / bin", 50, -2.5, 2.5);
                histManager->Fill1DHist(gen[i].DeltaR(leptons[j]),
                        "h1_GenMuDeltaR", ";#Delta_{R}(#mu_{gen},#mu_{reco});Muons / bin", 30, 0, 0.3);

                if (qxq == -1) {
                    histManager->Fill1DHist(gen[i].Pt(),
                            "h1_GenMuMisIdPt", "p_{T} of charge mis-id #mu;p_{T};Muons / 10 GeV", 29, 10., 300.);
                    histManager->Fill1DHist(gen[i].Eta(),
                            "h1_GenMuMisIdEta", "#eta of charge mis-id #mu;#eta;Muons / bin", 25, -2.5, 2.5);
                }
            }

            if (fabs(gen[i].GetPDGId()) == 13 && muMatched)
                histManager->FillProfile(gen[i].Pt(), 1, "p1_GenMuEffPt", ";p_{T,gen};Efficiency", 7, effBins);
        }
    }
}

<<<<<<< HEAD
//void fcncAnalyzer::FakePlots(vObj leptons)
//{
//
//    if (fakeables.size() >= 1) {
//        unsigned flCategory = GetHistCategory(2) - 10;
//        histManager->SetFileNumber(0);
//        histManager->SetDirectory(categoryNames[flCategory] + "/Fakes");
//        histManager->Fill1DHist(recoMET->DeltaPhi(fakeables[0].P2()),
//                "h1_MetFakeableDeltaPhi", "#Delta#phi(fakeable, MET);#Delta#phi(fakeable, MET);Entries / bin", 36, 0., TMath::Pi());
//
//        vector<TCJet> muFakeJets    = selector->GetSelectedJets("muFakes");
//        vector<TCJet> eleFakeJets   = selector->GetSelectedJets("eleFakes");
//
//        for (unsigned i = 0; i < muFakeJets.size(); ++i) {
//            histManager->Fill1DHist(muFakeJets[i].BDiscriminatorMap("CSV"),
//                    "h1_MatchedMuJetBDiscr", "matched #mu-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
//        }
//        for (unsigned i = 0; i < eleFakeJets.size(); ++i) {
//            histManager->Fill1DHist(eleFakeJets[i].BDiscriminatorMap("CSV"),
//                    "h1_MatchedEleJetBDiscr", "matched e-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
//        }
//    }
//}

=======
>>>>>>> parent of 17d9010... Adding new BDT weight files
void fcncAnalyzer::SetEventCategory(vObj leptons)
{
    evtCategory.reset();
    evtCategory.set(0);

    //!! Eta categories
    if (fabs(leptons[0].Eta()) < 1.442) 
        evtCategory.set(7);
    if (fabs(leptons[1].Eta()) < 1.442) 
        evtCategory.set(6);

    //!! Flavor categories
    if (leptons[0].Type() == "muon")
        evtCategory.set(11);
    if (leptons[1].Type() == "muon")
        evtCategory.set(10);

    //!! Charge categories
    if (leptons[0].Charge() == 1)
        evtCategory.set(15);
    if (leptons[1].Charge() == 1)
        evtCategory.set(14);

    //!! Dilepton OS category !!//
    if (leptons.size() == 2)
        if (leptons[0].Charge() != leptons[1].Charge())
            evtCategory.set(2);

    //!! Trilepton categories
    if (leptons.size() == 3) {
        evtCategory.set(3);
        if (fabs(leptons[2].Eta()) < 1.442) 
            evtCategory.set(5);
        if (leptons[2].Type() == "muon")
            evtCategory.set(9);
        if (leptons[2].Charge() == 1)
            evtCategory.set(13);

        if (fabs(leptons[0].Charge() + leptons[1].Charge() + leptons[2].Charge()) == 1) {

            //!! WH categories !!//
            bool isSet = false;
            for (unsigned i = 1; i < leptons.size(); ++i) {
                for (unsigned j = 0; j < i; ++j) {
                    if (leptons[i].Type() == leptons[j].Type()) {
                        if (
                                leptons[i].Charge()     == leptons[j].Charge()
                                && leptons[i].Type()    != leptons[3-i-j].Type()
                                && leptons[i].Charge()  != leptons[3-i-j].Charge()
                           ) {
                            evtCategory.set(16);
                            isSet = true;
                            break;

                        } else if (
                                leptons[i].Charge() != leptons[j].Charge()
                                ) {
                            evtCategory.set(17);
                            isSet = true;
                            break;
                        }
                        if (isSet) break;
                    }
                }
            }
        }
    }
}


int fcncAnalyzer::GetHistCategory(unsigned shift)
{

    // Returns an index that maps onto category names.  Categories are done in
    // groups of 4 bits so each shift switches the category type.  As of now,
    // the standard categories go as following (in bits),

    // 0 - 3:    set bit, empty, os/ss, dilepton/trilepton
    // 4 - 7:    eta categories
    // 8 - 11:   flavor categories
    // 12 - 15:  charge categores

    // Additionally there is OSSF and SSSF 3 lepton categories for syncing with
    // the WH analysis.  

    //  16: OSSF
    //  17: SSSF

    //unsigned lepCat     = (evtCategory.to_ulong() >> 2) & 0x3;
    unsigned lepCat     = evtCategory.test(2) + 2*evtCategory.test(3);
    unsigned category   = (evtCategory.to_ulong() >> 4*shift) & 0xF;
    unsigned sum = 0;

    for (int i = 0; i < 4; ++i) sum += (category >> i) & 1;

    unsigned histCategory = 0;

    if (shift < 3) 
        histCategory = 4 + 10*(shift - 1) + 3*lepCat + sum;
    else 
        if (evtCategory.test(16))
            histCategory = 14;
        else if (evtCategory.test(17))
            histCategory = 15;
        else
            histCategory = 0;

    return histCategory;
}

void fcncAnalyzer::FillYieldHists(string directory, float weight, unsigned cut)
{
    for (unsigned i = 0; i < N_CUTS; ++i) {
        histManager->SetFileNumber(i);

        histManager->SetDirectory(directory);
        histManager->SetWeight(weight);
        histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
        histManager->SetWeight(1);
        histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
        histManager->SetWeight(weight);
    }

}

void fcncAnalyzer::SetYields(unsigned cut)
{
    if (evtCategory.none()) {
        for (unsigned i = 0; i < N_CATEGORIES; ++i) {
            FillYieldHists(categoryNames[i] + "/" + subdir, evtWeight, cut);
        }
    } else {

        // inclusive
        FillYieldHists(categoryNames[0] + "/" + subdir, evtWeight, cut);
        // inclusive for lepton category
        unsigned lepCat = (evtCategory.to_ulong() >> 2) & 0x3;
        FillYieldHists(categoryNames[lepCat+1] + "/" + subdir, evtWeight, cut);
        // flavor category
        unsigned flCat = GetHistCategory(2) - 10;
        FillYieldHists(categoryNames[flCat] + "/" + subdir, evtWeight, cut);
        // WH category
        unsigned whCat = GetHistCategory(3);
        if (whCat != 0) {
            FillYieldHists(categoryNames[whCat] + "/" + subdir, evtWeight, cut);
        }
    }

    if (subdir == suffix) {
        ++eventCount[cut];
        eventCountWeighted[cut] += evtWeight;
    }

}

void fcncAnalyzer::SetEventVariables(vObj leptons, vector<TCJet> jets, vector<TCJet> bJets, TCMET met) 
{
    TLorentzVector sumP4;

    HT  = 0.;
    HTs = 0.;
    for (unsigned i = 0; i < bJets.size(); ++i) { 
        HT      += bJets[i].Pt();
        HTs     += bJets[i].Pt();

        sumP4   += bJets[i];
    }

    for (unsigned i = 0; i < jets.size(); ++i) {
        HT      += jets[i].Pt();
        HTs     += jets[i].Pt();

        sumP4   += jets[i];
    }

    if (evtCategory.test(3)) {
        flavorCat = 5 + (((evtCategory.to_ulong() >> 8) & 0xF) >> 1);
        chargeCat = 5 + (((evtCategory.to_ulong() >> 12) & 0xF) >> 1);
    } else {
        flavorCat = 1 + (((evtCategory.to_ulong() >> 8) & 0xF) >> 2);
        chargeCat = 1 + (((evtCategory.to_ulong() >> 12) & 0xF) >> 2);
    }

    if (leptons.size() == 3) 
        trileptonMass   = (leptons[0] + leptons[1] + leptons[2]).M();
    else if (leptons.size() == 2)
        trileptonMass   = -1.;


    // Reset OS variables for each event //
    zTagged     = false;
    ossfTagged  = false;
    dileptonMassOS = -1.;

    float zCandidateMass = 0.;
    for (unsigned i = 0; i < leptons.size(); ++i) {
        HTs     += leptons[i].Pt();
        sumP4   += leptons[i];

        for (unsigned j = leptons.size()-1; j > i; --j) {

            // Check for opposite-sign pair //
            if (leptons[i].Charge() != leptons[j].Charge()) {
                ossfTagged = true;

                // Is the pair mass consistent with the Z mass within a 30 GeV window?
                if (fabs((leptons[i] + leptons[j]).M() - 91.2) < 15 && leptons[i].Type() == leptons[j].Type()) {
                    zTagged = true;
                    zCandidateMass = (leptons[i] + leptons[j]).M();
                }

                // Select the pairing that is closest to the mass of the Z.
                if (zTagged) {
                    if (fabs(dileptonMassOS - 91.2) > fabs(zCandidateMass - 91.2)) {
                        dileptonP4      = leptons[i] + leptons[j];
                        lep1P4          = leptons[j];
                        lep2P4          = leptons[i];

                        dileptonMassOS  = dileptonP4.M();
                        dileptonDROS    = leptons[i].DeltaR(leptons[j]);

                        // If 3 leptons present, try to reconstruct a W from the unpaired lepton
                        if (leptons.size() == 3) {
                            lep3P4  = leptons[3 - (i + j)];
                            MT = CalculateTransMass(leptons[3 - (i + j)], met);
                        }
                    }
                } else if (!zTagged) { 
                    if ((leptons[i] + leptons[j]).M() > dileptonMassOS) { // Pick the highest mass OS pairing
                        dileptonP4      = leptons[i] + leptons[j];
                        lep1P4          = leptons[j];
                        lep2P4          = leptons[i];

                        dileptonMassOS  = dileptonP4.M();
                        dileptonDROS    = leptons[i].DeltaR(leptons[j]);

                        if (leptons.size() == 3) {
                            lep3P4  = leptons[3 - (i + j)];
                            MT = CalculateTransMass(leptons[3 - (i + j)], met);
                        }
                    }
                }
            }
        }
    }

    if (!ossfTagged && leptons.size() >= 2) {
        dileptonP4      = leptons[0] + leptons[1];
        lep1P4          = leptons[0];
        lep2P4          = leptons[1];

        dileptonMassOS  = dileptonP4.M();
        dileptonDROS    = leptons[0].DeltaR(leptons[1]);

        if (leptons.size() == 3) {
            lep3P4  = leptons[2];
            MT      = CalculateTransMass(leptons[2], met);
        }
    }

    MHT     = sumP4.Pt();
    MET     = met.Mod();
    metPhi  = met.Phi();

    METLD   = 0.00397*MET + 0.00265*MHT;
}

void fcncAnalyzer::FillLepMVA(vector<TCMuon> muons, vector<TCElectron> electrons, vector<TCJet> jets, TVector3 PV)
{
    for (unsigned i = 0; i < muons.size(); ++i) {

        sip3d       = 1.; // Update this once proper definition is known
        chPFIso     = muons[i].IsoMap("pfChargedHadronPt_R04");
        neuPFIso    = TMath::Max(0.0, (double)muons[i].IsoMap("pfPhotonEt_R04") + muons[i].IsoMap("pfNeutralHadronEt_R04")); // - TMath::Max(0.0, (double)rho25Factor*Selector::EffectiveArea(muons[i])));

        dz  = fabs(muons[i].Dz(&PV));
        dxy = fabs(muons[i].Dxy(&PV));

        if (jets.size() > 0) {
            TCJet closestJet; // = jets[0];
            float dRMin = 99.;

            for (unsigned j = 0; j < jets.size(); ++j) {
                if (muons[i].DeltaR(jets[j]) < dRMin) {
                    closestJet = jets[j];
                    dRMin = muons[i].DeltaR(jets[j]);
                }
            }

            drLepJet        = dRMin;
            ptRatioLepJet   = muons[i].Pt()/closestJet.Pt();
            btagLepJet      = closestJet.BDiscriminatorMap("CSV");

        } else {
            drLepJet        = -1.;
            ptRatioLepJet   = -1.;
            btagLepJet      = -999.;
        }

        muTree->Fill();
    }

    for (unsigned i = 0; i < electrons.size(); ++i) {

        sip3d       = 1; // Update this once proper definition is known
        chPFIso     = electrons[i].IsoMap("pfChIso_R04");
        neuPFIso    = TMath::Max(0.0, (double)electrons[i].IsoMap("pfPhoIso_R04") + electrons[i].IsoMap("pfNeuIso_R04")); // - TMath::Max(0.0, (double)rho25Factor*Selector::EffectiveArea(electrons[i])));

        mva         = electrons[i].IdMap("mva");
        missHits    = electrons[i].NumberOfLostPixelHits();

        if (jets.size() > 0) {
            TCJet closestJet;
            float dRMin = 99.;

            for (unsigned j = 0; j < jets.size(); ++j) {
                if (electrons[i].DeltaR(jets[j]) < dRMin) {
                    closestJet = jets[j];
                    dRMin = electrons[i].DeltaR(jets[j]);
                }
            }

            drLepJet        = dRMin;
            ptRatioLepJet   = electrons[i].Pt()/closestJet.Pt();
            btagLepJet      = closestJet.BDiscriminatorMap("CSV");
        }
        eleTree->Fill();
    }
}

void fcncAnalyzer::SetVarsMVA(vObj leptons, vector<TCJet> bJets, vector<TCJet> jets)
{
    bJetMult        = bJets.size();
    jetMult         = jets.size();

    // Jet multiplicities as floats for mva reader :( 
    f_flavorCat     = flavorCat;
    f_jetMult       = jets.size();
    f_bJetMult      = bJets.size();

    // lepton
    for (unsigned i = 0; i < leptons.size(); ++i) {  
        TCPhysObject* lep = new ((*selLeptons)[i]) TCPhysObject;
        lep->SetPxPyPzE(leptons[i].Px(), leptons[i].Py(), leptons[i].Pz(), leptons[i].E()); 
        lep->SetVtx(0., 0., 0.);
        lep->SetCharge(leptons[i].Charge());
        lep->SetType(leptons[i].Type());
    }

    // jets 
    for (unsigned i = 0; i < jets.size(); ++i) {  
        TCPhysObject* jet = new ((*selJets)[i]) TCPhysObject;
        jet->SetPxPyPzE(jets[i].Px(), jets[i].Py(), jets[i].Pz(), jets[i].E()); 
        jet->SetVtx(0., 0., 0.);
        jet->SetCharge(0);
        jet->SetType("jet");
    }

    for (unsigned i = 0; i < bJets.size(); ++i) {  
        TCPhysObject* jet = new ((*selJets)[i]) TCPhysObject;
        jet->SetPxPyPzE(bJets[i].Px(), bJets[i].Py(), bJets[i].Pz(), bJets[i].E()); 
        jet->SetVtx(0., 0., 0.);
        jet->SetCharge(0);
        jet->SetType("b-jet");
    }
}

float fcncAnalyzer::CalculateFourLeptonMass(vObj leptons) {

    TLorentzVector Z1, Z2;
    unsigned index1 = 0;
    unsigned index2 = 0;
    float minDeltaMass = 1e10;

    // Find Z1
    for (unsigned i = 1; i < leptons.size(); ++i) {
        for (unsigned j = 0; j < i; ++j) {
            if (
                    leptons[i].Type() == leptons[j].Type()
                    && leptons[i].Charge() != leptons[j].Charge()
                    //&& (leptons[i].Pt() > 20. || leptons[j].Pt() > 20.)
                    && fabs((leptons[i] + leptons[j]).M() - 80) < 40
                    && fabs((leptons[i] + leptons[j]).M() - 91.2) < minDeltaMass
               ) {
                Z1 = leptons[i] + leptons[j];
                index1 = i;
                index2 = j;
            }
        }
    }

    if (index1 == index2) 
        return 0;


    for (unsigned i = 1; i < leptons.size(); ++i) {
        if (i == index1 || i == index2) continue;
        for (unsigned j = 0; j < i; ++j) {
            if (j == index1 || j == index2) continue;

            if (
                    leptons[i].Type() == leptons[j].Type()
                    && leptons[i].Charge() != leptons[j].Charge()
                    //&& (leptons[i].Pt() > 20. || leptons[j].Pt() > 20.)
                    && ((leptons[i] + leptons[j]).M() > 4 && (leptons[i] + leptons[j]).M() < 120)
               )
                Z2 = leptons[i] + leptons[j];

        }
    }

    if (Z1.M() == 0 || Z2.M() == 0) 
        return 0;
    else
        return (Z1 + Z2).M();
}

bool fcncAnalyzer::CosmicMuonFilter(TCPhysObject muon1, TCPhysObject muon2)
{
    float dimuonAngle = muon1.Angle(muon2.Vect());
    if (TMath::Pi() - dimuonAngle  < 0.05)
        return true;
    else
        return false;
}

float fcncAnalyzer::CalculateTransMass(TCPhysObject lep, TCMET met)
{
    float M_T = sqrt(2*met.Mod()*lep.Pt()*(1 - cos(met.DeltaPhi(lep.P2()))));
    return M_T;
}

TLorentzVector fcncAnalyzer::CalculateNuP4(TLorentzVector lep, TCMET met)
{
    const float M_W = 80.4;

    float lepEt = (pow(lep.E(),2) - pow(lep.Pz(),2));
    float xTerm = (lep.Px()*met.Px() + lep.Px()*met.Px() + (pow(M_W,2))/2.);

    float pz = (lep.Pz()*xTerm + lep.E()*sqrt(pow(xTerm,2) - pow(met.Mod(),2)*lepEt))/lepEt;

    TLorentzVector nuP4(met.Px(), met.Py(), max(float(0.), pz), sqrt(pow(met.Px(),2) + pow(met.Py(),2) + pow(pz,2)));

    return nuP4;
}

