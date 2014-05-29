#define fakeAnalyzer_cxx

#include "fakeAnalyzer.h"

using namespace std;

/////////////////
//Analysis cuts//
/////////////////

const bool  doQCDDileptonCR = true;
const bool  doZPlusJetCR    = false;
const bool  doAntiIso3l     = false;
const bool  doPureLep       = false;
const bool  doSameSign      = true;

const bool  doGenMatching   = false;

const float jetPtCut[]  = {25., 15.};
const float muPtCut[]   = {10., 3.};
const float elePtCut[]  = {10., 10.};
const float phoPtCut[]  = {10., 10.};


unsigned  nMetBins      = 10;
unsigned  nPtBins       = 8;
float     metBins[]     = {0., 10., 20., 30., 40., 50., 60., 70., 80., 100., 300.}; 
float     ptBins[]      = {10., 15., 20., 25., 30., 35., 40., 45., 50.}; 
float     etaBinsMu[]   = {0., 1.5, 2.5};
float     etaBinsEle[]  = {0., 0.8, 1.479, 2.5};

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool BTagSortCondition(TCJet j1, TCJet j2) {return (j1.BDiscriminatorMap("CSV") > j2.BDiscriminatorMap("CSV"));} 

void fakeAnalyzer::Begin(TTree* tree) 
{
    // Job config
    TString     option  = GetOption();
    TObjArray   *args   = (TObjArray*)option.Tokenize(" ");

    suffix      = (string)((TObjString*)args->At(0))->GetString();
    selection   = (string)((TObjString*)args->At(1))->GetString();
    period      = (string)((TObjString*)args->At(2))->GetString();

    cout << suffix << " " << selection << " " << period << endl;

    // Get trigger names from jobTree
    vector<string>* triggerNames = 0;
    TFile   *inFile         = tree->GetCurrentFile();
    TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");

    jobTree->SetBranchAddress("triggerNames", &triggerNames);
    jobTree->GetEntry();

    // Initialize utilities and selectors here //
    selector        = new Selector(muPtCut, elePtCut, jetPtCut, phoPtCut);
    weighter        = new WeightUtils(suffix, period, selection, isRealData);
    triggerSelector = new TriggerSelector(selection, "2012", *triggerNames, true, true);

    // Add single lepton triggers for fake rates //
    vstring triggers;
    if (selection == "muon") {
        //triggers.push_back("HLT_Mu17_Mu8_v");
        //triggers.push_back("HLT_Mu17_TkMu8_v");
        triggers.push_back("HLT_Mu8_v");
        triggers.push_back("HLT_Mu17_v");
    } else if (selection == "electron") {
        //triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
        triggers.push_back("HLT_Ele8_CaloIdT_TrkIdVL_v");
        triggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
        triggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
        triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
        triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
    } else if (selection == "muEG") {
        //triggers.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
        //triggers.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    }
    //triggerSelector->AddTriggers(triggers);

    // Random numbers! //
    //rnGenerator = new TRandom3();

    // Initialize histograms //
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    histManager = new HistManager();
    histoFile   = new TFile("fakeHistograms.root", "RECREATE");

    histoFile->mkdir("inclusive", "inclusive");
    histoFile->GetDirectory("inclusive", "inclusive")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("AntiIso3l", "AntiIso3l");
    histoFile->GetDirectory("AntiIso3l", "AntiIso3l")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("QCD2l", "QCD2l");
    histoFile->GetDirectory("QCD2l", "QCD2l")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("ZPlusJet", "ZPlusJet");
    histoFile->GetDirectory("ZPlusJet", "ZPlusJet")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("PureLep", "PureLep");
    histoFile->GetDirectory("PureLep", "PureLep")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("SameSign", "SameSign");
    histoFile->GetDirectory("SameSign", "SameSign")->mkdir(suffix.c_str(), suffix.c_str());

    histManager->AddFile(histoFile);
    histManager->SetFileNumber(0);

    cout << endl;

    for (unsigned i = 0; i < 16; ++i) {
        eventCount[i] = 0;
    }
}

bool fakeAnalyzer::Process(Long64_t entry)
{  
    GetEntry(entry);
    selector->PurgeObjects();

    histManager->SetDirectory("inclusive/" + suffix);

    if (eventCount[1] == 0) {
        weighter->SetDataBit(isRealData);
        triggerSelector->SetDataBit(isRealData);
        selector->SetDataBit(isRealData);
    }

    ++eventCount[1];

    if (eventCount[1] % (int)1e4 == 0) cout << eventCount[1] << " analyzed!" << endl;

    bool triggerPass = false;
    triggerPass = triggerSelector->SelectTrigger(triggerStatus, hltPrescale);

    //cout << triggerPass << endl;
    if (!triggerPass) return kTRUE;
    ++eventCount[2];

    vstring passNames = triggerSelector->GetPassNames();
    if (passNames.size() == 0) 
        passNames.push_back("NULL");


    selector->PVSelector(primaryVtx);
    if (selector->GetSelectedPVs().size() < 1) 
        return kTRUE;
    else
        selectedVtx = selector->GetSelectedPVs()[0];
    ++eventCount[3];

    /////////////////////////////
    // Get gen level particles //
    /////////////////////////////

    vector<TCGenParticle> higgs, dubyas, Zeds, gElectrons, gMuons, gTaus, gLeptons;
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

        // Hack to split ttbar sample
        if (suffix == "ttbarHad" && gLeptons.size() > 1)
            return kTRUE;
        else if (suffix == "ttbarLep" && gLeptons.size() != 2)
            return kTRUE;
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

    // Get leptons for tag and probe
    vector<TCMuon>      muTags      = selector->GetSelectedMuons("QCD2l_CR_tag");
    vector<TCMuon>      muProbes    = selector->GetSelectedMuons("probe");
    vector<TCElectron>  eleProbes   = selector->GetSelectedElectrons("probe");

    // Get tight leptons for Z-tagging
    vector<TCMuon>      muons       = selector->GetSelectedMuons("tight");
    vector<TCElectron>  electrons   = selector->GetSelectedElectrons("tight");

    vObj leptons;
    leptons.insert(leptons.begin(), muons.begin(), muons.end());
    leptons.insert(leptons.begin(), electrons.begin(), electrons.end());

    // Get non-isolated leptons
    vector<TCMuon>     muonsNoIso     = selector->GetSelectedMuons("tight_id");
    vector<TCElectron> electronsNoIso = selector->GetSelectedElectrons("loose_id");
    sort(muonsNoIso.begin(), muonsNoIso.end(), P4SortCondition);
    sort(electronsNoIso.begin(), electronsNoIso.end(), P4SortCondition);

    // Get jets
    vector<TCJet> jets      = selector->GetSelectedJets("tight");
    vector<TCJet> bJetsL    = selector->GetSelectedJets("bJetsLoose");
    vector<TCJet> bJetsM    = selector->GetSelectedJets("bJetsMedium");

    //jets.insert(jets.end(), bJetsL.begin(), bJetsL.end());
    jets.insert(jets.end(), bJetsM.begin(), bJetsM.end());

    vector<TCJet> tagJets;
    tagJets.insert(tagJets.end(), bJetsM.begin(), bJetsM.end());
    tagJets.insert(tagJets.end(), bJetsL.begin(), bJetsL.end());

    // Order collections by pt
    sort(tagJets.begin(), tagJets.end(), P4SortCondition);
    sort(jets.begin(), jets.end(), P4SortCondition);
    sort(bJetsL.begin(), bJetsL.end(), BTagSortCondition);
    sort(bJetsM.begin(), bJetsM.end(), BTagSortCondition);
    sort(leptons.begin(), leptons.end(), P4SortCondition);


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    //                            //
    // Start fake rate estimation //
    // here...                    //
    //                            //
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

    // Prepare control regions for FR estimation...  
    // The control region is largely defined by the tag, but we'll first check
    // for any probes since they will be the same for both regions of interest

    if (
            doQCDDileptonCR 
            && (muProbes.size() >= 1 || eleProbes.size() >= 1)
            && muTags.size() == 1
       ) {

        histManager->SetDirectory("QCD2l/" + suffix);

        // For description of QCD dilepton control region, see section 7.4.1 of
        // ttH note (AN-13-159).  First thing is to find the tag lepton and the
        // probe lepton. For this control region, the tag is a muon that is
        // displaced from the PV and is anti-isolated.  The probe is a lepton
        // passing loose identification requirement without any isolation
        // requirement

        TCMuon      muProbe;
        TCElectron  eleProbe; 
        unsigned nMuProbes  = 0;
        unsigned nEleProbes = 0;

        tag = (TCPhysObject) muTags[0];
        vector<TCJet> cleanJets;
        // Make jet collection which does not include tag
        for (unsigned i = 0; i < jets.size(); ++i) {
            if (jets[i].DeltaR(tag) > 0.3) {
                cleanJets.push_back(jets[i]);
            }
        }

        // Make sure there is only one probe lepton and that it does not
        // overlap with the tag. Consider both muon and electron probes.

        bitset<2> isQCD2l; 
        if (muProbes.size() >= 1) {
            // Check for at least one probe muon that doesn't overlap with the
            // tag
            for (unsigned i = 0; i < muProbes.size(); ++i) {
                if (tag.DeltaR(muProbes[i]) > 0.5) {
                    muProbe = muProbes[i]; 
                    ++nMuProbes;
                }
            }

            if (nMuProbes == 1) 
                if (CheckQCD2lCR(tagJets, muProbe)) 
                    isQCD2l.set(0); 

        }

        if (eleProbes.size() >= 1) {
            // Remove electron overlap with tag muon and make sure there is
            // only one electron probe
            for (unsigned i = 0; i < eleProbes.size(); ++i) {
                if (tag.DeltaR(eleProbes[i]) > 0.5) {
                    eleProbe = eleProbes[i];
                    ++nEleProbes;
                }
            }

            if (nEleProbes == 1) 
                if (CheckQCD2lCR(tagJets, eleProbe)) 
                    isQCD2l.set(1);
        }

        bool singleProbe = true;
        if (isQCD2l.test(0) && isQCD2l.test(1)) {
            // Only allow events with muon and electron probes if they are
            // the same object
            if (eleProbe.DeltaR(muProbe) > 0.3)
                singleProbe = false;
        }

        if ((isQCD2l.test(0) || isQCD2l.test(1)) && singleProbe) {
            // Probe object is found and event is consistent with QCD 2l
            // control region requirements.

            // MC-truth matching for cleaning prompt lepton
            // contamination.  Turn off to study fakes from MC
            if (!isRealData && doGenMatching) { 
                if (isQCD2l.test(0))
                    if (!GenProbeMatcher(muProbe, gLeptons))
                        isQCD2l.reset(0);
                if (isQCD2l.test(1))
                    if (!GenProbeMatcher(eleProbe, gLeptons))
                        isQCD2l.reset(1);
            }
            
            if (isQCD2l.test(0)) {
                FillJetHists(muProbe, cleanJets);
                FillDenominatorHists(muProbe);

                if (muProbe.IdMap("IsoRel") < 0.12)
                    FillNumeratorHists(muProbe);
                else if (muProbe.IdMap("IsoRel") > 0.20)
                    FillClosureHists(muProbe);
            }

            if (isQCD2l.test(1)) {
                FillJetHists(eleProbe, cleanJets);
                FillDenominatorHists(eleProbe);

                if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                    FillNumeratorHists(eleProbe);
                else //if (!selector->ElectronMVA(&eleProbe) || eleProbe.IdMap("IsoRel") > 0.15)
                    FillClosureHists(eleProbe);
            }
        }
    } 

    if (
            doZPlusJetCR
            && leptons.size() >= 2
            && (eleProbes.size() >= 1 || muProbes.size() >= 1)
       ) {

        histManager->SetDirectory("ZPlusJet/" + suffix);

        // For description of Z plus jet control region, see section 7.4.1 of
        // ttH note (AN-13-159).  The main difference between this control
        // region and the dilepton control region is the replacement of
        // b-tagged muon tag with a dilepton tag.  The dilepton is required to
        // be consistent with a Z and a low MET requirement is applied to
        // remove contamination from real trilepton events originating from
        // WZ->3l+nu events.

        DoZTag(leptons, 15.);
        if (zTagged) {
            tag = TCPhysObject(ZP4, 0, "Zed");
            tag.SetIdMap("IsoRel", 0.);

            TCMuon      muProbe;
            TCElectron  eleProbe; 
            unsigned nMuProbes  = 0;
            unsigned nEleProbes = 0;

            vector<TCJet> cleanJets;
            // Make jet collection which does not include tag
            for (unsigned i = 0; i < jets.size(); ++i) {
                if (jets[i].DeltaR(lep1) > 0.3 && jets[i].DeltaR(lep2)) {
                    cleanJets.push_back(jets[i]);
                }
            }

            // Remove probes that overlap with leptons used in Z reconstruction
            for (unsigned i = 0; i < muProbes.size(); ++i) {
                if (lep1.DeltaR(muProbes[i]) > 0.5 && lep2.DeltaR(muProbes[i]) > 0.5) {
                    muProbe = muProbes[i];
                    ++nMuProbes;
                }
            }

            for (unsigned i = 0; i < eleProbes.size(); ++i) {
                if (lep1.DeltaR(eleProbes[i]) > 0.5 && lep2.DeltaR(eleProbes[i]) > 0.5) {
                    eleProbe = eleProbes[i];
                    ++nEleProbes;
                }
            }

            // Only allow events with muon and electron probes to proceed
            // if they are the same object
            bool singleProbe = true;
            if (nEleProbes == 1 && nMuProbes == 1) {
                if (eleProbe.DeltaR(muProbe) > 0.3)
                    singleProbe = false;
            }

            if ((nEleProbes == 1 || nMuProbes == 1) && singleProbe) {
                // Probe object is found and event is consistent with Z+jet 
                // control region requirements. Now fill histograms for
                // parameterizing fake rates by pt and eta.

                // MC-truth matching for cleaning prompt lepton
                // contamination.  Turn off to study fakes from MC
                if (!isRealData && doGenMatching) { 
                    if (nMuProbes == 1)
                        if (!GenProbeMatcher(muProbe, gLeptons))
                            nMuProbes = 0;
                    if (nEleProbes == 1)
                        if (!GenProbeMatcher(eleProbe, gLeptons))
                            nEleProbes = 0;
                }

                if (nEleProbes == 1) {
                    if (CheckZPlusJetCR(eleProbe)) {
                        FillJetHists(eleProbe, cleanJets);
                        FillDenominatorHists(eleProbe);

                        if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                            FillNumeratorHists(eleProbe);
                        else //if (!selector->ElectronMVA(&eleProbe) || eleProbe.IdMap("IsoRel") > 0.15)
                            FillClosureHists(eleProbe);
                    }                 
                } 

                if (nMuProbes == 1) {
                    if (CheckZPlusJetCR(muProbe)) {
                        FillJetHists(muProbe, cleanJets);
                        FillDenominatorHists(muProbe);

                        if (muProbe.IdMap("IsoRel") < 0.12)
                            FillNumeratorHists(muProbe);
                        else if (muProbe.IdMap("IsoRel") > 0.20)
                            FillClosureHists(muProbe);
                    } 
                } 
            }
        }
    }

    if (doAntiIso3l) {
        // Find 2 anti-isolated muons and remove them from the muonsNoIso collection
        vObj leptonsAntiIso;
        if (muonsNoIso.size() >= 2) {
            if (muonsNoIso[0].IdMap("IsoRel") > 0.4 && muonsNoIso[1].IdMap("IsoRel") > 0.4) { 
                leptonsAntiIso.push_back(muonsNoIso[0]);
                leptonsAntiIso.push_back(muonsNoIso[1]);
                muonsNoIso.erase(muonsNoIso.begin());
                muonsNoIso.erase(muonsNoIso.begin());
            }
        }

        histManager->SetDirectory("AntiIso3l/" + suffix);

        if (leptonsAntiIso.size() == 2 && (muonsNoIso.size() == 1 || electronsNoIso.size() >= 1)) {

            // A tag for this CR exists if the leading two leptons are
            // anti-isolated the probe is then the trailing (third in pt) lepton

            TCMuon      muProbe;
            TCElectron  eleProbe; 
            unsigned nMuProbes  = 0;
            unsigned nEleProbes = 0;

            tag = leptonsAntiIso[0]; // No clear how to define this for this case
            vector<TCJet> cleanJets;
            // Make jet collection which does not include tag
            for (unsigned i = 0; i < jets.size(); ++i) {
                if (jets[i].DeltaR(leptonsAntiIso[0]) > 0.3 && jets[i].DeltaR(leptonsAntiIso[1]) > 0.3) {
                    cleanJets.push_back(jets[i]);
                }
            }

            // Find probes.  Ensure that they don't overlap with tag leptons
            nMuProbes = 0;
            if (muonsNoIso.size() == 1) {
                float muISO = muonsNoIso[0].IdMap("IsoRel");
                if (
                        muISO < 1. && !(muISO > 0.12 && muISO < 0.2)
                        && muonsNoIso[0].DeltaR(leptonsAntiIso[0]) > 0.5
                        && muonsNoIso[0].DeltaR(leptonsAntiIso[1]) > 0.5
                   ) {
                    muProbe = muonsNoIso[0];
                    ++nMuProbes;
                }
            }

            nEleProbes = 0;
            for (unsigned i = 0; i < electronsNoIso.size(); ++i) {
                float eleISO = electronsNoIso[i].IdMap("IsoRel");
                if (
                        eleISO < 1. && !(eleISO > 0.15 && eleISO < 0.2)
                        && electronsNoIso[i].DeltaR(leptonsAntiIso[0]) > 0.5
                        && electronsNoIso[i].DeltaR(leptonsAntiIso[1]) > 0.5
                   ) {
                    eleProbe = electronsNoIso[i];
                    ++nEleProbes;
                }
            }

            bool singleProbe = true;
            if (nEleProbes == 1 && nMuProbes == 1) {
                if (eleProbe.DeltaR(muProbe) > 0.3)
                    singleProbe = false;
            }

            if ((nEleProbes == 1 || nMuProbes == 1) && singleProbe) {
                // Probe object is found and event is consistent with AntiIso3l 
                // control region requirements. Now fill histograms for
                // parameterizing fake rates by pt and eta.
                
                // MC-truth matching for cleaning prompt lepton
                // contamination.  Turn off to study fakes from MC
                if (!isRealData && doGenMatching) { 
                    if (nMuProbes == 1)
                        if (!GenProbeMatcher(muProbe, gLeptons))
                            nMuProbes = 0;
                    if (nEleProbes == 1)
                        if (!GenProbeMatcher(eleProbe, gLeptons))
                            nEleProbes = 0;
                }


                if (nEleProbes == 1) {
                    FillJetHists(eleProbe, cleanJets);
                    FillDenominatorHists(eleProbe);

                    if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                        FillNumeratorHists(eleProbe);
                    else //if (!selector->ElectronMVA(&eleProbe) || eleProbe.IdMap("IsoRel") > 0.15)
                        FillClosureHists(eleProbe);
                } 

                if (nMuProbes == 1) {
                    histManager->Fill1DHist((leptonsAntiIso[0] + leptonsAntiIso[1]).M(), 
                            "h1_Tag1Tag2Mass", "M_{tag1,tag2};M_{tag1,tag2};Entries / 4 GeV", 50, 0., 200.);
                    histManager->Fill1DHist(leptonsAntiIso[0].DeltaR(leptonsAntiIso[1]), 
                            "h1_Tag1Tag2DeltaR", "#Delta R(tag1, tag2);#Delta R(tag1, tag2);Entries", 50, 0., 5.);
                    histManager->Fill2DHist((leptonsAntiIso[0] + muProbe).M(), leptonsAntiIso[0].Eta(), 
                            "h2_Tag1MuProbeMassVsEta", "M_{tag1,probe} vs. #eta;M_{tag1,probe};#eta", 50, 0., 200., 4, -2.5, 2.5);
                    histManager->Fill2DHist(leptonsAntiIso[0].DeltaR(muProbe), leptonsAntiIso[0].Eta(), 
                            "h2_MuProbeTag1DeltaRVsEta", "#Delta R(tag1, probe) vs. #eta;#Delta R(tag1, probe);#eta", 50, 0., 5., 4, -2.5, 2.5);
                    histManager->Fill2DHist((leptonsAntiIso[1] + muProbe).M(), leptonsAntiIso[1].Eta(), 
                            "h2_Tag2MuProbeMassVsEta", "M_{tag2,probe} vs. #eta;M_{tag2,probe};#eta", 50, 0., 200., 4, -2.5, 2.5);
                    histManager->Fill2DHist(leptonsAntiIso[1].DeltaR(muProbe), leptonsAntiIso[1].Eta(), 
                            "h2_MuProbeTag2DeltaRVsEta", "#Delta R(tag2, probe) vs. #eta;#Delta R(tag2, probe);#eta", 50, 0., 5., 4, -2.5, 2.5);

                    FillJetHists(muProbe, cleanJets);
                    FillDenominatorHists(muProbe);

                    if (muProbe.IdMap("IsoRel") < 0.12)
                        FillNumeratorHists(muProbe);
                    else if (muProbe.IdMap("IsoRel") > 0.20)
                        FillClosureHists(muProbe);
                }
            } 
        }
    }

    if (doPureLep) {
        histManager->SetDirectory("PureLep/" + suffix);
        if (muonsNoIso.size() == 2) {

            if (
                    muonsNoIso[0].IdMap("IsoRel") < 0.12 
                    && fabs((muonsNoIso[0] + muonsNoIso[1]).M() - 91.2) < 15
                    && muonsNoIso[0].DeltaR(muonsNoIso[1]) > 0.5
                    && muonsNoIso[0].Charge() != muonsNoIso[1].Charge()
               ) {

                TCMuon muProbe = muonsNoIso[1];
                tag = muonsNoIso[0];

                // Make jet collection which does not include tag
                vector<TCJet> cleanJets;
                for (unsigned i = 0; i < jets.size(); ++i) {
                    if (jets[i].DeltaR(tag) > 0.3) {
                        cleanJets.push_back(jets[i]);
                    }
                }

                FillJetHists(muProbe, cleanJets);
                FillDenominatorHists(muProbe);
                if (muProbe.IdMap("IsoRel") < 0.12)
                    FillNumeratorHists(muProbe);
            }
        }     
    }

    if (doSameSign) {
        histManager->SetDirectory("SameSign/" + suffix);
        if (muonsNoIso.size() == 2) {

            if (
                    muonsNoIso[0].IdMap("IsoRel") < 0.12 
                    && muonsNoIso[0].DeltaR(muonsNoIso[1]) > 0.5
                    && muonsNoIso[0].Charge() == muonsNoIso[1].Charge()
               ) {

                TCMuon muProbe = muonsNoIso[1];
                tag = muonsNoIso[0];

                // Make jet collection which does not include tag
                vector<TCJet> cleanJets;
                for (unsigned i = 0; i < jets.size(); ++i) {
                    if (jets[i].DeltaR(tag) > 0.3) {
                        cleanJets.push_back(jets[i]);
                    }
                }

                FillJetHists(muProbe, cleanJets);
                FillDenominatorHists(muProbe);
                if (muProbe.IdMap("IsoRel") < 0.12)
                    FillNumeratorHists(muProbe);
            }
        }     
    }

    /*
       vObj triggerLeps;
       if (crType == "ZPlusJet") {
       triggerLeps = leptons;
       } else if (crType == "QCD2l") {
       triggerLeps.push_back(tag);
       if (nMuProbes == 1)
       triggerLeps.push_back(muProbe);
       else if (nEleProbes == 1)
       triggerLeps.push_back(eleProbe);
       } 

       if (!isRealData) {
       weighter->SetObjects(triggerLeps, jets, nPUVerticesTrue, passNames[0]);
       Float_t evtWeight = weighter->GetTotalWeight();
       histManager->SetWeight(evtWeight);
       } else
       histManager->SetWeight(1.);
       */

    return kTRUE;
}

void fakeAnalyzer::Terminate()
{
    histManager->SetDirectory("inclusive/" + suffix);

    histManager->SetWeight(eventCount[0]);
    histManager->Fill1DHist(1, "h1_YieldByCut", ";cut;Entries / bin", 2, 0.5, 2.5);
    histManager->SetWeight(eventCount[1]);
    histManager->Fill1DHist(2, "h1_YieldByCut", ";cut;Entries / bin", 2, 0.5, 2.5);

    cout << "\nFake rate estimation finished." << endl;

    // Clean-up and write to output files
    selector->Delete();
    triggerSelector->Delete();
    histManager->Delete();

    histoFile->Write();
    histoFile->Close();  
}

void fakeAnalyzer::DoZTag(vObj& leptons, float window)
{
    // Reset OS variables for each event of interest//
    zTagged      = false;
    dileptonMass = 0.;

    float zCandidateMass = 0.;
    for (unsigned i = 0; i < leptons.size(); ++i) {
        for (unsigned j = leptons.size()-1; j > i; --j) {

            //Find lepton pair that is consistent with a Z in some window
            if (
                    fabs((leptons[i] + leptons[j]).M() - 91.2) < window 
                    && leptons[i].Type() == leptons[j].Type()
                    && leptons[i].Charge() != leptons[j].Charge()
               ) {
                zTagged = true;
                zCandidateMass = (leptons[i] + leptons[j]).M();
            }

            // Select the pairing that is closest to the mass of the Z.
            if (zTagged && fabs(dileptonMass - 91.2) > fabs(zCandidateMass - 91.2)) {
                dileptonMass = zCandidateMass;
                lep1    = leptons[i];
                lep2    = leptons[j];
                ZP4     = (TLorentzVector)(leptons[i] + leptons[j]);

                if (leptons.size() == 3) {
                    lep3 = leptons[3-i-j];
                }
            }
        }
    }
}

void fakeAnalyzer::FillDenominatorHists(TCPhysObject& probe)
{
    string lepType = "";
    if (probe.Type() == "muon")
        lepType = "Mu";
    else if (probe.Type() == "electron")
        lepType = "Ele";

    histManager->Fill1DHist(tag.Pt(),
            "h1_" + lepType + "TagLepPt", "tag lepton p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
    histManager->Fill1DHist(tag.Eta(),
            "h1_" + lepType + "TagLepEta", "tag lepton #eta;#eta;Entries / bin", 25, -2.5, 2.5);
    histManager->Fill1DHist(tag.IdMap("IsoRel"),
            "h1_" + lepType + "TagIsoRel", "tag lepton ISO_{rel};ISO_{rel};Entries", 80, 0., 4.);
    histManager->Fill1DHist(fabs(tag.Dz(selectedVtx)),
            "h1_" + lepType + "TagDz", "tag  d_{z};d_{z};Entries", 50, 0., 2.);
    histManager->Fill1DHist(fabs(tag.Dxy(selectedVtx)),
            "h1_" + lepType + "TagDxy", "tag  d_{xy};d_{xy};Entries", 50, 0., 2.);

    histManager->Fill1DHist((tag + probe).M(),
            "h1_Tag" + lepType + "ProbeMass", "M_{tag,probe};M_{tag,probe};Entries / 4 GeV", 50, 0., 200.);
    histManager->Fill1DHist((tag + probe).M(),
            "h1_Tag" + lepType + "ProbeMassHiggs", "M_{tag,probe};M_{tag,probe};Entries / 2 GeV", 25, 100., 150.);
    histManager->Fill1DHist(tag.DeltaR(probe),
            "h1_Tag" + lepType + "ProbeDeltaR", "#Delta R(tag,probe);#Delta R(tag,probe);Entries / 3 GeV", 50, 0., 5.);

    // fake rate measurement plots
    histManager->Fill1DHist(probe.Pt(),
            "h1_" + lepType + "ProbeLepPt", "probe muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
    histManager->Fill1DHist(probe.Eta(),
            "h1_" + lepType + "ProbeLepEta", "probe muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);
    histManager->Fill1DHist(CalculateTransMass(probe, *recoMET),
            "h1_" + lepType + "ProbeTransverseMass", "MT muon probe;MT;Entries / bin", 75, 0., 150.);

    histManager->Fill2DHist(probe.IdMap("IsoRel"), probe.Pt(), 
            "h2_" + lepType + "ProbePtVsIso", "probe lepton p_{T} vs ISO_{rel};ISO_{rel};p_{T}", 60, 0., 6., 5, 10., 35.);

    if (probe.Pt() > 10. && probe.Pt() <= 15) {
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_" + lepType + "DenomIsoRelBin1", "probe muon IsoRel;IsoRel;Entries", 40, 0., 1.2);
    } else if (probe.Pt() > 15. && probe.Pt() <= 20) {
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_" + lepType + "DenomIsoRelBin2", "probe muon IsoRel;IsoRel;Entries", 40, 0., 1.2);
    }

    if (probe.Pt() < 50.) {
        if (probe.Type() == "muon") {
            histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                    "h1_" + lepType + "DenomEta", "probe muon #eta;#eta;Entries / bin", 2, etaBinsMu);
            histManager->Fill2DHistUnevenBins(probe.Pt(), fabs(probe.Eta()),
                    "h2_" + lepType + "Denom", "probe muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);
        } else if (probe.Type() == "electron") {
            histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                    "h1_" + lepType + "DenomEta", "probe electron #eta;#eta;Entries / bin", 3, etaBinsEle);
            histManager->Fill2DHistUnevenBins(probe.Pt(), fabs(probe.Eta()),
                    "h2_" + lepType + "Denom", "probe electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);
        }

        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_" + lepType + "DenomPt", "probe lepton p_{T};p_{T};Entries / bin", nPtBins, ptBins);

        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_" + lepType + "DenomIsoRel", "probe lepton IsoRel;IsoRel;Entries", 40, 0., 1.2);
        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_" + lepType + "DenomMet", "probe lepton Met;Met;Entries", nMetBins, metBins);
    }
}

void fakeAnalyzer::FillNumeratorHists(TCPhysObject& probe)
{
    string lepType = "";
    if (probe.Type() == "muon")
        lepType = "Mu";
    else if (probe.Type() == "electron")
        lepType = "Ele";

    histManager->Fill1DHist((tag + probe).M(),
            "h1_Tag" + lepType + "PassMass", "M_{tag,probe};M_{tag,probe};Entries / 4 GeV", 50, 0., 200.);
    histManager->Fill1DHist(tag.DeltaR(probe),
            "h1_Tag" + lepType + "PassDeltaR", "#Delta R(tag,probe);#Delta R(tag,probe);Entries / 3 GeV", 50, 0., 5.);

    histManager->Fill1DHist(probe.Pt(),
            "h1_" + lepType + "PassLepPt", "pass lepton p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
    histManager->Fill1DHist(probe.Eta(),
            "h1_" + lepType + "PassLepEta", "pass lepton #eta;#eta;Entries / bin", 25, -2.5, 2.5);

    if (probe.Type() == "muon") {
        if (probe.Pt() < 50.) {
            histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                    "h1_" + lepType + "NumerEta", "pass muon #eta;#eta;Entries", 2, etaBinsMu);
        }
        histManager->Fill2DHistUnevenBins(probe.Pt(), fabs(probe.Eta()),
                "h2_" + lepType + "Numer", "pass muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);
    } else if (probe.Type() == "electron") {
        if (probe.Pt() < 50.) {
            histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                    "h1_" + lepType + "NumerEta", "pass electron #eta;#eta;Entries", 3, etaBinsEle);
        }
        histManager->Fill2DHistUnevenBins(probe.Pt(), fabs(probe.Eta()),
                "h2_" + lepType + "Numer", "pass electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);
    }

    histManager->Fill1DHistUnevenBins(probe.Pt(),
            "h1_" + lepType + "NumerPt", "pass lepton p_{T};p_{T};Entries", nPtBins, ptBins);
    histManager->Fill1DHist(probe.IdMap("IsoRel"),
            "h1_" + lepType + "NumerIsoRel", "pass lepton IsoRel;IsoRel;Entries", 40, 0., 0.20);
    histManager->Fill1DHistUnevenBins(recoMET->Mod(),
            "h1_" + lepType + "NumerMet", "pass lepton Met;Met;Entries", nMetBins, metBins);
}

void fakeAnalyzer::FillClosureHists(TCPhysObject& probe)
{
    string lepType = "";
    if (probe.Type() == "muon")
        lepType = "Mu";
    else if (probe.Type() == "electron")
        lepType = "Ele";

    string categories[3] = {"QCD2l", "ZPlusJet", "AntiIso3l"};
    for (unsigned i = 0; i < 1; ++i) {
        string cat = categories[i];
        float fakeWeight = weighter->GetFakeWeight(probe, cat);
        histManager->SetWeight(fakeWeight);

        //cout << cat << ":" << fakeWeight << ":(" << probe.Type() << "; " << probe.Pt() << "; " << probe.Eta() << ")" << endl;

        histManager->Fill1DHist(probe.Pt(),
                "h1_" + lepType + "PtClosure_" + cat, "p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_" + lepType + "EtaClosure_" + cat, "#eta (closure);#eta;Entries / 3 GeV", 5, 0., 2.5);

        if (probe.Pt() < 50.) {
            if (probe.Type() == "muon") {
                histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                        "h1_" + lepType + "UnevenEtaClosure_" + cat, "muon #eta (closure);#eta;Entries / bin", 2, etaBinsMu);
            } else if (probe.Type() == "electron") {
                histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                        "h1_" + lepType + "UnevenEtaClosure_" + cat, "electron #eta (closure);#eta;Entries / bin", 3, etaBinsEle);
            }
        }

        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_" + lepType + "UnevenPtClosure_" + cat, "p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_" + lepType + "IsoRelClosure_" + cat, "pass lepton IsoRel;IsoRel;Entries", 20, 0., 1.0);
        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_" + lepType + "MetClosure_" + cat, "pass lepton IsoRel;IsoRel;Entries", nMetBins, metBins);
    }
    histManager->SetWeight(1.);
}

void fakeAnalyzer::FillJetHists(TCPhysObject& probe, vector<TCJet>& jets) 
{
    string lepType = "";
    if (probe.Type() == "muon")
        lepType = "Mu";
    else if (probe.Type() == "electron")
        lepType = "Ele";

    //Match jet to probe and remove it from collection of jets
    bool jetMatched = false;
    TCJet probeJet;
    for (unsigned i = 0; i < jets.size(); ++i) {
        if (jets[i].DeltaR(probe) < 0.3) {
            probeJet = jets[i];
            jetMatched = true;
            jets.erase(jets.begin() + i);
            --i;
        }
    }

    //cout << jets[jetIndex].BDiscriminatorMap("CSV") << endl;
    histManager->Fill1DHist(jets.size(),
            "h1_" + lepType + "JetMult", "jet multiplicity; N_{jets}; Entries / bin", 10, -0.5, 9.5);
    histManager->Fill2DHist(jets.size(), probe.IdMap("IsoRel"),
            "h2_" + lepType + "DenomIsoRelVsJetMult", "ISO_{rel} vs. jet multiplicity;N_{jets};ISO_{rel}", 6, -0.5, 5.5, 5, 0.15, 1.15);

    if (!isRealData && jetMatched) {
        unsigned jetFlavor = probeJet.JetFlavor();
        histManager->Fill1DHist(probeJet.BDiscriminatorMap("CSV"),
                "h1_Matched" + lepType + "JetBDiscr", "matched #mu-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
        histManager->Fill1DHist(jetFlavor,
                "h1_Matched" + lepType + "JetFlavor", "matched #mu-jet Flavor;Flavor;Entries / bin", 25, 0., 25);

        if (jetFlavor == 1 || jetFlavor == 2 || jetFlavor == 3 || jetFlavor == 21) {
            histManager->Fill1DHist(probe.IdMap("IsoRel"),
                    "h1_LightMatched" + lepType + "IsoRel", "matched lepton IsoRel;IsoRel;Entries / bin", 20, 0., 1.);
            histManager->Fill1DHist(probe.Pt(),
                    "h1_LightMatched" + lepType + "IsoRel", "matched lepton IsoRel;IsoRel;Entries / 5 GeV", 40, 10., 50.);
            histManager->Fill1DHist(fabs(probe.Eta()),
                    "h1_LightMatched" + lepType + "IsoRel", "matched lepton IsoRel;IsoRel;Entries / bin", 2, 0., 2.5);
        } else if (jetFlavor == 4 || jetFlavor == 5) {
            histManager->Fill1DHist(probe.IdMap("IsoRel"),
                    "h1_HeavyMatched" + lepType + "IsoRel", "matched lepton IsoRel;IsoRel;Entries / bin", 20, 0., 1.);
            histManager->Fill1DHist(probe.Pt(),
                    "h1_HeavyMatched" + lepType + "IsoRel", "matched lepton IsoRel;IsoRel;Entries / 5 GeV", 40, 10., 50.);
            histManager->Fill1DHist(fabs(probe.Eta()),
                    "h1_HeavyMatched" + lepType + "IsoRel", "matched lepton IsoRel;IsoRel;Entries / bin", 2, 0., 2.5);
        }
    }
}

bool fakeAnalyzer::CheckQCD2lCR(vector<TCJet>& tagJets, TCPhysObject& probe) 
{
    // Make sure event is consistent with bbbar production 
    // Match tag to a loose b-jet
    bool jetMatched = false;
    bool jetVeto    = false;
    for (unsigned i = 0; i < tagJets.size(); ++i) {
        //cout << tag.DeltaR(tagJets[i]) << "\t" << probe.DeltaR(tagJets[i]) << endl;
        if (tag.DeltaR(tagJets[i]) < 0.3) {
            jetMatched = true;
            continue;
        }
    }

    // Check tag/probe pair is back-to-back
    Float_t tpDeltaPhi  = tag.DeltaPhi(probe);
    Float_t tpBalance   = probe.Pt()/(tag.Pt()*(1 + tag.IdMap("IsoRel"))); 

    if (probe.Type() == "muon") {
        histManager->Fill1DHist(fabs(tpDeltaPhi),
                "h1_TagMuProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(tpBalance),
                "h1_TagMuProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);
    } else if (probe.Type() == "electron") {
        histManager->Fill1DHist(fabs(tpDeltaPhi),
                "h1_TagEleProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(tpBalance),
                "h1_TagEleProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);
    }

    if (fabs(tpDeltaPhi) > 2. && tpBalance < 1 && jetMatched) 
        return true;
    else
        return false;
}

bool fakeAnalyzer::CheckZPlusJetCR(TCPhysObject& probe)
{
    // Check tag/probe pair is back-to-back
    Float_t tpDeltaPhi  = tag.DeltaPhi(probe);
    Float_t tpBalance   = probe.Pt()/tag.Pt(); 

    histManager->SetDirectory("ZPlusJet/" + suffix);
    if (probe.Type() == "muon") {
        histManager->Fill1DHist(tag.M(),
                "h1_ZTagMassMu", "M_{ll};M_{ll} (GeV);Entries / GeV", 60, 60., 120.);
        histManager->Fill1DHist(fabs(tpDeltaPhi),
                "h1_TagMuProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(tpBalance),
                "h1_TagMuProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);
    } else if (probe.Type() == "electron") {
        histManager->Fill1DHist(tag.M(),
                "h1_ZTagMassEle", "M_{ll};M_{ll} (GeV);Entries / GeV", 60, 60., 120.);
        histManager->Fill1DHist(fabs(tpDeltaPhi),
                "h1_TagEleProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(tpBalance),
                "h1_TagEleProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);
    }

    if (recoMET->Mod() < 50)
        return true;
    else
        return false;
}

bool fakeAnalyzer::GenProbeMatcher(TCPhysObject& probe, vector<TCGenParticle>& gLeptons)
{
    bool genMatched = false;
    for (unsigned i = 0; i < gLeptons.size(); ++i) {
        if (probe.Type() == "muon") {
            histManager->Fill1DHist(probe.DeltaR(gLeptons[i]), 
                    "h1_GenLepMuProbeDeltaR", "#Delta R(gen, #mu);#Delta R(gen, #mu);Entries", 50, 0., 5.);
            histManager->Fill1DHist((probe.Pt() - gLeptons[i].Pt())/gLeptons[i].Pt(), 
                    "h1_GenLepMuProbePtBalance", "(p_{T,#mu} - p_{T,gen})/p_{T,gen};(p_{T,#mu} - p_{T,gen})/p_{T,gen};Entries", 150, -0.5, 1.);
        }
        if (probe.Type() == "electron") {
            histManager->Fill1DHist(probe.DeltaR(gLeptons[i]), 
                    "h1_GenLepEleProbeDeltaR", "#Delta R(gen, #ele);#Delta R(gen, #ele);Entries", 50, 0., 5.);
            histManager->Fill1DHist((probe.Pt() - gLeptons[i].Pt())/gLeptons[i].Pt(), 
                    "h1_GenLepEleProbePtBalance", "(p_{T,#ele} - p_{T,gen})/p_{T,gen};(p_{T,#ele} - p_{T,gen})/p_{T,gen};Entries", 150, -0.5, 1.);
        }

        if (gLeptons[i].DeltaR(probe) < 0.1) 
            genMatched = true;
    }

    if (genMatched) 
        return true;
    else {

        float leptonIso = probe.IdMap("IsoRel");
        if (probe.Type() == "muon" && leptonIso > 0.)
            histManager->Fill1DHist(leptonIso, "h1_FakeMuonIso", "fake muon ISO_{rel};ISO_{rel,#mu};Entries", 80, 0., 4.);
        if (probe.Type() == "electron" && leptonIso > 0.)
            histManager->Fill1DHist(leptonIso, "h1_FakeElectronIso", "fake electron ISO_{rel};ISO_{rel,#mu};Entries", 80, 0., 4.);

        return false;
    }
}

float fakeAnalyzer::CalculateTransMass(TCPhysObject& lep, TCMET& met)
{
    float M_T = sqrt(2*met.Mod()*lep.Pt()*(1 - cos(met.DeltaPhi(lep.P2()))));
    return M_T;
}

