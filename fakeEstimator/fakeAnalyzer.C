#define fakeAnalyzer_cxx

#include "fakeAnalyzer.h"

using namespace std;

/////////////////
//Analysis cuts//
/////////////////

const bool  doQCDDileptonCR = true;
const bool  doZPlusJetCR    = true;
const bool  doAntiIso3l     = true;

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

    histManager->AddFile(histoFile);
    histManager->SetFileNumber(0);

    cout << endl;

    for (unsigned i = 0; i < 2; ++i) {
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

    vstring passNames = triggerSelector->GetPassNames();
    if (passNames.size() == 0) 
        passNames.push_back("NULL");


    selector->PVSelector(primaryVtx);
    if (selector->GetSelectedPVs().size() < 1) 
        return kTRUE;
    else
        selectedVtx = selector->GetSelectedPVs()[0];

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
        if (suffix == "ttbarHad" && gLeptons.size() != 1)
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

    // Get jets
    jets = selector->GetSelectedJets("tight");
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

    bool isTP    = false;
    bool elPass  = false;
    bool muPass  = false;
    bool elFake  = false;
    bool muFake  = false;
    string crType  = "None";

    // Prepare control regions for FR estimation...  
    // The control region is largely defined by the tag, but we'll first check
    // for any probes since they will be the same for both regions of interest

    TCMuon      muProbe;
    TCElectron  eleProbe; 
    UInt_t nMuProbes    = muProbes.size();
    UInt_t nEleProbes   = eleProbes.size();

    // Veto events with excessive MET
    //if (recoMET->Mod() > 50.) return kTRUE;

    if (
            doQCDDileptonCR 
            && (nMuProbes == 1 || nEleProbes == 1)
            && muTags.size() == 1
       ) {

        // For description of QCD dilepton control region, see section 7.4.1 of
        // ttH note (AN-13-159).  First thing is to find the tag lepton and the
        // probe lepton. For this control region, the tag is a muon that is
        // displaced from the PV and is anti-isolated.  The probe is a lepton
        // passing loose identification requirement without any isolation
        // requirement

        tag = (TCPhysObject) muTags[0];

        // Make sure there is only one probe lepton and that it does not
        // overlap with the tag. Consider both muon and electron probes.

        bitset<2> isQCD2l; 
        if (nMuProbes == 1) {
            // Overlap between muon probe and tag ensure by selector
            muProbe = muProbes[0]; 
            if (CheckQCD2lCR(tagJets, muProbe)) 
                isQCD2l.set(0); 
        }

        if (nEleProbes >= 1) {
            // Remove electron overlap with tag muon and make sure there is
            // only one electron probe
            nEleProbes = 0;
            for (unsigned i = 0; i < eleProbes.size(); ++i) {
                TCElectron testProbe = eleProbes[i];
                if (tag.DeltaR(testProbe) > 0.1) {
                    eleProbe = testProbe;
                    ++nEleProbes;
                }
            }

            if (nEleProbes == 1) 
                if (CheckQCD2lCR(tagJets, eleProbe)) 
                    isQCD2l.set(1);
        }

        bool singleProbe = true;
        if (nMuProbes == 1 && nEleProbes == 1) {
            // Only allow events with muon and electron probes if they are
            // the same object
            if (eleProbe.DeltaR(muProbe) > 0.1)
                singleProbe = false;
        }

        if ((isQCD2l.test(0) || isQCD2l.test(1)) && singleProbe) {
            // Probe object is found and event is consistent with QCD 2l
            // control region requirements.
            isTP    = true;
            crType  = "QCD2l";

            // Apply additional tight selection requirement to probe leptons
            if (isQCD2l.test(1)) {
                if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                    elPass = true;
                else //if (!selector->ElectronMVA(&eleProbe) || eleProbe.IdMap("IsoRel") > 0.15)
                    elFake = true;
            }

            if (isQCD2l.test(0)) {
                if (muProbe.IdMap("IsoRel") < 0.12)
                    muPass = true;
                else //if (muProbe.IdMap("IsoRel") > 0.20)
                    muFake = true;
            }
        }
    } 

    if (
            doZPlusJetCR
            && leptons.size() >= 2
            && nEleProbes + nMuProbes > 0
       ) {
        // For description of Z plus jet control region, see section 7.4.1 of
        // ttH note (AN-13-159).  The main difference between this control
        // region and the dilepton control region is the replacement of
        // b-tagged muon tag with a dilepton tag.  The dilepton is required to
        // be consistent with a Z and a low MET requirement is applied to
        // remove contamination from real trilepton events originating from
        // WZ->3l+nu events.

        DoZTag(leptons);
        if (zTagged) {
            tag = TCPhysObject(ZP4, 0, "Zed");
            tag.SetIdMap("IsoRel", 0.);

            // Remove probes that overlap with leptons used in Z reconstruction
            nMuProbes = 0;
            for (unsigned i = 0; i < muProbes.size(); ++i) {
                TCMuon testProbe = muProbes[i];
                if (lep1.DeltaR(testProbe) > 0.1 && lep2.DeltaR(testProbe) > 0.1) {
                    muProbe = testProbe;
                    ++nMuProbes;
                }
            }

            nEleProbes = 0;
            for (unsigned i = 0; i < eleProbes.size(); ++i) {
                TCElectron testProbe = eleProbes[i];
                if (lep1.DeltaR(testProbe) > 0.1 && lep2.DeltaR(testProbe) > 0.1) {
                    eleProbe = testProbe;
                    ++nEleProbes;
                }
            }

            // Only allow events with muon and electron probes to proceed
            // if they are the same object
            bool singleProbe = true;
            if (nEleProbes == 1 && nMuProbes == 1) {
                if (eleProbe.DeltaR(muProbe) > 0.1)
                    singleProbe = false;
            }

            if ((nEleProbes == 1 || nMuProbes == 1) && singleProbe) {
                // Probe object is found and event is consistent with QCD 2l
                // control region requirements. Now fill histograms for
                // parameterizing fake rates by pt and eta.

                if (nEleProbes == 1) {
                    if (CheckZPlusJetCR(eleProbe)) {
                        isTP    = true;
                        crType  = "ZPlusJet";
                    }                 
                    if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                        elPass = true;
                    else //if (!selector->ElectronMVA(&eleProbe) || eleProbe.IdMap("IsoRel") > 0.15)
                        elFake = true;
                } 

                if (nMuProbes == 1) {
                    if (CheckZPlusJetCR(muProbe)) {
                        isTP    = true;
                        crType  = "ZPlusJet";
                    } 
                    if (muProbe.IdMap("IsoRel") < 0.12)
                        muPass = true;
                    else //if (muProbe.IdMap("IsoRel") > 0.20)
                        muFake = true;
                } 
            }
        }
    }

    if (doAntiIso3l && crType == "None") {
        vector<TCMuon>     muonsNoIso     = selector->GetSelectedMuons("tight_id");
        vector<TCElectron> electronsNoIso = selector->GetSelectedElectrons("loose_id");

        sort(muonsNoIso.begin(), muonsNoIso.end(), P4SortCondition);
        sort(electronsNoIso.begin(), electronsNoIso.end(), P4SortCondition);

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

        if (leptonsAntiIso.size() == 2 && (muonsNoIso.size() == 1 || electronsNoIso.size() >= 1)) {

            // A tag for this CR exists if the leading two leptons are
            // anti-isolated the probe is then the trailing (third in pt) lepton

            tag = leptonsAntiIso[0]; // No clear how to define this for this case

            // Find probes.  Ensure that electrons don't overlap with tag leptons
            nMuProbes = 0;
            if (muonsNoIso.size() == 1) {
                float muISO = muonsNoIso[0].IdMap("IsoRel");
                if (
                        true
                        //(leptonsAntiIso[0] + muonsNoIso[0]).M() > 12 
                        //(leptonsAntiIso[1] + muonsNoIso[0]).M() > 12 
                        //&& muISO < 0.6 && !(muISO > 0.12 && muISO < 0.2)
                   ) {
                    muProbe = muonsNoIso[0];
                    ++nMuProbes;
                }
            }

            nEleProbes = 0;
            for (unsigned i = 0; i < electronsNoIso.size(); ++i) {
                float eleISO = electronsNoIso[i].IdMap("IsoRel");
                if (
                        electronsNoIso[i].DeltaR(leptonsAntiIso[0]) > 0.3
                        && electronsNoIso[i].DeltaR(leptonsAntiIso[1]) > 0.3
                        //&& (leptonsAntiIso[0] + electronsNoIso[0]).M() > 20. 
                        //&& eleISO < 0.9 && !(eleISO > 0.15 && eleISO < 0.2)
                   ) {
                    eleProbe = electronsNoIso[i];
                    ++nEleProbes;
                }
            }

            bool singleProbe = true;
            //if (nEleProbes == 1 && nMuProbes == 1) {
            //    if (eleProbe.DeltaR(muProbe) > 0.1)
            //        singleProbe = false;
            //}

            if ((nEleProbes == 1 || nMuProbes == 1) && singleProbe) {
                // Probe object is found and event is consistent with QCD 2l
                // control region requirements. Now fill histograms for
                // parameterizing fake rates by pt and eta.

                isTP    = true;
                crType  = "AntiIso3l";
                jets    = CleanJetOverlap(leptonsAntiIso);

                if (nEleProbes == 1) {
                    if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                        elPass = true;
                    else //if (!selector->ElectronMVA(&eleProbe) || eleProbe.IdMap("IsoRel") > 0.15)
                        elFake = true;
                } 

                if (nMuProbes == 1) {
                    histManager->SetDirectory(crType + "/" + suffix);
                    histManager->Fill1DHist((leptonsAntiIso[0] + muProbe).M(), "h1_Tagr1MuProbeMass", "M_{tag1,probe};M_{tag1,probe};Entries / 4 GeV", 50, 0., 200.);
                    histManager->Fill1DHist(leptonsAntiIso[0].DeltaR(muProbe), "h1_MuProbeTag1DeltaR", "#Delta R(tag1, probe);#Delta R(tag1, probe);Enries", 50, 0., 5.);
                    histManager->Fill1DHist((leptonsAntiIso[1] + muProbe).M(), "h1_Tag2MuProbeMass", "M_{tag2,probe};M_{tag2,probe};Entries / 4 GeV", 50, 0., 200.);
                    histManager->Fill1DHist(leptonsAntiIso[1].DeltaR(muProbe), "h1_MuProbeTag2DeltaR", "#Delta R(tag2, probe);#Delta R(tag2, probe);Enries", 50, 0., 5.);

                    if (muProbe.IdMap("IsoRel") < 0.12)
                        muPass = true;
                    else //if (muProbe.IdMap("IsoRel") > 0.20)
                        muFake = true;
                }
            } 
        }
    }


    // Verify that a tag/probe pair is found is found and require that there is
    // only one tight lepton and no more than one b-jet.  Fill tag and
    // denominator/probe distributions
    if (!isTP) return kTRUE;

    // Use MC to estimate contamination.  Probes in MC should be required to come from a real prompt lepton
    if (!isRealData) {
        bool muMatched = false;
        bool elMatched = false;
        for (unsigned i = 0; i < gLeptons.size(); ++i) {

            if (nMuProbes == 1 && !muMatched) {
                if (gLeptons[i].DeltaR(muProbe) < 0.1) 
                    muMatched = true;
            }
            if (nEleProbes == 1 && !elMatched) {
                if (gLeptons[i].DeltaR(eleProbe) < 0.1) 
                    elMatched = true;
            }

            if (muMatched && elMatched) break;
        }

        if (!muMatched) 
            nMuProbes = 0; 
        if (!elMatched) 
            nEleProbes = 0; 

        if (!muMatched && !elMatched) return kTRUE;
    }

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

    //Remove probes and tag from jet collection
    for (unsigned i = 0; i < jets.size(); ++i) {
        if (crType == "QCD2l" && jets[i].DeltaR(tag) < 0.3) {
            jets.erase(jets.begin()+i);
            --i;
            continue;
        }
        if (nMuProbes == 1) {
            if (jets[i].DeltaR(muProbe) < 0.3) {
                jets.erase(jets.begin()+i);
                --i;
                continue;
            }
        }
        if (nEleProbes == 1) {
            if (jets[i].DeltaR(eleProbe) < 0.3) {
                jets.erase(jets.begin()+i);
                --i;
                continue;
            }
        }
    }

    histManager->SetDirectory(crType + "/" + suffix);
    histManager->Fill1DHist(leptons.size(),
            "h1_leptonMult", "lepton multiplicity; N_{leptons}; Entries / bin", 6, -0.5, 5.5);
    histManager->Fill1DHist(jets.size(),
            "h1_jetMult", "jet multiplicity; N_{jets}; Entries / bin", 10, -0.5, 9.5);
    histManager->Fill1DHist(bJetsM.size(),
            "h1_bJetMediumMult", "b-jet multiplicity (medium wp); N_{b-jet}; Entries / bin", 10, -0.5, 9.5);
    histManager->Fill1DHist(bJetsL.size(),
            "h1_bJetLooseMult", "b-jet multiplicity (loose wp); N_{b-jet}; Entries / bin", 10, -0.5, 9.5);


    // Do denominator and numerator histograms
    if (crType != "None") {
        if (nEleProbes == 1) {
            FillDenominatorHists(crType, eleProbe);
            FillJetFlavorHists(crType, eleProbe);
            if (elPass)
                FillNumeratorHists(crType, eleProbe);
            else if (elFake)
                FillClosureHists(crType, eleProbe);

        }
        if (nMuProbes == 1) {
            FillDenominatorHists(crType, muProbe);
            FillJetFlavorHists(crType, muProbe);
            if (muPass)
                FillNumeratorHists(crType, muProbe);
            else if (muFake)
                FillClosureHists(crType, muProbe);
        }
    } 

    return kTRUE;
}

void fakeAnalyzer::Terminate()
{
    histManager->SetDirectory("inclusive/" + suffix);

    histManager->SetWeight(eventCount[0]);
    histManager->Fill1DHist(1,
            "h1_YieldByCut", ";cut;Entries / bin", 2, 0.5, 2.5);

    histManager->SetWeight(eventCount[1]);
    histManager->Fill1DHist(2,
            "h1_YieldByCut", ";cut;Entries / bin", 2, 0.5, 2.5);

    cout << "\nFake rate estimation finished." << endl;

    // Clean-up and write to output files
    selector->Delete();
    triggerSelector->Delete();
    histManager->Delete();

    histoFile->Write();
    histoFile->Close();  
}

void fakeAnalyzer::DoZTag(vObj& leptons)
{
    // Reset OS variables for each event of interest//
    zTagged      = false;
    dileptonMass = 0.;

    float zCandidateMass = 0.;
    for (unsigned i = 0; i < leptons.size(); ++i) {
        for (unsigned j = leptons.size()-1; j > i; --j) {

            //Find lepton pair that is consistent with a Z in +- 15 GeV window
            if (
                    fabs((leptons[i] + leptons[j]).M() - 91.2) < 15 
                    && leptons[i].Type() == leptons[j].Type()
                    && leptons[i].Charge() != leptons[j].Charge()
               ) {
                zTagged = true;
                zCandidateMass = (leptons[i] + leptons[j]).M();
            }

            // Select the pairing that is closest to the mass of the Z.
            if (zTagged && fabs(dileptonMass - 91.2) > fabs(zCandidateMass - 91.2)) {
                dileptonMass = zCandidateMass;
                lep1 = leptons[i];
                lep2 = leptons[j];
                ZP4 = (TLorentzVector)(leptons[i] + leptons[j]);

                if (leptons.size() == 3) {
                    lep3 = leptons[3-i-j];
                }
            }
        }
    }
}

void fakeAnalyzer::FillDenominatorHists(string cat, TCPhysObject& probe)
{
    histManager->SetDirectory(cat + "/" + suffix);

    // Sanity check plots
    histManager->Fill1DHist(recoMET->Mod(), "h1_Met", "MET;MET;Entries / 4 GeV", 25, 0., 100.);
    if (probe.Type() == "electron") {
        // fake rate measurement plots
        histManager->Fill1DHist(tag.Pt(),
                "h1_EleTagLepPt", "tag lepton p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(tag.Eta(),
                "h1_EleTagLepEta", "tag lepton #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(tag.IdMap("IsoRel"),
                "h1_EleTagIsoRel", "tag lepton ISO_{rel};ISO_{rel};Entries", 80, 0., 4.);
        histManager->Fill1DHist(fabs(tag.Dz(selectedVtx)),
                "h1_EleTagDz", "tag  d_{z};d_{z};Entries", 50, 0., 2.);
        histManager->Fill1DHist(fabs(tag.Dxy(selectedVtx)),
                "h1_EleTagDxy", "tag  d_{xy};d_{xy};Entries", 50, 0., 2.);

        histManager->Fill1DHist((tag + probe).M(),
                "h1_TagEleProbeMass", "M_{tag,probe};M_{tag,probe};Entries / 4 GeV", 50, 0., 200.);
        histManager->Fill1DHist(tag.DeltaR(probe),
                "h1_TagEleProbeDeltaR", "#Delta R(tag,probe);#Delta R(tag,probe);Entries / 3 GeV", 50, 0., 5.);
        histManager->Fill2DHist(probe.IdMap("IsoRel"), probe.Pt(), 
                "h2_EleProbePtVsIso", "probe lepton p_{T} vs ISO_{rel};ISO_{rel};p_{T}", 60, 0., 6., 5, 10., 35.);
    }
    if (probe.Type() == "muon") {
        histManager->Fill1DHist(tag.Pt(),
                "h1_MuTagLepPt", "tag lepton p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(tag.Eta(),
                "h1_MuTagLepEta", "tag lepton #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(tag.IdMap("IsoRel"),
                "h1_MuTagIsoRel", "tag lepton ISO_{rel};ISO_{rel};Entries", 80, 0., 4.);
        histManager->Fill1DHist(fabs(tag.Dz(selectedVtx)),
                "h1_MuTagDz", "tag  d_{z};d_{z};Entries", 50, 0., 2.);
        histManager->Fill1DHist(fabs(tag.Dxy(selectedVtx)),
                "h1_MuTagDxy", "tag  d_{xy};d_{xy};Entries", 50, 0., 2.);
        histManager->Fill1DHist((tag + probe).M(),
                "h1_TagMuProbeMass", "M_{tag,probe};M_{tag,probe};Entries / 4 GeV", 50, 0., 200.);
        histManager->Fill1DHist(tag.DeltaR(probe),
                "h1_TagMuProbeDeltaR", "#Delta R(tag,probe);#Delta R(tag,probe);Entries / 3 GeV", 50, 0., 5.);

        histManager->Fill2DHist(probe.IdMap("IsoRel"), probe.Pt(), 
                "h2_MuProbePtVsIso", "probe lepton p_{T} vs ISO_{rel};ISO_{rel};p_{T}", 60, 0., 6., 5, 10., 35.);
    }

    string lepType = probe.Type();
    if (lepType == "muon") {
        // fake rate measurement plots
        histManager->Fill1DHist(probe.Pt(),
                "h1_MuProbeLepPt", "probe muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_MuProbeLepEta", "probe muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_MuDenomIsoRel", "probe muon IsoRel;IsoRel;Entries", 80, 0., 4.);
        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_MuDenomPt", "probe muon p_{T};p_{T};Entries / bin", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                "h1_MuDenomEta", "probe muon #eta;#eta;Entries / bin", 2, etaBinsMu);
        histManager->Fill2DHistUnevenBins(probe.Pt(), fabs(probe.Eta()),
                "h2_MuDenom", "probe muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);
        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_MuDenomMet", "probe muon Met;Met;Entries", nMetBins, metBins);
        histManager->Fill2DHist(jets.size(), probe.IdMap("IsoRel"),
                "h2_MuProbeIsoRelVsJetMult", "ISO_{rel} vs. jet multiplicity;N_{jets};ISO_{rel}", 6, -0.5, 5.5, 5, 0.15, 1.15);

    } else if (lepType == "electron") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_EleProbeLepPt", "probe electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_EleProbeLepEta", "probe electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_EleDenomIsoRel", "probe electron IsoRel;IsoRel;Entries", 80, 0., 4.);
        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_EleDenomPt", "probe electron p_{T};p_{T};Entries / bin", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(fabs(probe.Eta()),
                "h1_EleDenomEta", "probe electron #eta;#eta;Entries / bin", 3, etaBinsEle);
        histManager->Fill2DHistUnevenBins(probe.Pt(), fabs(probe.Eta()),
                "h2_EleDenom", "probe electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);
        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_EleDenomMet", "probe electron Met;Met;Entries", nMetBins, metBins);
        histManager->Fill2DHist(jets.size(), probe.IdMap("IsoRel"),
                "h2_EleProbeIsoRelVsJetMult", "ISO_{rel} vs. jet multiplicity;N_{jets};ISO_{rel}", 6, -0.5, 5.5, 5, 0.2, 1.2);
    }
}

void fakeAnalyzer::FillNumeratorHists(string cat, TCPhysObject& probeLep)
{
    histManager->SetDirectory(cat + "/" + suffix);

    string lepType = probeLep.Type();
    if (lepType == "muon") {
        // Apply tight muon selection to probe muon
        histManager->Fill1DHist(probeLep.Pt(),
                "h1_MuPassLepPt", "pass muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probeLep.Eta(),
                "h1_MuPassLepEta", "pass muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(probeLep.IdMap("IsoRel"),
                "h1_MuNumerIsoRel", "pass electron IsoRel;IsoRel;Entries", 40, 0., 0.20);

        histManager->Fill1DHistUnevenBins(probeLep.Pt(),
                "h1_MuNumerPt", "pass muon p_{T};p_{T};Entries", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(fabs(probeLep.Eta()),
                "h1_MuNumerEta", "pass muon #eta;#eta;Entries", 2, etaBinsMu);
        histManager->Fill2DHistUnevenBins(probeLep.Pt(), fabs(probeLep.Eta()),
                "h2_MuNumer", "pass muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_MuNumerMet", "pass muon Met;Met;Entries", nMetBins, metBins);
    } else if (lepType == "electron") {
        histManager->Fill1DHist(probeLep.Pt(),
                "h1_ElePassLepPt", "pass electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probeLep.Eta(),
                "h1_ElePassLepEta", "pass electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(probeLep.IdMap("IsoRel"),
                "h1_EleNumerIsoRel", "pass electron IsoRel;IsoRel;Entries", 40, 0., 0.20);

        histManager->Fill1DHistUnevenBins(probeLep.Pt(),
                "h1_EleNumerPt", "pass electron p_{T};p_{T};Entries / 3 GeV", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(fabs(probeLep.Eta()),
                "h1_EleNumerEta", "pass electron #eta;#eta;Entries / 3 GeV", 3, etaBinsEle);
        histManager->Fill2DHistUnevenBins(probeLep.Pt(), fabs(probeLep.Eta()),
                "h2_EleNumer", "pass electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_EleNumerMet", "pass electron Met;Met;Entries", nMetBins, metBins);
    }
}

void fakeAnalyzer::FillClosureHists(string category, TCPhysObject& probe)
{
    histManager->SetDirectory(category + "/" + suffix);

    string categories[3] = {"QCD2l", "ZPlusJet", "AntiIso3l"};
    for (unsigned i = 0; i < 3; ++i) {
        string cat = categories[i];
        float fakeWeight = weighter->GetFakeWeight(probe, cat);
        histManager->SetWeight(fakeWeight);

        //cout << cat << ":" << fakeWeight << ":(" << probe.Type() << "; " << probe.Pt() << "; " << probe.Eta() << ")" << endl;

        if (probe.Type() == "electron") {
            histManager->Fill1DHist(probe.Pt(),
                    "h1_ElePtClosure_" + cat, "electron p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
            histManager->Fill1DHist(fabs(probe.Eta()),
                    "h1_EleEtaClosure_" + cat, "electron #eta (closure);#eta;Entries / 3 GeV", 5, 0., 2.5);
            histManager->Fill1DHistUnevenBins(probe.Pt(),
                    "h1_EleUnevenPtClosure_" + cat, "electron p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
            histManager->Fill1DHist(probe.IdMap("IsoRel"),
                    "h1_EleIsoRelClosure_" + cat, "pass electron IsoRel;IsoRel;Entries", 20, 0., 1.0);
            histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                    "h1_EleMetClosure_" + cat, "pass electron IsoRel;IsoRel;Entries", nMetBins, metBins);
        } else if (probe.Type() == "muon") {
            histManager->Fill1DHist(probe.Pt(),
                    "h1_MuPtClosure_" + cat, "muon p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
            histManager->Fill1DHist(fabs(probe.Eta()),
                    "h1_MuEtaClosure_" + cat, "muon #eta (closure);#eta;Entries / 3 GeV", 5, 0., 2.5);
            histManager->Fill1DHistUnevenBins(probe.Pt(),
                    "h1_MuUnevenPtClosure_" + cat, "muon p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
            histManager->Fill1DHist(probe.IdMap("IsoRel"),
                    "h1_MuIsoRelClosure_" + cat, "pass muon IsoRel;IsoRel;Entries", 20, 0., 1.0);
            histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                    "h1_MuMetClosure_" + cat, "pass muon IsoRel;IsoRel;Entries", nMetBins, metBins);
        }
        histManager->SetWeight(1.);
    }
}

void fakeAnalyzer::FillJetFlavorHists(string cat, TCPhysObject& lepton) 
{
    histManager->SetDirectory(cat + "/" + suffix);

    unsigned jetIndex;
    bool jetMatched = false;
    for (unsigned i = 0; i < jets.size(); ++i) {
        if (lepton.DeltaR(jets[i]) < 0.25) {
            jetMatched = true;
            jetIndex = i;
            continue;
        }
    }

    if (jetMatched) {
        //cout << jets[jetIndex].BDiscriminatorMap("CSV") << endl;
        if (lepton.Type() == "muon")
            histManager->Fill1DHist(jets[jetIndex].BDiscriminatorMap("CSV"),
                    "h1_MatchedMuJetBDiscr", "matched #mu-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
        if (lepton.Type() == "electron")
            histManager->Fill1DHist(jets[jetIndex].BDiscriminatorMap("CSV"),
                    "h1_MatchedEleJetBDiscr", "matched e-jet b discriminator;CSV;Entries / bin", 50, -1., 1.5);
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

    histManager->Fill1DHist(fabs(tpDeltaPhi),
            "h1_TagProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
    histManager->Fill1DHist(fabs(tpBalance),
            "h1_TagProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);

    if (fabs(tpDeltaPhi) < 2. || tpBalance > 1 || !jetMatched) 
        return false;
    else
        return true;
}

bool fakeAnalyzer::CheckZPlusJetCR(TCPhysObject& probe)
{
    // Check tag/probe pair is back-to-back
    Float_t tpDeltaPhi  = tag.DeltaPhi(probe);
    Float_t tpBalance   = probe.Pt()/tag.Pt(); 

    histManager->Fill1DHist(fabs(tpDeltaPhi),
            "h1_TagProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
    histManager->Fill1DHist(fabs(tpBalance),
            "h1_TagProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);

    //if (fabs(tpDeltaPhi) < 2.5 || tpBalance > 1 || recoMET->Mod() > 50)
    //    return false;
    //else
    return true;
}

void fakeAnalyzer::GenMatcher(vObj& leptons, vector<TCGenParticle>& gLeptons)
{
    for (unsigned i = 0; i < leptons.size(); ++i) {
        //TCGenParticle matchedGen;
        bool matched = false; 
        for (unsigned j = 0; j < gLeptons.size(); ++j) {
            histManager->Fill1DHist(leptons[i].DeltaR(gLeptons[j]),
                    "h1_DeltaRGenLepton", "#Delta R (gen, reco);#Delta R (gen, reco);Entries", 100, 0., 1.);
            histManager->Fill1DHist((leptons[i].Pt() - gLeptons[j].Pt())/gLeptons[j].Pt(),
                    "h1_DeltaPtGenLepton", "#Delta p_{T} / p_{T,gen};#Delta p_{T} / p_{T,gen} ;Entries", 55, -1., 10.);
            if (
                    leptons[i].DeltaR(gLeptons[j]) < 0.3 || fabs(leptons[i].Pt() - gLeptons[j].Pt())/gLeptons[j].Pt() < 0.1
                    //|| (leptons[i].Type() == "muon" && fabs(gLeptons[j].GetPDGId()) == 13) 
                    //|| (leptons[i].Type() == "electron" && fabs(gLeptons[j].GetPDGId()) == 11)
               ) {
                //matchedGen = gLeptons[i];
                matched = true;
                break;
            }
        }

        if (!matched) {
            float leptonIso = leptons[i].IdMap("IsoRel");
            if (leptons[i].Type() == "muon" && leptonIso > 0.)
                histManager->Fill1DHist(leptonIso, "h1_FakeMuonIso", "fake muon ISO_{rel};ISO_{rel,#mu};Entries", 80, 0., 4.);
            if (leptons[i].Type() == "electron" && leptonIso > 0.)
                histManager->Fill1DHist(leptonIso, "h1_FakeElectronIso", "fake electron ISO_{rel};ISO_{rel,#mu};Entries", 80, 0., 4.);
        }
    }
}

vector<TCJet> fakeAnalyzer::CleanJetOverlap(vObj& objects)
{
    vector<TCJet> cleanJets;
    for (unsigned i = 0; i < jets.size(); ++i) {

        bool overlap = false;
        for (unsigned j = 0; j < objects.size(); ++j) {
            if (jets[i].DeltaR(objects[j]) < 0.3) {
                overlap = true;
                break;
            }
        }
        if (!overlap) 
            cleanJets.push_back(jets[i]);
    }
    return cleanJets;
}
