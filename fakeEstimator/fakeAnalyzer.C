#define fakeAnalyzer_cxx

#include "fakeAnalyzer.h"

using namespace std;

/////////////////
//Analysis cuts//
/////////////////

const bool  doQCDDileptonCR = true;
const bool  doZPlusJetCR    = true;
const bool  doFakeCR        = true;

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
    histoFile->mkdir("2l", "2l");
    histoFile->GetDirectory("2l", "2l")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("3l", "3l");
    histoFile->GetDirectory("3l", "3l")->mkdir(suffix.c_str(), suffix.c_str());

    histoFile->mkdir("QCD2l_inclusive", "QCD2l_inclusive");
    histoFile->GetDirectory("QCD2l_inclusive", "QCD2l_inclusive")->mkdir(suffix.c_str(), suffix.c_str());

    histoFile->mkdir("ZPlusJet_inclusive", "ZPlusJet_inclusive");
    histoFile->GetDirectory("ZPlusJet_inclusive", "ZPlusJet_inclusive")->mkdir(suffix.c_str(), suffix.c_str());

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
    vector<TCMuon>      muProbes    = selector->GetSelectedMuons("QCD2l_CR_probe");
    vector<TCElectron>  eleProbes   = selector->GetSelectedElectrons("QCD2l_CR_probe");

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
    string crType  = "none";

    // Prepare control regions for FR estimation...  
    // The control region is largely defined by the tag, but we'll first check
    // for any probes since they will be the same for both regions of interest

    TCMuon      muProbe;
    TCElectron  eleProbe; 
    UInt_t nMuProbes    = muProbes.size();
    UInt_t nEleProbes   = eleProbes.size();

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
            if (isQCD2l.test(0)) 
                if (muProbe.IdMap("IsoRel") < 0.12)
                    muPass = true;

            if (isQCD2l.test(1))
                if (selector->ElectronMVA(&eleProbe) && eleProbe.IdMap("IsoRel") < 0.15)
                    elPass = true;
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
                } 

                if (nMuProbes == 1) {
                    if (CheckZPlusJetCR(muProbe)) {
                        isTP    = true;
                        crType  = "ZPlusJet";
                    } 
                    if (muProbe.IdMap("IsoRel") < 0.12)
                        muPass = true;
                } 
            }
        }
    }

    if (doFakeCR) {
        vector<TCMuon>     muonsNoIso     = selector->GetSelectedMuons("tight_id");
        vector<TCElectron> electronsNoIso = selector->GetSelectedElectrons("loose_id");
        FakeCR(muonsNoIso, electronsNoIso, gLeptons);
    }

    // Verify that a tag/probe pair is found is found and require that there is
    // only one tight lepton and no more than one b-jet.  Fill tag and
    // denominator/probe distributions

    if (isTP) {
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

        //Remove probe and tag from jet collection
        for (unsigned i = 0; i < jets.size(); ++i) {
            if (jets[i].DeltaR(tag) < 0.3) {
                jets.erase(jets.begin()+i);
                continue;
            }
            if (nMuProbes == 1) {
                if (jets[i].DeltaR(muProbe) < 0.3) {
                    jets.erase(jets.begin()+i);
                    continue;
                }
            }
            if (nEleProbes == 1) {
                if (jets[i].DeltaR(eleProbe) < 0.3) {
                    jets.erase(jets.begin()+i);
                    continue;
                }
            }
        }

        histManager->SetDirectory(crType + "_inclusive/" + suffix);
        histManager->Fill1DHist(leptons.size(),
                "h1_leptonMult", "lepton multiplicity; N_{leptons}; Entries / bin", 6, -0.5, 5.5);
        histManager->Fill1DHist(jets.size(),
                "h1_jetMult", "jet multiplicity; N_{jets}; Entries / bin", 10, -0.5, 9.5);
        histManager->Fill1DHist(bJetsM.size(),
                "h1_bJetMediumMult", "b-jet multiplicity (medium wp); N_{b-jet}; Entries / bin", 10, -0.5, 9.5);
        histManager->Fill1DHist(bJetsL.size(),
                "h1_bJetLooseMult", "b-jet multiplicity (loose wp); N_{b-jet}; Entries / bin", 10, -0.5, 9.5);
    } else
        return kTRUE;

    if (recoMET->Mod() > 50.) return kTRUE;

    // Do denominator and numerator histograms
    if (crType == "QCD2l") {
        if (nEleProbes == 1) {
            FillDenominatorHists(crType + "_inclusive", eleProbe);
            FillJetFlavorHists(crType, eleProbe);

            if (elPass)
                FillNumeratorHists(crType + "_inclusive", eleProbe);
            else
                FillClosureHists(crType + "_inclusive", eleProbe);
        }
        if (nMuProbes == 1) {
            FillDenominatorHists(crType + "_inclusive", muProbe);
            FillJetFlavorHists(crType, muProbe);

            if (muPass)
                FillNumeratorHists(crType + "_inclusive", muProbe);
            else
                FillClosureHists(crType + "_inclusive", muProbe);
        }
    } else if (crType == "ZPlusJet") {
        if (nEleProbes == 1) {
            FillDenominatorHists(crType + "_inclusive", eleProbe);
            FillJetFlavorHists(crType, eleProbe);

            if (elPass)
                FillNumeratorHists(crType + "_inclusive", eleProbe);
            else
                FillClosureHists(crType + "_inclusive", eleProbe);
        }
        if (nMuProbes == 1) {
            FillDenominatorHists(crType + "_inclusive", muProbe);
            FillJetFlavorHists(crType, muProbe);
            if (muPass)
                FillNumeratorHists(crType + "_inclusive", muProbe);
            else
                FillClosureHists(crType + "_inclusive", muProbe);
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

            //Find lepton pair that is consistent with a Z in +- 10 GeV window
            if (
                    fabs((leptons[i] + leptons[j]).M() - 91.2) < 10 
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
        histManager->Fill1DHist((tag + probe).M(),
                "h1_TagEleProbeMass", "M_{tag,probe};M_{tag,probe};Entries / 3 GeV", 50, 0., 150.);
        histManager->Fill2DHist(probe.Pt(), probe.IdMap("IsoRel"),
                "h2_EleProbePtVsIso", "probe lepton p_{T} vs ISO_{rel};p_{T};ISO_{rel}", 50, 0., 150., 25, 0., 1.);
    }
    if (probe.Type() == "muon") {
        histManager->Fill1DHist((tag + probe).M(),
                "h1_TagMuProbeMass", "M_{tag,probe};M_{tag,probe};Entries / 3 GeV", 50, 0., 150.);
        histManager->Fill2DHist(probe.Pt(), probe.IdMap("IsoRel"),
                "h2_MuProbePtVsIso", "probe lepton p_{T} vs ISO_{rel};p_{T};ISO_{rel}", 50, 0., 150., 25, 0., 1.);
    }

    string lepType = probe.Type();
    if (lepType == "muon") {
        // fake rate measurement plots
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

        histManager->Fill1DHist(probe.Pt(),
                "h1_MuProbeLepPt", "probe muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_MuProbeLepEta", "probe muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_MuDenomIsoRel", "probe muon IsoRel;IsoRel;Entries", 80, 0., 4.);

        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_MuDenomPt", "probe muon p_{T};p_{T};Entries / bin", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(probe.Eta(),
                "h1_MuDenomEta", "probe muon #eta;#eta;Entries / bin", 2, etaBinsMu);
        histManager->Fill2DHistUnevenBins(probe.Pt(), probe.Eta(),
                "h2_MuDenom", "probe muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_MuDenomMet", "probe muon Met;Met;Entries", nMetBins, metBins);

        histManager->Fill2DHist(jets.size(), probe.IdMap("IsoRel"),
                "h2_MuProbeIsoRelVsJetMult", "ISO_{rel} vs. jet multiplicity;N_{jets};ISO_{rel}", 6, -0.5, 5.5, 5, 0.15, 1.15);

    } else if (lepType == "electron") {
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

        histManager->Fill1DHist(probe.Pt(),
                "h1_EleProbeLepPt", "probe electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_EleProbeLepEta", "probe electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);
        histManager->Fill1DHist(probe.IdMap("IsoRel"),
                "h1_EleDenomIsoRel", "probe electron IsoRel;IsoRel;Entries", 80, 0., 4.);

        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_EleDenomPt", "probe electron p_{T};p_{T};Entries / bin", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(probe.Eta(),
                "h1_EleDenomEta", "probe electron #eta;#eta;Entries / bin", 3, etaBinsEle);
        histManager->Fill2DHistUnevenBins(probe.Pt(), probe.Eta(),
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
        histManager->Fill1DHistUnevenBins(probeLep.Eta(),
                "h1_MuNumerEta", "pass muon #eta;#eta;Entries", 2, etaBinsMu);
        histManager->Fill2DHistUnevenBins(probeLep.Pt(), probeLep.Eta(),
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
        histManager->Fill1DHistUnevenBins(probeLep.Eta(),
                "h1_EleNumerEta", "pass electron #eta;#eta;Entries / 3 GeV", 3, etaBinsEle);
        histManager->Fill2DHistUnevenBins(probeLep.Pt(), probeLep.Eta(),
                "h2_EleNumer", "pass electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_EleNumerMet", "pass electron Met;Met;Entries", nMetBins, metBins);
    }
}

void fakeAnalyzer::FillClosureHists(string cat, TCPhysObject& probe)
{
    //histManager->SetWeight(weighter->GetFakeWeight(probe, cat));
    histManager->SetDirectory(cat + "/" + suffix);

    if (probe.Type() == "electron") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_ElePtClosure", "electron p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(fabs(probe.Eta()),
                "h1_EleEtaClosure", "electron #eta (closure);#eta;Entries / 3 GeV", 5, 0., 2.5);
        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_EleUnevenPtClosure", "electron p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
    } else if (probe.Type() == "muon") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_MuPtClosure", "muon p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(fabs(probe.Eta()),
                "h1_MuEtaClosure", "muon #eta (closure);#eta;Entries / 3 GeV", 5, 0., 2.5);
        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_MuUnevenPtClosure", "muon p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
    }
    histManager->SetWeight(1.);
}

void fakeAnalyzer::FillJetFlavorHists(string cat, TCPhysObject& lepton) 
{
    histManager->SetDirectory(cat + "_inclusive/" + suffix);

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

    if (fabs(tpDeltaPhi) < 2.5 || tpBalance > 1 || !jetMatched) 
        return false;

    // Correct for prompt lepton contamination for low pt probes; not really relevant in current analysis
    if (probe.Pt() < 10 && recoMET->Mod() > 15) 
        return false;

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

    if (fabs(tpDeltaPhi) < 2.5 || tpBalance > 1 || recoMET->Mod() > 50)
        return false;
    else
        return true;
}

void fakeAnalyzer::FakeCR(vector<TCMuon>& muons, vector<TCElectron>& electrons, vector<TCGenParticle>& gLeptons)
{
    vObj leptonsNoIso;
    leptonsNoIso.insert(leptonsNoIso.end(), muons.begin(), muons.end());
    leptonsNoIso.insert(leptonsNoIso.end(), electrons.begin(), electrons.end());
    sort(leptonsNoIso.begin(), leptonsNoIso.end(), P4SortCondition);

    if (leptonsNoIso.size() == 2)
        histManager->SetDirectory("2l/" + suffix);
    else if (leptonsNoIso.size() == 3)
        histManager->SetDirectory("3l/" + suffix);

    if (!isRealData) 
        GenMatcher(leptonsNoIso, gLeptons);

    if (leptonsNoIso.size() == 3 || leptonsNoIso.size() == 2) {
        vObj leptonsAntiIso, testLeptonsAntiIso, leptonsIso, testLeptonsIso;
        for (unsigned i = 0; i < leptonsNoIso.size(); ++i) {
            float leptonIso = leptonsNoIso[i].IdMap("IsoRel");

            if (leptonsNoIso[i].Type() == "muon")
                histManager->Fill1DHist(leptonIso, "h1_MuonIso_" + str(i+1), "muon ISO_{rel};ISO_{rel,#mu};Entries", 80, 0., 4.);
            else if (leptonsNoIso[i].Type() == "electron") 
                histManager->Fill1DHist(leptonIso, "h1_ElectronIso_" + str(i+1), "electron ISO_{rel};ISO_{rel,e};Entries", 80, 0., 4.);

            //cout << leptonsNoIso[i].Type() << ": " << leptonsNoIso[i].IsTriggered() << "\t";

            if (leptonsNoIso.size() == 3) {
                if (leptonIso > 0.4 && leptonsAntiIso.size() < 2) 
                    leptonsAntiIso.push_back(leptonsNoIso[i]);
                else if (leptonsAntiIso.size() == 2) 
                    testLeptonsAntiIso.push_back(leptonsNoIso[i]);

                if (
                        ((leptonsNoIso[i].Type() == "muon" && leptonIso < 0.12) 
                         || (leptonsNoIso[i].Type() == "electron" && leptonIso < 0.15)) 
                        && leptonsIso.size() < 2
                   )
                    leptonsIso.push_back(leptonsNoIso[i]);
                else if (leptonsIso.size() == 2) 
                    testLeptonsIso.push_back(leptonsNoIso[i]);

            }

            if (leptonsNoIso.size() == 2) {
                if (leptonsNoIso[0].Charge() == leptonsNoIso[1].Charge()) {
                    if (leptonIso > 0.4 && leptonsAntiIso.size() < 1) 
                        leptonsAntiIso.push_back(leptonsNoIso[i]);
                    else if (leptonsAntiIso.size() == 1) 
                        testLeptonsAntiIso.push_back(leptonsNoIso[i]);

                    if (
                            ((leptonsNoIso[i].Type() == "muon" && leptonIso < 0.12) 
                             || (leptonsNoIso[i].Type() == "electron" && leptonIso < 0.15)) 
                            && leptonsIso.size() < 1
                       ) 
                        leptonsIso.push_back(leptonsNoIso[i]);
                    else if (leptonsIso.size() == 1) 
                        testLeptonsIso.push_back(leptonsNoIso[i]);
                }
            }
        }
        //cout << endl;
        //cout << leptonsNoIso.size() << "\t" << leptonsAntiIso.size() << "\t" << testLeptons.size() << endl;

        if (leptonsAntiIso.size() >= 1 && testLeptonsAntiIso.size() == 1) {

            if (testLeptonsAntiIso[0].Type() == "muon") {
                histManager->Fill1DHist(testLeptonsAntiIso[0].IdMap("IsoRel"),
                        "h1_MuonIsoRel_AntiIso", "muon relative isolation;ISO_{rel};Entries", 80, 0., 4.);
                histManager->Fill1DHist(testLeptonsAntiIso[0].Pt(),
                        "h1_MuonPt_AntiIso", "muon p_{T};p_{T};Entries", 30, 0., 150.);

                if (testLeptonsAntiIso[0].IdMap("IsoRel") > 0.2 && testLeptonsAntiIso[0].IdMap("IsoRel") < 0.6) {
                    histManager->SetWeight(weighter->GetFakeWeight(testLeptonsAntiIso[0], "QCD2l"));
                    histManager->Fill1DHist(testLeptonsAntiIso[0].Pt(),
                            "h1_MuonPt_QCD2l_weight", "muon p_{T};p_{T};Entries", 30, 0., 150.);
                    histManager->Fill1DHist(fabs(testLeptonsAntiIso[0].Eta()),
                            "h1_MuonEta_QCD2l_weight", "muon p_{T};p_{T};Entries", 5, 0., 2.5);
                    histManager->Fill1DHist(recoMET->Mod(),
                            "h1_MuonMet_QCD2l_weight", ";MET;Entries", 20, 0., 100.);
                    histManager->Fill1DHist(jets.size(),
                            "h1_MuonJetMult_QCD2l_weight", ";N_{jets};Entries", 6, -0.5, 5.5);
                    histManager->SetWeight(1.);
                } else if (testLeptonsAntiIso[0].IdMap("IsoRel") < 0.12) {
                    histManager->Fill1DHist(testLeptonsAntiIso[0].Pt(),
                            "h1_MuonPt_QCD2l_tight", "muon p_{T};p_{T};Entries", 30, 0., 150.);
                    histManager->Fill1DHist(fabs(testLeptonsAntiIso[0].Eta()),
                            "h1_MuonEta_QCD2l_tight", "muon p_{T};p_{T};Entries", 5, 0., 2.5);
                    histManager->Fill1DHist(recoMET->Mod(),
                            "h1_MuonMet_QCD2l_tight", ";MET;Entries", 20, 0., 100.);
                    histManager->Fill1DHist(jets.size(),
                            "h1_MuonJetMult_QCD2l_tight", ";N_{jets};Entries", 6, -0.5, 5.5);
                }

            } else if (testLeptonsAntiIso[0].Type() == "electron") {
                histManager->Fill1DHist(testLeptonsAntiIso[0].IdMap("IsoRel"),
                        "h1_ElectronIsoRel_AntiIso", "electron relative isolation;ISO_{rel};Entries", 80, 0., 4.);
                histManager->Fill1DHist(testLeptonsAntiIso[0].Pt(),
                        "h1_ElectronPt_AntiIso", "electron p_{T};p_{T};Entries", 30, 0., 150.);

                if (testLeptonsAntiIso[0].IdMap("IsoRel") > 0.2) {
                    histManager->SetWeight(weighter->GetFakeWeight(testLeptonsAntiIso[0], "QCD2l"));
                    histManager->Fill1DHist(testLeptonsAntiIso[0].Pt(),
                            "h1_ElectronPt_QCD2l_weight", "electron p_{T};p_{T};Entries", 30, 0., 150.);
                    histManager->Fill1DHist(fabs(testLeptonsAntiIso[0].Eta()),
                            "h1_ElectronEta_QCD2l_weight", "electron #eta;#eta;Entries", 5, 0., 2.5);
                    histManager->Fill1DHist(recoMET->Mod(),
                            "h1_ElectronMet_QCD2l_weight", ";MET;Entries", 20, 0., 100.);
                    histManager->Fill1DHist(jets.size(),
                            "h1_ElectronJetMult_QCD2l_weight", ";N_{jets};Entries", 6, -0.5, 5.5);
                    histManager->SetWeight(1.);
                } else if (testLeptonsAntiIso[0].IdMap("IsoRel") < 0.15) {
                    histManager->Fill1DHist(testLeptonsAntiIso[0].Pt(),
                            "h1_ElectronPt_QCD2l_tight", "electron p_{T};p_{T};Entries", 30, 0., 150.);
                    histManager->Fill1DHist(fabs(testLeptonsAntiIso[0].Eta()),
                            "h1_ElectronEta_QCD2l_tight", "electron #eta;#eta;Entries", 5, 0., 2.5);
                    histManager->Fill1DHist(recoMET->Mod(),
                            "h1_ElectronMet_QCD2l_tight", ";MET;Entries", 20, 0., 100.);
                    histManager->Fill1DHist(jets.size(),
                            "h1_ElectronJetMult_QCD2l_tight", ";N_{jets};Entries", 6, -0.5, 5.5);

                }
            }
        }

        if (leptonsIso.size() >= 1 && testLeptonsIso.size() > 0) {
            for (unsigned i = 0; i < leptonsIso.size(); ++i) {
                histManager->Fill1DHist(leptonsIso[i].Pt(),
                        "h1_LeptonPt" + str(i+1) + "_Iso", "muon p_{T};p_{T};Entries", 30, 0., 150.);
            }

            if (testLeptonsIso[0].Type() == "muon") {
                histManager->Fill1DHist(testLeptonsIso[0].IdMap("IsoRel"),
                        "h1_MuonIsoRel_Iso", "muon relative isolation;ISO_{rel};Entries", 80, 0., 4.);
                histManager->Fill1DHist(testLeptonsIso[0].Pt(),
                        "h1_MuonPt_Iso", "muon p_{T};p_{T};Entries", 30, 0., 150.);
            } else if (testLeptonsIso[0].Type() == "electron") {
                histManager->Fill1DHist(testLeptonsIso[0].IdMap("IsoRel"),
                        "h1_ElectronIsoRel_Iso", "electron relative isolation;ISO_{rel};Entries", 80, 0., 4.);
                histManager->Fill1DHist(testLeptonsIso[0].Pt(),
                        "h1_ElectronPt_Iso", "electron p_{T};p_{T};Entries", 30, 0., 150.);
            }
        }
    }
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
