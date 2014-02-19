#define fakeAnalyzer_cxx

#include "fakeAnalyzer.h"

using namespace std;

/////////////////
//Analysis cuts//
/////////////////

const bool  doQCDDileptonCR = true;
const bool  doZPlusJetCR    = true;

const float jetPtCut[]        = {25., 15.};
const float muPtCut[]         = {10., 3.};
const float elePtCut[]        = {10., 10.};
const float phoPtCut[]        = {10., 10.};


unsigned  nMetBins      = 10;
unsigned  nPtBins       = 6;
float     metBins[]     = {0., 10., 20., 30., 40., 50., 60., 70., 80., 100., 300.}; 
float     ptBins[]      = {5., 10., 20., 30., 45., 60, 100.}; 
float     etaBinsMu[]   = {0., 1.5, 2.5};
float     etaBinsEle[]  = {0., 0.8, 1.479, 2.5};

// Do something about these: should just have one sort condition function
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
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
    triggerSelector = new TriggerSelector("mc", "2012", *triggerNames, true);

    // Add single lepton triggers for fake rates //
    vstring triggers;
    triggers.push_back("HLT_Mu8_v");
    triggers.push_back("HLT_Mu17_v");
    triggers.push_back("HLT_Ele8_CaloIdT_TrkIdVL_v");
    triggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    triggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
    triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
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

    histoFile->mkdir("QCD2l_inclusive", "QCD2l_inclusive");
    histoFile->GetDirectory("QCD2l_inclusive", "QCD2l_inclusive")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("QCD2l_low_met", "QCD2l_low_met");
    histoFile->GetDirectory("QCD2l_low_met", "QCD2l_low_met")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("QCD2l_high_met", "QCD2l_high_met");
    histoFile->GetDirectory("QCD2l_high_met", "QCD2l_high_met")->mkdir(suffix.c_str(), suffix.c_str());

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

    if (!triggerPass) return kTRUE;

    vstring passNames = triggerSelector->GetPassNames();
    if (passNames.size() == 0) passNames.push_back("NULL");

    selector->PVSelector(primaryVtx);
    if (selector->GetSelectedPVs().size() < 1) return kTRUE;

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
    vector<TCElectron>  olElectrons = selector->GetSelectedElectrons("tight_overlap");

    vObj leptons;
    leptons.insert(leptons.end(), muons.begin(), muons.end());
    leptons.insert(leptons.end(), electrons.begin(), electrons.end());

    // Get loose leptons
    vObj extraLeptons;
    vector<TCElectron>  looseElectrons  = selector->GetSelectedElectrons("loose");
    vector<TCMuon>      looseMuons      = selector->GetSelectedMuons("loose");
    extraLeptons.insert(extraLeptons.end(), looseElectrons.begin(), looseElectrons.end());
    extraLeptons.insert(extraLeptons.end(), looseMuons.begin(), looseMuons.end());

    // Get jets
    vector<TCJet> jets      = selector->GetSelectedJets("tight");
    vector<TCJet> bJetsL    = selector->GetSelectedJets("bJetsLoose");
    vector<TCJet> bJetsM    = selector->GetSelectedJets("bJetsMedium");
    vector<TCJet> fwdJets   = selector->GetSelectedJets("forward");
    vector<TCJet> muJets    = selector->GetSelectedJets("muJets");
    vector<TCJet> eleJets   = selector->GetSelectedJets("eleJets");

    jets.insert(jets.end(), fwdJets.begin(), fwdJets.end());
    //jets.insert(jets.end(), bJetsL.begin(), bJetsL.end());
    jets.insert(jets.end(), bJetsM.begin(), bJetsM.end());
    jets.insert(jets.end(), muJets.begin(), muJets.end());
    jets.insert(jets.end(), eleJets.begin(), eleJets.end());

    vector<TCJet> tagJets;
    tagJets.insert(tagJets.end(), bJetsM.begin(), bJetsM.end());
    tagJets.insert(tagJets.end(), bJetsL.begin(), bJetsL.end());

    // Order collections by pt
    sort(extraLeptons.begin(), extraLeptons.end(), P4SortCondition);
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

    isTP    = false;
    crType  = "none";

    // Prepare control regions for FR estimation...  
    // The control region is largely defined by the tag, but we'll first check
    // for any probes since they will be the same for both regions of interest

    vector<TCMuon>      muProbes    = selector->GetSelectedMuons("QCD2l_CR_probe");
    vector<TCElectron>  eleProbes   = selector->GetSelectedElectrons("QCD2l_CR_probe");

    TCPhysObject   muProbe, muPass; 
    TCPhysObject   eleProbe, elePass; 

    UInt_t nMuProbes    = muProbes.size();
    UInt_t nEleProbes   = eleProbes.size();

    if (doQCDDileptonCR && leptons.size() < 2) {
        // For description of QCD dilepton control region, see section 7.4.1 of
        // ttH note (AN-13-159).  First thing is to find the tag lepton and the
        // probe lepton. For this control region, the tag is a muon that is
        // displaced from the PV and is anti-isolated.  The probe is a lepton
        // passing loose identification requirement without any isolation
        // requirement

        UInt_t nTags = selector->GetSelectedMuons("QCD2l_CR_tag").size();
        if (nTags == 1) { 
            tag = (TCPhysObject)selector->GetSelectedMuons("QCD2l_CR_tag")[0];

            // Make sure there is only one probe lepton and that it does not overlap with the tag 
            // Consider probes that are both muons and electrons
            bool isQCD2l = false;
            if (nMuProbes == 1) {
                // Overlap between muon probe and tag ensure by selector
                muProbe = (TCPhysObject)muProbes[0];
                isQCD2l = CheckQCD2lCR(tagJets, muProbe); 
            }

            if (nEleProbes >= 1) {
                nEleProbes = 0;
                for (unsigned i = 0; i < eleProbes.size(); ++i) {
                    TCPhysObject testProbe = (TCPhysObject)eleProbes[i];
                    if (tag.DeltaR(testProbe) > 0.1) {
                        eleProbe = testProbe;
                        ++nEleProbes;
                    }
                }
                if (nEleProbes == 1) 
                    isQCD2l = CheckQCD2lCR(tagJets, eleProbe); 
            }

            bool singleProbe = true;
            if (nEleProbes == 1 && nMuProbes == 1) {
                // Only allow events with muon and electron probes 
                // if they are the same object
                if (eleProbe.DeltaR(muProbe) > 0.1)
                    singleProbe = false;
            }


            if (isQCD2l && (nEleProbes == 1 || nMuProbes == 1) && singleProbe) {
                // Probe object is found and event is consistent with QCD 2l
                // control region requirements. Now fill histograms for
                // parameterizing fake rates by pt and eta.

                isTP    = true;
                crType  = "QCD2l";

            }
        } 
    } 

    if (doZPlusJetCR && leptons.size() >= 2) {
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

            nMuProbes = 0;
            for (unsigned i = 0; i < muProbes.size(); ++i) {
                TCPhysObject testProbe = (TCPhysObject)muProbes[i];
                if (
                        lep1.DeltaR(testProbe) > 0.3 
                        && lep2.DeltaR(testProbe) > 0.3
                        && (tag + testProbe).M() > 12
                   ) {
                    muProbe = testProbe;
                    ++nMuProbes;
                }
            }

            nEleProbes = 0;
            for (unsigned i = 0; i < eleProbes.size(); ++i) {
                TCPhysObject testProbe = (TCPhysObject)eleProbes[i];
                if (
                        lep1.DeltaR(testProbe) > 0.3 
                        && lep2.DeltaR(testProbe) > 0.3
                        && (tag + testProbe).M() > 12
                   ) {
                    eleProbe = testProbe;
                    ++nEleProbes;
                }
            }

            bool singleProbe = true;
            if (nEleProbes == 1 && nMuProbes == 1) {
                // Only allow events with muon and electron probes to proceed
                // if they are the same object
                if (eleProbe.DeltaR(muProbe) > 0.1)
                    singleProbe = false;
            }

            if ((nEleProbes == 1 || nMuProbes == 1) && singleProbe) {
                if (CheckZPlusJetCR()) {

                    // Probe object is found and event is consistent with QCD 2l
                    // control region requirements. Now fill histograms for
                    // parameterizing fake rates by pt and eta.

                    isTP    = true;
                    crType  = "ZPlusJet";
                }
            }
        }
    }


    // Verify that a tag/probe pair is found is found and require that there is
    // only one tight lepton and no more than one b-jet.  Fill tag and
    // denominator/probe distributions

    if (isTP) {

        vObj triggerLeps;
        if (crType == "ZPlusJet") {
            triggerLeps = leptons;
        } /*else if (crType == "QCD2l") {
            triggerLeps.push_back(tag);
            triggerLeps.push_back(probe);
            }*/ 

        if (!isRealData) {
            weighter->SetObjects(triggerLeps, jets, nPUVerticesTrue, passNames[0]);
            Float_t evtWeight = weighter->GetTotalWeight();
            histManager->SetWeight(evtWeight);
        } else
            histManager->SetWeight(1.);

        histManager->SetDirectory(crType + "_inclusive/" + suffix);
        histManager->Fill1DHist(leptons.size(),
                "h1_leptonMult", "lepton multiplicity; N_{leptons}; Entries / bin", 6, -0.5, 5.5);
        histManager->Fill1DHist(jets.size(),
                "h1_jetMult", "jet multiplicity; N_{jets}; Entries / bin", 10, -0.5, 9.5);
        histManager->Fill1DHist(bJetsM.size(),
                "h1_bJetMult", "b-jet multiplicity; N_{b-jet}; Entries / bin", 10, -0.5, 9.5);
    } else
        return kTRUE;

    if (crType == "QCD2l" && leptons.size() <= 1) {


        if (nEleProbes == 1) {
            if (recoMET->Mod() < 50) {
                FillDenominatorHists(crType + "_inclusive", eleProbe);
                FillJetFlavorHists(crType, eleProbe, jets);
            }

            if (recoMET->Mod() < 20)
                FillDenominatorHists(crType + "_low_met", eleProbe);
            else if (recoMET->Mod() > 45 && recoMET->Mod() < 80)
                FillDenominatorHists(crType + "_high_met", eleProbe);
        }

        if (nMuProbes == 1) {
            if (recoMET->Mod() < 50){
                FillDenominatorHists(crType + "_inclusive", muProbe);
                FillJetFlavorHists(crType, muProbe, jets);
            }

            if (recoMET->Mod() < 20)
                FillDenominatorHists(crType + "_low_met", muProbe);
            else if (recoMET->Mod() > 45 && recoMET->Mod() < 80)
                FillDenominatorHists(crType + "_high_met", muProbe);
        }
    } else if (crType == "ZPlusJet" && leptons.size() <= 3) {
        if (nEleProbes == 1) {
            FillJetFlavorHists(crType, eleProbe, jets);
            FillDenominatorHists(crType + "_inclusive", eleProbe);
        }
        if (nMuProbes == 1) {
            FillJetFlavorHists(crType, muProbe, jets);
            FillDenominatorHists(crType + "_inclusive", muProbe);
        }
    } else
        return kTRUE;


    // Match probe leptons to tight leptons
    bool eleMatched = false;
    bool muMatched  = false;
    if (crType == "QCD2l" && leptons.size() == 1) {
        if (nEleProbes == 1) {
            if (leptons[0].Type() == "muon" && olElectrons.size() == 1)
                leptons[0] = (TCPhysObject)olElectrons[0];

            if (eleProbe.DeltaR(leptons[0]) < 0.05 && eleProbe.Type() == leptons[0].Type()) {
                elePass = leptons[0];
                eleMatched = true;
            }
        }

        if (nMuProbes == 1) {
            if (muProbe.DeltaR(leptons[0]) < 0.05 && muProbe.Type() == leptons[0].Type()) {
                muPass = leptons[0];
                muMatched = true;
            }
        }
    } else if (crType == "ZPlusJet" && leptons.size() == 3) {
        if (nEleProbes == 1) {
            if (lep3.Type() == "muon" && olElectrons.size() == 1)
                if (lep3.DeltaR(olElectrons[0]) < 0.1)
                    lep3 = (TCPhysObject)olElectrons[0];

            if (eleProbe.DeltaR(lep3) < 0.05 && eleProbe.Type() == lep3.Type()) {
                elePass = lep3;
                eleMatched = true;
            }
        }

        if (nMuProbes == 1) {
            if (muProbe.DeltaR(lep3) < 0.05 && muProbe.Type() == lep3.Type()) {
                muPass = lep3;
                muMatched = true;
            }
        }
    } 


    if (eleMatched) {
        if (recoMET->Mod() < 50)
            FillNumeratorHists(crType + "_inclusive", elePass);

        if (crType == "QCD2l") {
            if (recoMET->Mod() < 20)
                FillNumeratorHists(crType + "_low_met", elePass);
            else if (recoMET->Mod() > 45 && recoMET->Mod() < 80)
                FillNumeratorHists(crType + "_high_met", elePass);
        }
    } else if (nEleProbes == 1)
        FillClosureHists(crType, eleProbe);

    if (muMatched) {
        if (recoMET->Mod() < 50)
            FillNumeratorHists(crType + "_inclusive", muPass);

        if (crType == "QCD2l") {
            if (recoMET->Mod() < 20)
                FillNumeratorHists(crType + "_low_met", muPass);
            else if (recoMET->Mod() > 45 && recoMET->Mod() < 80)
                FillNumeratorHists(crType + "_high_met", muPass);
        }
    } else  if (nMuProbes == 1)
        FillClosureHists(crType, muProbe);


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
    dileptonMass = -1.;

    float zCandidateMass = 0.;
    for (unsigned i = 0; i < leptons.size(); ++i) {
        for (unsigned j = leptons.size()-1; j > i; --j) {

            //Find lepton pair that is consisten with a Z in +- 10 GeV window
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
    histManager->Fill1DHist(recoMET->Mod(),
            "h1_Met", "MET;MET;Entries / 4 GeV", 25, 0., 100.);
    histManager->Fill2DHist(probe.Pt(), probe.IsoMap("IsoRel"),
            "h2_ProbeLepPtVsIso", "probe lepton p_{T} vs ISO_{rel};p_{T};ISO_{rel}", 50, 0., 150., 25, 0., 1.);

    // fake rate measurement plots
    histManager->Fill1DHist(tag.Pt(),
            "h1_TagLepPt", "tag lepton p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
    histManager->Fill1DHist(tag.Eta(),
            "h1_TagLepEta", "tag lepton #eta;#eta;Entries / bin", 25, -2.5, 2.5);

    float probePt = probe.Pt();
    //if (probe.Pt() < 100.)
    //    probePt = probe.Pt();
    //else (probe.Pt() > 100.)
    //    probePt = 99.;

    string lepType = probe.Type();
    if (lepType == "muon") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_MuProbeLepPt", "probe muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_MuProbeLepEta", "probe muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(probePt,
                "h1_MuDenomPt", "probe muon p_{T};p_{T};Entries / bin", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(probe.Eta(),
                "h1_MuDenomEta", "probe muon #eta;#eta;Entries / bin", 2, etaBinsMu);
        histManager->Fill2DHistUnevenBins(probePt, probe.Eta(),
                "h2_MuDenom", "probe muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_MuDenomMet", "probe muon Met;Met;Entries", nMetBins, metBins);

    } else if (lepType == "electron") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_EleProbeLepPt", "probe electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probe.Eta(),
                "h1_EleProbeLepEta", "probe electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(probePt,
                "h1_EleDenomPt", "probe electron p_{T};p_{T};Entries / bin", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(probe.Eta(),
                "h1_EleDenomEta", "probe electron #eta;#eta;Entries / bin", 3, etaBinsEle);
        histManager->Fill2DHistUnevenBins(probePt, probe.Eta(),
                "h2_EleDenom", "probe electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_EleDenomMet", "probe electron Met;Met;Entries", nMetBins, metBins);
    }
}

void fakeAnalyzer::FillNumeratorHists(string cat, TCPhysObject& passLep)
{
    histManager->SetDirectory(cat + "/" + suffix);

    float passPt = passLep.Pt();
    //if (passLep.Pt() < 100.)
    //    passPt = passLep.Pt();
    //else
    //    passPt = 99.;

    string lepType = passLep.Type();
    if (lepType == "muon") {
        histManager->Fill1DHist(passLep.Pt(),
                "h1_MuPassLepPt", "pass muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(passLep.Eta(),
                "h1_MuPassLepEta", "pass muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(passPt,
                "h1_MuNumerPt", "pass muon p_{T};p_{T};Entries", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(passLep.Eta(),
                "h1_MuNumerEta", "pass muon #eta;#eta;Entries", 2, etaBinsMu);
        histManager->Fill2DHistUnevenBins(passPt, passLep.Eta(),
                "h2_MuNumer", "pass muon p_{T};p_{T};#eta", nPtBins, ptBins, 2, etaBinsMu);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_MuNumerMet", "pass muon Met;Met;Entries", nMetBins, metBins);

    } else if (lepType == "electron") {
        histManager->Fill1DHist(passLep.Pt(),
                "h1_ElePassLepPt", "pass electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(passLep.Eta(),
                "h1_ElePassLepEta", "pass electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(passPt,
                "h1_EleNumerPt", "pass electron p_{T};p_{T};Entries / 3 GeV", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(passLep.Eta(),
                "h1_EleNumerEta", "pass electron #eta;#eta;Entries / 3 GeV", 3, etaBinsEle);
        histManager->Fill2DHistUnevenBins(passPt, passLep.Eta(),
                "h2_EleNumer", "pass electron p_{T};p_{T};#eta", nPtBins, ptBins, 3, etaBinsEle);

        histManager->Fill1DHistUnevenBins(recoMET->Mod(),
                "h1_EleNumerMet", "pass electron Met;Met;Entries", nMetBins, metBins);
    }
}

void fakeAnalyzer::FillClosureHists(string cat, TCPhysObject& probe)
{
    histManager->SetWeight(weighter->GetFakeWeight(probe, cat));
    histManager->SetDirectory(cat + "_inclusive/" + suffix);

    //cout << weighter->GetFakeWeight(tmpObj) << endl;

    if (probe.Type() == "electron") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_ElePtClosure", "electron p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_EleUnevenPtClosure", "electron p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
    } else if (probe.Type() == "muon") {
        histManager->Fill1DHist(probe.Pt(),
                "h1_MuPtClosure", "muon p_{T} (closure);p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHistUnevenBins(probe.Pt(),
                "h1_MuUnevenPtClosure", "muon p_{T} (closure);p_{T};Entries / Bin", nPtBins, ptBins);
    }
    histManager->SetWeight(1.);
}

void fakeAnalyzer::FillJetFlavorHists(string cat, TCPhysObject& lepton, vector<TCJet>& jets) 
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
    // Make sure event is consistent with bbbar production //

    // Match the tag to a loose b-jet
    bool jetMatched = false;
    //bool jetVeto    = false;
    for (unsigned i = 0; i < tagJets.size(); ++i) {
        //cout << tag.DeltaR(tagJets[i]) << "\t" << probe.DeltaR(tagJets[i]) << endl;
        if (tag.DeltaR(tagJets[i]) < 0.3) {
            jetMatched = true;
            continue;
        }
    }

    if (!jetMatched) 
        return false;

    // Check tag/probe pair is back-to-back
    Float_t tpDeltaPhi  = tag.DeltaPhi(probe);
    Float_t tpBalance   = probe.Pt()/(tag.Pt()*(1 + tag.IsoMap("IsoRel"))); 

    histManager->Fill1DHist(fabs(tpDeltaPhi),
            "h1_TagProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
    histManager->Fill1DHist(fabs(tpBalance),
            "h1_TagProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, 0., 4.);

    if (fabs(tpDeltaPhi) < 2.5 || tpBalance > 1) 
        return false;

    // Correct for prompt lepton contamination for low pt probes
    if (probe.Pt() < 10 && recoMET->Mod() > 15) 
        return false;

    return true;
}

bool fakeAnalyzer::CheckZPlusJetCR()
{
    if (recoMET->Mod() > 30)
        return false;
    else
        return true;
}
