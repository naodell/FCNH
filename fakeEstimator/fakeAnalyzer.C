#define fakeAnalyzer_cxx

#include "fakeAnalyzer.h"

using namespace std;

/////////////////
//Analysis cuts//
/////////////////

const bool  doQCDDileptonCR = true;
const bool  doGenPrint      = false;

const float jetPtCut[]        = {25., 15.};
const float muPtCut[]         = {5., 3.};
const float elePtCut[]        = {5., 10.};
const float phoPtCut[]        = {10., 10.};

unsigned  nPtBins     = 6;
unsigned  nEtaBins    = 2;
float     ptBins[]    = {5., 10., 20., 35., 50., 70., 100.}; 
float     etaBins[]   = {0., 1.5, 2.5};

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
    histoFile->mkdir("low_met", "low_met");
    histoFile->GetDirectory("low_met", "low_met")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("high_met", "high_met");
    histoFile->GetDirectory("high_met", "high_met")->mkdir(suffix.c_str(), suffix.c_str());

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
        triggerSelector->SetDataBit(isRealData);
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

    // Get jets
    vector<TCJet> jets      = selector->GetSelectedJets("tight");
    vector<TCJet> bJetsL    = selector->GetSelectedJets("bJetsLoose");
    vector<TCJet> bJetsM    = selector->GetSelectedJets("bJetsMedium");
    vector<TCJet> fwdJets   = selector->GetSelectedJets("forward");
    vector<TCJet> muJets    = selector->GetSelectedJets("muJets");
    vector<TCJet> eleJets   = selector->GetSelectedJets("eleJets");

    jets.insert(jets.end(), fwdJets.begin(), fwdJets.end());
    jets.insert(jets.end(), bJetsL.begin(), bJetsL.end());

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

    isTP = false;

    // Prepare control regions for FR estimation...  

    if (doQCDDileptonCR) {
        // For description of QCD dilepton control region, see section 7.4.1 of
        // ttH note (AN-13-159).
        // First thing is to find the tag lepton and the probe lepton. For this
        // control region, the tag is a muon that is displaced from the PV and
        // is anti-isolated.  The probe is a lepton passing loose
        // identification requirement without any isolation requirement

        UInt_t nTags    = selector->GetSelectedMuons("QCD2l_CR_tag").size();
        UInt_t nProbes  = selector->GetSelectedMuons("QCD2l_CR_probe").size();

        if (nTags == 1 && nProbes == 1) {

            tagLep    = (TCPhysObject)selector->GetSelectedMuons("QCD2l_CR_tag")[0];
            probeLep  = (TCPhysObject)selector->GetSelectedMuons("QCD2l_CR_probe")[0];

            // Next we  make sure event is consistent with bbbar production
            //probes.push_back(selector->GetSelectedElectrons("QCD2l_CR_probe")); // <-- Add these

            // Match the tag to a loose b-jet
            bool jetMatched = false;
            bool jetVeto    = false;
            for (unsigned i = 0; i < tagJets.size(); ++i) {
                //cout << tagLep.DeltaR(tagJets[i]) << "\t" << probeLep.DeltaR(tagJets[i]) << endl;
                if (tagLep.DeltaR(tagJets[i]) < 0.3) {
                    jetMatched = true;
                    continue;
                }

                if (
                        tagJets.size() > 1 
                        && (probeLep.DeltaR(tagJets[i]) < 0.7 || tagLep.DeltaR(tagJets[i]) < 0.7 )
                   )
                    jetVeto = true;
            }

            if (!jetMatched) return kTRUE;

            // Check tag/probe pair is back-to-back
            Float_t tpDeltaPhi  = tagLep.DeltaPhi(probeLep);
            Float_t tpBalance   = probeLep.Pt()/(tagLep.Pt()*(1 + tagLep.IsoMap("IsoRel"))); 

            histManager->Fill1DHist(fabs(tpDeltaPhi),
                    "h1_TagProbeDeltaPhi", "#Delta #phi (tag,probe);#Delta #phi (tag,probe);Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(tpBalance),
                    "h1_TagProbePtBalance", "balance (tag,probe);balance (tag,probe);Entries / bin", 40, -2., 2.);

            if (fabs(tpDeltaPhi) < 2.5 || tpBalance > 1) return kTRUE;

            // Correct for prompt lepton contamination for low pt probes
            if (probeLep.Pt() < 10 && recoMET->Mod() > 15) return kTRUE;

            // Fakeable object is found!!! Now fill histograms for parameterizing fake rates by pt and eta
            isTP = true;
            FillDenominatorHists("inclusive");

            if (recoMET->Mod() < 20)
                FillDenominatorHists("low_met");
            else if (recoMET->Mod() > 45 && recoMET->Mod() < 80)
                FillDenominatorHists("high_met");
        } else
            return kTRUE;
    } else
        return kTRUE;

    // Verify that a tag/probe pair is found is found
    if (!isTP) return kTRUE;

    // Require that there is only one tight lepton
    if (leptons.size() != 1)
        return kTRUE;

    // Match probe lepton to tight leptons
    bool matched = false;
    if (probeLep.DeltaR(leptons[0]) < 0.01) {
        passLep = leptons[0];
        matched = true;
    }

    histManager->SetDirectory("inclusive/" + suffix);
    histManager->Fill1DHist(leptons.size(),
            "h1_leptonMult", "lepton multiplicity; N_{leptons}; Entries / bin", 6, -0.5, 5.5);
    histManager->Fill1DHist(jets.size(),
            "h1_jetMult", "jet multiplicity; N_{jets}; Entries / bin", 10, -0.5, 9.5);
    histManager->Fill1DHist(bJetsM.size(),
            "h1_bJetMult", "b-jet multiplicity; N_{b-jet}; Entries / bin", 10, -0.5, 9.5);

    if (matched) {
        FillNumeratorHists("inclusive");
        if (recoMET->Mod() < 20)
            FillNumeratorHists("low_met");
        else if (recoMET->Mod() > 45 && recoMET->Mod() < 80)
            FillNumeratorHists("high_met");
    } else
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


void fakeAnalyzer::DoZTag(vObj leptons)
{
    // Reset OS variables for each event //
    zTagged     = false;
    ossfTagged  = false;
    dileptonMassOS = -1.;

    //cout << leptons.size() << ":\t";

    float zCandidateMass = 0.;
    for (unsigned i = 0; i < leptons.size(); ++i) {
        for (unsigned j = leptons.size()-1; j > i; --j) {

            // Check for opposite-sign pair //
            if (leptons[i].Charge() != leptons[j].Charge()) {
                ossfTagged = true;

                // Is the pair mass consistent with the Z mass within a 20 GeV window?
                if (
                        fabs((leptons[i] + leptons[j]).M() - 91.2) < 10 
                        && leptons[i].Type() == leptons[j].Type()
                   ) {
                    zTagged = true;
                    zCandidateMass = (leptons[i] + leptons[j]).M();
                }
            }
        }
    }
}

void fakeAnalyzer::FillDenominatorHists(string cat)
{
    histManager->SetDirectory(cat + "/" + suffix);

    // Sanity check plots
    histManager->Fill1DHist(recoMET->Mod(),
            "h1_Met", "MET;MET;Entries / 4 GeV", 25, 0., 100.);

    // fake rate measurement plots
    histManager->Fill1DHist(tagLep.Pt(),
            "h1_TagLepPt", "tag lepton p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
    histManager->Fill1DHist(tagLep.Eta(),
            "h1_TagLepEta", "tag lepton #eta;#eta;Entries / bin", 25, -2.5, 2.5);

    string lepType = probeLep.Type();
    if (lepType == "muon") {
        histManager->Fill1DHist(probeLep.Pt(),
                "h1_MuProbeLepPt", "probe muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probeLep.Eta(),
                "h1_MuProbeLepEta", "probe muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(probeLep.Pt(),
                "h1_MuDenomPt", "probe muon p_{T};p_{T};Entries / 3 GeV", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(probeLep.Eta(),
                "h1_MuDenomEta", "probe muon #eta;#eta;Entries / 3 GeV", nEtaBins, etaBins);
        histManager->Fill2DHistUnevenBins(probeLep.Pt(), probeLep.Eta(),
                "h2_MuDenom", "probe muon p_{T};p_{T};#eta", nPtBins, ptBins, nEtaBins, etaBins);

    }   else if (lepType == "electron") {
        histManager->Fill1DHist(probeLep.Pt(),
                "h1_EleProbeLepPt", "probe electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(probeLep.Eta(),
                "h1_EleProbeLepEta", "probe electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(probeLep.Pt(),
                "h1_EleDenomPt", "probe electron p_{T};p_{T};Entries / 3 GeV", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(probeLep.Eta(),
                "h1_EleDenomEta", "probe electron #eta;#eta;Entries / 3 GeV", nEtaBins, etaBins);
        histManager->Fill2DHistUnevenBins(probeLep.Pt(), probeLep.Eta(),
                "h2_EleDenom", "probe electron p_{T};p_{T};#eta", nPtBins, ptBins, nEtaBins, etaBins);
    }
}

void fakeAnalyzer::FillNumeratorHists(string cat)
{
    histManager->SetDirectory(cat + "/" + suffix);

    string lepType = passLep.Type();
    if (lepType == "muon") {
        histManager->Fill1DHist(passLep.Pt(),
                "h1_MuPassLepPt", "pass muon p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(passLep.Eta(),
                "h1_MuPassLepEta", "pass muon #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(passLep.Pt(),
                "h1_MuNumerPt", "pass muon p_{T};p_{T};Entries / 3 GeV", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(passLep.Eta(),
                "h1_MuNumerEta", "pass muon #eta;#eta;Entries / 3 GeV", nEtaBins, etaBins);
        histManager->Fill2DHistUnevenBins(passLep.Pt(), passLep.Eta(),
                "h2_MuNumer", "pass muon p_{T};p_{T};#eta", nPtBins, ptBins, nEtaBins, etaBins);

    }   else if (lepType == "electron") {
        histManager->Fill1DHist(passLep.Pt(),
                "h1_ElePassLepPt", "pass electron p_{T};p_{T};Entries / 3 GeV", 50, 0., 150);
        histManager->Fill1DHist(passLep.Eta(),
                "h1_ElePassLepEta", "pass electron #eta;#eta;Entries / bin", 25, -2.5, 2.5);

        histManager->Fill1DHistUnevenBins(passLep.Pt(),
                "h1_EleNumerPt", "pass electron p_{T};p_{T};Entries / 3 GeV", nPtBins, ptBins);
        histManager->Fill1DHistUnevenBins(passLep.Eta(),
                "h1_EleNumerEta", "pass electron #eta;#eta;Entries / 3 GeV", nEtaBins, etaBins);
        histManager->Fill2DHistUnevenBins(passLep.Pt(), passLep.Eta(),
                "h2_EleNumer", "pass electron p_{T};p_{T};#eta", nPtBins, ptBins, nEtaBins, etaBins);
    }
}
