#define fakeAnalyzer_cxx

#include "fakeAnalyzer.h"

using namespace std;

/////////////////
//Analysis cuts//
/////////////////

const string suffix = "SUFFIX";

const bool  doGen = false;

const float jetPtCut[]        = {30., 15.};
const float muPtCut[]         = {10., 3.};
const float elePtCut[]        = {10., 10.};
const float phoPtCut[]        = {10., 10.};
const float leptonPtCut[]     = {20., 20.};
const float metCut[]          = {40., 30.};
const float bJetVeto          = 1e9;

float ptBins[]   = {10., 15., 20., 25., 30., 35., 50., 80., 120., 250.};
float etaBins[]  = {0., 1., 1.479, 2., 2.5};    

// Do something about these: should just have one sort condition function
bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 

void fakeAnalyzer::Begin(TTree* tree) 
{
    // Get trigger names from jobTree
    vector<string>* triggerNames = 0;
    TFile   *inFile         = tree->GetCurrentFile();
    TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");

    jobTree->SetBranchAddress("triggerNames", &triggerNames);
    jobTree->GetEntry();

    // Initialize utilities and selectors here //
    selector        = new Selector(muPtCut, elePtCut, jetPtCut, phoPtCut);
    triggerSelector = new TriggerSelector("fakes", "2012", *triggerNames, false);

    // Add single lepton triggers for fake rates //
    vstring triggers;
    triggers.push_back("HLT_Mu8_v");
    triggers.push_back("HLT_Mu17_v");
    triggers.push_back("HLT_Ele8_CaloIdT_TrkIdVL_v");
    triggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    triggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
    triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
    triggerSelector->AddTriggers(triggers);

    // Random numbers! //
    //rnGenerator = new TRandom3();

    // Initialize histograms //
    TString option = GetOption();
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    histManager = new HistManager();

    histoFile   = new TFile("histos/fakeHistograms.root", "RECREATE");

    histoFile->mkdir("inclusive", "inclusive");
    histoFile->GetDirectory("inclusive", "inclusive")->mkdir(suffix.c_str(), suffix.c_str());

    histoFile->mkdir("ele_v1", "ele_v1");
    histoFile->GetDirectory("ele_v1", "ele_v1")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("ele_v2", "ele_v2");
    histoFile->GetDirectory("ele_v2", "ele_v2")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("ele_v3", "ele_v3");
    histoFile->GetDirectory("ele_v3", "ele_v3")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("ele_v4", "ele_v4");
    histoFile->GetDirectory("ele_v4", "ele_v4")->mkdir(suffix.c_str(), suffix.c_str());

    histoFile->mkdir("mu_v1", "mu_v1");
    histoFile->GetDirectory("mu_v1", "mu_v1")->mkdir(suffix.c_str(), suffix.c_str());
    histoFile->mkdir("mu_v2", "mu_v2");
    histoFile->GetDirectory("mu_v2", "mu_v2")->mkdir(suffix.c_str(), suffix.c_str());

    histManager->AddFile(histoFile);
    histManager->SetFileNumber(0);

    cout << endl;
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

    if (!isRealData) {
        selector->GenParticleSelector(genParticles, 25, "Higgs");
        selector->GenParticleSelector(genParticles, 11, "electrons");
        selector->GenParticleSelector(genParticles, 24, "Dubya");
        selector->GenParticleSelector(genParticles, 23, "Zeds");
        selector->GenParticleSelector(genParticles, 13, "muons");
        selector->GenParticleSelector(genParticles, 15, "taus");

        higgs       = selector->GetSelectedGenParticles("Higgs");
        dubyas      = selector->GetSelectedGenParticles("Dubya");
        Zeds        = selector->GetSelectedGenParticles("Zeds");
        gElectrons  = selector->GetSelectedGenParticles("electrons");
        gMuons      = selector->GetSelectedGenParticles("muons");
        gTaus       = selector->GetSelectedGenParticles("taus");

        gLeptons.insert(gLeptons.end(), gElectrons.begin(), gElectrons.end());
        gLeptons.insert(gLeptons.end(), gMuons.begin(), gMuons.end());
        gLeptons.insert(gLeptons.end(), gTaus.begin(), gTaus.end());
    }

    if (doGen) {
        for (unsigned i = 0; i < higgs.size(); ++i)
            if (higgs[i].GetStatus() == 3) 
                cout << higgs[i].GetStatus() << ", " << higgs[i].M() << ", " << higgs[i].Mother() << ", higgs" << endl;
        for (unsigned i = 0; i < dubyas.size(); ++i)
            if (dubyas[i].GetStatus() == 3) 
                cout << "\t status = " << dubyas[i].GetStatus() <<  ", mass = " << dubyas[i].M() << ", pt = " << dubyas[i].Pt() << ", mother = " << dubyas[i].Mother() << ", dubyas" << endl;
        for (unsigned i = 0; i < Zeds.size(); ++i)
            if (Zeds[i].GetStatus() == 3) 
                cout << "\t status = " << Zeds[i].GetStatus() << ", mass = " << Zeds[i].M() << ", pt = " << Zeds[i].Pt() << ", mother = " << Zeds[i].Mother() << ", Zeds" << endl;
        for (unsigned i = 0; i < gElectrons.size(); ++i)
            if (gElectrons[i].GetStatus() == 3) 
                cout << "\t\t" << gElectrons[i].GetStatus() << ", " << gElectrons[i].Pt() << ", " << gElectrons[i].Mother() << ", electrons" << endl;
        for (unsigned i = 0; i < gMuons.size(); ++i)
            if (gMuons[i].GetStatus() == 3) 
                cout << "\t\t" << gMuons[i].GetStatus() << ", " << gMuons[i].Pt() << ", " << gMuons[i].Mother() << ", muons" << endl;
        for (unsigned i = 0; i < gTaus.size(); ++i)
            if (gTaus[i].GetStatus() == 3) 
                cout << "\t\t" << gTaus[i].GetStatus() << ", " << gTaus[i].Pt() << ", " << gTaus[i].Mother() << ", taus" << endl;

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
    vector<TCJet> allJets;
    vector<TCJet> jets      = selector->GetSelectedJets("tight");
    vector<TCJet> bJets     = selector->GetSelectedJets("bJets");
    vector<TCJet> fwdJets   = selector->GetSelectedJets("forward");
    jets.insert(jets.end(), fwdJets.begin(), fwdJets.end());
    allJets.insert(allJets.end(), jets.begin(), jets.end());
    allJets.insert(allJets.end(), bJets.begin(), bJets.end());

    // Order collections by pt
    sort(extraLeptons.begin(), extraLeptons.end(), P4SortCondition);
    sort(allJets.begin(), allJets.end(), P4SortCondition);
    sort(jets.begin(), jets.end(), P4SortCondition);
    sort(bJets.begin(), bJets.end(), P4SortCondition);
    sort(leptons.begin(), leptons.end(), P4SortCondition);


    histManager->Fill1DHist(leptons.size(),
            "h1_leptonMult", "lepton multiplicity; N_{leptons}; Entries / bin", 6, -0.5, 5.5);
    histManager->Fill1DHist(jets.size(),
            "h1_jetMult", "jet multiplicity; N_{jets}; Entries / bin", 10, -0.5, 9.5);
    histManager->Fill1DHist(bJets.size(),
            "h1_bJetMult", "b-jet multiplicity; N_{b-jet}; Entries / bin", 10, -0.5, 9.5);



    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    //                            //
    // Start fake rate estimation //
    // here...                    //
    //                            //
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!//


    // Prepare sample for FR estimation.  Require only one reconstructed
    // lepton, at least one jet, and MET < 30 (20) as in HWW analysis

    /*bool zTag = false;
    pair<unsigned, unsigned> lepIndex;

    for (unsigned i = 1; i < leptons.size(); ++i) {
        for (unsigned j = 0; j < i; ++j) {
            if (selector->IsZCandidate(&leptons[i], &leptons[j], 15.)) {
                zTag = true;
                lepIndex.first  = i;
                lepIndex.second = j;
            }
        }
    }*/

    if (
            (recoMuons->GetSize() + recoElectrons->GetSize()) > 1
            || recoMET->Mod() > 20 
            || allJets.size() == 0
            ) return kTRUE;

    // Do numerator/denominator counting for estimate. 
    // electrons first..

    for (unsigned i = 0; i < 4; ++i) { 

        string index = str(i+1);
        vector<TCElectron> eleDenom = selector->GetSelectedElectrons("denom_v" + index);
        unsigned nMatched   = 0;
        unsigned nDenom     = 0;

        // Select lead jet to mold fake pt spectrum

        if (allJets[0].Pt() < 35) break;

        histManager->SetDirectory("ele_v" + index + "/" + suffix);

        for (unsigned j = 0; j < eleDenom.size(); ++j) {

            if (
                    eleDenom[j].DeltaR(allJets[0]) < 1.
                    //|| eleDenom[j].DeltaR(leptons[lepIndex.first]) == 0. 
                    //|| eleDenom[j].DeltaR(leptons[lepIndex.second]) == 0.
                    ) continue;

            //if (j == 0) {
            //    histManager->Fill1DHist((leptons[lepIndex.first] + leptons[lepIndex.second]).M(), "h1_DileptonMass", "dilepton mass;M_{ll};Entries / 2 GeV", 30, 60, 120);
            //}

            histManager->Fill1DHistUnevenBins(eleDenom[j].Pt(), "h1_DenomPt", "electron fakeable p_{T};p_{T};Entries / 10 GeV", 9, ptBins);
            histManager->Fill1DHistUnevenBins(fabs(eleDenom[j].Eta()), "h1_DenomEta", "electron fakeable #eta;#eta;Entries / bin", 4, etaBins);
            histManager->Fill2DHistUnevenBins(fabs(eleDenom[j].Eta()), eleDenom[j].Pt(), "h2_DenomPtVsEta", "electron fakeable ;#eta;p_{T}", 4, etaBins, 9, ptBins);

            ++nDenom;

            for (unsigned k = 0; k < leptons.size(); ++k) {

                if (leptons[k].Type() != "electron" || eleDenom[j].DeltaR(leptons[k]) != 0.) continue;

                histManager->Fill1DHistUnevenBins(eleDenom[j].Pt(), "h1_NumerPt", "electron fakeable p_{T};p_{T};Entries / 10 GeV", 9, ptBins);
                histManager->Fill1DHistUnevenBins(fabs(eleDenom[j].Eta()), "h1_NumerEta", "electron fakeable #eta;#eta;Entries / bin", 4, etaBins);
                histManager->Fill2DHistUnevenBins(fabs(eleDenom[j].Eta()), eleDenom[j].Pt(), "h2_NumerPtVsEta", "electron fakeable;#eta;p_{T}", 4, etaBins, 9, ptBins);

                ++nMatched; 
            }

            histManager->Fill1DHist(nDenom, "h1_DenomMult", "electron denominator;N_{denom};Entries / bin", 4, -0.5, 3.5);
        }
    }


    // Now the muons...

    for (unsigned i = 0; i < 2; ++i) { 

        // Apply jet conditions based on HWW analysis
        // Might not be necessary, but...
        if (allJets[0].Pt() < 15) break;

        string index = str(i+1);
        vector<TCMuon> muDenom = selector->GetSelectedMuons("denom_v" + index);
        unsigned nMatched   = 0;
        unsigned nDenom     = 0;

        histManager->SetDirectory("mu_v" + index + "/" + suffix);

        for (unsigned j = 0; j < muDenom.size(); ++j) {

            if (
                    muDenom[j].DeltaR(allJets[0]) < 1.
                    //|| muDenom[j].DeltaR(leptons[lepIndex.first]) == 0. 
                    //|| muDenom[j].DeltaR(leptons[lepIndex.second]) == 0.
                    ) continue;

            //if (j == 0) {
            //    histManager->Fill1DHist((leptons[lepIndex.first] + leptons[lepIndex.second]).M(), "h1_DileptonMass", "dilepton mass;M_{ll};Entries / 2 GeV", 30, 60, 120);
            //}

            histManager->Fill1DHistUnevenBins(muDenom[j].Pt(), "h1_DenomPt", "muon fakeable p_{T};p_{T};Entries / 10 GeV", 9, ptBins);
            histManager->Fill1DHistUnevenBins(fabs(muDenom[j].Eta()), "h1_DenomEta", "muon fakeable #eta;#eta;Entries / bin", 4, etaBins);
            histManager->Fill2DHistUnevenBins(fabs(muDenom[j].Eta()), muDenom[j].Pt(), "h2_DenomPtVsEta", "muon fakeable ;#eta;p_{T}", 4, etaBins, 9, ptBins);

            ++nDenom;

            for (unsigned k = 0; k < leptons.size(); ++k) {

                if (leptons[k].Type() != "muon" || muDenom[j].DeltaR(leptons[k]) != 0.) continue;

                histManager->Fill1DHistUnevenBins(muDenom[j].Pt(), "h1_NumerPt", "muon fakeable p_{T};p_{T};Entries / 10 GeV", 9, ptBins);
                histManager->Fill1DHistUnevenBins(fabs(muDenom[j].Eta()), "h1_NumerEta", "muon fakeable #eta;#eta;Entries / bin", 4, etaBins);
                histManager->Fill2DHistUnevenBins(fabs(muDenom[j].Eta()), muDenom[j].Pt(), "h2_NumerPtVsEta", "muon fakeable ;#eta;p_{T}", 4, etaBins, 9, ptBins);

                ++nMatched; 
            }

            histManager->Fill1DHist(nDenom, "h1_DenomMult", "muon denominator;N_{denom};Entries / bin", 4, -0.5, 3.5);
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


