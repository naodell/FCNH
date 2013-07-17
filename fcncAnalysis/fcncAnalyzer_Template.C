#define fcncAnalyzer_cxx

#include "fcncAnalyzer.h"

using namespace std;

/////////////////////////////
//Specify parameters here. //
/////////////////////////////


const string  suffix        = "SUFFIX";
const string  selection     = "SELECTION";
const string  period        = "PERIOD";
const bool    doPrintout    = false;
const bool    doGenPrint    = false;
const bool    doWH          = false;
const bool    doMVA         = true;
const bool    doTrees       = false;


/////////////////
//Analysis cuts//
/////////////////


const float   jetPtCut[]        = {30., 15.};
const float   muPtCut[]         = {10., 3.};
const float   elePtCut[]        = {10., 10.};
const float   phoPtCut[]        = {10., 10.};
const float   leptonPtCut[]     = {20., 10.};
const float   metCut[]          = {40., 30.};
const float   bJetVeto          = 1e9;

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
bool BTagSortCondition(TCJet& j1, TCJet& j2) {return (j1.BDiscriminatorMap("CSV") > j2.BDiscriminatorMap("CSV"));} 


void fcncAnalyzer::Begin(TTree* tree) 
{
    // Get trigger names from jobTree
    vector<string>* triggerNames = 0;
    TFile   *inFile         = tree->GetCurrentFile();
    TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");

    jobTree->SetBranchAddress("triggerNames", &triggerNames);
    jobTree->GetEntry();


    // Initialize utilities and selectors here //
    selector        = new Selector(muPtCut, elePtCut, jetPtCut, phoPtCut);
    weighter        = new WeightUtils(suffix, period, selection, isRealData);
    triggerSelector = new TriggerSelector(selection, period, *triggerNames, true);

    // Random numbers! //
    //rnGenerator = new TRandom3();

    // Initialize histograms //
    TString option = GetOption();
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    histManager  = new HistManager();

    for (unsigned iCut = 0; iCut < N_CUTS; ++iCut) {
        string index = str(iCut+1);

        histoFile[iCut] = new TFile(("histos/fcncHistograms_cut" + index + ".root").c_str(), "RECREATE");
        histManager->AddFile(histoFile[iCut]);
        histManager->SetFileNumber(iCut);

        histoFile[iCut]->mkdir("TESTS", "TESTS");
        histoFile[iCut]->GetDirectory("TESTS")->mkdir(suffix.c_str());

        for (unsigned i = 0; i < N_CATEGORIES; ++i) { 

            histoFile[iCut]->mkdir(categoryNames[i].c_str(), categoryNames[i].c_str());
            histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(suffix.c_str(), suffix.c_str());
        }

    }

    // Initialize pass tree for MVA input //
    passTree = new TTree(("passTree_" + suffix).c_str(), "Tree for input into MVA");

    passTree->Branch("met", &MET, "met/F");
    passTree->Branch("HT", &H_T, "HT/F");
    passTree->Branch("MT", &M_T, "MT/F");
    passTree->Branch("TrileptonMass", &TrileptonMass, "TrileptonMass/F");
    passTree->Branch("DileptonMassOS", &DileptonMassOS, "DileptonMassOS/F");
    passTree->Branch("DileptonDROS", &DileptonDROS, "DileptonDROS/F");
    passTree->Branch("AvgBDiscr", &AvgBDiscr, "AvgBDiscr/F");
    passTree->Branch("flavorCat", &flavorCat, "flavorCat/I");
    passTree->Branch("JetMult", &JetMult, "JetMult/I");
    passTree->Branch("BJetMult", &BJetMult, "BJetMult/I");
    passTree->Branch("weights", &weights, "weights/F");

    if (doMVA) {
        mvaReader = new TMVA::Reader("!Color:!Silent");

        mvaReader->AddVariable("met", &MET);
        mvaReader->AddVariable("HT", &H_T);
        mvaReader->AddVariable("MT", &M_T);
        mvaReader->AddVariable("TrileptonMass", &TrileptonMass);
        mvaReader->AddVariable("DileptonMassOS", &DileptonMassOS);
        mvaReader->AddVariable("DileptonDROS", &DileptonDROS);
        mvaReader->AddVariable("AvgBDiscr", &AvgBDiscr);
        mvaReader->AddVariable("flavorCat", &f_flavorCat);
        mvaReader->AddVariable("JetMult", &f_JetMult);
        mvaReader->AddVariable("BJetMult", &f_BJetMult);

        mvaReader->BookMVA("test", "../data/weights/TMVAClassification_BDT.weights.xml");
    }
}

bool fcncAnalyzer::Process(Long64_t entry)
{  
    GetEntry(entry);
    selector->PurgeObjects();

    histManager->SetFileNumber(0);
    histManager->SetDirectory(categoryNames[0] + "/" + suffix);

    if (eventCount[1] == 0) {
        weighter->SetDataBit(isRealData);
        triggerSelector->SetDataBit(isRealData);
    }

    float evtWeight = 1.;
    std::bitset<18> evtCategory;
    SetYields(1, evtCategory, evtWeight);

    if (eventCount[1] % (int)1e4 == 0) cout << eventCount[5] << " events passed of " << eventCount[1] << " checked!" << endl;


    //////////////////
    //Trigger status//
    //////////////////


    histManager->SetFileNumber(0);
    histManager->SetDirectory(categoryNames[0] + "/" + suffix);

    for(int i = 0; i < 64; ++i) {
        if (triggerStatus & ULong64_t(0x1) << i) {
            histManager->Fill1DHist(i+1, "h1_TriggerStatus", "Triggers", 64, 0.5, 64.5);
        }
    } 

    bool triggerPass = false;
    triggerPass = triggerSelector->SelectTrigger(triggerStatus, hltPrescale);

    // Double electron workaround.  Gets rid of hopelessly prescaled events fo July 20-26, 2011
    //if (selection == "electron" && (runNumber > 171200 && runNumber < 171600)) return kTRUE;

    if (!triggerPass) return kTRUE;
    SetYields(2, evtCategory, evtWeight);

    vstring passNames = triggerSelector->GetPassNames();

    if (passNames.size() == 0) passNames.push_back("NULL");


    ////////////////////////////
    //Check the event vertices//
    ////////////////////////////


    if (!isRealData) {
        histManager->SetFileNumber(0);
        histManager->SetDirectory(categoryNames[0] + "/" + suffix);

        histManager->Fill1DHist(nPUVertices,
                "h1_SimVertexMult", "Multiplicity of simulated vertices", 500, 0., 100.);
        histManager->Fill1DHist(nPUVerticesTrue,
                "h1_SimVertexMultTrue", "True simulated PU", 500, 0., 100.);
    }

    selector->PVSelector(primaryVtx);

    if (selector->GetSelectedPVs().size() < 1) return kTRUE;
    SetYields(3, evtCategory, evtWeight);

    TVector3 selectedVtx = *selector->GetSelectedPVs()[0];


    //////////////////
    // Data quality //
    //////////////////


    if (isRealData && (
                NoiseFilters_isCSCTightHalo
                //|| NoiseFilters_isNoiseHcalHBHE 
                //|| NoiseFilters_isScraping 
                )
       ) return kTRUE;

    SetYields(4, evtCategory, evtWeight);

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
    vector<TCJet> allJets;
    vector<TCJet> jets      = selector->GetSelectedJets("tight");
    vector<TCJet> bJets     = selector->GetSelectedJets("bJets");
    vector<TCJet> fwdJets   = selector->GetSelectedJets("forward");

    jets.insert(jets.end(), fwdJets.begin(), fwdJets.end());
    allJets.insert(allJets.end(), jets.begin(), jets.end());
    allJets.insert(allJets.end(), bJets.begin(), bJets.end());

    // Order collections by pt
    sort(extraLeptons.begin(), extraLeptons.end(), P4SortCondition);
    sort(jets.begin(), jets.end(), P4SortCondition);
    sort(bJets.begin(), bJets.end(), P4SortCondition);
    sort(leptons.begin(), leptons.end(), P4SortCondition);


    //!!!!!!!!!!!!!!!!!!!!!!!!//
    //                        //
    //  Event categorization  //
    //  and weighting....     //
    //                        //
    //!!!!!!!!!!!!!!!!!!!!!!!!//


    histManager->SetFileNumber(0);
    histManager->SetDirectory(categoryNames[0] + "/" + suffix);
    histManager->Fill1DHist(leptons.size(), "h1_LeptonMult", "lepton multiplicity;N_{leptons};Events / bin", 6, -0.5, 5.5);


    if (!isRealData) {
        histManager->SetFileNumber(0);

        histManager->SetDirectory("inclusive/" + suffix);
        GenPlots(gLeptons, leptons);
    }

    if (leptons.size() == 2) {
        //!!! Dilepton selection !!!//
        if ( leptons[0].Pt() < leptonPtCut[0] || leptons[1].Pt() < leptonPtCut[1]) return kTRUE;
    } else if (leptons.size() == 3) {
        //!!! Trilepton selection !!!//
        if ( 
                leptons[0].Pt() < leptonPtCut[0] || leptons[1].Pt() < leptonPtCut[1] || leptons[2].Pt() < leptonPtCut[1]
                || fabs(leptons[0].Charge() + leptons[1].Charge() + leptons[2].Charge()) != 1) 
            return kTRUE;
    } else if (leptons.size() == 4) {
        //!!! Tetralepton selection !!!//
        if (leptons[0].Pt() > leptonPtCut[0]) {

            float mass4L = CalculateFourLeptonMass(leptons);

            if (mass4L != 0)
                histManager->Fill1DHist(mass4L,
                        "h1_4lMass", "4l Mass;M_{4l};Entries / 4 GeV", 50, 5., 250.);
        }
        return kTRUE;
    } else
        return kTRUE;

    //!! low mass resonance rejection !!//
    bool lowMassOS  = false;
    bool highMassOS = false;
    bool isCosmics  = false;

    for (unsigned i = 1; i < leptons.size(); ++i) {
        for (unsigned j = 0; j < i; ++j) {
            if (
                    leptons[i].Type() == leptons[j].Type()
                    && leptons[i].Charge() != leptons[j].Charge()
               ) {
                if ((leptons[i] + leptons[j]).M() < 12)
                    lowMassOS = true;

                if ((leptons[i] + leptons[j]).M() > 85)
                    highMassOS = true;

                if (isRealData && leptons[i].Type() == "muon" && leptons[j].Type() == "muon")
                    if (CosmicMuonFilter(leptons[i], leptons[j]))
                        isCosmics = true;
            }
        }
    }

    if (lowMassOS || isCosmics) return kTRUE;

    GetEventCategory(leptons, evtCategory);
    weighter->SetObjects(leptons, jets, nPUVerticesTrue, passNames[0]);
    evtWeight *= weighter->GetTotalWeight();
    histManager->SetWeight(evtWeight);

    // Fill MVA ntuples //
    if (leptons.size() == 3) {
        SetVarsMVA(leptons, bJets, jets, *recoMET, evtCategory, evtWeight);
        if (doTrees) passTree->Fill();
    }


    if (!isRealData) {
        unsigned cat = 1 + ((evtCategory.to_ulong() >> 2) & 0x3);
        histManager->SetDirectory(categoryNames[cat] + "/" + suffix);
        GenPlots(gLeptons, leptons);

        cat = GetHistCategory(evtCategory, 2) - 10;
        histManager->SetDirectory(categoryNames[cat] + "/" + suffix);
        GenPlots(gLeptons, leptons);
    }


    MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 0);
    SetYields(5, evtCategory, evtWeight);

    if (leptons.size() == 3 && doMVA) {
        float mvaValue = mvaReader->EvaluateMVA("test");

        histManager->SetFileNumber(4);
        histManager->SetDirectory("3l_inclusive/" + suffix);
        histManager->Fill1DHist(mvaValue, "h1_BDT", "BDT value;Entries / bin;BDT", 36, -1., 0.2);

        if (mvaValue > -0.3) {
            SetYields(15, evtCategory, evtWeight);
            MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 4);
        }
        if (mvaValue > -0.2) {
            SetYields(14, evtCategory, evtWeight);
        }
        if (mvaValue > -0.1) {
            SetYields(13, evtCategory, evtWeight);
        }
    }



    //!!!!!!!!!!!!!!!!!!!!!!//
    //                      //
    //  Analysis selection  //
    //  Cut n' Count!!!     //
    //                      //
    //!!!!!!!!!!!!!!!!!!!!!!//


    // For WH sync
    if (doWH) {

        //!! MET cut !!//
        if (evtCategory.test(16)) {
            if (recoMET->Mod() < metCut[1]) 
                return kTRUE;
        } else if (evtCategory.test(17)) {
            if (recoMET->Mod() < metCut[0]) 
                return kTRUE;
        }

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 1);
        SetYields(6, evtCategory, evtWeight);

        //!! Z-veto !!//
        if (leptons.size() == 3 && suffix != "fakes") {
            for (unsigned i = 1; i < leptons.size(); ++i) {
                for (unsigned j = 0; j < i; ++j) {
                    if (selector->IsZCandidate(&leptons[i], &leptons[j], 25.)) return kTRUE;
                }
            }
        }

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 2);
        SetYields(7, evtCategory, evtWeight);

        //!! select on top decay !!//
        if (bJets.size() != 0) return kTRUE;
        if (jets.size() > 0)
            if (jets[0].Pt() > 40)
                return kTRUE;

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 3);
        SetYields(8, evtCategory, evtWeight);

        //!! Delta R dilepton cut !!//
        if (leptons.size() == 3) {
            for (unsigned i = 1; i < leptons.size(); ++i) {
                for (unsigned j = 0; j < i; ++j) {
                    if (leptons[i].Charge() != leptons[j].Charge()) {
                        if (
                                (leptons[i] + leptons[j]).M() > 100.
                                || leptons[i].DeltaR(leptons[j]) > 2.
                           ) 
                            return kTRUE;
                    }
                }
            }
        }

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 4);
        SetYields(9, evtCategory, evtWeight);
    } else {

        //!! Z-veto !!//
        if (leptons.size() == 3) {
            float trileptonMass = (leptons[0] + leptons[1] + leptons[2]).M();
            for (unsigned i = 1; i < leptons.size(); ++i) {
                for (unsigned j = 0; j < i; ++j) {
                    if ( 
                            (selector->IsZCandidate(&leptons[i], &leptons[j], 7.5) && trileptonMass > 100)
                            || (leptons[i].Type() == leptons[j].Type() && leptons[i].Charge() != leptons[j].Charge()
                                && (leptons[i] + leptons[j]).M() > 40 && fabs(trileptonMass - 90.) < 7.5)
                       ) return kTRUE;
                }
            }
        } else if (leptons.size() == 2 && leptons[0].Charge() == leptons[1].Charge()) {
            if (selector->IsZCandidate(&leptons[0], &leptons[1], 10.) ) 
                return kTRUE;
        }

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 1);
        SetYields(6, evtCategory, evtWeight);

        //!! MET+HT cut !!//

        // Calculate HT 
        float HT = 0.;
        for (unsigned i = 0; i < bJets.size(); ++i) HT += bJets[i].Pt();
        for (unsigned i = 0; i < jets.size(); ++i) HT += jets[i].Pt();

        if (leptons.size() == 2){
            if (leptons[0].Charge() == leptons[1].Charge()) 
                if (recoMET->Mod() < metCut[1]) 
                    return kTRUE;
        } else {
            //if (recoMET->Mod() < metCut[0] && HT < 75.) 
            if (
                    (pow(1.5*recoMET->Mod(),2) + pow(HT,2)) < pow(150.,2)
                    || HT < 50 || recoMET->Mod() < 20
               ) 
                return kTRUE;
        }

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 2);
        SetYields(7, evtCategory, evtWeight);

        //!! Require at least one b-jet !!//
        if (bJets.size() == 0) return kTRUE;
        //if (jets.size() == 0) return kTRUE;

        MakePlots(leptons, jets, bJets, *recoMET, selectedVtx, evtWeight, evtCategory, 3);
        SetYields(8, evtCategory, evtWeight);
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

    if (doWH) {
        // WH sync selection
        cout<<"| MET > 30 (40):                     |\t" << eventCount[6]  << "\t|\t" << eventCountWeighted[6] << "\t|"<<endl;
        cout<<"| Z veto:                            |\t" << eventCount[7]  << "\t|\t" << eventCountWeighted[7] << "\t|"<<endl;
        cout<<"| Top veto:                          |\t" << eventCount[8]  << "\t|\t" << eventCountWeighted[8] << "\t|"<<endl;
        cout<<"| Dilepton and Delta R cuts:         |\t" << eventCount[9]  << "\t|\t" << eventCountWeighted[9] << "\t|"<<endl;
    } else {

        // FCNH selection
        cout<<"| Z veto:                            |\t" << eventCount[6]  << "\t|\t" << eventCountWeighted[6] << "\t|"<<endl;
        cout<<"| MET > 40 and HT > 75:              |\t" << eventCount[7]  << "\t|\t" << eventCountWeighted[7] << "\t|"<<endl;
        cout<<"| At least one b-jet:                |\t" << eventCount[8]  << "\t|\t" << eventCountWeighted[8] << "\t|"<<endl;
        cout<<"| BDT > -0.3:                        |\t" << eventCount[15]  << "\t|\t" << eventCountWeighted[15] << "\t|"<<endl;

    }

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


void fcncAnalyzer::GetEventCategory(vObj leptons, bitset<18>& evtCategory)
{
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

        if (fabs(leptons[0].Charge()*leptons[0].Charge()*leptons[0].Charge()) == 1) {
            //!! WH categories !!//
            bool isSet = false;
            for (unsigned j = 1; j < leptons.size(); ++j) {
                for (unsigned k = 0; k < j; ++k) {
                    if (leptons[j].Type() == leptons[k].Type()) {
                        if (
                                leptons[j].Type()       != leptons[3-k-j].Type()
                                && leptons[j].Charge()  != leptons[3-k-j].Charge()
                                && leptons[j].Charge()  == leptons[k].Charge()
                           ) {
                            evtCategory.set(16);
                            isSet = true;
                            break;

                        } else if (
                                leptons[j].Charge() != leptons[k].Charge()
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


int fcncAnalyzer::GetHistCategory(bitset<18> evtCategory, unsigned shift)
{
    /*
        Returns an index that maps onto category names.  Categories are done
        in groups of 4 bits so each shift switches the category type.  As of now, 
        the standard categories go as following (in bits),

        0 - 3:    set bit, os/ss, number of leptons
        4 - 7:    eta categories
        8 - 11:   flavor categories
        12 - 15:  charge categores

        Additionally there is OSSF and SSSF 3 lepton categories for syncing with
        the WH analysis.  

        16: OSSF
        17: SSSF
     */

    unsigned lepCat     = (evtCategory.to_ulong() >> 2) & 0x3;
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


void fcncAnalyzer::SetYields(unsigned cut, bitset<18> evtCategory, float evtWeight)
{
    histManager->SetWeight(evtWeight);

    for (unsigned i = 0; i < N_CUTS; ++i) {
        histManager->SetFileNumber(i);

        if (evtCategory.none()) {
            for (unsigned j = 0; j < N_CATEGORIES; ++j) {
                histManager->SetDirectory(categoryNames[j] + "/" + suffix);
                histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
                histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
            }
        } else {
            // inclusive
            histManager->SetDirectory(categoryNames[0] + "/" + suffix);
            histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
            histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);

            // inclusive for lepton category
            unsigned lepCat = (evtCategory.to_ulong() >> 2) & 0x3;
            histManager->SetDirectory(categoryNames[lepCat+1] + "/" + suffix);
            histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
            histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);

            /*
            // eta category
            unsigned etaCat = GetHistCategory(evtCategory, 1);
            histManager->SetDirectory(categoryNames[etaCat] + "/" + suffix);
            histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
            histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
             */

            // flavor category
            unsigned flCat = GetHistCategory(evtCategory, 2) - 10;
            histManager->SetDirectory(categoryNames[flCat] + "/" + suffix);
            histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
            histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);

            // WH category
            unsigned whCat = GetHistCategory(evtCategory, 3);
            if (whCat != 0) {
                histManager->SetDirectory(categoryNames[whCat] + "/" + suffix);
                histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
                histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 16, 0.5, 16.5);
            }
        }
    }

    ++eventCount[cut];
    eventCountWeighted[cut] += evtWeight;
}


void fcncAnalyzer::MakePlots(vObj leptons, vector<TCJet> jets, vector<TCJet> bJets, TCMET met, TVector3 PV, float evtWeight, bitset<18> evtCategory, unsigned cutLevel)
{

    // Fill histograms that correspond to event category
    // and histograms that are inclusive for 2 and 3 lepton
    // selections separately

    histManager->SetFileNumber(cutLevel);

    //cout << "\n" << endl;

    for (int i = 0; i < 4; ++i) { // fixing 3l selection, removing categories
        unsigned histCategory = 0;
        unsigned flCategory, chCategory;

        switch (i) {
            case 1:
                // inclusive lepton categories
                histCategory = 1 + ((evtCategory.to_ulong() >> 2) & 0x3);
                break;
            case 2:
                // flavor categories
                histCategory = GetHistCategory(evtCategory.to_ulong(), 2) - 10;
                break;
            case 3:
                // WH categories
                histCategory = GetHistCategory(evtCategory.to_ulong(), 3);
                break;
            case 4:
                // eta categories
                histCategory = GetHistCategory(evtCategory.to_ulong(), 2);
                break;
            default:
                histCategory = 0;
        }

        if (i != 0 && histCategory == 0) continue;

        //cout << i << ", " << histCategory << " |\t";

        if (evtCategory.test(3)) {
            chCategory = 5 + (((evtCategory.to_ulong() >> 12) & 0xF) >> 1);
            flCategory = 5 + (((evtCategory.to_ulong() >> 8) & 0xF) >> 1);
        } else {
            chCategory = 1 + (((evtCategory.to_ulong() >> 12) & 0xF) >> 2);
            flCategory = 1 + (((evtCategory.to_ulong() >> 8) & 0xF) >> 2);
        }

        /*
           if (cutLevel == 4 && leptons.size() == 3) { // && i == 3) {
           cout << "\n" << categoryNames[histCategory] << ", " << cutLevel << ", " << met.Mod() << endl;;
        //cout << "\n" << categoryNames[histCategory] << ", " << leptons.size() << endl;;
        for (unsigned j = 0; j < leptons.size(); ++j) 
        cout << setw(10) << left << "\t" << leptons[j].Charge() << "\t" << leptons[j].Type() << endl; //"\t" << leptons[j].Pt() << endl;
        }
         */

        histManager->SetDirectory(categoryNames[histCategory] + "/" + suffix);
        histManager->Fill1DHist(chCategory, "h1_LeptonCharge", "lepton charge categories;charge category;Entries / bin", 12, 0.5, 12.5);
        histManager->Fill1DHist(flCategory, "h1_LeptonFlavor", "lepton flavor categories;flavor category;Entries / bin", 12, 0.5, 12.5);
        histManager->Fill2DHist(chCategory, flCategory, "h2_LepChargeVsFlavor", "lepton flavor vs. charge;charge category;flavor category", 12, 0.5, 12.5, 12, 0.5, 12.5);

        LeptonPlots(leptons, met, jets, bJets, PV, histCategory);
        MetPlots(met, leptons, histCategory);
        JetPlots(jets, bJets, histCategory);
        DileptonPlots2D(leptons, histCategory);


        //////////
        // misc //
        //////////


        histManager->SetDirectory(categoryNames[histCategory] + "/" + suffix);

        histManager->Fill1DHist(evtWeight,
                "h1_EventWeight", "event weight", 100, 0., 3.);
        histManager->Fill1DHist(primaryVtx->GetSize(),
                "h1_PvMult", "Multiplicity of PVs", 51, -0.5, 50.);
        //histManager->Fill1DHist(primaryVtx[0].Z(),
        //        "h1_PvZ", "z_{PV};z_{PV};Entries / bin" 50, 0.5, 50.5);

        histManager->SetWeight(1);
        histManager->Fill1DHist(primaryVtx->GetSize(),
                "h1_PvMultUnweighted", "Multiplicity of PVs", 51, -0.5, 50.);
        histManager->SetWeight(evtWeight);
}
}

void fcncAnalyzer::LeptonPlots(vObj leptons, TCMET met, vector<TCJet> jets, vector<TCJet> bJets, TVector3 PV, unsigned evtCategory)
{
    histManager->SetDirectory(categoryNames[evtCategory] + "/" + suffix);

    float HTs = 0.;
    float HT = 0.;
    unsigned centralCount = 0;
    TLorentzVector sumP4(met.Px(), met.Py(), 0, met.Mod());

    if (leptons.size() == 3) {
        TLorentzVector trileptonP4 = leptons[0] + leptons[1] + leptons[2];
        histManager->Fill1DHist(trileptonP4.M(),
                "h1_TrileptonMass", "M_{lll};M_{lll};Entries / 5 GeV", 58, 10., 300.);
        histManager->Fill1DHist(trileptonP4.Pt(),
                "h1_TrileptonPt", "p_{T,3l};p_{T,3l};Entries / 5 GeV", 39, 10., 400.);

        // Check Trilepton mass for different run ranges
        unsigned runBin = 1;
        float runs[] = {190456, 190782, 190949, 193621, 193833, 196531, 198022, 198913, 198934, 203746, 203768, 208686}; 

        for (unsigned i = 0; i < 11; ++i) {
            if (runNumber > runs[i] && runNumber < runs[i+1]) {
                histManager->Fill2DHist(runBin, trileptonP4.M(),
                        "h1_TrileptonMassVsRunNumber", "M_{lll} vs. run number;run number;M_{lll}", 11, 0, 11, 58, 10., 300.);
                break;
            }
            ++runBin;
        }
        
    }

    for (unsigned i = 0; i < leptons.size(); ++i) {
        string index = str(i+1);

        // Balance variables
        HTs     += leptons[i].Pt();
        sumP4   += leptons[i];

        histManager->Fill1DHist(leptons[i].Pt(),
                "h1_Lepton" + index + "Pt", "p_{T} lepton " + index + " ;p_{T}^{l" + index + "};Entries / 5 GeV", 100, 0., 500.);
        histManager->Fill1DHist(leptons[i].Eta(),
                "h1_Lepton" + index + "Eta", "#eta lepton " + index + ";#eta^{l" + index + "};Entries / bin", 50, -2.5, 2.5);
        histManager->Fill1DHist(leptons[i].Phi(),
                "h1_Lepton" + index + "Phi", "#phi lepton " + index + ";#phi^{l" + index + "};Entries / bin", 36, -TMath::Pi(), TMath::Pi());


        histManager->Fill1DHist(leptons[i].Dxy(&PV), 
                "h1_Lepton" + index + " dxy", "d_{xy} leptons " + index + ";d_{xy};Entries / bin", 100., -.03, 0.03);
        histManager->Fill1DHist(leptons[i].Dz(&PV), 
                "h1_Lepton" + index + " dz", "d_{z} leptons " + index + ";d_{z};Entries / bin", 100., -0.3, 0.3);

        if (leptons[i].Type() == "electron") {
            histManager->Fill1DHist(leptons[i].Pt(),
                    "h1_ElectronPt", "p_{T} Electron;p_{T,e};Entries / 5 GeV", 100, 0., 500.);
            histManager->Fill1DHist(leptons[i].Eta(),
                    "h1_ElectronEta", "#eta Electron;#eta_{e};Entries / bin", 50, -2.5, 2.5);
            histManager->Fill1DHist(leptons[i].Phi(),
                    "h1_ElectronPhi", "#phi Electron;#phi_{e};Entries / bin", 36, -TMath::Pi(), TMath::Pi());
        } else if (leptons[i].Type() == "muon") {
            histManager->Fill1DHist(leptons[i].Pt(),
                    "h1_MuonPt", "p_{T} muon;p_{T,#mu};Entries / 5 GeV", 100, 0., 500.);
            histManager->Fill1DHist(leptons[i].Eta(),
                    "h1_MuonEta", "#eta muon;#eta_{#mu};Entries / bin", 50, -2.5, 2.5);
            histManager->Fill1DHist(leptons[i].Phi(),
                    "h1_MuonPhi", "#phi muon;#phi_{#mu};Entries / bin", 36, -TMath::Pi(), TMath::Pi());
        }

        if (fabs(leptons[i].Eta()) < 1.) 
            ++centralCount;

        if (bJets.size() > 0) {
            histManager->Fill1DHist(fabs(leptons[i].DeltaPhi(bJets[0])),
                    "h1_Lepton" + index + "BJetDeltaPhi", "#Delta #phi;#Delta #phi_{l" + index + ",b};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(leptons[i].Eta() - bJets[0].Eta()),
                    "h1_Lepton" + index + "BJetDeltaEta", "#Delta #eta; #eta_{l" + index + ",b};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(leptons[0].DeltaR(bJets[0])),
                    "h1_Lepton" + index + "BJetDeltaR", "#Delta R;#Delta R_{l" + index + ",b};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(leptons[i].Pt() - bJets[0].Pt())/(leptons[i].Pt() + bJets[0].Pt()),
                    "h1_Lepton" + index + "BJetDeltaPt", "#Delta p_{T,l" + index + "b}/#Sigma p_{T,l" + index + "b};#Delta p_{T,l" + index + "b}/#Sigma p_{T,l" + index + "b};Entries / bin", 100, 0., 1.);
        }
        for (unsigned j = 0; j < i; ++j) {
            string jndex = str(j+1);

            histManager->Fill1DHist((leptons[i] + leptons[j]).M(),
                    "h1_DileptonMass" + index + jndex, "dilepton M_{" + index + jndex + "};M_{" + index + jndex + "};Entries / 5 GeV", 80, 0., 400.);
            histManager->Fill1DHist((leptons[i] + leptons[j]).Mt(),
                    "h1_DileptonTransMass" + index + jndex, "dilepton M_{T," + index + jndex + "};M_{T," + index + jndex + "};Entries / 5 GeV", 100, 0., 500.);
            histManager->Fill1DHist((leptons[i] + leptons[j]).Pt(),
                    "h1_DileptonQt" + index + jndex, "dilepton q_{T," + index + jndex + "};q_{T}^{" + index + jndex + "};Entries / 5 GeV", 100, 0., 500.);

            histManager->Fill1DHist(fabs(leptons[i].DeltaPhi(leptons[j])),
                    "h1_DileptonDeltaPhi" + index + jndex, "dilepton #Delta #phi_{" + index + jndex + "};#Delta #phi_{" + index + jndex + "};Entries / bin", 36, 0., TMath::Pi());
            histManager->Fill1DHist(fabs(leptons[i].Eta() - leptons[j].Eta()),
                    "h1_DileptonDeltaEta" + index + jndex, "dilepton #Delta #eta_{" + index + jndex + "};#Delta #eta_{" + index + jndex + "};Entries / bin", 60, 0., 6.);
            histManager->Fill1DHist(fabs(leptons[i].DeltaR(leptons[j])),
                    "h1_DileptonDeltaR" + index + jndex, "dilepton #Delta R_{" + index + jndex + "};#Delta R_{" + index + jndex + "};Entries / bin", 70, 0., 7.);
            histManager->Fill1DHist(fabs(leptons[i].Pt() - leptons[j].Pt())/(leptons[i].Pt() + leptons[j].Pt()),
                    "h1_DileptonDeltaPt" + index + jndex, "dilepton #Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};#Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};Entries / bin", 100, 0., 1.);
            histManager->Fill1DHist(fabs((leptons[i] + leptons[j]).Pt())/(leptons[i].Pt() + leptons[j].Pt()),
                    "h1_DileptonBalance" + index + jndex, "dilepton #Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};#Delta p_{T" + index + jndex + "}/#Sigma p_{T" + index + jndex + "};Entries / bin", 100, 0., 1.);

            // find the OSSF pair and make physics
            if (
                    leptons.size() == 3 
                    && leptons[i].Charge()  != leptons[j].Charge()
                    && leptons[i].Type()    == leptons[j].Type()
               ) {

                TCPhysObject osPair = TCPhysObject(leptons[i] + leptons[j], 0);
                TCPhysObject upLep  = leptons[3 - (i + j)];

                // OS dilepton pair variables
                histManager->Fill1DHist(osPair.M(),
                        "h1_DileptonOSMass", "OS dilepton M;M_{OS};Entries / 4 GeV", 100, 0., 400.);
                histManager->Fill1DHist(osPair.Mt(),
                        "h1_DileptonOSTransMass", "OS dilepton MT;MT_{OS};Entries / 5 GeV", 100, 0., 500.);
                histManager->Fill1DHist(osPair.Pt(),
                        "h1_DileptonOSQt", "dilepton q_{T,OS};q_{T}^{OS};Entries / 5 GeV", 100, 0., 500.);
                histManager->Fill1DHist(fabs(leptons[i].DeltaPhi(leptons[j])),
                        "h1_DileptonOSDeltaPhi", "dilepton #Delta #phi_{OS};#Delta #phi_{OS};Entries / bin", 36, 0., TMath::Pi());
                histManager->Fill1DHist(fabs(leptons[i].Eta() - leptons[j].Eta()),
                        "h1_DileptonOSDeltaEta", "dilepton #Delta #eta_{OS};#Delta #eta_{OS};Entries / bin", 60, 0., 6.);
                histManager->Fill1DHist(fabs(leptons[i].DeltaR(leptons[j])),
                        "h1_DileptonOSDeltaR", "dilepton #Delta R_{OS};#Delta R_{OS};Entries / bin", 70, 0., 7.);
                histManager->Fill1DHist(fabs(leptons[i].Pt() - leptons[j].Pt())/(leptons[i].Pt() + leptons[j].Pt()),
                        "h1_DileptonOSDeltaPt", "dilepton #Delta p_{T, OS}/#Sigma p_{T, OS};#Delta p_{T, OS}/#Sigma p_{T, OS};Entries / bin", 100, 0., 1.);
                histManager->Fill1DHist(osPair.Pt()/(leptons[i].Pt() + leptons[j].Pt()),
                        "h1_DileptonOSBalance", "dilepton #Delta p_{T, OS}/#Sigma p_{T, OS};#Delta p_{T, OS}/#Sigma p_{T, OS};Entries / bin", 100, 0., 1.);

                histManager->Fill1DHist(osPair.DeltaR(leptons[3-i-j]), 
                        "h1_DileptonLepDeltaR", "#Delta R(ll,l);#Delta R(ll,l);Entries / bin", 70, 0., 7.);
                histManager->Fill1DHist(fabs(osPair.DeltaPhi(leptons[3-i-j])), 
                        "h1_DileptonLepDeltaPhi", "#Delta #phi(ll,l);#Delta #phi(ll,l);Entries / bin", 36, 0., TMath::Pi());
                histManager->Fill1DHist(fabs(osPair.Eta() - leptons[3-i-j].Eta()), 
                        "h1_DileptonLepDeltaEta", "#Delta #eta(ll,l);#Delta #eta(ll,l);Entries / bin", 60, 0., 6.);


                histManager->Fill1DHist(CalculateTransMass(upLep, met),
                        "h1_Top2TransMass", "MT of t->W;MT;Entries / 5 GeV", 80, 0., 400.);

                if (fabs(91.2 - (leptons[i] + leptons[j]).M()) < 15) {
                    histManager->Fill1DHist(fabs(TCPhysObject(leptons[i] + leptons[j], 0).P2().DeltaPhi(met + TCPhysObject(leptons[3 - (i + j)], 0).P2())),
                            "h1_DeltaPhiWZ", "#Delta#phi(W,Z);#Delta#phi(W,Z);Entries / bin", 36, 0, TMath::Pi());
                    histManager->Fill1DHist((leptons[i] + leptons[j]).M(), "h1_DileptonMassWZ", "M_{ll};M_{ll};Entries 2 GeV", 25, 66., 116.);
                }

                if (jets.size() > 0 && bJets.size() > 0) {

                    TCPhysObject top1 = TCPhysObject(leptons[i] + leptons[j] + jets[0], 0);
                    TCPhysObject top2 = TCPhysObject(upLep + bJets[0], 0);

                    histManager->Fill1DHist(top1.M(),
                            "h1_Top1Mass", "top 1 M;M^{t1};Entries / 5 GeV", 40, 0., 400.);
                    histManager->Fill1DHist(top1.Mt(),
                            "h1_Top1TransMass", "top 1 M_{t};M_{t}^{t1};Entries / 5 GeV", 50, 0., 500.);
                    histManager->Fill1DHist(top1.Pt(),
                            "h1_Top1Pt", "t_{1} p_{T};p_{T}^{t1};Entries / 5 GeV", 40, 0., 400.);

                    histManager->Fill1DHist(fabs(top1.P2().DeltaPhi(met)),
                            "h1_DeltaPhiTop1Met", "#Delta #phi(MET, t_{1});#Delta #phi(MET, t_{1});Entries / bin", 36, 0, TMath::Pi());
                    histManager->Fill1DHist(fabs(top2.P2().DeltaPhi(met)),
                            "h1_DeltaPhiTop2Met", "#Delta #phi(MET, t_{2});#Delta #phi(MET, t_{2});Entries / bin", 36, 0, TMath::Pi());

                    histManager->Fill1DHist(top2.Pt(),
                            "h1_Top2Pt", "t_{2} p_{T};p_{T}^{t2};Entries / 5 GeV", 80, 0., 400.);
                    histManager->Fill1DHist(top1.P2().DeltaPhi(met + top2.P2()),
                            "h1_DeltaPhiTop1Top2Met", "#Delta #phi(MET + t_{2}, t_{1});#Delta #phi(MET, t_{1});Entries / bin", 36, 0, TMath::Pi());
                }
            }
        }
    }

    if (leptons.size() == 2) {
        if (leptons[0].Charge() == leptons[1].Charge()){
            histManager->Fill1DHist(leptons[0].Charge()*int(jets.size()),
                    "h1_JetMultCharge", "q_{ll}*N_{jets} for ss dileptons;q_{ll}*N_{jets};Entries / bin", 21, -10.5, 10.5);
        }
    }

    for (unsigned i = 0; i < bJets.size(); ++i) {
        HT      += bJets[i].Pt();
        HTs     += bJets[i].Pt();
        sumP4   += bJets[i];

        if (fabs(bJets[i].Eta()) < 1.)
            ++centralCount;
    }
    for (unsigned i = 0; i < jets.size(); ++i) {
        HT      += jets[i].Pt();
        HTs     += jets[i].Pt();
        sumP4   += jets[i];

        if (fabs(jets[i].Eta()) < 1.)
            ++centralCount;
    }

    histManager->Fill1DHist(centralCount,
            "h1_Centrality", "Event Centrality;centrality;Entries / bin", 10, -0.5, 9.5);
    histManager->Fill1DHist(HT,
            "h1_HT", "H_{T};H_{T};Entries / 10 GeV", 100, 0., 2000.);
    histManager->Fill1DHist(HTs + met.Mod(),
            "h1_HTs", "HT;HT;Entries / 10 GeV", 100, 0., 2000.);
    histManager->Fill1DHist(sumP4.Pt()/HTs,
            "h1_EventBalance", "#Sigma p_{T}/HT_{s};#Sigma p_{T}/HT_{s};Entries / bin", 50, 0., 1.); 

    histManager->Fill2DHist(HT, met.Mod(),
            "h2_metVsHt", "MET vs HT;HT;MET", 50, 0., 1000., 35, 0., 350.); 
    histManager->Fill2DHist(sqrt(HT), met.Mod(),
            "h2_metVsSqrtHt", "MET vs #sqrt{HT};#sqrt{HT};MET", 50, 0., 40., 35, 0., 350.); 
}


void fcncAnalyzer::MetPlots(TCMET met, vObj leptons, unsigned evtCategory)
{
    histManager->SetDirectory(categoryNames[evtCategory] + "/" + suffix);

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


void fcncAnalyzer::JetPlots(vector<TCJet> jets, vector<TCJet> bJets, unsigned evtCategory)
{
    histManager->SetDirectory(categoryNames[evtCategory] + "/" + suffix);

    histManager->Fill1DHist(jets.size(),
            "h1_JetMult", "Multiplicity of jets;N_{jets};Entries / bin", 11, -0.5, 10.5);
    histManager->Fill1DHist(bJets.size(),
            "h1_BJetMult", "Multiplicity of b-jets;N_{b-jets};Entries / bin", 6, -0.5, 5.5);

    float avgBDiscrJet = 0;
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

        avgBDiscrJet += jets[i].BDiscriminatorMap("CSV")/jets.size();
    }

    float avgBDiscrBJet = 0;
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

        avgBDiscrBJet += bJets[i].BDiscriminatorMap("CSV")/bJets.size();
    }

    histManager->Fill1DHist(avgBDiscrBJet,
            "h1_BJetAvgBDiscriminator", "Average b-jet b-discriminator;#Sigma CSV_{b-jet};Entries", 70, 0.5, 1.2);
    histManager->Fill1DHist(avgBDiscrJet,
            "h1_JetAvgBDiscriminator", "Average jet b-discriminator;#Sigma CSV_{jet};Entries", 46., -4., 0.6);

    histManager->Fill2DHist(avgBDiscrJet, avgBDiscrBJet,
            "h2_BJetVsJetAvgBDiscriminator", "Average b-discriminator;#Sigma CSV_{jet};#Sigma CSV_{b-jet}", 23, -4., 0.6, 35, 0.5, 1.2);
}


void fcncAnalyzer::DileptonPlots2D(vObj leptons, unsigned evtCategory)
{
    histManager->SetDirectory(categoryNames[evtCategory] + "/" + suffix);

    for (unsigned i = 1; i < leptons.size(); ++i) {
        string index = str(i+1);

        for (unsigned j = 0; j < i; ++j) {
            unsigned k = 3-i-j;

            string jndex = str(j+1);
            string kndex = str(k+1);

            /*
               histManager->Fill2DHist(leptons[i].Pt(), (leptons[i] + leptons[j]).M(),
               "h2_DileptonM" + index + jndex + "VsLepPt" + index, "M_{" + index + "," + jndex + "} vs p^{l" + index + "}_{T};p^{l" + index + "}_{T};M_{" + index + "," + jndex + "}", 30, 0., 300., 60, 0., 300.);
               histManager->Fill2DHist((leptons[i] + leptons[j]).Pt(), (leptons[i] + leptons[j]).M(),
               "h2_DileptonM" + index + jndex + "VsQt" + index + jndex, "M_{" + index + "," + jndex + "} vs Q^{" + index + jndex + "}_{T};Q^{" + index + jndex + "}_{T};M_{" + index + "," + jndex + "}", 98, 0., 500., 60, 0., 300.);
               histManager->Fill2DHist(fabs(leptons[i].DeltaPhi(leptons[j])), (leptons[i] + leptons[j]).M(),
               "h2_DileptonM" + index + jndex + "VsDeltaPhi" + index + jndex, "M_{"+index + "," + jndex + "} vs #Delta #phi;#Delta #phi_{" + index + "," + jndex + "};M_{" + index + "," + jndex + "}", 36, 0., TMath::Pi(), 98, 0., 500.);
               histManager->Fill2DHist(leptons[i].DeltaR(leptons[j]), (leptons[i] + leptons[j]).M(),
               "h2_DileptonM" + index + jndex + "VsDeltaR" + index + jndex, "M_{"+index + "," + jndex + "} vs #Delta R;#Delta R_{" + index + "," + jndex + "};M_{" + index + "," + jndex + "}", 70, 0., 7., 98, 0., 500.);
             */

            // Dalitz Plots //
            histManager->Fill2DHist((leptons[i] + leptons[j]).M(), (leptons[j] + leptons[k]).M(),
                    "h2_DileptonM" + jndex + kndex + "VsM" + index + jndex, "M_{" + jndex + kndex + "} vs M_{" + index + jndex + "};M_{" + jndex + kndex + "};M_{" + index + jndex + "}", 30, 0., 300., 30, 0., 300.);
            histManager->Fill2DHist(pow((leptons[i] + leptons[j]).M(), 2), pow((leptons[j] + leptons[k]).M(), 2),
                    "h2_DalitzM" + jndex + kndex + "VsM" + index + jndex, "M^{2}_{" + jndex + kndex + "} vs M^{2}_{" + index + jndex + "};M^{2}_{" + jndex + kndex + "};M^{2}_{" + index + jndex + "}", 30, 0., 90000., 30, 0., 90000.);

            if (
                    leptons.size() == 3 
                    && leptons[i].Charge()  != leptons[j].Charge()
                    && leptons[i].Type()    == leptons[j].Type()
               ) {

                /*if (evtCategory == 0) {
                  cout << i << "\t" << j << "\t" << eventNumber << endl;
                  cout << leptons[i].Charge() << "\t" << leptons[j].Charge() << endl;
                  cout << leptons[i].Type() << "\t" << leptons[j].Type() << endl;
                  cout << (leptons[i] + leptons[j]).M() << endl;
                  }*/

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

        //cout << fabs(gen[i].GetPDGId()) << ", " << leptons[j].Type() << ", " << leptons[j].DeltaR(gen[i]) << endl;

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

                //cout << gen[i].Charge() << ", " << leptons[j].Charge() << endl;
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

                //cout << gen[i].Charge() << ", " << leptons[j].Charge() << endl;
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


float fcncAnalyzer::CalculateTransMass(TCPhysObject lep, TCMET met)
{
    float MT = sqrt(2*met.Mod()*lep.Pt()*(1 - cos(met.DeltaPhi(lep.P2()))));
    return MT;
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

bool fcncAnalyzer::CosmicMuonFilter(TCPhysObject muon1, TCPhysObject muon2)
{
    float dimuonAngle = muon1.Angle(muon2.Vect());
    if (TMath::Pi() - dimuonAngle  < 0.05)
        return true;
    else
        return false;
}

void fcncAnalyzer::SetVarsMVA(vObj leptons, vector<TCJet> bJets, vector<TCJet> jets, TCMET met, std::bitset<18> evtCategory, float evtWeight)
{
    float HT = 0.;
    float avgBDiscrJet = 0;
    float avgBDiscrBJet = 0;

    for (unsigned i = 0; i < bJets.size(); ++i) { 
        HT += bJets[i].Pt();
        avgBDiscrBJet += bJets[i].BDiscriminatorMap("CSV")/bJets.size();
    }

    for (unsigned i = 0; i < jets.size(); ++i) {
        HT += jets[i].Pt();
        avgBDiscrJet += jets[i].BDiscriminatorMap("CSV")/jets.size();
    }

    unsigned flCategory;
    if (evtCategory.test(3)) {
        flCategory = 5 + (((evtCategory.to_ulong() >> 8) & 0xF) >> 1);
    } else {
        flCategory = 1 + (((evtCategory.to_ulong() >> 8) & 0xF) >> 2);
    }


    flavorCat       = flCategory;
    MET             = met.Mod();
    H_T             = HT;
    TrileptonMass   = (leptons[0] + leptons[1] + leptons[2]).M();
    AvgBDiscr       = avgBDiscrBJet + avgBDiscrJet;
    BJetMult        = bJets.size();
    JetMult         = jets.size();
    weights         = evtWeight;

    // Jet multiplicities as floats for mva reader :( 
    f_flavorCat     = flCategory;
    f_JetMult       = jets.size();
    f_BJetMult      = bJets.size();

    for (unsigned i = 1; i < leptons.size(); ++i) {
        for (unsigned j = 0; j < i; ++j) {

            if (leptons[i].Type() == leptons[j].Type() && leptons[i].Charge() != leptons[j].Charge()) {
                DileptonMassOS  = (leptons[i] + leptons[j]).M();
                DileptonDROS    = leptons[i].DeltaR(leptons[j]);
                M_T             = CalculateTransMass(leptons[3 - (i + j)], met);
            }
        }
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

                //cout << index1 << "\t" << index2 << "\t\t" << i << "\t" << j << endl;
        }
    }

    if (Z1.M() == 0 || Z2.M() == 0) 
        return 0;
    else
        return (Z1 + Z2).M();
}

//float fcncAnalyzer::CalculateChi2Mass(vObj leptons, vObj jets, vObj bJets, TCMET met)
//{
//    float chi2 = 1e10;
//
//    return chi2;
//}

