#define analyzer_cxx
#include "analyzer.h"

using namespace std;

/////////////////////////////
//Specify parameters here. //
/////////////////////////////

const bool      doPrintout  = false;
const bool      doGenPrint  = false;

// MVA switches
const bool      doMVACut    = true;

// Data-driven BG estimation switches
bool doQFlips = true;
bool doFakes  = true;


/////////////////
//Analysis cuts//
/////////////////


const float   leptonPtCut[]     = {20., 10.};
const float   metCut[]          = {30., 0.};
const float   htCut[]           = {13., 14.};

bool P4SortCondition(const TLorentzVector& p1, const TLorentzVector& p2) {return (p1.Pt() > p2.Pt());} 
bool BTagSortCondition(TCJet j1, TCJet j2) {return (j1.BDiscriminatorMap("CSV") > j2.BDiscriminatorMap("CSV"));} 

void analyzer::Begin(TTree * /*tree*/)
{

    // Job config
    TString option  = GetOption();
    TObjArray *args = (TObjArray*)option.Tokenize(" ");
    string suffix   = (string)((TObjString*)args->At(0))->GetString();
    selection       = (string)((TObjString*)args->At(1))->GetString();
    period          = (string)((TObjString*)args->At(2))->GetString();

    cout << suffix << " " << selection << " " << period << endl;

    weighter        = new WeightUtils(suffix, period, selection, isRealData);
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

        for (unsigned i = 0; i < N_CATEGORIES; ++i) { 

            histoFile[iCut]->mkdir(categoryNames[i].c_str(), categoryNames[i].c_str());
            histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(suffix.c_str(), suffix.c_str());

            if (doQFlips && (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "TEST")) 
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir("QFlips", "QFlips");
            else
                doQFlips = false;

            // Samples for fake bg cleanup
            string fakeMC_2l[] = {"ZJets_M-50", "ZJets_M-10To50", "ttbar", "WWJets2L2Nu", "ZZJets2L2Nu", "ZZJets2L2Q"};
            string fakeMC_3l[] = {"WZJets3LNu", "ZZ4mu", "ZZ4e", "ZZ4tau", "ZZ2e2mu", "ZZ2mu2tau", "ZZ2e2tau", "ttZ", "ttW", "WWW", "WWZ", "WZZ", "ZZZ"};

            if (doFakes && (suffix == "DATA_ELECTRON" || suffix == "DATA_MUEG" || suffix == "DATA_MUON")) {
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir("Fakes_e", "Fakes_e");
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir("Fakes_mu", "Fakes_mu");
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir("Fakes_ll", "Fakes_ll");
            } else if (categoryNames[i].substr(0,2) != "os" || suffix.substr(0, 4) != "FCNC") {
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(("Fakes_e_"+suffix).c_str(), ("Fakes_e_"+suffix).c_str());
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(("Fakes_mu_"+suffix).c_str(), ("Fakes_mu_"+suffix).c_str());
                histoFile[iCut]->GetDirectory(categoryNames[i].c_str())->mkdir(("Fakes_ll_"+suffix).c_str(), ("Fakes_ll_"+suffix).c_str());
            }
        }
    }

    if (doMVACut) {
        //string mva3lCats[4] = {"eee", "eemu", "emumu", "mumumu"};
        string mva3lCats[4] = {"inclusive"};
        for (unsigned i = 0; i < 1; ++i) {
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

            mvaReader->BookMVA("test", ("../data/weights/20131217_205058/TMVAClassification_3l_" + mva3lCats[i] + "_BDTG.weights.xml").c_str());
            mva3lReader.push_back(mvaReader);
        }

        string mvaSSCats[3] = {"inclusive"};
        for (unsigned i = 0; i < 1; ++i) {
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

            mvaReader->BookMVA("test", ("../data/weights/20131217_205058/TMVAClassification_SS_" + mvaSSCats[i] + "_BDTG.weights.xml").c_str());
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

Bool_t analyzer::Process(Long64_t entry)
{
    GetEntry(entry);

    return kTRUE;
}

void analyzer::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    histoFile->Write();
    histoFile->Close();
}

bool analyzer::AnalysisSelection(vObj leptons, vector<TCJet> jets, vector<TCJet> bJetsM, vector<TCJet> bJetsL, TVector3 PV, string histDir)
{
    subdir = histDir;
    //cout << subdir << endl;
    //
    //for (unsigned i = 0; i < leptons.size(); ++i)
    //    cout << leptons[i].Type() << "\t" << leptons[i].IsFake() << ", " ;
    //cout << flavorCat << endl;

    // ZZ control region //
    if (leptons.size() == 4) {
        if ( bJetsM.size() == 0) {
            Make4lPlots(leptons, *recoMET);
            SetYields(13);
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
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 5);
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
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 6);
        SetYields(11);
    }

    // Z+fake control region //
    if (
            zTagged 
            && leptons.size() == 3 
            && MET < 30
       ) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 7);
        SetYields(12);
    }


    //!!!!!!!!!!!!!!//
    // Do MVA Trees //
    //!!!!!!!!!!!!!!//


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

    if ( // Problematic region for same-sign
            leptons.size() == 2 
            && leptons[0].Type() == "muon" 
            && leptons[1].Type() == "muon"
            && (bJetsM.size() + jets.size()) == 0
            && (leptons[0] + leptons[1]).M() < 30
            && recoMET->Mod() < 40
            //&& leptons[0].Pt() < 50
            //&& fabs(leptons[0].Eta() - leptons[1].Eta()) <  1.
       ) 
        return true;

    //    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 1);
    //    SetYields(6);

    //!! Z-veto !!//
    if (
            leptons.size() == 2 
            && leptons[0].Charge() == leptons[1].Charge()
            && leptons[0].Type() == "electron"
            && leptons[1].Type() == "electron"
            && fabs((leptons[0] + leptons[1]).M() - 91.2) < 15. 
       ) 
        return true;
    else if (leptons.size() == 3 && (zTagged || (dileptonMassOS > 50 && fabs(trileptonMass - 91.2) < 7.5))) 
        return true;

    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 1);
    SetYields(6);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    // 0-jet CR //
    if (jets.size() + bJetsM.size() == 0) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 8);
        SetYields(14);
    }

    // 1-jet CR //
    if (jets.size() + bJetsM.size() == 1) {
        MakePlots(leptons, jets, bJetsM, *recoMET, PV, 9);
        SetYields(15);
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

    //!! Require at least two jets !!//
    if (bJetsM.size() + jets.size() <= 1) return true;
    MakePlots(leptons, jets, bJetsM, *recoMET, PV, 2);
    SetYields(7);


    //!! Do mva selection !!//
    if (doMVACut) {
        float mvaValue = -99.;
        float mvaCut = -99.;
        if (leptons.size() == 3) {
            mvaValue = mva3lReader[0]->EvaluateMVA("test");
            mvaCut   = -0.6489;

            //if (flavorCat == 5) {
            //    mvaValue = mva3lReader[0]->EvaluateMVA("test");
            //    mvaCut   = -0.1578;
            //} else if (flavorCat == 6 || flavorCat == 7 || flavorCat == 9) {
            //    mvaValue = mva3lReader[1]->EvaluateMVA("test");
            //    mvaCut   = -0.0289;
            //} else if (flavorCat == 8 || flavorCat == 10 || flavorCat == 11) {
            //    mvaValue = mva3lReader[2]->EvaluateMVA("test");
            //    mvaCut   = -0.0854;
            //} else if (flavorCat == 12) {
            //    mvaValue = mva3lReader[3]->EvaluateMVA("test");
            //    mvaCut   = -0.1532;
            //}

            histManager->SetFileNumber(1);
            histManager->SetDirectory("3l_inclusive/" + subdir);
            histManager->Fill1DHist(mvaValue, "h1_BDT", "BDT value;Entries / bin;BDT", 36, -1., 0.2);

        } else if (leptons.size() == 2 && leptons[0].Charge() == leptons[1].Charge()) {
            mvaValue = mvaSSReader[0]->EvaluateMVA("test");
            mvaCut   = -0.8628;

            //if (flavorCat == 1) {
            //    mvaValue = mvaSSReader[0]->EvaluateMVA("test");
            //    mvaCut   = 0.1390;
            //} else if (flavorCat == 2 || flavorCat == 3) {
            //    mvaValue = mvaSSReader[1]->EvaluateMVA("test");
            //    mvaCut   = -0.2249;
            //} else if (flavorCat == 4) {
            //    mvaValue = mvaSSReader[2]->EvaluateMVA("test");
            //    mvaCut   = -0.499;
            //}

            histManager->SetFileNumber(1);
            histManager->SetDirectory("ss_inclusive/" + subdir);
            histManager->Fill1DHist(mvaValue, "h1_BDT", "BDT value;Entries / bin;BDT", 36, -1., 0.2);
        }

        if (mvaValue > mvaCut) {
            MakePlots(leptons, jets, bJetsM, *recoMET, PV, 4);
            SetYields(9);
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

    return true;
}

void analyzer::MakePlots(vObj leptons, vector<TCJet> jets, vector<TCJet> bJets, TCMET met, TVector3 PV, unsigned cutLevel)
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
                histCategory = 1 + ((evtCategory.to_ulong() & 0xF) >> 2);
                break;
            case 2:
                // eta categories
                histCategory = GetHistCategory(1) + 16;
                break;
            case 3:
                // flavor categories
                histCategory = GetHistCategory(2) + 4;
                break;
            case 4:
                // WH categories
                histCategory = GetHistCategory(3);
                break;
            default:
                histCategory = 0;
        }


        if ((i != 0 && histCategory == 0) || histCategory >= N_CATEGORIES) continue;

        //if (cutLevel == 0) {
        //    cout << categoryNames[histCategory] << ":\t" << leptons[0].Type() << ", " << leptons[0].Eta() << ";\t" << leptons[1].Type() << ", " << leptons[1].Eta() <<  "|\t" << (bJets.size() + jets.size()) << endl;
        //}

        histManager->SetDirectory(categoryNames[histCategory] + "/" + subdir);

        LeptonPlots(leptons, jets, bJets, PV);
        MetPlots(met, leptons);
        JetPlots(jets, bJets);
        DileptonPlots2D(leptons);
        FakePlots(leptons, jets, bJets, PV);
        MiscPlots();

        histManager->SetWeight(1);
        histManager->Fill1DHist(primaryVtx->GetSize(),
                "h1_PvMultUnweighted", "Multiplicity of PVs", 51, -0.5, 50.);
        histManager->SetWeight(evtWeight);
    }
}

void analyzer::LeptonPlots(vObj leptons, vector<TCJet> jets, vector<TCJet> bJets, TVector3 PV)
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
        histManager->Fill1DHist(fabs(lep1.DeltaPhi(lep2)),
                "h1_DileptonOSDeltaPhi", "dilepton #Delta #phi_{OS};#Delta #phi_{OS};Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(lep1.Eta() - lep2.Eta()),
                "h1_DileptonOSDeltaEta", "dilepton #Delta #eta_{OS};#Delta #eta_{OS};Entries / bin", 60, 0., 6.);
        histManager->Fill1DHist(fabs(lep2.DeltaR(lep1)),
                "h1_DileptonOSDeltaR", "dilepton #Delta R_{OS};#Delta R_{OS};Entries / bin", 70, 0., 7.);
        histManager->Fill1DHist(fabs(lep1.Pt() - lep2.Pt())/(lep1.Pt() + lep2.Pt()),
                "h1_DileptonOSDeltaPt", "dilepton #Delta p_{T, OS}/#Sigma p_{T, OS};#Delta p_{T, OS}/#Sigma p_{T, OS};Entries / bin", 50, 0., 1.);
        histManager->Fill1DHist(dileptonP4.Pt()/(lep1.Pt() + lep2.Pt()),
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

        histManager->Fill1DHist(dileptonP4.DeltaR(lep3), 
                "h1_DileptonLepDeltaR", "#Delta R(OS,l3);#Delta R(ll,l);Entries / bin", 70, 0., 7.);
        histManager->Fill1DHist(fabs(dileptonP4.DeltaPhi(lep3)), 
                "h1_DileptonLepDeltaPhi", "#Delta #phi(OS,l3);#Delta #phi(ll,l);Entries / bin", 36, 0., TMath::Pi());
        histManager->Fill1DHist(fabs(dileptonP4.Eta() - lep3.Eta()), 
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


void analyzer::MetPlots(TCMET met, vObj leptons)
{
    histManager->Fill1DHist(met.Mod(), 
            "h1_Met", "MET;MET;Entries / 10 GeV", 35, 0., 350.);
    histManager->Fill1DHist(met.Phi() - TMath::Pi(),
            "h1_MetPhi", "#phi MET;#phi;Entries / 0.087 rad", 36, -TMath::Pi(), TMath::Pi());
    histManager->Fill1DHist(met.SumEt(),
            "h1_MetSumEt", "#Sigma E_{T} of MET;#Sigma E_{T};Entries / 20 GeV", 75, 0., 2600.);
    //histManager->Fill1DHist(met.Significance(),
    //        "h1_MetSig", "MET/#sigma_{MET};MET/#sigma_{MET};Entries", 50, 0., 10.);
    //histManager->Fill1DHist(met.Mod()/met.Significance(),
    //        "h1_MetOverMetSig", "MET/#sigma_{MET};MET/#sigma_{MET};Entries", 50, 0., 10.);

    //cout << met.Significance() << endl;

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


void analyzer::JetPlots(vector<TCJet> jets, vector<TCJet> bJets)
{
    float ptBins[] = {15, 30, 40, 50, 60, 70, 90, 120, 160, 200, 1000};

    histManager->Fill1DHist(jets.size(),
            "h1_JetMult", "Multiplicity of jets;N_{jets};Entries / bin", 11, -0.5, 10.5);
    histManager->Fill1DHist(bJets.size(),
            "h1_BJetMult", "Multiplicity of b-jets;N_{b-jets};Entries / bin", 6, -0.5, 5.5);
    histManager->Fill1DHist(jets.size() + bJets.size(),
            "h1_AllJetMult", "Multiplicity of all jets;N_{jets};Entries / bin", 11, -0.5, 10.5);

    if (jets.size() > 1) 
        histManager->Fill1DHist((jets[0] + jets[1]).M(),
                "h1_DijetMass", "M_{jj};M_{jj};Entries / 20 GeV", 50, 0, 1000);

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


void analyzer::DileptonPlots2D(vObj leptons)
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

void analyzer::MiscPlots()
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

    // Check 3l fakes
    if (flavorCat > 4) {
        unsigned lepIndex = 1;
        //if (lep1.Type() == "electron" && lep1.Type() == lep2.Type())
        //    lepIndex = 1;
        //else if (lep1.Type() == "muon" && lep1.Type() == lep2.Type())
        //    lepIndex = 3;
        //else if (lep1.Type() != lep2.Type())
        //    lepIndex = 5;

        if (lep3.Type() == "electron")
            lepIndex += 0;
        else if (lep3.Type() == "muon")
            lepIndex += 1;

        if (subdir.substr(0,5) == "Fakes") {
            if (!lep1.IsFake() && !lep2.IsFake() && lep3.IsFake()) 
                histManager->Fill1DHist(lepIndex, "h1_FakeCategory", "fake category;fake cat;Entries", 2, 0.5, 2.5);
        } else {
            histManager->Fill1DHist(lepIndex, "h1_FakeCategory", "fake category;fake cat;Entries", 2, 0.5, 2.5);
        }
    }
}

void analyzer::MakeQMisIDPlots(vObj electrons)
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

void analyzer::Make4lPlots(vObj leptons, TCMET met)//, vector<TCJet> jets, vector<TCJet> bJets) 
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
void analyzer::GenPlots(vector<TCGenParticle> gen, vObj leptons)
{

    float    effBins[]  = {10., 20., 30., 45., 65., 90., 140., 200., 400.};

    if (gen.size() == 3) {
        histManager->Fill1DHist((gen[0] + gen[1] + gen[2]).M(),
                "h1_GenTrileptonMass", "gen trilepton mass;M_{lll};Entries / 5 GeV", 60, 0., 300.);
    }

    unsigned higgsLepCount  = 0;
    //unsigned topLepCount    = 0;

    for (unsigned i = 0; i < gen.size(); ++i) {

        /*if (gen.size() == 3) {
            if (gen[i].Mother()->Mother()->GetPDGId() == 25) {
                string index = str(higgsLepCount + 1);
                histManager->Fill1DHist(gen[i].Pt(),
                        "h1_GenHiggsLeptonPt" + index, "Higgs lepton " + index + " p_{T};p_{T," + index + "};Entries / 5 GeV", 40, 0., 200.);
                histManager->Fill1DHist(gen[i].Eta(),
                        "h1_GenHiggsLeptonEta" + index, "Higgs lepton " + index + " #eta;#eta_{" + index + "};Entries / bin", 50, -5., 5.);

                ++higgsLepCount;

            }

            if (fabs(gen[i].Mother()->Mother()->GetPDGId()) == 6) {
                histManager->Fill1DHist(gen[i].Pt(),
                        "h1_TopLeptonPt", "top lepton p_{T};p_{T};Entries / 5 GeV", 40, 0., 200.);
                histManager->Fill1DHist(gen[i].Eta(),
                        "h1_TopLeptonEta", "top lepton #eta;#eta;Entries / bin", 50, -5., 5.);
            }
        }*/



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

int analyzer::GetHistCategory(unsigned shift)
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
        histCategory = 3*lepCat + sum;
    else 
        if (evtCategory.test(16))
            histCategory = 14;
        else if (evtCategory.test(17))
            histCategory = 15;
        else
            histCategory = 0;

    return histCategory;
}

string analyzer::GetFakeCategory(vObj fakeables)
{
    string cat = "";
    if (fakeables.size() == 1) {
        if (fakeables[0].Type() == "electron") 
            cat = "e";
        if (fakeables[0].Type() == "muon") 
            cat = "mu";
    } 

    if (fakeables.size() == 2) {
        cat = "ll";
        //if (fakeables[0].Type() == "electron" && fakeables[0].Type() == "electron") 
        //    cat = "e";
        //if (fakeables[0].Type() == "muon" && fakeables[0].Type() == "electron") 
        //    cat = "emu";
        //if (fakeables[0].Type() == "muon" && fakeables[0].Type() == "muon") 
        //    cat = "mumu";
    }
    return cat;
}

void analyzer::FillYieldHists(string directory, float weight, unsigned cut)
{
    for (unsigned i = 0; i < N_CUTS; ++i) {
        histManager->SetFileNumber(i);
        histManager->SetDirectory(directory);
        histManager->SetWeight(weight);

        histManager->Fill1DHist(cut+1, "h1_YieldByCut", "Weighted number of events passing cuts by cut; cut; Entries", 24, 0.5, 24.5);
        histManager->SetWeight(1);
        histManager->Fill1DHist(cut+1, "h1_YieldByCutRaw", "Raw number of events passing cuts by cut; cut; Entries", 24, 0.5, 24.5);
        histManager->SetWeight(weight);
    }
}

void analyzer::SetYields(unsigned cut)
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
        unsigned flCat = GetHistCategory(2) + 4;
        FillYieldHists(categoryNames[flCat] + "/" + subdir, evtWeight, cut);
        // WH category
        unsigned whCat = GetHistCategory(3);
        if (whCat != 0) {
            FillYieldHists(categoryNames[whCat] + "/" + subdir, evtWeight, cut);
        }
        // geometric category
        unsigned geoCat = GetHistCategory(1) + 16;
        if (geoCat < N_CATEGORIES) {
            FillYieldHists(categoryNames[geoCat] + "/" + subdir, evtWeight, cut);
        }
    }

    if (subdir == suffix) {
        ++eventCount[cut];
        eventCountWeighted[cut] += evtWeight;
    }
}

float analyzer::CalculateTransMass(TCPhysObject lep, TCMET met)
{
    float M_T = sqrt(2*met.Mod()*lep.Pt()*(1 - cos(met.DeltaPhi(lep.P2()))));
    return M_T;
}

TLorentzVector analyzer::CalculateNuP4(TLorentzVector lep, TCMET met)
{
    const float M_W = 80.4;

    float lepEt = (pow(lep.E(),2) - pow(lep.Pz(),2));
    float xTerm = (lep.Px()*met.Px() + lep.Px()*met.Px() + (pow(M_W,2))/2.);

    float pz = (lep.Pz()*xTerm + lep.E()*sqrt(pow(xTerm,2) - pow(met.Mod(),2)*lepEt))/lepEt;

    TLorentzVector nuP4(met.Px(), met.Py(), max(float(0.), pz), sqrt(pow(met.Px(),2) + pow(met.Py(),2) + pow(pz,2)));

    return nuP4;
}
