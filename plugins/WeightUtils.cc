#include "WeightUtils.h"

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData)
{
    _sampleName = sampleName;
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    Initialize();

    // Muon reco efficiencies
    TFile* f_muRecoSF2012_TRIG  = new TFile("../data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root"); // Trigger scale factors
    TFile* f_muRecoSF2012_ID    = new TFile("../data/MuonEfficiencies_Run2012ReReco_53X.root", "OPEN"); // ID scale factors
    TFile* f_muRecoSF2012_ISO   = new TFile("../data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "OPEN"); // ISO scale factors

    _muSF2012_ID[0] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta<0.9");
    _muSF2012_ID[1] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
    _muSF2012_ID[2] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");
    _muSF2012_ID[3] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4");

    _muSF2012_ISO[0] = (TGraphErrors*)f_muRecoSF2012_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9");
    _muSF2012_ISO[1] = (TGraphErrors*)f_muRecoSF2012_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2");
    _muSF2012_ISO[2] = (TGraphErrors*)f_muRecoSF2012_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1");
    _muSF2012_ISO[3] = (TGraphErrors*)f_muRecoSF2012_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4");

    for (unsigned i = 0; i < 4; ++i) {
        Int_t nBins = _muSF2012_ID[i]->GetN();
        Double_t* xPointsID = _muSF2012_ID[i]->GetX();
        Double_t* yPointsID = new Double_t[nBins];

        for (unsigned j = 0; j < nBins; ++j) {
            yPointsID[j] = _muSF2012_ID[i]->GetErrorY(j);
        }
        _muSF2012_ID_err[i] = new TGraph(nBins, xPointsID, yPointsID);

        Double_t* xPointsISO = _muSF2012_ISO[i]->GetX();
        Double_t* yPointsISO = new Double_t[nBins];

        for (unsigned j = 0; j < nBins; ++j) {
            yPointsISO[j] = _muSF2012_ISO[i]->GetErrorY(j);
        }
        _muSF2012_ISO_err[i] = new TGraph(nBins, xPointsISO, yPointsISO);
    }

    h2_MuTriggerSFs[0] = (TH2D*)f_muRecoSF2012_TRIG->Get("DATA_over_MC_Mu17Mu8_OR_Mu17TkMu8_Tight_Mu1_10To20_&_Mu2_20ToInfty_with_STAT_uncrt");
    h2_MuTriggerSFs[1] = (TH2D*)f_muRecoSF2012_TRIG->Get("DATA_over_MC_Mu17Mu8_OR_Mu17TkMu8_Tight_Mu1_20ToInfty_&_Mu2_20ToInfty_with_STAT_uncrt");

    // Electron reco (MVA) efficiencies
    TFile* f_elRecoFile2012 = new TFile("../data/electrons_scale_factors.root", "OPEN");
    h2_EleMVASF = (TH2D*)f_elRecoFile2012->Get("electronsDATAMCratio_FO_ID_ISO");

    // PU weights
    TFile* f_puFile = new TFile("../data/puReweight.root", "OPEN");
    puReweight["2011"]  = (TH1D*)f_puFile->Get("pileupWeights");
    puReweight["2012"]  = (TH1D*)f_puFile->Get("pileupWeights");

    // weights for fake background
    TFile* f_fakeFile = new TFile("../data/fakeRates.root", "OPEN");
    //g_MuonFakesPtB["QCD2l"]         = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l/g_MuonFake_1");
    //g_MuonFakesPtE["QCD2l"]         = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l/g_MuonFake_2");
    //g_ElectronFakesPtB["QCD2l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l/g_ElectronFake_1");
    //g_ElectronFakesPtG["QCD2l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l/g_ElectronFake_2");
    //g_ElectronFakesPtE["QCD2l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l/g_ElectronFake_3");

    //g_MuonFakesPtB["ZPlusJet"]      = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet/g_MuonFake_1");
    //g_MuonFakesPtE["ZPlusJet"]      = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet/g_MuonFake_2");
    //g_ElectronFakesPtB["ZPlusJet"]  = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet/g_ElectronFake_1");
    //g_ElectronFakesPtG["ZPlusJet"]  = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet/g_ElectronFake_2");
    //g_ElectronFakesPtE["ZPlusJet"]  = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet/g_ElectronFake_3");

    //g_MuonFakesPtB["AntiIso3l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("AntiIso3l/g_MuonFake_1");
    //g_MuonFakesPtE["AntiIso3l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("AntiIso3l/g_MuonFake_2");
    //g_ElectronFakesPtB["AntiIso3l"] = (TGraphAsymmErrors*)f_fakeFile->Get("AntiIso3l/g_ElectronFake_1");
    //g_ElectronFakesPtG["AntiIso3l"] = (TGraphAsymmErrors*)f_fakeFile->Get("AntiIso3l/g_ElectronFake_2");
    //g_ElectronFakesPtE["AntiIso3l"] = (TGraphAsymmErrors*)f_fakeFile->Get("AntiIso3l/g_ElectronFake_3");

    h2_MuonFakes["QCD2l"]           = (TH2D*)f_fakeFile->Get("QCD2l/h2_MuonFake");
    h2_MuonFakes["ZPlusJet"]        = (TH2D*)f_fakeFile->Get("ZPlusJet/h2_MuonFake");
    h2_MuonFakes["AntiIso3l"]       = (TH2D*)f_fakeFile->Get("AntiIso3l/h2_MuonFake");
    h2_MuonFakes["Combined"]        = (TH2D*)f_fakeFile->Get("Combined/h2_MuonFake");

    h2_ElectronFakes["QCD2l"]       = (TH2D*)f_fakeFile->Get("QCD2l/h2_ElectronFake");
    h2_ElectronFakes["ZPlusJet"]    = (TH2D*)f_fakeFile->Get("ZPlusJet/h2_ElectronFake");
    h2_ElectronFakes["AntiIso3l"]   = (TH2D*)f_fakeFile->Get("AntiIso3l/h2_ElectronFake");
    h2_ElectronFakes["Combined"]    = (TH2D*)f_fakeFile->Get("Combined/h2_ElectronFake");

    // Weights for charge flip background
    TFile* f_misQFile = new TFile("../data/electronQMisID.root", "OPEN");
    h2_DielectronMisQ = (TH2D*)f_misQFile->Get("inclusive/h2_DielectronMisQ");
    g_QFlipBB_Low   = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQLowJet_BB");
    g_QFlipBE_Low   = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQLowJet_BE");
    g_QFlipEE_Low   = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQLowJet_EE");
    g_QFlipBB_High  = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQHighJet_BB");
    g_QFlipBE_High  = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQHighJet_BE");
    g_QFlipEE_High  = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQHighJet_EE");
    g_QFlipBB = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQ_BB");
    g_QFlipBE = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQ_BE");
    g_QFlipEE = (TGraph*)f_misQFile->Get("inclusive/g_DielectronMisQ_EE");

    // Weight files for AIC background
    TFile* f_aicFile = new TFile("../data/AIC.root", "OPEN");
    g_AIC["mumumu"] = (TGraph*)f_aicFile->Get("inclusive/g_mumumu");
    g_AIC["emumu"]  = (TGraph*)f_aicFile->Get("inclusive/g_emumu");
    g_AIC["eemu"]   = (TGraph*)f_aicFile->Get("inclusive/g_eemu");
    g_AIC["eee"]    = (TGraph*)f_aicFile->Get("inclusive/g_eee");
}

void WeightUtils::Initialize()
{
    _puWeight = 1.;
    _zzWeight = 1.;
    _vbfWeight = 1.;
    _recoWeight = 1.;
    _triggerWeight = 1.;
}

void WeightUtils::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}

void WeightUtils::SetDataPeriod(string dataPeriod)
{
    _dataPeriod = dataPeriod;
}

void WeightUtils::SetSampleName(string sampleName)
{
    _sampleName = sampleName;
}

void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
}

void WeightUtils::SetPassTrigger(string passTrig)
{
    _passTrig = passTrig;
}

void WeightUtils::SetObjects(vector<TCPhysObject>& leptons, vector<TCJet>& jets, float nPU, string passTrig)
{
    _leptons    = leptons;
    _jets       = jets;
    _passTrig   = passTrig;
    _nPU        = nPU;
}

float WeightUtils::GetTotalWeight()
{
    Initialize();
    float weight = 1.;

    if (!_isRealData) {
        weight *= PUWeight();
        weight *= RecoWeight();
    } 

    return weight;
}

float WeightUtils::PUWeight()
{
    if (puReweight.count(_dataPeriod) != 0)
        _puWeight = puReweight[_dataPeriod]->GetBinContent(puReweight[_dataPeriod]->FindBin(_nPU)); 
    else
        _puWeight = 1.;

    return _puWeight;
}

float WeightUtils::RecoWeight()
{
    _triggerWeight = 1.;
    _recoWeight    = 1.;

    // selection efficiencies
    for (unsigned i = 0; i < _leptons.size(); ++i) {
        TCPhysObject lep = _leptons[i];
        if (lep.Type() == "muon") 
            _recoWeight    *= GetMuEff(lep);
        if (lep.Type() == "electron") 
            _recoWeight    *= GetElectronEff(lep);
    }

    // trigger efficiencies
    if (_leptons.size() == 2) {
        if (_leptons[0].Type() == "muon" && _leptons[1].Type() == "muon") {
            _triggerWeight = GetMuTriggerEff(_leptons[0], _leptons[1]);
            //cout << _leptons[0].Eta() << ", (" << _leptons[1].Eta() << ", " << _leptons[1].Pt() << ")\t" << _triggerWeight << endl;
        }
        if (_leptons[0].Type() == "electron" && _leptons[1].Type() == "electron")  {
            _triggerWeight = GetEleTriggerEff(_leptons[0], _leptons[1]);
        }
    } else {
        _triggerWeight = 1.;
    }

    return _triggerWeight*_recoWeight;
}

float WeightUtils::VBFHiggsWeight(float genMass, int higgsMass)
{
    _vbfWeight = float(higgsMass)/genMass;
    return _vbfWeight;
}

float WeightUtils::GetMuTriggerEff(TLorentzVector& lep1, TLorentzVector& lep2) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    float weight = 1.;

    if (lep2.Pt() < 20.) {
        unsigned etaBin1, etaBin2;
        if (fabs(lep1.Eta()) > 0. && fabs(lep1.Eta()) <= 1.2) {
            etaBin1 = 1;
        } else {
            etaBin1 = 2;
        }

        if (fabs(lep2.Eta()) > 0. && fabs(lep2.Eta()) <= 1.2) {
            etaBin2 = 1;
        } else {
            etaBin2 = 2;
        }

        weight = h2_MuTriggerSFs[0]->GetBinContent(etaBin2, etaBin1);

    } else {
        unsigned etaBin1, etaBin2;
        for (int i = 0; i < 4; ++i) {
            if (fabs(lep1.Eta()) > binningEta[i] && fabs(lep1.Eta()) <= binningEta[i+1]) {
                etaBin1 = i+1;
                break;
            }
        }
        for (int i = 0; i < 4; ++i) {
            if (fabs(lep2.Eta()) > binningEta[i] && fabs(lep2.Eta()) <= binningEta[i+1]) {
                etaBin2 = i+1;
                break;
            }
        }

        if (fabs(etaBin2) < fabs(etaBin1)) {
            weight = h2_MuTriggerSFs[1]->GetBinContent(etaBin1, etaBin2);
        } else {
            weight = h2_MuTriggerSFs[1]->GetBinContent(etaBin2, etaBin1);
        }
    }
    return weight;
}

float WeightUtils::GetEleTriggerEff(TLorentzVector& lep1, TLorentzVector& lep2) const
{
    int etaBin[]  = {0,0};
    int ptBin[]   = {0,0};

    float ptBins1[]  = {10., 15., 20., 30., 40., 50., 9999.};
    float ptBins2[]  = {20., 30., 40., 50., 9999.};

    for (int i = 0; i < 4; ++i) {
        if (lep1.Pt() > ptBins1[i] && lep1.Pt() <= ptBins1[i+1]) {
            ptBin[0] = i;
            break;
        }
    }

    if (fabs(lep1.Eta()) < 1.4442) 
        etaBin[0] = 0;
    else if (fabs(lep1.Eta()) < 1.566)
        etaBin[0] = 1;
    else
        etaBin[0] = 2;


    for (int i = 0; i < 4; ++i) {
        if (lep2.Pt() > ptBins2[i] && lep2.Pt() <= ptBins2[i+1]) {
            ptBin[1] = i;
            break;
        }
    }

    if (fabs(lep2.Eta()) < 1.4442) 
        etaBin[1] = 0;
    else if (fabs(lep2.Eta()) < 1.566)
        etaBin[1] = 1;

    // scale factors from Brian
    float _HLTEl17El8_8Leg2012[3][6] = {
        //10<pt<15 15<pt<20 20<pt<30 30<pt<40 40<pt<50 50<pt
        { 0.9085 , 0.9750 , 0.9869 , 0.9908 , 0.9914 , 0.9929 }, // |eta| < 1.4442, trailing
        { 0.9493 , 1.0459 , 1.0033 , 0.9929 , 0.9940 , 0.9944 }, // 1.4442 > |eta| > 1.566
        { 0.9010 , 0.9715 , 0.9966 , 0.9954 , 0.9977 , 0.9979 } // |eta| > 1.566, trailing
    };


    float _HLTEl17El8_17Leg2012[3][4] = {
        //20<pt<30 30<pt<40 40<pt<60 50<pt
        { 0.9863 , 0.9910 , 0.9920 , 0.9933 }, // |eta| < 1.4442, leading
        { 0.9664 , 0.9645 , 0.9752 , 0.9868 }, // 1.4442 > |eta| > 1.566
        { 0.9892 , 0.9965 , 0.9991 , 0.9998 } // |eta| > 1.566, leading
    };

    return _HLTEl17El8_8Leg2012[etaBin[0]][ptBin[0]]*_HLTEl17El8_17Leg2012[etaBin[1]][ptBin[1]];
}

float WeightUtils::GetElectronEff(TLorentzVector& lep) const
{
    float weight = 1.;
    float errStat, errSys;
    if (lep.Pt() < 200.) {
        weight  = h2_EleMVASF->GetBinContent(h2_EleMVASF->FindBin(fabs(lep.Eta()), lep.Pt()));
        errStat = h2_EleMVASF->GetBinError(h2_EleMVASF->FindBin(fabs(lep.Eta()), lep.Pt()));

    } else {
        weight = h2_EleMVASF->GetBinContent(h2_EleMVASF->FindBin(fabs(lep.Eta()), 199.));
        errStat = h2_EleMVASF->GetBinError(h2_EleMVASF->FindBin(fabs(lep.Eta()), 199.));
    }

    errSys  = 0.02*weight;
    float error = sqrt(errStat*errStat + errSys*errSys);
    //weight += error;
    //weight -= error;
    
    //cout << weight << "+/-" << errSys << "+/-" << errStat << " :: " << error << endl;

    return weight;
}

float WeightUtils::GetMuEff(TLorentzVector& lep) const
{
    int etaBin = 0;
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    float weight = 1.;

    for (int i = 0; i < 4; ++i) {
        if (fabs(lep.Eta()) > binningEta[i] && fabs(lep.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    if (lep.Pt() < 300.) {

        float idWeight  = _muSF2012_ID[etaBin]->Eval(lep.Pt());
        float isoWeight = _muSF2012_ISO[etaBin]->Eval(lep.Pt());

        float sysErrID  = 0.005*idWeight;
        float statErrID = _muSF2012_ID_err[etaBin]->Eval(lep.Pt());
        float errID     = sqrt(sysErrID*sysErrID + statErrID*statErrID);

        float sysErrISO  = 0.002*isoWeight;
        float statErrISO = _muSF2012_ISO_err[etaBin]->Eval(lep.Pt());
        float errISO     = sqrt(sysErrISO*sysErrISO + statErrISO*statErrISO);

        //cout << idWeight << ", " << isoWeight << ": " << sysErrID << ", " << statErrID << " :: " << sysErrISO << ", " << statErrISO << endl;

        weight = idWeight*isoWeight;
        //weight = idWeight*(1. + errID)*isoWeight*(1 + errISO);
        //weight = idWeight*(1. - errID)*isoWeight*(1 - errISO);
    } else {
        weight = 1;
    }

    return weight;
}

float WeightUtils::GetFakeWeight(TCPhysObject& fakeable, string controlRegion)
{
    float fakeWeight    = 1.;
    float fakeRate      = 0.;

    unsigned  iPt = 0;
    unsigned  nPtBins = 8;
    float     ptBins[] = {10., 15., 20., 25., 30., 35., 40., 45., 50.}; 

    if (fakeable.Pt() > 40.) {
        iPt = 6;
    } else {
        for (unsigned j = 0; j < nPtBins; ++j) {
            if (fakeable.Pt() > ptBins[j] && fakeable.Pt() < ptBins[j + 1]) {
                iPt = j+1;
                break;
            }
        }
    }

    if (fakeable.Type() == "muon") {

        float fakeablePt = fakeable.Pt();
        //if (fakeable.Pt() < 35) 
        //    fakeablePt = fakeable.Pt();
        //else
        //    fakeablePt = 35;

        if (fabs(fakeable.Eta()) < 1.5) {
            //fakeRate  = g_MuonFakesPtB[controlRegion]->Eval(fakeablePt);
            fakeRate  = h2_MuonFakes[controlRegion]->GetBinContent(iPt, 1);
        } else if (fabs(fakeable.Eta()) >= 1.5) {
            //fakeRate  = g_MuonFakesPtE[controlRegion]->Eval(fakeablePt);
            fakeRate  = h2_MuonFakes[controlRegion]->GetBinContent(iPt, 2);
        }

    } else if (fakeable.Type() == "electron") {

        float fakeablePt = fakeable.Pt();
        if (fakeable.Pt() < 35) 
            fakeablePt = fakeable.Pt();
        else
            fakeablePt = 35;

        if (fabs(fakeable.Eta()) < 0.8) {
            //fakeRate  = g_ElectronFakesPtB[controlRegion]->Eval(fakeablePt);
            fakeRate  = h2_ElectronFakes[controlRegion]->GetBinContent(iPt, 1);
        } else if (fabs(fakeable.Eta()) >= 0.8 && fabs(fakeable.Eta()) < 1.479) {
            //fakeRate  = g_ElectronFakesPtG[controlRegion]->Eval(fakeablePt);
            fakeRate  = h2_ElectronFakes[controlRegion]->GetBinContent(iPt, 2);
        } else if (fabs(fakeable.Eta()) >= 1.479) {
            //fakeRate  = g_ElectronFakesPtE[controlRegion]->Eval(fakeablePt);
            fakeRate  = h2_ElectronFakes[controlRegion]->GetBinContent(iPt, 3);
        }
    }
    //cout << fakeRate << endl;
    if (fakeRate < 0. || fakeRate >= 0.5)
        return 0.;
    else 
        return fakeRate / (1 - fakeRate);
}


float WeightUtils::GetFakeUncertainty(TCPhysObject& fakeable, string controlRegion) 
{
    float fakeError   = 0.;
    unsigned iPt = 0;
    unsigned  nPtBins = 5;
    float     ptBins[] = {10., 15., 20., 25., 30., 35., 40.}; 

    if (fakeable.Pt() > ptBins[nPtBins]) {
        iPt = nPtBins;
    } else {
        for (unsigned j = 0; j < nPtBins; ++j) {
            if (fakeable.Pt() > ptBins[j] && fakeable.Pt() < ptBins[j + 1]) {
                iPt = j+1;
                break;
            }
        }
    }

    float fakeablePt;
    if (fakeable.Type() == "muon") {

        fakeablePt = fakeable.Pt();
        if (fakeable.Pt() < 35)
            fakeablePt = fakeable.Pt();
        else
            fakeablePt = 35;

        if (fabs(fakeable.Eta()) < 1.5) {
            fakeError = g_MuonFakesPtB[controlRegion]->GetErrorY(iPt);
        } else if (fabs(fakeable.Eta()) >= 1.5) {
            fakeError = g_MuonFakesPtE[controlRegion]->GetErrorY(iPt);
        }
    } else if (fakeable.Type() == "electron") {

        fakeablePt = fakeable.Pt();
        if (fakeable.Pt() < 35)
            fakeablePt = fakeable.Pt();
        else
            fakeablePt = 35;


        if (fabs(fakeable.Eta()) < 0.8) {
            fakeError = g_ElectronFakesPtB[controlRegion]->GetErrorY(iPt);
        } else if (fabs(fakeable.Eta()) >= 0.8 && fabs(fakeable.Eta()) < 1.479) {
            fakeError = g_ElectronFakesPtG[controlRegion]->GetErrorY(iPt);
        } else if (fabs(fakeable.Eta()) >= 1.479) {
            fakeError = g_ElectronFakesPtE[controlRegion]->GetErrorY(iPt);
        }
    }

    return fakeError;
}

float WeightUtils::GetQFlipWeight(unsigned nJets, string weightType)
{
    // Set iEta bins for leading and trailing electrons
    float weight = 0.;

    if (weightType == "2D") {
        unsigned iEta1, iPt1, iEta2, iPt2;
        // Set iEta bins for leading and trailing leptons
        if (fabs(_leptons[0].Eta()) < 0.8)
            iEta1 = 0;
        else if (fabs(_leptons[0].Eta()) >= 0.8 && fabs(_leptons[0].Eta()) < 1.479)
            iEta1 = 1;
        else if (fabs(_leptons[0].Eta()) >= 1.479 && fabs(_leptons[0].Eta()) < 2.1)
            iEta1 = 2;

        if (fabs(_leptons[1].Eta()) < 0.8)
            iEta2 = 0;
        else if (fabs(_leptons[1].Eta()) >= 0.8 && fabs(_leptons[1].Eta()) < 1.479)
            iEta2 = 1;
        else if (fabs(_leptons[1].Eta()) >= 1.479 && fabs(_leptons[1].Eta()) < 2.1)
            iEta2 = 2;

        // Set iPt bins for leading and trailing leptons
        if (_leptons[0].Pt() >= 10. && _leptons[0].Pt() < 25.)
            iPt1 = 1;
        else if (_leptons[0].Pt() >= 25. && _leptons[0].Pt() < 40.)
            iPt1 = 2;
        else if (_leptons[0].Pt() >= 40. && _leptons[0].Pt() < 60.)
            iPt1 = 3;
        else if (_leptons[0].Pt() >= 60.)
            iPt1 = 4;

        if (_leptons[1].Pt() >= 10. && _leptons[1].Pt() < 25.)
            iPt2 = 1;
        else if (_leptons[1].Pt() >= 25. && _leptons[1].Pt() < 40.)
            iPt2 = 2;
        else if (_leptons[1].Pt() >= 40. && _leptons[1].Pt() < 60.)
            iPt2 = 3;
        else if (_leptons[1].Pt() >= 60.)
            iPt2 = 4;

        weight = h2_DielectronMisQ->GetBinContent(4*iEta1 + iPt1, 4*iEta2 + iPt2); 

    } else if (weightType == "fit") {
        for (unsigned i = 0; i < _leptons.size(); ++i) {
            if (_leptons[i].Type() != "electron") continue;

            float electronPt    = _leptons[i].Pt();
            float electronEta   = _leptons[i].Eta();
            //if (fabs(_leptons[i].Eta()) < 0.8) {
            //    weight += g_QFlipBB->Eval(electronPt);
            //} else if (fabs(_leptons[i].Eta()) >= 0.8 && fabs(_leptons[i].Eta()) < 1.479) {
            //    weight += g_QFlipBE->Eval(electronPt);
            //} else if (fabs(_leptons[i].Eta()) >= 1.479)
            //    weight += g_QFlipEE->Eval(electronPt);

            if (nJets < 2) {
                if (fabs(_leptons[i].Eta()) < 0.8)
                    weight += g_QFlipBB_Low->Eval(electronPt);
                else if (fabs(electronEta) >= 0.8 && fabs(electronEta) < 1.479)
                    weight += g_QFlipBE_Low->Eval(electronPt);
                else if (fabs(electronEta) >= 1.479)
                    weight += g_QFlipEE_Low->Eval(electronPt);
            } else if (nJets >= 2) {
                if (fabs(electronEta) < 0.8)
                    weight += g_QFlipBB_High->Eval(electronPt);
                else if (fabs(electronEta) >= 0.8 && fabs(electronEta) < 1.479)
                    weight += g_QFlipBE_High->Eval(electronPt);
                else if (fabs(electronEta) >= 1.479)
                    weight += g_QFlipEE_High->Eval(electronPt);
            }
            //if (i == 0 && electronPt < 25.) weight *= 0.75;
        }

        // correction for jet multiplicity
        float jet_corrections[] = {1., 1.2};
        if (nJets == 0)
            weight *= jet_corrections[0];
        else if (nJets >= 1)
            weight *= jet_corrections[1];
    }

    //cout << weight << endl;
    return weight;
}

float WeightUtils::GetAICWeight(const TCPhoton& photon, const string& type)
{
    float aicWeight = 1.;
    float photonPt;

    if (photon.Pt() < 45.)
        photonPt = photon.Pt();
    else 
        photonPt = 45.;

    aicWeight = g_AIC[type]->Eval(photonPt);

    return aicWeight;
}
