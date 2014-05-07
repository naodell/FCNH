#include "WeightUtils.h"

WeightUtils::WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData)
{
    _sampleName = sampleName;
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    Initialize();

    // Muon reco efficiencies
    //TFile* f_muRecoSF2012 = new TFile("../data/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "OPEN"); 
    TFile* f_muRecoSF2012 = new TFile("../data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "OPEN"); 

    _muSF2012[0] = (TGraphErrors*)f_muRecoSF2012->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9");
    _muSF2012[1] = (TGraphErrors*)f_muRecoSF2012->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2");
    _muSF2012[2] = (TGraphErrors*)f_muRecoSF2012->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1");
    _muSF2012[3] = (TGraphErrors*)f_muRecoSF2012->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4");

    // Electron reco (MVA) efficiencies
    TFile* f_elRecoFile2012 = new TFile("../data/CombinedMethod_ScaleFactors_IdIsoSip.root", "OPEN");
    h2_EleMVASF = (TH2D*)f_elRecoFile2012->Get("h_electronScaleFactor_IdIsoSip");

    // PU weights
    TFile* f_puFile = new TFile("../data/puReweight.root", "OPEN");
    puReweight["2011"]  = (TH1D*)f_puFile->Get("pileupWeights");
    puReweight["2012"]  = (TH1D*)f_puFile->Get("pileupWeights");

    // weights for fake background
    TFile* f_fakeFile = new TFile("../data/fakeRates.root", "OPEN");
    g_MuonFakesPtB["QCD2l"]         = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l_inclusive/g_MuonFake_1");
    g_MuonFakesPtE["QCD2l"]         = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l_inclusive/g_MuonFake_2");
    g_ElectronFakesPtB["QCD2l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l_inclusive/g_ElectronFake_1");
    g_ElectronFakesPtG["QCD2l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l_inclusive/g_ElectronFake_2");
    g_ElectronFakesPtE["QCD2l"]     = (TGraphAsymmErrors*)f_fakeFile->Get("QCD2l_inclusive/g_ElectronFake_3");

    g_MuonFakesPtB["ZPlusJet"]      = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet_inclusive/g_MuonFake_1");
    g_MuonFakesPtE["ZPlusJet"]      = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet_inclusive/g_MuonFake_2");
    g_ElectronFakesPtB["ZPlusJet"]  = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet_inclusive/g_ElectronFake_1");
    g_ElectronFakesPtG["ZPlusJet"]  = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet_inclusive/g_ElectronFake_2");
    g_ElectronFakesPtE["ZPlusJet"]  = (TGraphAsymmErrors*)f_fakeFile->Get("ZPlusJet_inclusive/g_ElectronFake_3");

    // Weights for charge flip background
    TFile* f_misQFile = new TFile("../data/electronQMisID.root", "OPEN");
    h2_DielectronMisQ = (TH2D*)f_misQFile->Get("inclusive/h2_DielectronMisQ");
    g_QFlipB = (TGraph*)f_misQFile->Get("inclusive/g_QFlipB");
    g_QFlipE = (TGraph*)f_misQFile->Get("inclusive/g_QFlipE");

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
        if (_leptons[0].Type() == "muon" && _leptons[1].Type() == "muon") 
            _triggerWeight = GetMuTriggerEff(_leptons[0], _leptons[1]);
        if (_leptons[0].Type() == "electron" && _leptons[1].Type() == "electron") 
            _triggerWeight = GetEleTriggerEff(_leptons[0], _leptons[1]);
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
    float _DataEff_HLTMu17Mu8_8Leg2012[9][3] = {
        //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
        {0.967172 , 0.930476 , 0.916004}, // 10<pt<20
        {0.967066 , 0.93713 , 0.920873}, // 20<pt<30
        {0.965728 , 0.935549 , 0.916849}, // 30<pt<35
        {0.965991 , 0.932407 , 0.915271}, // 35<pt<40
        {0.965568 , 0.935851 , 0.918777}, // 40<pt<50
        {0.964775 , 0.937591 , 0.917744}, // 50<pt<60
        {0.96494 , 0.933094 , 0.918446}, // 60<pt<90
        {0.960397 , 0.939106 , 0.909972}, // 90<pt<140
        {0.961868 , 0.92921 , 0.937057} // 140<pt<500
    };
    float _MCEff_HLTMu17Mu8_8Leg2012[9][3] = {
        //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
        {0.969444 , 0.925885 , 0.921075} , // 10<pt<20
        {0.976136 , 0.945697 , 0.927715} , // 20<pt<30
        {0.976453 , 0.948453 , 0.926418} , // 30<pt<35
        {0.975895 , 0.944976 , 0.925758} , // 35<pt<40
        {0.976162 , 0.946528 , 0.928904} , // 40<pt<50
        {0.975507 , 0.950925 , 0.931956} , // 50<pt<60
        {0.976274 , 0.951396 , 0.926831} , // 60<pt<90
        {0.972801 , 0.951352 , 0.932385} , // 90<pt<140
        {0.971679 , 0.973072 , 0.939368} // 140<pt<500
    };
    float _DataEff_HLTMu17Mu8_17Leg2012[9][3] = {
        //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
        {0.609746 , 0.496018 , 0.428991} , // 10<pt<20
        {0.964723 , 0.923791 , 0.892096} , // 20<pt<30
        {0.964065 , 0.924091 , 0.896823} , // 30<pt<35
        {0.964584 , 0.923641 , 0.898797} , // 35<pt<40
        {0.964363 , 0.928434 , 0.90573} , // 40<pt<50
        {0.963617 , 0.930997 , 0.907169} , // 50<pt<60
        {0.963878 , 0.925745 , 0.908756} , // 60<pt<90
        {0.960051 , 0.935225 , 0.901006} , // 90<pt<140
        {0.959271 , 0.92921 , 0.937057} // 140<pt<500
    };
    float _MCEff_HLTMu17Mu8_17Leg2012[9][3] = {
        //|eta|<0.9 , 0.9<|eta|<1.2 , 1.2<|eta|<2.4
        {0.617508 , 0.488784 , 0.428354} , // 10<pt<20
        {0.97418 , 0.935211 , 0.893312} , // 20<pt<30
        {0.975246 , 0.93891 , 0.903676} , // 30<pt<35
        {0.974711 , 0.93787 , 0.907107} , // 35<pt<40
        {0.975291 , 0.939777 , 0.915754} , // 40<pt<50
        {0.974371 , 0.94515 , 0.920956} , // 50<pt<60
        {0.975252 , 0.946933 , 0.917094} , // 60<pt<90
        {0.972801 , 0.945771 , 0.92517} , // 90<pt<140
        {0.971679 , 0.973072 , 0.931013} // 140<pt<500
    };

    float _HLTMu17Mu8_2012[3][10] = {
        //10<pt<20 , 20<pt<25 , 25<pt<30 , 30<pt<35 , 35<pt<40 , 40<pt<50 , 50<pt<60 , 60<pt<90 , 90<pt<140 , 140<pt<500
        {0.991 , 0.989 , 0.989 , 0.988 , 0.989 , 0.988 , 0.987 , 0.989 , 0.990 , 0.982}, // |eta| < 0.9
        {1.00 , 0.985 , 0.986 , 0.983 , 0.984 , 0.984 , 0.986 , 0.983 , 0.982 , 0.964}, // 0.9 < |eta| < 1.2
        {1.01 , 1.00 , 0.994 , 0.991 , 0.987 , 0.986 , 0.985 , 0.984 , 0.972 , 1.01} // |eta| > 1.2
    };

    float muTrigSF1 = 1.0;
    float muTrigSF2 = 1.0;

    float muTrigDataA8  = 1.0;
    float muTrigDataA17 = 1.0;
    float muTrigDataB8  = 1.0;
    float muTrigDataB17 = 1.0;

    float muTrigMCA8    = 1.0;
    float muTrigMCA17   = 1.0;
    float muTrigMCB8    = 1.0;
    float muTrigMCB17   = 1.0;

    // 2012 use the 2D arrays

    int ptBinV1 = 0;
    int ptBinV2 = 0;
    int etaBin = 0;
    float binningPtV1[] = {10., 20., 25., 30., 35., 40., 50., 60., 90., 140., 500.};
    float binningPtV2[] = {10., 20., 30., 35., 40., 50., 60., 90., 140., 500.};

    if (fabs(lep1.Eta()) < 0.9) {
        etaBin = 0;
    }else if (fabs(lep1.Eta()) < 1.2){
        etaBin = 1;
    }else{
        etaBin = 2;
    }
    for (int i = 0; i < 10; ++i) {
        if (lep1.Pt() >= binningPtV1[i] && lep1.Pt() < binningPtV1[i+1]) {
            ptBinV1 = i;
            break;
        }
    }
    for (int i = 0; i < 9; ++i) {
        if (lep1.Pt() >= binningPtV2[i] && lep1.Pt() < binningPtV2[i+1]) {
            ptBinV2 = i;
            break;
        }
    }
    muTrigSF1 = _HLTMu17Mu8_2012[etaBin][ptBinV1];
    muTrigDataA17 = _DataEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigDataA8 = _DataEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];
    muTrigMCA17 = _MCEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigMCA8 = _MCEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];

    if (fabs(lep2.Eta()) < 0.9) {
        etaBin = 0;
    }else if (fabs(lep2.Eta()) < 1.2){
    }else{
        etaBin = 2;
    }
    for (int i = 0; i < 10; ++i) {
        if (lep2.Pt() >= binningPtV1[i] && lep2.Pt() < binningPtV1[i+1]) {
            ptBinV1 = i;
            break;
        }
    }
    for (int i = 0; i < 9; ++i) {
        if (lep2.Pt() >= binningPtV2[i] && lep2.Pt() < binningPtV2[i+1]) {
            ptBinV2 = i;
            break;
        }
    }

    muTrigSF2 = _HLTMu17Mu8_2012[etaBin][ptBinV1];
    muTrigDataB17 = _DataEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigDataB8 = _DataEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];
    muTrigMCB17 = _MCEff_HLTMu17Mu8_17Leg2012[etaBin][ptBinV2];
    muTrigMCB8  = _MCEff_HLTMu17Mu8_8Leg2012[etaBin][ptBinV2];

    //return muTrigSF1*muTrigSF2;
    return (muTrigDataA8*muTrigDataB17 + muTrigDataA17*muTrigDataB8 - muTrigDataA17*muTrigDataB17)/
        (muTrigMCA8*muTrigMCB17 + muTrigMCA17*muTrigMCB8 - muTrigMCA17*muTrigMCB17);
}

float WeightUtils::GetEleTriggerEff(TLorentzVector& lep1, TLorentzVector& lep2) const
{
    int etaBin[]  = {0,0};
    int ptBin[]   = {0,0};
    float weight = 1.;
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
    else
        etaBin[1] = 2;

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

    if (lep.Pt() < 200) 
        weight = h2_EleMVASF->GetBinContent(h2_EleMVASF->FindBin(lep.Pt(), lep.Eta()));
    else
        weight = h2_EleMVASF->GetBinContent(h2_EleMVASF->FindBin(199, lep.Eta()));

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

    if (lep.Pt() < 500.)
        weight = _muSF2012[etaBin]->Eval(lep.Pt());
    else
        weight = 1;

    return weight;
}

float WeightUtils::GetFakeWeight(TCPhysObject& fakeable, string controlRegion)
{
    float fakeWeight    = 1.;
    float fakeRate      = 0.;

    if (fakeable.Type() == "muon") {
        
        float fakeablePt = fakeable.Pt();
        if (fakeable.Pt() < 35) 
            fakeablePt = fakeable.Pt();
        else
            fakeablePt = 35;

        if (fabs(fakeable.Eta()) < 1.5) {
            fakeRate  = g_MuonFakesPtB[controlRegion]->Eval(fakeablePt);
        } else if (fabs(fakeable.Eta()) >= 1.5) {
            fakeRate  = g_MuonFakesPtE[controlRegion]->Eval(fakeablePt);
        }
    } else if (fakeable.Type() == "electron") {

        float fakeablePt = fakeable.Pt();
        if (fakeable.Pt() < 50) 
            fakeablePt = fakeable.Pt();
        else
            fakeablePt = 50;

        if (fabs(fakeable.Eta()) < 0.8) {
            fakeRate  = g_ElectronFakesPtB[controlRegion]->Eval(fakeablePt);
        } else if (fabs(fakeable.Eta()) >= 0.8 && fabs(fakeable.Eta()) < 1.479) {
            fakeRate  = g_ElectronFakesPtG[controlRegion]->Eval(fakeablePt);
        } else if (fabs(fakeable.Eta()) >= 1.479) {
            fakeRate  = g_ElectronFakesPtE[controlRegion]->Eval(fakeablePt);
        }
    }
    fakeWeight = fakeRate / (1 - fakeRate);

    return fakeWeight;
}

float WeightUtils::GetQFlipWeight()
{
    float weight = 0.;

    // Set iEta bins for leading and trailing electrons
    for (unsigned i = 0; i < _leptons.size(); ++i) {
        if (_leptons[i].Type() != "electron") continue;

        float electronPt;
        if (_leptons[i].Pt() < 125.)
            electronPt = _leptons[i].Pt();
        else
            electronPt = 125.;

        if (fabs(_leptons[0].Eta()) < 1.5)
            weight += g_QFlipB->Eval(electronPt);
        else 
            weight += g_QFlipE->Eval(electronPt);
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

    //cout << type << ": " << photonPt << ", " << aicWeight << endl;

    return aicWeight;
}

float WeightUtils::GetFakeUncertainty() const
{
    return _fakeWeightErr;
}
