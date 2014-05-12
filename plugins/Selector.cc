#include "Selector.h" 

Selector::Selector() 
{
}

Selector::~Selector()
{
    delete electronMVA;
}

Selector::Selector(const float* muPtCuts, const float* elePtCuts, const float* jetPtCuts, const float* phoPtCuts) 
{

    _muPtCuts   = muPtCuts;
    _elePtCuts  = elePtCuts;
    _jetPtCuts  = jetPtCuts;
    _phoPtCuts  = phoPtCuts;
    _vtxIndex   = 0;
    _isRealData = false;

    // b-tag mc efficiencies
    TFile* f_bEff = new TFile("../data/bEff_ttbar_2012.root");
    _misTagEff  = (TGraphAsymmErrors*)f_bEff->Get("g_MistagEff");
    _bTagEff    = (TGraphAsymmErrors*)f_bEff->Get("g_bTagEff");
    _cTagEff    = (TGraphAsymmErrors*)f_bEff->Get("g_cTagEff");

    muCorrector = new rochcor2012(229);

    rnGen = new TRandom3(1337);

}

void Selector::PurgeObjects()
{
    _vtxIndex = 0;

    if (_selMuons.begin() != _selMuons.end()) 
        _selMuons.clear();
    if (_selElectrons.begin() != _selElectrons.end()) 
        _selElectrons.clear();
    if (_selPhotons.begin() != _selPhotons.end()) 
        _selPhotons.clear();
    if (_selJets.begin() != _selJets.end()) 
        _selJets.clear();
    if (_selVertices.begin() != _selVertices.end())
        _selVertices.clear();
    if (_selGenParticles.begin() != _selGenParticles.end())
        _selGenParticles.clear();
    if (_selGenJets.begin() != _selGenJets.end())
        _selGenJets.clear();
}

void Selector::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}

void Selector::SetRho(float rho)
{
    _rho = rho;
}


vector<TVector3*> Selector::GetSelectedPVs()
{
    return _selVertices;
}

vector<TCMuon> Selector::GetSelectedMuons(string key) 
{
    return _selMuons[key];
}

vector<TCElectron> Selector::GetSelectedElectrons(string key)  
{
    return _selElectrons[key];
}

vector<TCPhoton> Selector::GetSelectedPhotons(string key)  
{
    return _selPhotons[key];
}

vector<TCJet> Selector::GetSelectedJets(string key)
{
    return _selJets[key];
}

vector<TCGenParticle> Selector::GetSelectedGenParticles(string key)
{
    return _selGenParticles[key];
}

vector<TCGenJet> Selector::GetSelectedGenJets()
{
    return _selGenJets;
}

bool Selector::IsZCandidate(TCPhysObject& cand1, TCPhysObject& cand2, float window)
{
    bool isZCandidate = false;

    if (cand1.Charge() != cand2.Charge()
            && cand1.Type() == cand2.Type()
            && fabs(91.2 - (cand1 + cand2).M()) < window) 
        isZCandidate = true ; 

    return isZCandidate;
}

float* Selector::PhotonEffectiveArea(TCPhysObject* photon)
{
    //Get effective area components
    float EAPho[7][3] = {
        {0.012, 0.030, 0.148}, // eta < 1.0
        {0.010, 0.057, 0.130}, // 1.0 < eta < 1.479
        {0.014, 0.039, 0.112}, // 1.479 < eta < 2.0
        {0.012, 0.015, 0.216}, // 2.0 < eta < 2.2
        {0.016, 0.024, 0.262}, // 2.2 < eta < 2.3
        {0.020, 0.039, 0.260}, // 2.3 < eta < 2.4
        {0.012, 0.072, 0.266} // 2.5 < eta
    };

    float eta = fabs(photon->Eta());
    if (eta < 1.)
        return EAPho[0];
    else if (eta > 1. && eta < 1.479)
        return EAPho[1];
    else if (eta > 1.479 && eta < 2.)
        return EAPho[2];
    else if (eta > 2. && eta < 2.2)
        return EAPho[3];
    else if (eta > 2.2 && eta < 2.3)
        return EAPho[4];
    else if (eta > 2.3 && eta < 2.4)
        return EAPho[5];
    else
        return EAPho[6];
}

float Selector::LeptonEffectiveArea(TCPhysObject* lepton) 
{

    if (lepton->Type() == "electron") {
        float eta = fabs(lepton->Eta());
        if (eta < 1.0)          
            return 0.674;
        else if  (eta < 1.5) 
            return 0.565;
        else if  (eta < 2.0) 
            return 0.442;
        else if  (eta < 2.2) 
            return 0.515;
        else if  (eta < 2.3) 
            return 0.821;
        else if  (eta < 2.4) 
            return 0.660;
        else                               
            return 0.00;
    }
}


//////////////////////
// primary vertices //
//////////////////////

void Selector::PVSelector(TClonesArray* pv)
{
    float bestNDof = 0;

    for (int i = 0; i < pv->GetSize(); ++i) {
        TCPrimaryVtx* pVtx = (TCPrimaryVtx*) pv->At(i);

        if (
                !pVtx->IsFake() 
                && pVtx->NDof()         > 4.
                && fabs(pVtx->z())      <= 24.
                && fabs(pVtx->Perp())   <= 2.
           ) {
            _selVertices.push_back(pVtx);

            if (pVtx->NDof() > bestNDof) {
                bestNDof    = pVtx->NDof();
                _vtxIndex   = _selVertices.size() - 1;
            }
        }
    }
}


///////////
// muons //
///////////


bool Selector::MuonTightID(TCMuon* muon)
{
    bool pass = false;
    if (
            //(muon->IsTRK() && muon->IsGLB() && muon->IsPF())
            //&& muon->PtError()/muon->Pt() < 0.1
            (muon->IsGLB() && muon->IsPF())
            && muon->NormalizedChi2() < 10
            && muon->NumberOfValidMuonHits()  > 0
            && muon->NumberOfMatchedStations() > 1
            && ((
                    fabs(muon->Eta()) < 1.5 
                    && fabs(muon->Dz(_selVertices[0]))  < 0.05
                    && fabs(muon->Dxy(_selVertices[0])) < 0.015
                ) || (
                    fabs(muon->Eta()) > 1.5 
                    && fabs(muon->Dz(_selVertices[0]))  < 0.05
                    && fabs(muon->Dxy(_selVertices[0])) < 0.015 // Should probably reduce this to 0.005
                    ))
            && muon->NumberOfValidPixelHits() > 0
            && muon->TrackLayersWithMeasurement() > 5
       ) pass = true;

    return pass;
}

bool Selector::MuonLooseID(TCMuon* muon)
{
    bool pass = false;
    if (
            muon->NumberOfValidTrackerHits() > 10
            //&& muon->TMLastStationAngTight() <-- Add this to ntuples
            && muon->IsTRK()
            && fabs(muon->Dxy(_selVertices[0])) < 0.2
            && fabs(muon->Dz(_selVertices[0]))  < 0.2 
       ) pass = true;

    return pass;
}


void Selector::MuonSelector(TClonesArray* muons) 
{
    for (int i = 0; i < muons->GetSize(); ++ i) {
        TCMuon* thisMuon = (TCMuon*) muons->At(i);    
        thisMuon->SetType("muon");

        if (fabs(thisMuon->Eta()) > 2.4) continue;	

        // momentum scale corrections (Rochestor corrections)
        TLorentzVector tmpP4 = *thisMuon;

        float muPtErr = 1.;
        if (_isRealData) {
            muCorrector->momcor_data(tmpP4, (float)thisMuon->Charge(), 0, muPtErr);
        } else {
            muCorrector->momcor_mc(tmpP4, (float)thisMuon->Charge(), 0, muPtErr);
        }
        thisMuon->SetPtEtaPhiM(tmpP4.Pt(), tmpP4.Eta(), tmpP4.Phi(), tmpP4.M());

        // isolation
        float muISO = (thisMuon->PfIsoChargedHad() + max(0.,(double)thisMuon->PfIsoNeutral()+ thisMuon->PfIsoPhoton() - 0.5*thisMuon->PfIsoPU()))/thisMuon->Pt();
        thisMuon->SetIdMap("IsoRel", muISO);

        // pt cuts, identification, and isolation
        if (thisMuon->Pt() > _muPtCuts[0]) {
            // QCD dilepton control region tag and probe
            if (
                    _selMuons["QCD2l_CR_tag"].size() == 0
                    && thisMuon->Dz(_selVertices[0]) > 0.2
                    && thisMuon->Dxy(_selVertices[0]) > 0.2 // Replacement for SIP3D inverted cut -- needs to be tuned
                    && muISO > 0.2
               ) {
                _selMuons["QCD2l_CR_tag"].push_back(*thisMuon);
            } else if (
                    MuonTightID(thisMuon) 
                    //&& muISO < 0.6 
                    //&& !(muISO > 0.12 && muISO < 0.2)
                    )
                _selMuons["probe"].push_back(*thisMuon);

            // analysis lepton selection
            if (MuonTightID(thisMuon)) { 
                _selMuons["tight_id"].push_back(*thisMuon);
                if (muISO < 0.12) {
                    _selMuons["tight"].push_back(*thisMuon);
                } else if (muISO > 0.2 && muISO < 0.6) {
                    thisMuon->SetFake(true);
                    _selMuons["fakeable"].push_back(*thisMuon);
                }
            }
        } else if (thisMuon->Pt() > _muPtCuts[1]) 
            if (
                    MuonLooseID(thisMuon)
                    && ((thisMuon->Pt() > 20 && muISO > 0.2) || (thisMuon->Pt() < 20 && muISO < 0.25))
               )
                _selMuons["loose"].push_back(*thisMuon);
    }
}



///////////////
// Electrons //
///////////////


bool Selector::ElectronMVA(TCElectron* electron)
{
    bool pass = false;
    double mvaValue = electron->MvaID_Old();
    electron->SetIdMap("mva", mvaValue);

    if (fabs(electron->Eta()) < 0.8) {
        if (electron->Pt() > 20 && mvaValue > 0.94)
            pass = true;
        else if (electron->Pt() < 20 && mvaValue > 0.)
            pass = true;
    } else if (fabs(electron->Eta()) > 0.8 && fabs(electron->Eta()) < 1.48) {
        if (electron->Pt() > 20 && mvaValue > 0.85)
            pass = true;
        else if (electron->Pt() < 20 && mvaValue > 0.1)
            pass = true;
    } else if (fabs(electron->Eta()) > 1.48 && fabs(electron->Eta()) < 2.5) {
        if (electron->Pt() > 20 && mvaValue > 0.92)
            pass = true;
        else if (electron->Pt() < 20 && mvaValue > 0.62)
            pass = true;
    }

    return pass;
}

bool Selector::ElectronTightID(TCElectron* electron)
{
    bool pass = false;
    if (
            ((fabs(electron->Eta()) < 1.442     
              //&& electron->PtError()/electron->Pt()  < 0.1
              && electron->SigmaIEtaIEta()          < 0.01  
              && fabs(electron->SCDeltaPhi())       < 0.06  
              && fabs(electron->SCDeltaEta())       < 0.004 
              && electron->HadOverEm()              < 0.12      
             ) ||
             (fabs(electron->Eta()) >  1.556  
              //&& electron->PtError()/electron->Pt()  < 0.1
              && electron->SigmaIEtaIEta()          < 0.03  
              && fabs(electron->SCDeltaPhi())       < 0.03  
              && fabs(electron->SCDeltaEta())       < 0.007 
              && electron->HadOverEm()              < 0.1
             ))
            && electron->IdMap("fabsEPDiff")        < 0.05
            && fabs(electron->Dxy(_selVertices[0])) < 0.02
            && fabs(electron->Dz(_selVertices[0]))  < 0.1
            && electron->PassConversionVeto()
       ) pass = true;

    return pass;
}

bool Selector::ElectronLooseID(TCElectron* electron)
{
    bool pass = false;
    if (
            ((fabs(electron->Eta()) < 1.442     
              && electron->SigmaIEtaIEta()          < 0.01  
              && fabs(electron->SCDeltaPhi())       < 0.15  
              && fabs(electron->SCDeltaEta())       < 0.007 
              && electron->HadOverEm()              < 0.12      
             ) ||
             (fabs(electron->Eta()) >  1.556  
              && electron->SigmaIEtaIEta()           < 0.03  
              && fabs(electron->SCDeltaPhi())        < 0.1   
              && fabs(electron->SCDeltaEta())        < 0.009 
              && electron->HadOverEm()               < 0.1
             ))
            && fabs(electron->Dxy(_selVertices[0])) < 0.2
            && fabs(electron->Dz(_selVertices[0]))  < 0.2
            && electron->PassConversionVeto()
       ) pass = true;

    return pass;
}

bool Selector::ElectronMVAPreSel(TCElectron* electron)
{
    bool pass = false;
    if (
            electron->IdMap("gsf_numberOfLostHits") == 0
            && (electron->IdMap("dr03TkSumPt")) / electron->Pt() < 0.2
            && (electron->IdMap("dr03EcalRecHitSumEt")) /electron->Pt() < 0.2
            && (electron->IdMap("dr03HcalTowerSumEt")) / electron->Pt() < 0.2
       ) {

        if (fabs(electron->Eta()) < 1.479) {
            if (electron->SigmaIEtaIEta()< 0.014 && electron->IdMap("hadronicOverEm") < 0.15)
                pass = true;
        } else { //endcap
            if (electron->SigmaIEtaIEta()< 0.035 && electron->IdMap("hadronicOverEm") < 0.10) 
                pass = true;
        }
    }

    return pass;
}

void Selector::ElectronSelector(TClonesArray* electrons) 
{

    for (int i = 0; i <  electrons->GetSize(); ++i) {
        TCElectron* thisElec = (TCElectron*) electrons->At(i);    

        thisElec->SetType("electron");

        // Use regression corrected Pt
        thisElec->SetPtEtaPhiE(thisElec->RegressionMomCombP4().Pt(),
                thisElec->RegressionMomCombP4().Eta(),
                thisElec->RegressionMomCombP4().Phi(),
                thisElec->RegressionMomCombP4().E());

        bool muOverlap = false;
        for (unsigned j = 0; j < _selMuons["tight"].size(); ++j) 
            if (thisElec->DeltaR(_selMuons["tight"][j]) < 0.1) 
                muOverlap = true;

        // electron preselection
        if ( 
                thisElec->Pt() < _elePtCuts[0] 
                || fabs(thisElec->Eta()) > 2.5 
                || !thisElec->PassConversionVeto() 
                || fabs(thisElec->Dz(_selVertices[0])) > 0.05 
                || fabs(thisElec->Dxy(_selVertices[0])) > 0.015
           ) continue;

        float pfPhoIso_corr = ElectronPhoIsoHack(*thisElec);
        float eleISO = (thisElec->PfIsoCharged() + max(0.,(double)thisElec->PfIsoNeutral() 
                    + pfPhoIso_corr - _rho*thisElec->EffArea()))/thisElec->Pt();
        float eleISO_uncorr = (thisElec->PfIsoCharged() + max(0.,(double)thisElec->PfIsoNeutral() 
                    + thisElec->PfIsoPhoton() - _rho*thisElec->EffArea()))/thisElec->Pt();

        thisElec->SetIdMap("IsoRel", eleISO);
        thisElec->SetIdMap("IsoRel_uncorr", eleISO_uncorr);
        thisElec->SetIdMap("pfPhoIso_corr", pfPhoIso_corr);

        //eleISO = eleISO_uncorr;
        if (ElectronLooseID(thisElec)) {
            _selElectrons["loose_id"].push_back(*thisElec);

            if (true || (eleISO < 0.9 && !(eleISO > 0.15 && eleISO < 0.2))) {
                _selElectrons["probe"].push_back(*thisElec);
                if (eleISO > 0.2 && !ElectronMVA(thisElec)){
                    thisElec->SetFake(true);
                    _selElectrons["fakeable"].push_back(*thisElec);
                }
            }
        }

        // analysis electrons
        if (ElectronMVAPreSel(thisElec)) {
            if (ElectronMVA(thisElec)) {
                if (!muOverlap) {
                    _selElectrons["tight_id"].push_back(*thisElec);

                    if (eleISO < 0.15) {
                        _selElectrons["tight"].push_back(*thisElec);			
                    }                 
                } 

                if (eleISO < 0.15 && muOverlap) {
                    _selElectrons["tight_overlap"].push_back(*thisElec);
                }
            }
        } 
    }
}


/////////////
// Photons //
/////////////

bool Selector::PhotonTightID(TCPhoton* photon)
{
    bool pass = false;
    if (
            ((fabs(photon->Eta()) < 1.442     
              && photon->ConversionVeto()           
              && photon->SigmaIEtaIEta()            < 0.01  
              && photon->HadOverEm()                < 0.05      
             ) ||
             (fabs(photon->Eta()) >  1.556  
              && photon->SigmaIEtaIEta()            < 0.028 
              && photon->HadOverEm()                < 0.065
             ))
            && photon->ConversionVeto()           
            //&& photon->R9()                       > 0.94
            //&& fabs(photon->Dxy(_selVertices[0])) < 0.02
            //&& fabs(photon->Dz(_selVertices[0]))  < 0.1
       ) pass = true;

    return pass;
}

bool Selector::PhotonIsolation(TCPhoton* photon)
{
    float* EA = PhotonEffectiveArea(photon);

    float chIsoCor = photon->PfIsoCharged() - _rho*EA[0];
    float nhIsoCor = photon->PfIsoNeutral() - _rho*EA[1];
    float phIsoCor = photon->PfIsoPhoton()  - _rho*EA[2];

    float eta = photon->Eta();
    if (
            eta < 1.442
            && chIsoCor < 1.5
            && nhIsoCor < 1.0 + 0.04*photon->Pt()
            && phIsoCor < 0.7 + 0.005*photon->Pt()
       )
        return true;
    else
        return false;

    if (
            eta > 1.566
            && chIsoCor < 1.2
            && nhIsoCor < 1.5 + 0.04*photon->Pt()
            && phIsoCor < 1.0 + 0.005*photon->Pt()
       )
        return true;
    else
        return false;
}

void Selector::PhotonSelector(TClonesArray* photons) 
{
    for (int i = 0; i <  photons->GetSize(); ++i) {
        TCPhoton* thisPho = (TCPhoton*) photons->At(i);    

        thisPho->SetType("photon");

        // check for electron overlap
        bool elOverlap = false;
        bool muOverlap = false;
        for (unsigned j = 0; j < _selElectrons["tight"].size(); ++j) 
            if (thisPho->DeltaR(_selElectrons["tight"][j]) < 0.1) 
                elOverlap = true;
        for (unsigned j = 0; j < _selMuons["tight"].size(); ++j) 
            if (thisPho->DeltaR(_selMuons["tight"][j]) < 0.1) 
                muOverlap = true;

        // photon preselection
        if (thisPho->Pt() < _phoPtCuts[0] || fabs(thisPho->Eta()) > 2.4) continue;

        bool passIso = PhotonIsolation(thisPho);
        thisPho->SetIdMap("IsoRel", 0.01);

        _selPhotons["noCuts"].push_back(*thisPho);

        // analysis photons
        if (
                PhotonTightID(thisPho) 
                && !elOverlap 
                && !muOverlap
           ) {
            _selPhotons["tight_noIso"].push_back(*thisPho);
            if (passIso) 
                _selPhotons["tight"].push_back(*thisPho);			
        } 
    }
}

//////////
// Jets //
//////////


void Selector::JetSelector(TClonesArray* jets) 
{
    for (int i = 0; i < jets->GetSize(); ++i) {
        TCJet* thisJet = (TCJet*) jets->At(i);     

        // Prevent lepton overlap //
        std::bitset<4> overlap;
        for (int j = 0; j < (int)_selMuons["tight"].size(); ++j) 
            if (thisJet->DeltaR(_selMuons["tight"][j]) < 0.3) overlap.set(0);

        for (int j = 0; j < (int)_selElectrons["tight"].size(); ++j) 
            if (thisJet->DeltaR(_selElectrons["tight"][j]) < 0.3) overlap.set(1);

        for (int j = 0; j < (int)_selMuons["fakeable"].size(); ++j) 
            if (thisJet->DeltaR(_selMuons["fakeable"][j]) < 0.3) overlap.set(2);

        for (int j = 0; j < (int)_selElectrons["fakeable"].size(); ++j) 
            if (thisJet->DeltaR(_selElectrons["fakeable"][j]) < 0.3) overlap.set(3);

        // Apply JER corrections; maybe better to do in the analysis code...
        TCJet *corJet = this->JERCorrections(thisJet);

        if (fabs(corJet->Eta()) < 2.4) {
            if (
                    corJet->Pt() > _jetPtCuts[0]
                    && corJet->NumConstit()  > 1
                    && corJet->NeuHadFrac()  < 0.99
                    && corJet->NeuEmFrac()   < 0.99
                    && corJet->ChHadFrac()   > 0.
                    && corJet->NumChPart()   > 0.
                    && corJet->ChEmFrac()    < 0.99
               ) {

                if (overlap[0]) 
                    _selJets["muJets"].push_back(*corJet);
                else if (overlap[1]) 
                    _selJets["eleJets"].push_back(*corJet);
                else {
                    if (BTagModifier(corJet, "CSVM")) {
                        _selJets["bJetsMedium"].push_back(*corJet);

                        if (!overlap[2] && !overlap[3])
                            _selJets["bJetsMedium_NoFakes"].push_back(*corJet);
                        else if (overlap[2])
                            _selJets["muFakes"].push_back(*corJet);
                        else if (overlap[3])
                            _selJets["eleFakes"].push_back(*corJet);

                    } else if ( true
                            //corJet.VtxNTracks() > 0
                            //&& corJet.VtxSumPtFrac() > 0. 
                            //&& ((int)corJet.VtxSumPtIndex() == 1)
                            ) {
                        _selJets["tight"].push_back(*corJet);

                        if (!overlap[2] && !overlap[3])
                            _selJets["tight_NoFakes"].push_back(*corJet);
                        else if (overlap[2])
                            _selJets["muFakes"].push_back(*corJet);
                        else if (overlap[3])
                            _selJets["eleFakes"].push_back(*corJet);
                    }

                    if ( corJet->BDiscriminatorMap("CSV") > 0.244 && corJet->BDiscriminatorMap("CSV") < 0.679) {
                        _selJets["bJetsLoose"].push_back(*corJet);

                        if (!overlap[2] && !overlap[3])
                            _selJets["bJetsLoose_NoFakes"].push_back(*corJet);
                    }
                }
            }
        } else if (fabs(corJet->Eta()) < 4.7) {
            if (
                    corJet->Pt() > _jetPtCuts[0]
                    && corJet->NumConstit() > 1
                    && corJet->NeuHadFrac() < 0.99
                    && corJet->NeuEmFrac() < 0.99
               ) { 
                if (overlap[0]) 
                    _selJets["muJets"].push_back(*corJet);
                else if (overlap[1]) 
                    _selJets["eleJets"].push_back(*corJet);
                else {
                    _selJets["forward"].push_back(*corJet); 

                    if (!overlap[2] && !overlap[3])
                        _selJets["forward_NoFakes"].push_back(*corJet);
                    else if (overlap[2])
                        _selJets["muFakes"].push_back(*corJet);
                    else if (overlap[3])
                        _selJets["eleFakes"].push_back(*corJet);
                }
            }
        }
    }
}

TCJet* Selector::JERCorrections(TCJet *inJet)
{

    float sfEta[]     = {1.052, 1.057, 1.096, 1.134, 1.288};
    float etaBins[]   = {0., 0.5, 1.1, 1.7, 2.3, 5.0};

    // Match jet to generator level jet
    TCJet*   jet = inJet;
    float    matchedJetPt   = 0.;
    unsigned count          = 0;

    for (unsigned i = 0; i < _selGenJets.size(); ++i) {
        if (inJet->DeltaR(_selGenJets[i]) < 0.15) {
            //cout << _selGenJets[i].Pt() << "\t" << inJet->Pt() << ", " << inJet->DeltaR(_selGenJets[i]) << "\t\t";
            matchedJetPt += _selGenJets[i].Pt();
            ++count;
        }
    }


    if (count > 0) {
        for (unsigned i = 0; i < 5; ++i) {
            if (fabs(inJet->Eta()) > etaBins[i] && fabs(inJet->Eta()) < etaBins[i+1]) {

                if ((matchedJetPt - inJet->Pt())/inJet->Pt() < 0.3*inJet->Pt()) {
                    double corJetPt = matchedJetPt + sfEta[i]*(inJet->Pt() - matchedJetPt);
                    jet->SetPtEtaPhiE(corJetPt, inJet->Eta(), inJet->Phi(), inJet->E());

                    //cout << "\t" << count << ", " << matchedJetPt << ", " << inJet->Pt() << ", " << corJetPt << endl;

                }
            }
        }
    }  

    return jet;
}


bool Selector::BTagModifier(TCJet* jet, string bTag)
{
    float jetPt     = jet->Pt();
    float jetEta    = jet->Eta();
    int   jetFlavor = jet->JetFlavor();
    bool  isBTagged = false;

    // Get b-tagging efficiencies scale factors for jet depending on it's pt 
    float bTagSF       = 1.;
    if (bTag == "CSVL") {
        if (jet->BDiscriminatorMap("CSV") > 0.244) isBTagged = true;
        bTagSF = 0.981149*((1.+(-0.000713295*jetPt))/(1.+(-0.000703264*jetPt)));
    } else if (bTag == "CSVM") {
        if (jet->BDiscriminatorMap("CSV") > 0.679) isBTagged = true;
        bTagSF = 0.726981*((1.+(0.253238*jetPt))/(1.+(0.188389*jetPt)));
    } else if (bTag == "CSVT") {
        if (jet->BDiscriminatorMap("CSV") > 1.) isBTagged = true;
        bTagSF = 0.869965*((1.+(0.0335062*jetPt))/(1.+(0.0304598*jetPt)));
    }

    // Get mistag scale factors dependent on jet pt and eta
    float bMistagSF = 0.;
    if( bTag == "CSVM") {
        if (fabs(jetEta) < 0.8) 
            bMistagSF = 1.07541 + 0.00231827*jetPt - 4.74249e-06*pow(jetPt,2) + 2.70862e-09*pow(jetPt, 3);
        if (fabs(jetEta) > 0.8 && fabs(jetEta) < 1.6) 
            bMistagSF = 1.05613 + 0.00114031*jetPt - 2.56066e-06*pow(jetPt,2) + 1.67792e-09*pow(jetPt,3);
        if (fabs(jetEta) > 1.6 && fabs(jetEta) < 2.4) 
            bMistagSF = 1.05625 + 0.00487231*jetPt - 2.22792e-06*pow(jetPt,2) + 1.70262e-09*pow(jetPt,3);
    } else if( bTag == "CSVL") {
        if (fabs(jetEta) < 0.5) 
            bMistagSF = 1.01177 + 0.00231827*jetPt - 4.74249e-06*pow(jetPt,2) + 2.70862e-09*pow(jetPt, 3); // Need to update non-zeroth order terms
        if (fabs(jetEta) > 0.5 && fabs(jetEta) < 1.) 
            bMistagSF = 0.97596 + 0.00114031*jetPt - 2.56066e-06*pow(jetPt,2) + 1.67792e-09*pow(jetPt,3);
        if (fabs(jetEta) > 1. && fabs(jetEta) < 1.5) 
            bMistagSF = 0.93821 + 0.00114031*jetPt - 2.56066e-06*pow(jetPt,2) + 1.67792e-09*pow(jetPt,3);
        if (fabs(jetEta) > 1.5 && fabs(jetEta) < 2.4) 
            bMistagSF = 1.00022 + 0.00487231*jetPt - 2.22792e-06*pow(jetPt,2) + 1.70262e-09*pow(jetPt,3);
    }

    // Upgrade or downgrade jet
    float rNumber = rnGen->Uniform(1.);

    if (abs(jetFlavor) == 5 || abs(jetFlavor) == 4) {
        float bTagEff;
        if (abs(jetFlavor) == 4) 
            bTagEff = _cTagEff->Eval(jetPt);
        else if (abs(jetFlavor) == 5) 
            bTagEff = _bTagEff->Eval(jetPt);

        if(bTagSF > 1){  // use this if SF>1
            if (!isBTagged) {
                //upgrade to tagged
                float mistagPercent = (1.0 - bTagSF) / (1.0 - (bTagSF/bTagEff) );
                if(rNumber < mistagPercent ) isBTagged = true;
            }
        } else if (bTagSF < 1) {
            //downgrade tagged to untagged
            if( isBTagged && rNumber > bTagSF ) isBTagged = false;
        }

    } else if (abs(jetFlavor) > 0) {
        float mistagEff = _misTagEff->Eval(jetPt);
        if(bMistagSF > 1){  // use this if SF>1
            if (!isBTagged) {
                //upgrade to tagged
                float mistagPercent = (1.0 - bMistagSF) / (1.0 - (bMistagSF/mistagEff));
                if(rNumber < mistagPercent ) isBTagged = true;
            }
        } else if (bMistagSF < 1) {
            //downgrade tagged to untagged
            if( isBTagged && rNumber > bMistagSF ) isBTagged = false;
        }
    }

    return isBTagged;
}


/////////////////////////
// Generator particles //
/////////////////////////


void Selector::GenParticleSelector(TClonesArray* gen, unsigned pdgId, unsigned status, string type)
{
    for (int i = 0; i < gen->GetSize(); ++i) {
        TCGenParticle* iGen = (TCGenParticle*) gen->At(i);

        if (fabs(iGen->GetPDGId()) == pdgId && iGen->GetStatus() == status) {
            //cout << iGen->Pt() << ", " << type << endl;
            iGen->SetType(type);
            _selGenParticles[type].push_back(*iGen);
        }
    }
}

void Selector::GenJetSelector(TClonesArray* genJets)
{
    for (int i = 0; i < genJets->GetSize(); ++i) {
        TCGenJet* iGenJet = (TCGenJet*) genJets->At(i);

        // Specify some cuts based on the type here (maybe...?)

        _selGenJets.push_back(*iGenJet);
    }
}


/////////////////////////////
// Electron isolation hack //
/////////////////////////////

float Selector::ElectronPhoIsoHack(TCElectron& electron)
{
    float pfPhoIso = electron.PfIsoPhoton();
    if (fabs(pfPhoIso - electron.Pt())/electron.Pt() < 0.2) {
        for (unsigned i = 0; i < _selPhotons["noCuts"].size(); ++i) {
            TCPhoton* photon = &_selPhotons["noCuts"][i];
            if (electron.DeltaR(*photon) < 0.15 && fabs(electron.Pt()/photon->Pt() - 1) < 0.3) {
                pfPhoIso = max(0., double(pfPhoIso - photon->Pt()));
            }
        }
    }
    return pfPhoIso;
}
