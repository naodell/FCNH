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
    _vtxIndex   = 0;
    _isRealData = false;

    // b-tag mc efficiencies
    TFile* f_bEff = new TFile("../data/bEff_ttbar_2012.root");
    _misTagEff  = (TGraphAsymmErrors*)f_bEff->Get("g_MistagEff");
    _bTagEff    = (TGraphAsymmErrors*)f_bEff->Get("g_bTagEff");

    // Initialize electron MVA tool 
    std::vector<std::string> WeightsMVA;
    WeightsMVA.push_back("../data/weights/Electrons_BDTG_TrigV0_Cat1.weights.xml");
    WeightsMVA.push_back("../data/weights/Electrons_BDTG_TrigV0_Cat2.weights.xml");
    WeightsMVA.push_back("../data/weights/Electrons_BDTG_TrigV0_Cat3.weights.xml");
    WeightsMVA.push_back("../data/weights/Electrons_BDTG_TrigV0_Cat4.weights.xml");
    WeightsMVA.push_back("../data/weights/Electrons_BDTG_TrigV0_Cat5.weights.xml");
    WeightsMVA.push_back("../data/weights/Electrons_BDTG_TrigV0_Cat6.weights.xml");

    electronMVA = new EGammaMvaEleEstimator();
    electronMVA->initialize("BDT", EGammaMvaEleEstimator::kTrig, true, WeightsMVA);

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

bool Selector::IsZCandidate(TCPhysObject* cand1, TCPhysObject* cand2, float window)
{
    bool isZCandidate = false;

    if (cand1->Charge() != cand2->Charge()
            && cand1->Type() == cand2->Type()
            && fabs(91.2 - (*cand1 + *cand2).M()) < window) 
        isZCandidate = true ; 

    return isZCandidate;
}

float Selector::EffectiveArea(TCPhysObject *lepton) 
{
    float area = 0;

    if (fabs(lepton->Eta())<1.0)          
        area = 0.674;
    else if  (fabs(lepton->Eta())<1.5) 
        area = 0.565;
    else if  (fabs(lepton->Eta())<2.0) 
        area = 0.442;
    else if  (fabs(lepton->Eta())<2.2) 
        area = 0.515;
    else if  (fabs(lepton->Eta())<2.3) 
        area = 0.821;
    else if  (fabs(lepton->Eta())<2.4) 
        area = 0.660;
    else                               
        area = 0.00;

    return area;
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
            (muon->IsTRK() && muon->IsGLB() && muon->IsPF())
            && muon->NormalizedChi2()  < 10
            && muon->NumberOfValidMuonHits()  > 0
            && muon->NumberOfMatchedStations() > 1
            && muon->NumberOfValidPixelHits() > 0
            && muon->TrackLayersWithMeasurement() > 5
            && fabs(muon->Dz(_selVertices[0]))  < 0.05 
            && fabs(muon->Dxy(_selVertices[0])) < 0.01

            //&& muon->NumberOfMatches() > 1
            //&& muon->NumberOfValidTrackerHits() > 10 // Possibly not valid in 2012
            //&& muon->PtError()/muon->Pt() < 0.1

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

    //cout << "Muons(" << muons->GetSize() << "): ";

    for (int i = 0; i < muons->GetSize(); ++ i) {
        TCMuon* thisMuon = (TCMuon*) muons->At(i);    

        if (fabs(thisMuon->Eta()) > 2.4) continue;	

        thisMuon->SetType("muon");

        // momentum scale corrections (Rochestor corrections)
        TLorentzVector tmpP4 = *thisMuon;
        float muPtErr = 1.;
        if (_isRealData) {
            muCorrector->momcor_data(tmpP4, (float)thisMuon->Charge(), 0, muPtErr);
        } else {
            muCorrector->momcor_data(tmpP4, (float)thisMuon->Charge(), 0, muPtErr);
        }

        thisMuon->SetPtEtaPhiM(tmpP4.Pt(), tmpP4.Eta(), tmpP4.Phi(), tmpP4.M());

        // isolation
        float muISO = 0.;
        muISO = (
                thisMuon->IsoMap("pfChargedHadronPt_R04") + TMath::Max(0.0, (double)thisMuon->IsoMap("pfPhotonEt_R04") 
                    + thisMuon->IsoMap("pfNeutralHadronEt_R04") - TMath::Max(0.0, (double)_rho*EffectiveArea(thisMuon)))
                )/thisMuon->Pt();

        thisMuon->SetIsoMap("IsoRel", muISO);
        //cout << "(" << thisMuon->Pt() << ", " << thisMuon->Eta() << "),\t";

        // pt cuts, identification, and isolation
        if (thisMuon->Pt() > _muPtCuts[0]) {

            // QCD dilepton control region tag and probe
            if (
                    sqrt(pow(thisMuon->Dz(_selVertices[0]), 2) + pow(thisMuon->Dxy(_selVertices[0]), 2)) > 1. // Replacement for SIP3D inverted cut -- needs to be tuned
                    && muISO > 0.2
                    ) 
                _selMuons["QCD2l_CR_tag"].push_back(*thisMuon);
            else if (
                    thisMuon->IsPF()
                    ) 
                _selMuons["QCD2l_CR_probe"].push_back(*thisMuon);

            // analysis lepton selection
            if (MuonTightID(thisMuon) && muISO < 0.12) 
                _selMuons["tight"].push_back(*thisMuon);
            else if (thisMuon->IsPF())
                _selMuons["fakeable"].push_back(*thisMuon);

        } else if ( thisMuon->Pt() > _muPtCuts[1]  ) 
                if (MuonLooseID(thisMuon)
                //&& (muISO > 0.1 && thisMuon->Pt() > 20)
                //&& (muISO < 0.15 && thisMuon->Pt() < 20)
                   )
            _selMuons["loose"].push_back(*thisMuon);
    }

    //cout << endl;
}



///////////////
// Electrons //
///////////////


bool Selector::ElectronMVA(TCElectron* electron)
{
    bool pass = false;

    double mvaValue = electronMVA->mvaValue( electron->IdMap("fbrem"),
            electron->IdMap("kfChi2"), electron->IdMap("kfNLayers"),
            electron->IdMap("gsfChi2"), electron->IdMap("dEta"),
            electron->IdMap("dPhi"), electron->IdMap("dEtaAtCalo"),
            electron->SigmaIEtaIEta(), electron->IdMap("SigmaIPhiIPhi"),
            electron->IdMap("SCEtaWidth"), electron->IdMap("SCPhiWidth"),
            electron->IdMap("ome1x5oe5x5"), electron->IdMap("R9"),
            electron->HadOverEm(), electron->IdMap("EoP"),
            electron->IdMap("ooemoopV1"), electron->IdMap("eopOut"),
            electron->IdMap("preShowerORaw"), electron->IdMap("d0"),
            electron->IdMap("ip3d"), electron->SCEta(), electron->Pt(), false);                

    electron->SetIdMap("mva", mvaValue);

    if (fabs(electron->Dz(_selVertices[0]))  < 0.1 && fabs(electron->Dxy(_selVertices[0])) < 0.015) {
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
    }

    return pass;
}

bool Selector::ElectronTightID(TCElectron* electron)
{
    bool pass = false;
    if (
            ((fabs(electron->Eta()) < 1.442     
              //&& electron->PtError()/electron->Pt()  < 0.1
              && electron->SigmaIEtaIEta()           < 0.01  
              && fabs(electron->DphiSuperCluster())  < 0.06  
              && fabs(electron->DetaSuperCluster())  < 0.004 
              && electron->HadOverEm()               < 0.12      
             ) ||
             (fabs(electron->Eta()) >  1.556  
              //&& electron->PtError()/electron->Pt()  < 0.1
              && electron->SigmaIEtaIEta()           < 0.03  
              && fabs(electron->DphiSuperCluster())  < 0.03  
              && fabs(electron->DetaSuperCluster())  < 0.007 
              && electron->HadOverEm()               < 0.1
             ))
            && electron->IdMap("fabsEPDiff")       < 0.05
            && fabs(electron->Dxy(_selVertices[0])) < 0.02
            && fabs(electron->Dz(_selVertices[0]))  < 0.1
            && electron->ConversionVeto()
       ) pass = true;

    return pass;
}

bool Selector::ElectronLooseID(TCElectron* electron)
{
    bool pass = false;
    if (
            ((fabs(electron->Eta()) < 1.442     
              && electron->SigmaIEtaIEta()           < 0.01  
              && fabs(electron->DphiSuperCluster())  < 0.15  
              && fabs(electron->DetaSuperCluster())  < 0.007 
              && electron->HadOverEm()               < 0.12      
             ) ||
             (fabs(electron->Eta()) >  1.556  
              && electron->SigmaIEtaIEta()           < 0.03  
              && fabs(electron->DphiSuperCluster())  < 0.1   
              && fabs(electron->DetaSuperCluster())  < 0.009 
              && electron->HadOverEm()               < 0.1
             ))
            && fabs(electron->Dxy(_selVertices[0])) < 0.2
            && fabs(electron->Dz(_selVertices[0]))  < 0.2
            && electron->ConversionVeto()
       ) pass = true;

    return pass;
}

void Selector::ElectronSelector(TClonesArray* electrons) 
{
    //cout << "Electrons (" << electrons->GetSize() << "): ";

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

        //cout << "(" << thisElec->Pt() << ", " << thisElec->Eta() << "),\t";

        // electron preselection
        if ( thisElec->Pt() < _elePtCuts[0] || fabs(thisElec->Eta()) > 2.5 ) continue;

        float eleISO = (thisElec->IsoMap("pfChIso_R04") + max(0., (double)(thisElec->IsoMap("pfPhoIso_R04") 
                        + thisElec->IsoMap("pfNeuIso_R04") - _rho*thisElec->IsoMap("EffArea_R04"))))/thisElec->Pt(); 

        thisElec->SetIsoMap("IsoRel", eleISO);

        // analysis electrons
        //if (ElectronTightID(thisElec)) 

        if (thisElec->IdMap("preSelPassV1")) {

            _selElectrons["QCD2l_CR_probe"].push_back(*thisElec);

            if (ElectronMVA(thisElec) && !muOverlap)
                _selElectrons["premva"].push_back(*thisElec);

            if (ElectronMVA(thisElec) && eleISO < 0.15) 
                if (!muOverlap)
                    _selElectrons["tight"].push_back(*thisElec);			
                else
                    _selElectrons["tight_overlap"].push_back(*thisElec);			
            else if (!muOverlap)
                _selElectrons["fakeable"].push_back(*thisElec);

        } else if (
                ElectronLooseID(thisElec)
                && (thisElec->Pt() > 20 && eleISO > 0.20)
                && !muOverlap
                ) _selElectrons["loose"].push_back(*thisElec);
    }
    //cout << endl;
}


/////////////
// Photons //
/////////////

bool Selector::PhotonTightID(TCPhoton* photon)
{
    bool pass = false;
    if (
            ((fabs(photon->Eta()) < 1.442     
              && photon->SigmaIEtaIEta()            < 0.01  
              && photon->HadOverEm()                < 0.05      
             ) ||
             (fabs(photon->Eta()) >  1.556  
              && photon->SigmaIEtaIEta()            < 0.028 
              && photon->HadOverEm()                < 0.065
             ))
            && photon->R9()                       > 0.94
            //&& fabs(photon->Dxy(_selVertices[0])) < 0.02
            //&& fabs(photon->Dz(_selVertices[0]))  < 0.1
            //&& photon->ConversionVeto()
       ) pass = true;

    return pass;
}

void Selector::PhotonSelector(TClonesArray* photons) 
{

    for (int i = 0; i <  photons->GetSize(); ++i) {
        TCPhoton* thisPho = (TCPhoton*) photons->At(i);    

        thisPho->SetType("photon");

        // photon preselection
        if (thisPho->Pt() < _phoPtCuts[0] || fabs(thisPho->Eta()) > 2.5 ) continue;


        float phoISO = (thisPho->IsoMap("chIso03") 
                + max(0., (double)(thisPho->IsoMap("nhIso03") 
                        + thisPho->IsoMap("phIso03") 
                        - _rho*EffectiveArea(thisPho))))/thisPho->Pt(); 


        // analysis photons
        if (PhotonTightID(thisPho) && phoISO < 0.15) {
            _selPhotons["tight"].push_back(*thisPho);			

        } else {
            _selPhotons["loose"].push_back(*thisPho);

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
            if (thisJet->DeltaR(_selMuons["tight"][j]) < 0.5) overlap.set(0);

        for (int j = 0; j < (int)_selElectrons["tight"].size(); ++j) 
            if (thisJet->DeltaR(_selElectrons["tight"][j]) < 0.5) overlap.set(1);

        for (int j = 0; j < (int)_selMuons["fakeable"].size(); ++j) 
            if (thisJet->DeltaR(_selMuons["fakeable"][j]) < 0.5) overlap.set(2);

        for (int j = 0; j < (int)_selElectrons["fakeable"].size(); ++j) 
            if (thisJet->DeltaR(_selElectrons["fakeable"][j]) < 0.5) overlap.set(3);

        // Apply JER corrections; maybe better to do in the analysis code...
        TCJet corJet = this->JERCorrections(thisJet);

        if (fabs(corJet.Eta()) < 2.4) {
            if (
                    corJet.Pt() > _jetPtCuts[0]
                    && corJet.NumConstit() > 1
                    && corJet.NeuHadFrac() < 0.99
                    && corJet.NeuEmFrac() < 0.99
                    && corJet.ChHadFrac() > 0.
                    && corJet.NumChPart() > 0.
                    && corJet.ChEmFrac() < 0.99
               ) {

                if (overlap[0]) 
                    _selJets["muJets"].push_back(corJet);
                else if (overlap[1]) 
                    _selJets["eleJets"].push_back(corJet);
                else {
                    if (BTagModifier(corJet, "CSVM")) {
                        _selJets["bJetsMedium"].push_back(corJet);

                        if (!overlap[2] && !overlap[3])
                            _selJets["bJetsMedium_NoFakes"].push_back(corJet);

                    } else if (BTagModifier(corJet, "CSVL")) {
                        _selJets["bJetsLoose"].push_back(corJet);
                        if (!overlap[2] && !overlap[3])
                            _selJets["bJetsLoose_NoFakes"].push_back(corJet);

                    } else if (
                            corJet.VtxNTracks() > 0
                            && corJet.VtxSumPtFrac() > 0. 
                            && ((int)corJet.VtxSumPtIndex() == 1)
                            ) { 
                        _selJets["tight"].push_back(corJet);

                        if (!overlap[2] && !overlap[3])
                            _selJets["tight_NoFakes"].push_back(corJet);
                    }
                }
            }
        } else if (fabs(corJet.Eta()) < 3.5) {
            if (
                    corJet.Pt() > _jetPtCuts[0]
                    && corJet.NumConstit() > 1
                    && corJet.NeuHadFrac() < 0.99
                    && corJet.NeuEmFrac() < 0.99
               ) { 
                    if (overlap[0]) 
                        _selJets["muJets"].push_back(corJet);
                    else if (overlap[1]) 
                        _selJets["eleJets"].push_back(corJet);
                    else {
                        _selJets["forward"].push_back(corJet); 
                        _selJets["tight_NoFakes"].push_back(corJet);
                    }
            }
        }
    }
}

TCJet Selector::JERCorrections(TCJet *inJet)
{

    float sfEta[]     = {1.052, 1.057, 1.096, 1.134, 1.288};
    float etaBins[]   = {0., 0.5, 1.1, 1.7, 2.3, 5.0};

    // Match jet to generator level jet
    TCJet    jet = *inJet;
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
                    jet.SetPtEtaPhiE(corJetPt, inJet->Eta(), inJet->Phi(), inJet->E());

                    //cout << "\t" << count << ", " << matchedJetPt << ", " << inJet->Pt() << ", " << corJetPt << endl;

                }
            }
        }
    }  

    return jet;
}


bool Selector::BTagModifier(TCJet jet, string bTag)
{
    float jetPt     = jet.Pt();
    float jetEta    = jet.Eta();
    int   jetFlavor = jet.JetFlavor();
    bool  isBTagged = false;

    // Get b-tagging efficiencies scale factors for jet depending on it's pt 
    float bTagSF       = 1.;
    if (bTag == "CSVL") {
        if (jet.BDiscriminatorMap("CSV") > 0.244) isBTagged = true;
        bTagSF = 0.981149*((1.+(-0.000713295*jetPt))/(1.+(-0.000703264*jetPt)));
    } else if (bTag == "CSVM") {
        if (jet.BDiscriminatorMap("CSV") > 0.679) isBTagged = true;
        bTagSF = 0.726981*((1.+(0.253238*jetPt))/(1.+(0.188389*jetPt)));
    } else if (bTag == "CSVT") {
        if (jet.BDiscriminatorMap("CSV") > 1.) isBTagged = true;
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
        float bTagEff   = _bTagEff->Eval(jetPt);
        if (abs(jetFlavor) == 4) bTagSF = bTagSF/5.;
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

        if (fabs(iGen->GetPDGId()) == pdgId and iGen->GetStatus() == status) 
            _selGenParticles[type].push_back(*iGen);
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
