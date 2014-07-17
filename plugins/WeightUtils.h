/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

// c++ libraries
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// boost libraries
//#include <boost/array.hpp>

// custom libraries
#include "../interface/TCPhysObject.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCJet.h"
#include "../interface/TCGenJet.h"

using namespace std;

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSampleName(string);
        void    SetSelection(string);
        void    SetPassTrigger(string); 
        void    SetObjects(vector<TCPhysObject>&, vector<TCJet>&, float, string);

        float   PUWeight();
        float   RecoWeight();
        float   VBFHiggsWeight(float, int);
        float   GetTotalWeight();
        float   GetFakeWeight(TCPhysObject&, string);
        float   GetFakeUncertainty(TCPhysObject&, string);
        float   GetQFlipWeight(unsigned);
        float   GetAICWeight(const TCPhoton&, const string&);

        // lepton reco efficiencies
        float GetElectronEff(TLorentzVector&) const;
        float GetMuEff(TLorentzVector&) const; 

        // lepton trigger efficiencies
        float GetMuTriggerEff(TLorentzVector&, TLorentzVector&) const;
        float GetEleTriggerEff(TLorentzVector&, TLorentzVector&) const;

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        // input objects
        vector<TCPhysObject>    _leptons;
        vector<TCJet>           _jets;
        string                  _passTrig;
        float                   _nPU;

        //sources
        map<string, TH1D*>  puReweight;

        TGraphErrors *_muSF2012_ID[4];
        TGraphErrors *_muSF2012_ISO[4];

        TH2D    *h2_MuTriggerSFs[2]; // Good for Mu17_Mu8 or Mu17_TkMu8
        TH2D    *h2_EleMVASF;
        TH2D    *h2_DielectronMisQ;

        TGraph  *g_QFlipBB_Low, *g_QFlipEE_Low, *g_QFlipBB_High, *g_QFlipEE_High;
        TGraph  *g_QFlipBB, *g_QFlipBE, *g_QFlipEE;
        map<string, TGraph*>  g_AIC;

        map<string, TGraphAsymmErrors*> g_MuonFakesPtB;
        map<string, TGraphAsymmErrors*> g_MuonFakesPtE;
        map<string, TGraphAsymmErrors*> g_ElectronFakesPtB;
        map<string, TGraphAsymmErrors*> g_ElectronFakesPtG;
        map<string, TGraphAsymmErrors*> g_ElectronFakesPtE;

        map<string, TH2D*> h2_MuonFakes;
        map<string, TH2D*> h2_ElectronFakes;

        map<string, TH1D*> h1_MuonFakes;
        map<string, TH1D*> h1_ElectronFakes;

        //weights
        float _puWeight;
        float _zzWeight;
        float _vbfWeight;
        float _recoWeight;
        float _triggerWeight;

        // error on weights
        float _puWeightErr;
        float _zzWeightErr;
        float _vbfWeightErr;
        float _recoWeightErr;
        float _triggerWeightErr;
        float _qFlipWeightErr;
        float _fakeWeightErr;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
