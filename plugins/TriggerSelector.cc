#include "TriggerSelector.h"

TriggerSelector::TriggerSelector() {}

TriggerSelector::~TriggerSelector() {}

TriggerSelector::TriggerSelector(string type, string dataPeriod, vstring triggerNames, bool triggerDefaults)
{
    _type           = type;
    _dataPeriod     = dataPeriod;
    _triggerNames   = triggerNames;
    _isRealData     = false;

    cout << "\ndata type: " << _type << "\t data period: " << dataPeriod << endl;

    if (triggerDefaults) {
        TriggerDefaults();
        SetSelectedBits();
    }
}

int TriggerSelector::GetEventPrescale() const
{
    return _eventPrescale;
}

vstring TriggerSelector::GetPassNames() const
{
    return _passNames;
}

void TriggerSelector::SetDataBit(bool b) 
{
    _isRealData = b;
}

void TriggerSelector::AddTriggers(vstring trigs)
{
    _triggers.insert(_triggers.begin(), trigs.begin(), trigs.end());
    SetSelectedBits();
}


bool TriggerSelector::CheckPrescales(ULong64_t trigs, UInt_t* hltPrescales)
{
    bool isPrescaled = true;

    for (unsigned i = 0; i < 64; ++i) {
        if (trigs & ULong64_t(0x1) << i) {
            if (hltPrescales[i] == 1) {
                isPrescaled = false;
                break;
            }
        }
    }

    return isPrescaled; 
}

void TriggerSelector::TriggerDefaults()
{
    /*
       Sets triggers for analysis that mixes dilepton
       stream data.  During selection, priority is given
       to sf streams with dimuons taking highest precedent. */

    // double muon triggers
    _triggers.push_back("HLT_Mu17_Mu8_v");
    _triggers.push_back("HLT_Mu17_TkMu8_v");
    _triggers.push_back("HLT_Mu22_TkMu8_v");
    _triggers.push_back("HLT_Mu22_TkMu22_v");

    // double electron triggers
    _triggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

    // muEG triggers
    _triggers.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    _triggers.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

    if (!_isRealData) {
        _triggers.push_back("HLT_Mu17_Ele8_CaloIdL_v");
        _triggers.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v");
        _triggers.push_back("HLT_Mu8_Ele17_CaloIdL_v");
        _triggers.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v");
        _triggers.push_back("HLT_Mu13_Mu8_v");
        _triggers.push_back("HLT_DoubleMu7_v");
        _triggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v"); 
    }

}

void TriggerSelector::SetSelectedBits()
{
    /* Convert list of trigger names to a corresponding
       bit array by comparing to stored trigger names. */

    unsigned count = 0;
    _passTriggers = 0x0;

    cout << endl;
    for (vstring::const_iterator iName = _triggerNames.begin(); iName != _triggerNames.end(); ++iName) {
        for (vstring::const_iterator iTrig = _triggers.begin(); iTrig != _triggers.end(); ++iTrig) {

            if (*iTrig == *iName) {
                cout << *iTrig << endl;
                _passTriggers |= ULong64_t(0x1) << count;
                break;
            }
        }
        ++count;
    }

    cout << endl;

    //cout << hex << _passTriggers << dec << endl;
}

bool TriggerSelector::CheckOverlap() {

    /* Depending on data stream, this method will
       return true if the event should be vetoed, i.e.,
       if a double muon trigger fires for an event from 
       the MuEG stream. */

    //cout << "Checking for overlap...";

    string muonTrigs[]  = {"HLT_Mu17_Mu8_v", "HLT_Mu17_TkMu8_v", "HLT_Mu22_TkMu8_v", "HLT_Mu22_TkMu22_v"};
    string eleTrigs[]   = {"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"};

    bool overlap = false;

    if (_type == "muEG" || _type == "electron" ) {
        for (unsigned i = 0; i < 4; ++i) {
            if (find(_passNames.begin(), _passNames.end(), muonTrigs[i]) != _passNames.end())
                overlap = true;
        }
    }

    if (_type == "muEG") {
        for (unsigned i = 0; i < 1; ++i) {
            if (find(_passNames.begin(), _passNames.end(), eleTrigs[i]) != _passNames.end())
                overlap = true;
        }
    }

    //if (overlap)
    //    cout << "overlap found" << endl;
    //else
    //    cout << "no overlap found" << endl;

    return overlap;
}

bool TriggerSelector::SelectTrigger(ULong64_t triggerStatus, UInt_t* hltPrescales)
{
    _eventPrescale  = -1;
    _eventPass      = false; 
    ULong64_t trigs = triggerStatus & _passTriggers;

    if (_passNames.size() > 0) {
        _passNames.clear();
    }

    if (trigs) _eventPass = true;

    //int  prescale   = 1;

    if (_eventPass) {

    //cout << hex << triggerStatus << "\t" << _passTriggers << dec << endl;

        for (unsigned i = 0; i < _triggerNames.size(); ++i) {
            if (trigs & ULong64_t(0x1) << i) {
                _passNames.push_back(_triggerNames[i]);
            }
        }

        if (_isRealData) { 
            //if (CheckOverlap() || CheckPrescales(trigs, hltPrescales)) _eventPass = false;
            if (CheckOverlap()) _eventPass = false;
        }
    }

    if (_type == "mc") {
        _eventPass = true;
        _passNames.push_back("NADA");
    }

    return _eventPass;
}
