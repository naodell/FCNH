#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

const Int_t N_CATEGORIES = 2;
const unsigned nToys = 10;
const string categoryNames[] = {
    "3l_inclusive",
    //"3l_OSSF",
    //"3l_SSSF",
    //"3l_eee",
    //"3l_eemu",
    //"3l_emumu", 
    //"3l_mumumu", 
    "ss_inclusive"
    //"ss_mumu", 
    //"ss_ee"
    //"ss_emu", 
    //"3l_BDT",
    //"3l_BDT_OSSF",
    //"3l_BDT_SSSF"
};

const Double_t cs_ttbar     = 252.;
const Double_t br_tWb       = 1.;
const Double_t br_tHj       = 0.01;
const Double_t br_HWW       = 0.215;
const Double_t br_HZZ       = 0.215;
const Double_t br_HTauTau   = 0.063;
const Double_t br_Wlnu      = 0.308;
const Double_t br_Zll       = 0.1;
const Double_t br_Znunu     = 0.2;
const Double_t br_Zjj       = 0.7;
const Double_t br_total     = cs_ttbar*(2*br_tWb*br_tHj)*br_Wlnu*(br_HWW*(br_Wlnu*br_Wlnu) + 2*br_Wlnu*(1-br_Wlnu) + br_HZZ*(br_Zll*br_Zll + 2*br_Zll*br_Znunu + 2*br_Zll*br_Zjj) + br_HTauTau);

void simpleLimit()
{

    Double_t lumi       = 19.71e3;                   // luminosity of gathered data
    Double_t lumiErr    = 0.02;                     // luminosity of gathered data
    Double_t sigInit    = 6500;                     // initial number of signal events accounting for cross-section
    //Double_t sigInitErr = 24.99;                  // error on initial number of events
    //Double_t sigInit    = br_total*lumi;            // initial number of signal events accounting for cross-section
    Double_t sigInitErr = 0.1*sqrt(sigInit);            // error on initial number of events

    Double_t sigList[]  = {9.11, 49.84};           // Final yields for signal
    Double_t sSigList[] = {1.15, 3.28};           // Errors on final yields of signal
    Double_t bckList[]  = {33.41, 408.14};          // Final yields for background
    Double_t sbckList[] = {3.71, 16.01};           // Errors on final yields of background

    cout << sigInit << endl;

    Int_t    nObs[]     = {29, 488};

    for (unsigned i = 0; i < N_CATEGORIES; ++i) {
        //if (categoryNames[i] != "ss_inclusive") continue;

        Double_t eff    = sigList[i]/sigInit;

        //Double_t eff    = (sigList[i]*(0.0074/br_total))/sigInit;
        //Double_t seff   = 0.2 * eff; // <--???

        //Double_t eff    = sigList[i]/sigInit;
        Double_t seff   = weighted_binomial_error(sigList[i], sSigList[i], sigInit, sigInitErr);

        set_limit(lumi, lumiErr*lumi, eff, seff, bckList[i], sbckList[i], nObs[i], nToys, i); 
    }
}

Double_t weighted_binomial_error(Double_t nFin, Double_t nFinErr, Double_t nInit, Double_t nInitErr)
{
    float eff = nFin/nInit;
    float sigma_eff = sqrt(((1 - eff)*nFinErr)**2 + nInit*(1 - eff)*eff**2)/nInit;
    return sigma_eff;
}

void set_limit(Double_t ilum, Double_t slum, 
        Double_t eff, Double_t seff, 
        Double_t bck, Double_t sbck, 
        Int_t nObs, Int_t nToys,
        Int_t iCat) 
{

    // Why C++!? Why!?!?
    string category = categoryNames[iCat];
    string outName = "limits_" + category + ".txt";

    ofstream outFile(outName.c_str(), ofstream::binary);

    outFile << "\n\t Category " << category << "\n\n";

    cout << ilum << ", " << slum << ", " << eff << ", " << seff << ", " << bck << ", " << sbck << endl;

    Double_t limit0   = roostats_cl95(ilum,slum,eff,seff,bck,sbck,nObs);
    LimitResult limit = GetExpectedLimit(ilum, slum, eff, seff, bck, sbck, 10, "bayesian");

    Double_t exp_limit = limit.GetExpectedLimit();
    Double_t exp_up    = limit.GetOneSigmaHighRange();
    Double_t exp_down  = limit.GetOneSigmaLowRange();
    Double_t exp_2up   = limit.GetTwoSigmaHighRange();
    Double_t exp_2down = limit.GetTwoSigmaLowRange();

    outFile << "\nexp_limit = " << exp_limit << "\n";
    outFile << "\t1-sigma range: " << exp_up << "\t" << exp_down << "\n";
    outFile << "\t2-sigma range: " << exp_2up << "\t" << exp_2down << "\n";
    outFile << "\nobs_limit = " << limit0 << "\n";


    // very crude check
    Double_t nEvUL = 1.95 * sqrt(bck);
    Double_t csUL = nEvUL/(eff*ilum);
    outFile << "\ncrude check: " << csUL << "\n";

    //
    // convert to an upper limit on the branching ratio
    //

    Double_t brUL_exp = exp_limit/(br_total/br_tHj);
    Double_t brUL_obs = limit0/(br_total/br_tHj);

    outFile << "\nConvert to upper limit on BR(t->Hq) :\n";
    outFile << "\texpected:\t" << brUL_exp << "\n";
    outFile << "\tobserved:\t" << brUL_obs << "\n";

    outFile << "\n\t1-sigma on expected: " << brUL_exp*exp_up/exp_limit << "\t" << brUL_exp*exp_down/exp_limit << "\n";
    outFile << "\t2-sigma on expected: " << brUL_exp*exp_2up/exp_limit << "\t" << brUL_exp*exp_2down/exp_limit << "\n";

    outFile.close();
}

