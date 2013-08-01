#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

const Int_t N_CATEGORIES = 3;
const unsigned nToys = 10;
const string categoryNames[] = {
    "3l_inclusive",
    "3l_OSSF",
    "3l_SSSF"
    //"3l_mumumu", 
    //"3l_emumu", 
    //"3l_eemu",
    //"3l_eee",
    //"ss_inclusive"
    //"ss_mumu", 
    //"ss_emu", 
    //"ss_ee"
};


void simpleLimit()
{

    Double_t lumi       = 19.7e3;                   // luminosity of gathered data
    Double_t lumiErr    = 0.04;                     // luminosity of gathered data
    Double_t sigInit    = 291.6;                    // initial number of signal events accounting for cross-section
    Double_t sigInitErr = 24.99;                    // error on initial number of events

    Double_t sigList[]  = {2.76, 1.98, 0.78};           // Final yields for signal
    Double_t sSigList[] = {0.07, 0.06, 0.04};           // Errors on final yields of signal
    Double_t bckList[]  = {71.54, 51.18, 20.35};          // Final yields for background
    Double_t sbckList[] = {6.36, 5.38, 3.41};           // Errors on final yields of background

    Int_t    nObs[]     = {71, 47, 24};


    for (unsigned i = 0; i < N_CATEGORIES; ++i) {
        Double_t eff    = sigList[i]/sigInit;
        Double_t seff   = weighted_binomial_error(sigList[i], sSigList[i], sigInit, sigInitErr);

        set_limit(lumi, lumiErr*lumi, eff, seff, bckList[i], sbckList[i], nObs[i], nToys, i+1); 
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
        Int_t category) 
{

    // Why C++!? Why!?!?
    stringstream conv;
    conv << category;
    string outName = "test_limits_" + conv.str() + ".txt";

    ofstream outFile(outName.c_str(), ofstream::binary);

    outFile << "\n\t Category " << category << "\n\n";

    Double_t limit0   = roostats_cl95(ilum,slum,eff,seff,bck,sbck,nObs);
    LimitResult limit = roostats_clm(ilum,slum,eff,seff,bck,sbck,nToys);

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

    Double_t brWL  = 0.1080;   // pm 0.0009
    Double_t brHWW = 0.216*0.105;
    Double_t csTOP = 227.;    // pm ???
    //  Double_t csTOP = 167.5;    // pm ???

    Double_t topBR = 2.*csTOP*2.*brWL*brHWW;

    Double_t brUL_exp = exp_limit/topBR;
    Double_t brUL_obs = limit0/topBR;

    outFile << "\nConvert to upper limit on BR(t->Hq) :\n";
    outFile << "\texpected:\t" << brUL_exp << "\n";
    outFile << "\tobserved:\t" << brUL_obs << "\n";

    outFile << "\n\t1-sigma on expected: " << brUL_exp*exp_up/exp_limit << "\t" << brUL_exp*exp_down/exp_limit << "\n";
    outFile << "\t2-sigma on expected: " << brUL_exp*exp_2up/exp_limit << "\t" << brUL_exp*exp_2down/exp_limit << "\n";

    outFile.close();
}

