//=====================================================================
//                                                                    
//      roostats_twobin.C                                     
//
//      version 1.0-rc1
//                                                                    
// Set up the model for a combination of two coupled 
// 'counting experiments' in a signal and a background-dominated 
// control regions, with systematic uncertainties                                      
//                                                                    
// Authors: Gena Kukartsev, Fedor Ratnikov
//                                                                    
// May 2011                                                           
//                                                                    
//=====================================================================
//
// The original specifications from Fedor Ratnicov
//
//- we have two channels counting experiment: "signal" region and
//  "control" region.
//- in the signal region there are two background contribution:
//  - independently estimated b_signal_1 with its uncertainty
//  - b_signal_2 extrapolated from the control region
//- background contributions to control region b_control and
//  another contribution to the signal region b_signal_2 are
//  related by tau-factor:   b_signal_2 = tau * b_control
//- tau is independently estimated with its uncertainty
//- signal contributes to both signal (stronger contribution) and
//  control (weaker contribution) areas, corresponding signal
//  yields yield_signal and yield_background are independently
//  estimated with their uncertainties
//- there is luminosity uncertainty affecting both expected signal
//  contributions in correlated way.
//
//Thus, in the model we have:
//- parameter of interest: signal "s"
//- model parameters: b_signal_1, tau,  yield_signal,
//  yield_background, lumi, all with corresponding uncertainties
//- model variable: b_control, which is unknown a priori
//- observables: n_signal, n_control
//- counting experiment likelihood is then
//    Poisson(n_signal,
//            s*yield_signal*lumi+b_signal1+tau*b_control)*
//    Poisson(n_control, s*yield_control*lumi+b_control)
//
//In total, one poi, one variable, 5 nuisances, two measurements
//
//
//////////////////////////////////
//
// More math-oriented definition of the model
//
// a,b          - signal and control regions
//
// na, nb       - observables
// Na, Nb       - observed values
// xsec         - POI
// bg_b         - data-driven nuisance parameter
// xsec_bg_a    - constant (could be a nuisance parameter if uncertainty is not negligible)
// lumi, eff_a, eff_b, tau                           - nuisance parameters
// nom_lumi, nom_eff_a, nom_eff_b, nom_tau, nom_bg_a - best estimated values from auxillary measurements
// d_lumi, d_eff_a, d_eff_b, d_tau, d_bg_a           - relative uncertainties from auxillary measurements
//
// mu_a, mu_b   - Poisson means for a, b
// mu_a = sig_a + bg_a + tau*bg_b // 4 + 5 + 0.5*10
// mu_b = sig_b + bg_b            // 1 + 10
// sig_a = lumi*xsec*eff_a        // 100*0.2*0.2
// sig_b = lumi*xsec*eff_b        // 100*0.2*0.05
// bg_a  = lumi*xsec_bg_a         // 100*0.05
// bg_b                           // 10
//
// P(Na|mu_a)*P(Nb|mu_b)*PRODUCT[L(nom_X|X,d_X)]
// 
// (product over likelihoods of auxillary measurements for all nuisance parameters

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooArgSet.h"

#include "RooStats/RooStatsUtils.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"

using namespace RooFit;
using namespace RooStats;

// declarations
void TwoBinInstructional( void );
void SetConstants(RooWorkspace * w, RooStats::ModelConfig * mc);
void SetConstant(const RooArgSet * vars, Bool_t value );

// implementation
void TwoBinInstructional( void ){
  
  // let's time this example
  TStopwatch t;
  t.Start();

  // set RooFit random seed for reproducible results
  RooRandom::randomGenerator()->SetSeed(4357);

  // make model
  RooWorkspace * pWs = new RooWorkspace("ws");

  // derived from data
  pWs->factory("xsec[0.2,0,2]"); // POI
  pWs->factory("bg_b[10,0,50]");    // data driven nuisance

  // predefined nuisances
  pWs->factory("lumi[100,0,1000]");
  pWs->factory("eff_a[0.2,0,1]");
  pWs->factory("eff_b[0.05,0,1]");
  pWs->factory("tau[0,1]");
  pWs->factory("xsec_bg_a[0.05]"); // constant
  pWs->var("xsec_bg_a")->setConstant(1);

  // channel a (signal): lumi*xsec*eff_a + lumi*bg_a + tau*bg_b
  pWs->factory("prod::sig_a(lumi,xsec,eff_a)");
  pWs->factory("prod::bg_a(lumi,xsec_bg_a)");
  pWs->factory("prod::tau_bg_b(tau, bg_b)");
  pWs->factory("Poisson::pdf_a(na[14,0,100],sum::mu_a(sig_a,bg_a,tau_bg_b))");

  // channel b (control): lumi*xsec*eff_b + bg_b
  pWs->factory("prod::sig_b(lumi,xsec,eff_b)");
  pWs->factory("Poisson::pdf_b(nb[11,0,100],sum::mu_b(sig_b,bg_b))");

  // nuisance constraint terms (systematics)
  pWs->factory("Lognormal::l_lumi(lumi,nom_lumi[100,0,1000],sum::kappa_lumi(1,d_lumi[0.1]))");
  pWs->factory("Lognormal::l_eff_a(eff_a,nom_eff_a[0.20,0,1],sum::kappa_eff_a(1,d_eff_a[0.05]))");
  pWs->factory("Lognormal::l_eff_b(eff_b,nom_eff_b[0.05,0,1],sum::kappa_eff_b(1,d_eff_b[0.05]))");
  pWs->factory("Lognormal::l_tau(tau,nom_tau[0.50,0,1],sum::kappa_tau(1,d_tau[0.05]))");
  //pWs->factory("Lognormal::l_bg_a(bg_a,nom_bg_a[0.05,0,1],sum::kappa_bg_a(1,d_bg_a[0.10]))");

  // complete model PDF
  pWs->factory("PROD::model(pdf_a,pdf_b,l_lumi,l_eff_a,l_eff_b,l_tau)");

  // Now create sets of variables. Note that we could use the factory to
  // create sets but in that case many of the sets would be duplicated
  // when the ModelConfig objects are imported into the workspace. So,
  // we create the sets outside the workspace, and only the needed ones
  // will be automatically imported by ModelConfigs

  // observables
  RooArgSet obs(*pWs->var("na"), *pWs->var("nb"), "obs");

  // global observables
  RooArgSet globalObs(*pWs->var("nom_lumi"), *pWs->var("nom_eff_a"), *pWs->var("nom_eff_b"), 
		      *pWs->var("nom_tau"),
		      "global_obs");

  // parameters of interest
  RooArgSet poi(*pWs->var("xsec"), "poi");

  // nuisance parameters
  RooArgSet nuis(*pWs->var("lumi"), *pWs->var("eff_a"), *pWs->var("eff_b"), *pWs->var("tau"), "nuis");

  // priors (for Bayesian calculation)
  pWs->factory("Uniform::prior_xsec(xsec)"); // for parameter of interest
  pWs->factory("Uniform::prior_bg_b(bg_b)"); // for data driven nuisance parameter
  pWs->factory("PROD::prior(prior_xsec,prior_bg_b)"); // total prior

  // create data
  pWs->var("na")->setVal(14);
  pWs->var("nb")->setVal(11);
  RooDataSet * pData = new RooDataSet("data","",obs);
  pData->add(obs);
  pWs->import(*pData);
  //pData->Print();

  // signal+background model
  ModelConfig * pSbModel = new ModelConfig("SbModel");
  pSbModel->SetWorkspace(*pWs);
  pSbModel->SetPdf(*pWs->pdf("model"));
  pSbModel->SetPriorPdf(*pWs->pdf("prior"));
  pSbModel->SetParametersOfInterest(poi);
  pSbModel->SetNuisanceParameters(nuis);
  pSbModel->SetObservables(obs);
  pSbModel->SetGlobalObservables(globalObs);

  // set all but obs, poi and nuisance to const
  SetConstants(pWs, pSbModel);
  pWs->import(*pSbModel);


  // background-only model
  // use the same PDF as s+b, with xsec=0
  // POI value under the background hypothesis
  Double_t poiValueForBModel = 0.0;
  ModelConfig* pBModel = new ModelConfig(*(RooStats::ModelConfig *)pWs->obj("SbModel"));
  pBModel->SetName("BModel");
  pBModel->SetWorkspace(*pWs);
  pWs->import(*pBModel);


  // find global maximum with the signal+background model
  // with conditional MLEs for nuisance parameters
  // and save the parameter point snapshot in the Workspace
  //  - safer to keep a default name because some RooStats calculators
  //    will anticipate it
  RooAbsReal * pNll = pSbModel->GetPdf()->createNLL(*pData);
  RooAbsReal * pProfile = pNll->createProfile(RooArgSet());
  pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
  RooArgSet * pPoiAndNuisance = new RooArgSet();
  if(pSbModel->GetNuisanceParameters())
    pPoiAndNuisance->add(*pSbModel->GetNuisanceParameters());
  pPoiAndNuisance->add(*pSbModel->GetParametersOfInterest());
  cout << "\nWill save these parameter points that correspond to the fit to data" << endl;
  pPoiAndNuisance->Print("v");
  pSbModel->SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;

  // Find a parameter point for generating pseudo-data
  // with the background-only data.
  // Save the parameter point snapshot in the Workspace
  pNll = pBModel->GetPdf()->createNLL(*pData);
  pProfile = pNll->createProfile(poi);
  ((RooRealVar *)poi.first())->setVal(poiValueForBModel);
  pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
  pPoiAndNuisance = new RooArgSet();
  if(pBModel->GetNuisanceParameters())
    pPoiAndNuisance->add(*pBModel->GetNuisanceParameters());
  pPoiAndNuisance->add(*pBModel->GetParametersOfInterest());
  cout << "\nShould use these parameter points to generate pseudo data for bkg only" << endl;
  pPoiAndNuisance->Print("v");
  pBModel->SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;

  // inspect workspace
  pWs->Print();

  // save workspace to file
  pWs->writeToFile("ws_twobin.root");

  // clean up
  delete pWs;
  delete pData;
  delete pSbModel;
  delete pBModel;

} // ----- end of tutorial ----------------------------------------



// helper functions

void SetConstants(RooWorkspace * pWs, RooStats::ModelConfig * pMc){
  //
  // Fix all variables in the PDF except observables, POI and
  // nuisance parameters. Note that global observables are fixed.
  // If you need global observables floated, you have to set them
  // to float separately.
  //

  pMc->SetWorkspace(*pWs);

  RooAbsPdf * pPdf = pMc->GetPdf(); // we do not own this

  RooArgSet * pVars = pPdf->getVariables(); // we do own this

  RooArgSet * pFloated = new RooArgSet(*pMc->GetObservables());
  pFloated->add(*pMc->GetParametersOfInterest());
  pFloated->add(*pMc->GetNuisanceParameters());

  TIterator * pIter = pVars->createIterator(); // we do own this

  for(TObject * pObj = pIter->Next(); pObj; pObj = pIter->Next() ){
    std::string _name = pObj->GetName();
    RooRealVar * pFloatedObj = (RooRealVar *)pFloated->find(_name.c_str());
    if (pFloatedObj){
      ((RooRealVar *)pObj)->setConstant(kFALSE);
    }
    else{
      ((RooRealVar *)pObj)->setConstant(kTRUE);
    }
    //pObj->Print();
  }

  delete pIter;
  delete pVars;
  delete pFloated;

  return;
}



void SetConstant(const RooArgSet * vars, Bool_t value ){
  //
  // Set the constant attribute for all vars in the set
  //

  TIterator * pIter = vars->createIterator(); // we do own this

  for(TObject * pObj = pIter->Next(); pObj; pObj = pIter->Next() ){
    ((RooRealVar *)pObj)->setConstant(value);
    //pObj->Print();
  }

  delete pIter;

  return;
}
