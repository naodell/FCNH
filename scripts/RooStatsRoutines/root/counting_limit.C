//////////////////////////////////////////////////////////////////
//
// counting_limit.C
//
// This is an example how to use roostats_cl95.C
// to compute a limit on the signal cross section
// for a counting experiment
//
// How to change various options from default values
//
// Gena Kukartsev, May 2012
//
// Usage:
//        - start root
//        - compile and run this routine
//
//          [1] .L counting_limit.C+
//          [2] limit("bayesian")
//          [3] limit("mcmc")
//          [4] limit("cls")
//

#include <iostream>

#include "roostats_cl95.C"


void limit( std::string method = "bayesian" ){
  //
  // compute limit in a counting experiment
  //

  double ilum = 1.0;
  double slum = 0.1;
  double  eff = 1.0;
  double seff = 0.1;
  double  bkg = 1.0;
  double sbkg = 0.0;

  int    data = 2;

  // optional: set some parameters
  SetParameter("Optimize", false);
  SetParameter("CorrelatedLumiSyst", false);
  SetParameter("MakePlot", true);
  SetParameter("GaussianStatistics", false);

  SetParameter("NClsSteps", 10);
  SetParameter("NToys", 1000);
  SetParameter("CalculatorType", 0); // 0 for frequentist
  SetParameter("TestStatType", 3); // LHC-style 1-sided profile likelihood
  SetParameter("Verbosity", 3);
  SetParameter("RandomSeed", 0);

  SetParameter("ConfidenceLevel", 0.95);
  SetParameter("NToysRatio", 2.0);

  // bayesian limit
  if (method.find("bayesian")!=std::string::npos){
    double bayes_limit = GetBayesianLimit(ilum, slum,
					  eff, seff,
					  bkg, sbkg,
					  data);
    
    std::cout << "observed 95% CL bayesian limit = " << bayes_limit << std::endl;


    // expected limit
    // will work with any method, example shown for Bayesian
    LimitResult bayes_expected = GetExpectedLimit(ilum, slum,
						  eff, seff,
						  bkg, sbkg,
						  100,
						  "bayesian");
    
    std::cout << " expected limit (median) " << bayes_expected.GetExpectedLimit() << std::endl;
    std::cout << " expected limit (-1 sig) " << bayes_expected.GetOneSigmaLowRange() << std::endl;
    std::cout << " expected limit (+1 sig) " << bayes_expected.GetOneSigmaHighRange() << std::endl;
    std::cout << " expected limit (-2 sig) " << bayes_expected.GetTwoSigmaLowRange() << std::endl;
    std::cout << " expected limit (+2 sig) " << bayes_expected.GetTwoSigmaHighRange() << std::endl;
  }

  
  // bayesian limit (MCMC)
  if (method.find("mcmc")!=std::string::npos){
    double bayes_limit = GetBayesianLimit(ilum, slum,
					  eff, seff,
					  bkg, sbkg,
					  data,
					  "mcmc");
    
    std::cout << "observed 95% CL bayesian limit = " << bayes_limit << std::endl;
  }

  
  // CLs limit
  if (method.find("cls")!=std::string::npos){
    LimitResult cls_result = GetClsLimit(ilum, slum,
					 eff, seff,
					 bkg, sbkg,
					 data);
    
    double cls_limit = cls_result.GetObservedLimit();
    
    std::cout << "observed 95% CL limit = " << cls_limit << std::endl;
  }



  return;
}

