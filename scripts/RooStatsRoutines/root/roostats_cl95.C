static const char* desc =
"=====================================================================\n"
"|                                                                    \n"
"|\033[1m        roostats_cl95.C                               \033[0m\n"
"|                                                                    \n"
"| Standard c++ routine for 95% C.L. limit calculation                \n"
"| for cross section in a 'counting experiment'                       \n"
"|                                                                    \n"
"|\033[1m Gena Kukartsev                                       \033[0m\n"
"|\033[1m Lorenzo Moneta (CLs core)                            \033[0m\n"
"|\033[1m Michael Segala (Feldman-Cousins)                     \033[0m\n"
"|\033[1m Stefan Schmitz, Gregory Schott                       \033[0m\n"
"|                                                                    \n"
"| For more information, see                                          \n"
"|                                                                    \n"
"|       https://twiki.cern.ch/twiki/bin/view/CMS/RooStatsCl95        \n"
"|                                                                    \n"
"=====================================================================\n"
"                                                                     \n";
//
//
//Prerequisites:                                                       
//                ROOT version 5.32.00 or higher                       
//                                                                     
//                                                                     
//                                                                     
//The code should be compiled in ROOT:                                 
//                                                                     
//root -l                                                              
//                                                                     
//.L roostats_cl95.C+                                                  
//                                                                     
//Usage:                                                               
// Double_t             limit = roostats_cl95(ilum, slum, eff, seff, bck, sbck, n, gauss = false, nuisanceModel, method, plotFileName, seed); 
// LimitResult expected_limit = roostats_clm(ilum, slum, eff, seff, bck, sbck, ntoys, nuisanceModel, method, seed); 
// Double_t     average_limit = roostats_cla(ilum, slum, eff, seff, bck, sbck, nuisanceModel, method, seed); 
//                                                                     
// LimitResult limit = roostats_limit(ilum, slum, eff, seff, bck, sbck, n, gauss = false, nuisanceModel, method, plotFileName, seed); 
// Double_t obs_limit = limit.GetObservedLimit();                      
// Double_t exp_limit = limit.GetExpectedLimit();                      
// Double_t exp_up    = limit.GetOneSigmaHighRange();                  
// Double_t exp_down  = limit.GetOneSigmaLowRange();                   
// Double_t exp_2up   = limit.GetTwoSigmaHighRange();                  
// Double_t exp_2down = limit.GetTwoSigmaLowRange();                   
//                                                                     
//Inputs:                                                              
//       ilum          - Nominal integrated luminosity (pb^-1)         
//       slum          - Absolute error on the integrated luminosity   
//       eff           - Nominal value of the efficiency times         
//                       acceptance (in range 0 to 1)                  
//       seff          - Absolute error on the efficiency times        
//                       acceptance                                    
//       bck           - Nominal value of the background estimate      
//       sbck          - Absolute error on the background              
//       n             - Number of observed events (not used for the   
//                       expected limit)                               
//       ntoys         - Number of pseudoexperiments to perform for    
//                       expected limit calculation)                   
//       gauss         - if true, use Gaussian statistics for signal   
//                       instead of Poisson; automatically false       
//                       for n = 0.                                    
//                       Always false for expected limit calculations  
//       nuisanceModel - distribution function used in integration over
//                       nuisance parameters:                          
//                       0 - Gaussian (default), 1 - lognormal,        
//                       2 - gamma;                                    
//                       (automatically 0 when gauss == true)          
//       method        - method of statistical inference:              
//                       \"bayesian\"  - Bayesian with numeric         
//                                       integration (default),        
//                       \"mcmc\"      - another implementation of     
//                                       Bayesian, not optimized,      
//                                       to be used for cross checks   
//                                       only!                         
//                       \"cls\"       - CLs observed limit. We suggest
//                                       using the dedicated interface 
//                                       roostats_limit() instead        
//                       \"fc\"        - Feldman Cousins with numeric  
//                                     integration,                    
//                       \"workspace\" - only create workspace and save
//                                     to file, no interval calculation
//       plotFileName  - file name for the control plot to be created  
//                       file name extension will define the format,   
//                       <plot_cl95.pdf> is the default value,         
//                       specify empty string if you do not want       
//                       the plot to be created (saves time)           
//       seed          - seed for random number generation,            
//                       specify 0 for unique irreproducible seed      
//                                                                     
//                                                                     
//The statistics model in this routine: the routine addresses the task 
//of a Bayesian evaluation of limits for a one-bin counting experiment 
//with systematic uncertainties on luminosity and efficiency for the   
//signal and a global uncertainty on the expected background (implying 
//no correlated error on the luminosity for signal and  background,    
//which will not be suitable for all use cases!). The observable is the
//measured number of events.                                           
//                                                                     
//For more details see                                                 
//        https://twiki.cern.ch/twiki/bin/view/CMS/RooStatsCl95        
//                                                                     
//\033[1m       Note!                                           \033[0m
//If you are running nonstandard ROOT environment, e.g. in CMSSW,      
//you need to make sure that the RooFit and RooStats header files      
//can be found since they might be in a nonstandard location.          
//                                                                     
//For CMSSW_4_2_0_pre8 and later, add the following line to your       
//rootlogon.C:                                                         
//      gSystem -> SetIncludePath( \"-I$ROOFITSYS/include\" );


#include <algorithm>
#include <limits>

#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SequentialProposal.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

// FIXME: remove unnecessary headers
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"


// FIXME: remove namespaces
using namespace RooFit;
using namespace RooStats;
using namespace std;

class LimitResult;


// -----> User interface


LimitResult
GetClsLimit( Double_t ilum, Double_t slum,
	     Double_t eff, Double_t seff,
	     Double_t bck, Double_t sbck,
	     Int_t n );



double
GetBayesianLimit( Double_t ilum, Double_t slum,
		  Double_t eff, Double_t seff,
		  Double_t bck, Double_t sbck,
		  Int_t n,
		  std::string option = "" );



LimitResult 
GetExpectedLimit( Double_t ilum, Double_t slum,
		  Double_t eff, Double_t seff,
		  Double_t bck, Double_t sbck,
		  Int_t nit = 200,
		  std::string method = "bayesian" );


double
roostats_zscore( Double_t ilum, Double_t slum,
		 Double_t eff, Double_t seff,
		 Double_t bck, Double_t sbck,
		 Int_t n,
		 Bool_t gauss = false,
		 Int_t nuisanceModel = 1,
		 std::string method = "bayesian",
		 std::string plotFileName = "plot_zscore.pdf",
		 UInt_t seed = 12345 );



void SetParameter(const char * name, const char * value);
void SetParameter(const char * name, bool value);
void SetParameter(const char * name, int value);
void SetParameter(const char * name, double value);

bool   GetParameter(const char * name, bool);
int    GetParameter(const char * name, int);
double GetParameter(const char * name, double);




// legacy interfaces

Double_t
roostats_cl95( Double_t ilum, Double_t slum,
	       Double_t eff, Double_t seff,
	       Double_t bck, Double_t sbck,
	       Int_t n,
	       Bool_t gauss = kFALSE,
	       Int_t nuisanceModel = 1,
	       std::string method = "bayesian",
	       std::string plotFileName = "plot_cl95.pdf",
	       UInt_t seed = 12345,
	       LimitResult * pLimitResult = 0);

LimitResult 
roostats_clm(  Double_t ilum, Double_t slum,
	       Double_t eff, Double_t seff,
	       Double_t bck, Double_t sbck,
	       Int_t nit = 200, Int_t nuisanceModel = 1,
	       std::string method = "bayesian",
	       UInt_t seed = 12345);

LimitResult 
roostats_limit(Double_t ilum, Double_t slum,
	       Double_t eff, Double_t seff,
	       Double_t bck, Double_t sbck,
	       Int_t n,
	       Bool_t gauss,
	       Int_t nuisanceModel,
	       std::string method,
	       std::string plotFileName,
	       UInt_t seed);

Double_t roostats_cla(Double_t ilum, Double_t slum,
		      Double_t eff, Double_t seff,
		      Double_t bck, Double_t sbck,
		      Int_t nuisanceModel = 0,
		      std::string method = "bayesian",
		      UInt_t seed = 12345);




// ---> implementation below --------------------------------------------

//double confidenceLevel = 0.95;
//bool includeLumiSystIntoBkg = true;
//int verbosity = 3;
//int calculatorType = 0;    // 0-freq, 1-hybrid, 2-asymptotic
//int testStatisticType = 3; // 1-sided PL
//int numberOfClsScanPoints = 50;

class LimitResult{

  friend class CL95Calc;
  
public:
  LimitResult():
    _observed_limit(0),
    _observed_limit_error(0),
    _expected_limit(0),
    _low68(0),
    _high68(0),
    _low95(0),
    _high95(0),
    _cover68(0),
    _cover95(0){};

  // copy constructor
  LimitResult(const LimitResult & other):
    _observed_limit(other._observed_limit),
    _observed_limit_error(other._observed_limit_error),
    _expected_limit(other._expected_limit),
    _low68(other._low68),
    _high68(other._high68),
    _low95(other._low95),
    _high95(other._high95),
    _cover68(other._cover68),
    _cover95(other._cover95){}

  ~LimitResult(){};

  Double_t GetObservedLimit(){return _observed_limit;};
  Double_t GetObservedLimitError(){return _observed_limit_error;};
  Double_t GetExpectedLimit(){return _expected_limit;};

  Double_t GetOneSigmaLowRange(){return _low68;};
  Double_t GetOneSigmaHighRange(){return _high68;};
  Double_t GetOneSigmaCoverage(){return _cover68;};

  Double_t GetTwoSigmaLowRange(){return _low95;};
  Double_t GetTwoSigmaHighRange(){return _high95;};
  Double_t GetTwoSigmaCoverage(){return _cover95;};

  void Clear( void ){
    _observed_limit = 0;
    _observed_limit_error = 0;
    _expected_limit = 0;
    _low68 = 0;
    _high68 = 0;
    _low95 = 0;
    _high95 = 0;
    _cover68 = 0;
    _cover95 = 0;
  }

private:
  Double_t _observed_limit;
  Double_t _observed_limit_error;
  Double_t _expected_limit;
  Double_t _low68;
  Double_t _high68;
  Double_t _low95;
  Double_t _high95;
  Double_t _cover68;
  Double_t _cover95;

  void SetObservedLimit(Double_t limit){_observed_limit = limit;};
  void SetObservedLimitError(Double_t error){_observed_limit_error = error;};
  void SetExpectedLimit(Double_t limit){_expected_limit = limit;};
  void SetOneSigmaLowRange(Double_t band){_low68 = band;};
  void SetOneSigmaHighRange(Double_t band){_high68 = band;};
  void SetTwoSigmaLowRange(Double_t band){_low95 = band;};
  void SetTwoSigmaHighRange(Double_t band){_high95 = band;};

};



class CL95Calc{

public:

  // no public constructors - this is a singleton class
  // use static method CL95Calc::GetInstance() to get the instance
  ~CL95Calc();

  static CL95Calc * GetInstance(void){
    if (!mspInstance) mspInstance = new CL95Calc();
    return mspInstance;
  }

  RooWorkspace * MakeWorkspace(Double_t ilum, Double_t slum,
			       Double_t eff, Double_t seff,
			       Double_t bck, Double_t sbck,
			       Int_t n,
			       Bool_t gauss,
			       Int_t nuisanceModel);

  RooWorkspace *          getWorkspace(){ return pWs;}
  RooStats::ModelConfig * GetModelConfig( std::string mcName = "SbModel" );
  RooAbsData *            GetData( void ){ return data; }

  RooStats::HypoTestResult * 
  GetHypoTest( std::string hypoName, std::string plotName );

  LikelihoodInterval *    GetPlrInterval( double conf_level );

  RooAbsData * makeData(Int_t n);

  Double_t cl95(std::string method = "bayesian", LimitResult * result = 0);

  Double_t cla( Double_t ilum, Double_t slum,
		Double_t eff, Double_t seff,
		Double_t bck, Double_t sbck,
		Int_t nuisanceModel,
		std::string method );
  
  LimitResult clm(Double_t ilum, Double_t slum,
		  Double_t eff, Double_t seff,
		  Double_t bck, Double_t sbck,
		  Int_t nit = 200, Int_t nuisanceModel = 0,
		  std::string method = "bayesian");
  
  int makePlot( std::string method,
		std::string plotFileName = "plot_cl95.pdf" );

  Double_t FC_calc(int Nbins, float conf_int, float ULprecision, bool UseAdaptiveSampling = true, bool CreateConfidenceBelt = true);

  void SetSeed(UInt_t seed);
  void SetParameter(const char * name, const char * value);
  void SetParameter(const char * name, bool value);
  void SetParameter(const char * name, int value);
  void SetParameter(const char * name, double value);

  bool   GetParameter(const char * name, bool);
  int    GetParameter(const char * name, int);
  double GetParameter(const char * name, double);


private:

  CL95Calc();
  CL95Calc(const CL95Calc &); // stop default
  //CL95Calc( UInt_t seed );

  // methods
  void Init( UInt_t seed ); //  to be called by constructor

  int CheckInputs(Double_t ilum, Double_t slum,
		  Double_t eff, Double_t seff,
		  Double_t bck, Double_t sbck,
		  Int_t n,
		  Int_t nuisanceModel);

  void PrintMethodInfo( std::string method );

  int CreateSystTerm( std::string varName,
		      double value,
		      double error,
		      int nuisanceModel,
		      std::string extraVar = "" );

  Double_t GetRandom( std::string pdf, std::string var );
  Long64_t LowBoundarySearch(std::vector<Double_t> * cdf, Double_t value);
  Long64_t HighBoundarySearch(std::vector<Double_t> * cdf, Double_t value);
  MCMCInterval * GetMcmcInterval(double conf_level,
				 int n_iter,
				 int n_burn,
				 double left_side_tail_fraction,
				 int n_bins);
  void makeMcmcPosteriorPlot( std::string filename );
  double printMcmcUpperLimit( std::string filename = "" );
  
  // CLs limit calculator
  HypoTestInverterResult * 
  RunHypoTestInverter(RooWorkspace * w = 0,
		  const char * modelSBName = "ModelConfig",
		  const char * modelBName = "",
		  const char * dataName = "obsData",                 
		  int mCalculatorType = 0,
		  int testStatType = 0, 
		  bool useCLs = true ,  
		  int npoints = 6,   
		  double poimin = 0,  
		  double poimax = 5, 
		  int ntoys=1000,
		  bool useNumberCounting = false,
		  const char * nuisPriorName = 0);

  Double_t RoundUpperBound(Double_t bound);

  // output verbosity
  enum eVerbosity { mERROR = 0, mWARNING = 1, mINFO = 2 };

  // data members
  RooWorkspace * pWs;
  RooStats::ModelConfig SbModel;
  RooStats::ModelConfig BModel;
  RooAbsData * data;
  BayesianCalculator * bcalc;
  RooStats::SimpleInterval * sInt;
  double nsig_rel_err;
  double nbkg_rel_err;
  Int_t _nuisance_model;

  // attributes
  bool mbHaveSigErr;
  bool mbHaveLumiErr;
  bool mbHaveEffErr;
  bool mbHaveBkgErr;
  //bool mbGaussianStatistics;

  // settable parameters
  bool mbOptimize;
  bool mbCorrelatedLumiSyst;
  bool mbMakePlot;
  bool mbGaussianStatistics;
  int mNClsSteps;
  int mNToys;
  int mCalculatorType;
  int mTestStatType;
  int mVerbosity;
  int mRandomSeed;
  double mNToysRatio;
  double mConfidenceLevel;

  // for Bayesian MCMC calculation
  MCMCInterval * mcInt;
  
  // for Feldman-Cousins Calculator
  FeldmanCousins * fcCalc;

  // for profile likelihood ratio
  LikelihoodInterval * pPlrInt;

  // random numbers
  TRandom3 mRandom;

  // expected limits
  Double_t _expected_limit;
  Double_t _low68;
  Double_t _high68;
  Double_t _low95;
  Double_t _high95;
  
  // pointer to class instance
  static CL95Calc * mspInstance;
};


// static pointer to the instance
CL95Calc * CL95Calc::mspInstance = 0;



// default constructor
CL95Calc::CL95Calc(){
  Init(0);
}



void CL95Calc::Init(UInt_t seed){
  pWs = 0;
  data = 0;

  sInt = 0;
  bcalc = 0;
  mcInt = 0;
  fcCalc = 0;
  pPlrInt = 0;
  SbModel.SetName("SbModel");
  SbModel.SetTitle("ModelConfig for roostats_cl95");

  nsig_rel_err = -1.0; // default non-initialized value
  nbkg_rel_err = -1.0; // default non-initialized value

  // set random seed
  SetSeed(seed);

  // default Gaussian nuisance model
  _nuisance_model = 0;

  // set default attributes
  mbHaveSigErr = false;
  mbHaveLumiErr = false;
  mbHaveEffErr = false;
  mbHaveBkgErr = false;
  mbGaussianStatistics = false;

  // set settable parameter defaults
  mbOptimize = false;
  mbCorrelatedLumiSyst = false;
  mbMakePlot = false;
  mbGaussianStatistics = false;
  mNClsSteps = 10;
  mNToys = 1000;
  mNToysRatio = 2.0;
  mCalculatorType = 0;
  mTestStatType = 3;
  mVerbosity = 3;
  mRandomSeed = 12345;
  mConfidenceLevel = 0.95;

  return;
}


CL95Calc::~CL95Calc(){
  delete pWs;
  delete data;
  delete sInt;
  delete bcalc;
  delete mcInt;
  delete fcCalc;
  delete pPlrInt;
}



void CL95Calc::SetSeed( UInt_t seed ){
  //
  // Set random seed. If 0, set unique random.
  //
  std::string _legend = "[CL95Calc::SetSeed]: ";

  mRandomSeed = seed;

  if (seed == 0){
    mRandom.SetSeed();
    UInt_t _seed = mRandom.GetSeed();
    UInt_t _pid = gSystem->GetPid();
    std::cout << _legend << "random seed: " << _seed << std::endl;
    std::cout << _legend << "process ID: " << _pid << std::endl;
    _seed = 31*_seed+_pid;
    std::cout << _legend << "new random seed (31*seed+pid): " << _seed << std::endl;
    mRandom.SetSeed(_seed);
    
    // set RooFit random seed (it has a private copy)
    RooRandom::randomGenerator()->SetSeed(_seed);
  }
  else{
    std::cout << _legend 
	      << "random seed: " << seed 
	      << std::endl;
    mRandom.SetSeed(seed);
    
    // set RooFit random seed (it has a private copy)
    RooRandom::randomGenerator()->SetSeed(seed);
  }

  return;
}



void
CL95Calc::SetParameter(const char * name, bool value){
   //
   // set boolean parameters
   //

   std::string s_name(name);

   if (s_name.find("Optimize") != std::string::npos) mbOptimize = value;
   else if (s_name.find("CorrelatedLumiSyst") != std::string::npos) mbCorrelatedLumiSyst = value;
   else if (s_name.find("MakePlot") != std::string::npos) mbMakePlot = value;
   else if (s_name.find("GaussianStatistics") != std::string::npos) mbGaussianStatistics = value;
   else{
     Info("SetParameter","Unknown parameter %s, ignored", name);
   }
   return;
}



void
CL95Calc::SetParameter(const char * name, int value){
   //
   // set integer parameters
   //

   std::string s_name(name);

   if (s_name.find("NClsSteps") != std::string::npos) mNClsSteps = value;
   else if (s_name.find("NToys") != std::string::npos) mNToys = value;
   else if (s_name.find("CalculatorType") != std::string::npos) mCalculatorType = value;
   else if (s_name.find("TestStatType") != std::string::npos) mTestStatType = value;
   else if (s_name.find("Verbosity") != std::string::npos) mVerbosity = value;
   else if (s_name.find("RandomSeed") != std::string::npos) SetSeed(value);
   else{
     Info("SetParameter","Unknown parameter %s, ignored", name);
   }

   return;
}



void
CL95Calc::SetParameter(const char * name, double value){
   //
   // set double precision parameters
   //

   std::string s_name(name);

   if (s_name.find("NToysRatio") != std::string::npos) mNToysRatio = value;
   else if (s_name.find("ConfidenceLevel") != std::string::npos) mConfidenceLevel = value;
   else{
     Info("SetParameter","Unknown parameter %s, ignored", name);
   }

   return;
}



void
CL95Calc::SetParameter(const char * name, const char * value){
   //
   // set string parameters
   //

   std::string s_name(name);

   // silence warning
   (void)value;

   if (0) {}
   else{
     Info("SetParameter","Unknown parameter %s, ignored", name);
   }

   return;
}



bool
CL95Calc::GetParameter(const char * name, bool){
   //
   // get boolean parameters
   //

   std::string s_name(name);

   if (s_name.find("Optimize") != std::string::npos) return mbOptimize;
   else if (s_name.find("CorrelatedLumiSyst") != std::string::npos) return mbCorrelatedLumiSyst;
   else if (s_name.find("MakePlot") != std::string::npos) return mbMakePlot;
   else if (s_name.find("GaussianStatistics") != std::string::npos) return mbGaussianStatistics;
   else{
     Info("GetParameter","Unknown parameter %s, ignored", name);
     return -1.0;
   }
}



int
CL95Calc::GetParameter(const char * name, int){
   //
   // get integer parameters
   //

   std::string s_name(name);

   if (s_name.find("NClsSteps") != std::string::npos) return mNClsSteps;
   else if (s_name.find("NToys") != std::string::npos) return mNToys;
   else if (s_name.find("CalculatorType") != std::string::npos) return mCalculatorType;
   else if (s_name.find("TestStatType") != std::string::npos) return mTestStatType;
   else if (s_name.find("RandomSeed") != std::string::npos) return mRandomSeed;
   else if (s_name.find("Verbosity") != std::string::npos) return mVerbosity;
   else{
     Info("GetParameter","Unknown parameter %s, ignored", name);
     return -1;
   }
}



double
CL95Calc::GetParameter(const char * name, double){
   //
   // get double precision parameters
   //

   std::string s_name(name);

   if (s_name.find("NToysRatio") != std::string::npos) return mNToysRatio;
   else if (s_name.find("ConfidenceLevel") != std::string::npos) return mConfidenceLevel;
   else{
     Info("GetParameter","Unknown parameter %s, ignored", name);
     return -1.0;
   }
}



int CL95Calc::CheckInputs(Double_t ilum, Double_t slum,
			  Double_t eff, Double_t seff,
			  Double_t bck, Double_t sbck,
			  Int_t n,
			  Int_t nuisanceModel){
  //
  // Check inputs, return ok(0), warning(1) or error(-1)
  //

  std::string _legend = "[CL95Calc::CheckInputs]: ";

  // silence unused var warnings
  (void)bck;
  //(void)n;
  (void)nuisanceModel;

  int _status = 0;

  if (eff <= 0.0){
    _status = -1;
    if (mVerbosity >= mERROR){
      std::cout << _legend 
		<< "ERROR: efficiency is zero!"
		<< std::endl;
    }
  }
  
  if (ilum <= 0.0){
    _status = -1;
    if (mVerbosity >= mERROR){
      std::cout << _legend 
		<< "ERROR: integrated luminosity is zero!"
		<< std::endl;
    }
  }
  
  if (slum > 0.0){
    mbHaveLumiErr = true;
    if (mVerbosity >= mINFO){
      std::cout << _legend 
		<< "nonzero uncertainty for integrated luminosity"
		<< std::endl;
    }
  }
  
  if (seff > 0.0){
    mbHaveEffErr = true;
    if (mVerbosity >= mINFO){
      std::cout << _legend 
		<< "nonzero uncertainty for acceptance*efficiency"
		<< std::endl;
    }
  }
  
  if (sbck > 0.0){
    mbHaveBkgErr = true;
    if (mVerbosity >= mINFO){
      std::cout << _legend 
		<< "nonzero uncertainty for estimated background"
		<< std::endl;
    }
  }
  
  if (n < 0){
    std::cout << "Negative observed number of events specified, exiting" << std::endl;
    return -1.0;
  }

  if (n == 0) mbGaussianStatistics = kFALSE;

  if (mbGaussianStatistics){
    nuisanceModel = 0;
    std::cout << _legend << "Gaussian statistics used" << endl;
  }
  else{
    std::cout << _legend << "Poisson statistics used" << endl;
  }

  return _status;
}



int CL95Calc::CreateSystTerm(std::string varName,
			     double value,
			     double error,
			     int nuisanceModel,
			     std::string extraVar){
  //
  // create a variable, its nuisance parameter,
  // and the corresponding global observable and
  //  the constraint term,
  // according to the requested model: gauss/lognormal/gamma
  //
  // no validation here: all is expected to be correct here
  //
  // If extraVar is specified, its uncertainty gets multiplied
  // in, additionally to error
  //

  // silence warnings about unused variables
  (void)nuisanceModel;

  std::string sFactory;

  if ( error > 0.0 ){
    sFactory = varName+"_nom[1.0]";
    pWs->factory( sFactory.c_str() );
    sFactory = varName+"_nom";
    pWs->var( sFactory.c_str() )->setVal(value);
    pWs->var( sFactory.c_str() )->setConstant(true);
    
    // Lognormal
    sFactory = varName+"_kappa[1.0]";
    pWs->factory( sFactory.c_str() );
    sFactory = varName+"_kappa";
    pWs->var( sFactory.c_str() )->setVal(1.0+error/value);
    pWs->var( sFactory.c_str() )->setConstant(true);
    
    sFactory = "expr::alpha_"+varName+"('pow("+varName
      +"_kappa,beta_"+varName+")',"+varName+"_kappa,beta_"
      +varName+"[0,-5,5])";
    pWs->factory( sFactory.c_str() );
    
    if ( extraVar.size() > 0 ){
      sFactory = "prod::"+varName+"("+varName+"_nom,alpha_"
	+varName+",alpha_"+extraVar+")";
    }
    else{
      sFactory = "prod::"+varName+"("+varName+"_nom,alpha_"
	+varName+")";
    }
    pWs->factory( sFactory.c_str() );
    
    sFactory = "Gaussian::constr_"+varName+"(beta_"+varName
      +",glob_"+varName+"[0,-5,5],1)";
    pWs->factory( sFactory.c_str() );

    // set global observable to const
    sFactory = "glob_" + varName;
    pWs->factory( sFactory.c_str() );
  }
  else if ( extraVar.size() > 0 ){
    sFactory = varName+"_nom[1.0]";
    pWs->factory( sFactory.c_str() );
    sFactory = varName+"_nom";
    pWs->var( sFactory.c_str() )->setVal(value);
    pWs->var( sFactory.c_str() )->setConstant(true);
    
    sFactory = "prod::"+varName+"("+varName+"_nom,alpha_"
      +extraVar+")";
    pWs->factory( sFactory.c_str() );
  }
  else{
    sFactory = varName+"[1.0]";
    pWs->factory( sFactory.c_str() );
    pWs->var( varName.c_str() )->setVal(value);
    pWs->var( varName.c_str() )->setConstant(true);
  }

  return 0;
}



RooWorkspace * 
CL95Calc::MakeWorkspace(Double_t ilum, Double_t slum,
			Double_t eff, Double_t seff,
			Double_t bck, Double_t sbck,
			Int_t n,
			Bool_t gauss,
			Int_t nuisanceModel){
  
  std::string _legend = "[CL95Calc::MakeWorkspace]: ";

  delete pWs;
  pWs = new RooWorkspace("ws");

  mbGaussianStatistics = gauss;

  int _input_status = CheckInputs( ilum, slum,
				   eff, seff,
				   bck, sbck,
				   n,
				   nuisanceModel );

  if ( _input_status < 0 ){
    std::cout << _legend << "critical error in inputs, cannot continue."
	      << std::endl;
    std::exit(-1);
  }

  _nuisance_model = nuisanceModel;
    
  // observable: number of events
  pWs->factory( "n[0]" );

  // integrated luminosity with systematics
  CreateSystTerm("lumi", ilum, slum, nuisanceModel);

  // cross section - parameter of interest
  pWs->factory( "xsec[0,0,1]" );

  // selection efficiency * acceptance with systematics
  CreateSystTerm("efficiency", eff, seff, nuisanceModel);

  // signal yield
  pWs->factory( "prod::nsig(lumi,xsec,efficiency)" );

  // background yield with systematics
  if (mbHaveLumiErr && mbCorrelatedLumiSyst){ // extra uncertainty from lumi
    CreateSystTerm("nbkg", bck, sbck, nuisanceModel, "lumi");
  }
  else{
    CreateSystTerm("nbkg", bck, sbck, nuisanceModel);
  }
  //pWs->Print();

  // core model:
  pWs->factory("sum::yield(nsig,nbkg)");
  if (mbGaussianStatistics){
    // Poisson probability with mean signal+bkg
    std::cout << "[CL95Calc]: creating Gaussian probability as core model..." << std::endl;
    pWs->factory( "Gaussian::model_core(n,yield,expr('sqrt(yield)',yield))" );
  }
  else{
    // Poisson probability with mean signal+bkg
    std::cout << "[CL95Calc]: creating Poisson probability as core model..." << std::endl;
    pWs->factory( "Poisson::model_core(n,yield)" );
  }

  // compose the full model including the constraint terms
  if (mbHaveLumiErr || mbHaveEffErr || mbHaveBkgErr){
    std::string sFactory = "PROD::model(model_core";

    if (mbHaveLumiErr){
      sFactory += ",constr_lumi";
    }

    if (mbHaveEffErr){
      sFactory += ",constr_efficiency";
    }

    if (mbHaveBkgErr){
      sFactory += ",constr_nbkg";
    }

    sFactory += ")";
    std::cout << sFactory << std::endl;
    pWs->factory( sFactory.c_str() );
  }
  else{
    pWs->pdf("model_core")->SetName("model");
  }

  // flat prior for the parameter of interest
  pWs->factory( "Uniform::prior(xsec)" );  

  // floating parameters ranges
  // crude estimates! Need to know data to do better
  pWs->var("n")        ->setRange( 0.0, bck+(5.0*sbck)+10.0*((double)n+1.0)); // ad-hoc range for obs
  //pWs->var("xsec")     ->setRange( 0.0, 15.0*(1.0+nsig_rel_err)/ilum/eff ); // ad-hoc range for POI
  Double_t xsec_upper_bound = 4.0*(std::max(3.0,n-bck)+sqrt(n)+sbck)/ilum/eff;  // ad-hoc range for POI
  xsec_upper_bound = RoundUpperBound(xsec_upper_bound);
  pWs->var("xsec")     ->setRange( 0.0, xsec_upper_bound );

  // observables
  RooArgSet obs(*pWs->var("n"), "obs");

  // global observables
  RooArgSet globalObs("global_obs");
  if (mbHaveLumiErr) pWs->var("glob_lumi")->setConstant(true);
  if (mbHaveEffErr)  pWs->var("glob_efficiency")->setConstant(true);
  if (mbHaveBkgErr)  pWs->var("glob_nbkg")->setConstant(true);
  if (mbHaveLumiErr) globalObs.add( *pWs->var("glob_lumi") );
  if (mbHaveEffErr)  globalObs.add( *pWs->var("glob_efficiency") );
  if (mbHaveBkgErr)  globalObs.add( *pWs->var("glob_nbkg") );

  // parameters of interest
  RooArgSet poi(*pWs->var("xsec"), "poi");

  // nuisance parameters
  RooArgSet nuis("nuis");
  if (mbHaveLumiErr) nuis.add( *pWs->var("beta_lumi") );
  if (mbHaveEffErr)  nuis.add( *pWs->var("beta_efficiency") );
  if (mbHaveBkgErr)  nuis.add( *pWs->var("beta_nbkg") );

  // setup the S+B model
  SbModel.SetWorkspace(*pWs);
  SbModel.SetPdf(*(pWs->pdf("model")));
  SbModel.SetParametersOfInterest(poi);
  SbModel.SetPriorPdf(*(pWs->pdf("prior")));
  SbModel.SetNuisanceParameters(nuis);
  SbModel.SetObservables(obs);
  SbModel.SetGlobalObservables(globalObs);

  // will import the model config once the snapshot is saved

  // background-only model
  // use the same PDF as s+b, with xsec=0
  // (poi zero value will be set in the snapshot)
  //BModel = *(RooStats::ModelConfig *)pWs->obj("SbModel");
  BModel = SbModel;
  BModel.SetName("BModel");
  BModel.SetWorkspace(*pWs);

  // create data
  pWs->var("n")         ->setVal(n);
  delete data;
  data = new RooDataSet("data","",*(SbModel.GetObservables()));
  data->add( *(SbModel.GetObservables()));
  data->SetName("observed_data");
  pWs->import(*data);

  // make RooFit quiet
  // cash the current message level first
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // Now set up parameter snapshots for the S+B and B models

  // find global maximum with the signal+background model
  // with conditional MLEs for nuisance parameters
  // and save the parameter point snapshot in the Workspace
  //  - safer to keep a default name because some RooStats calculators
  //    will anticipate it
  RooAbsReal * pNll = SbModel.GetPdf()->createNLL(*data);
  //RooAbsReal * pProfile = pNll->createProfile(RooArgSet());
  RooAbsReal * pProfile = pNll->createProfile( globalObs );
  pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
  RooArgSet * pPoiAndNuisance = new RooArgSet("poiAndNuisance");
  if(SbModel.GetNuisanceParameters())
    pPoiAndNuisance->add(*SbModel.GetNuisanceParameters());
  pPoiAndNuisance->add(*SbModel.GetParametersOfInterest());
  std::cout << "\nWill save these parameter points that correspond to the fit to data" << std::endl;
  pPoiAndNuisance->Print("v");
  SbModel.SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;

  // Find a parameter point for generating pseudo-data
  // with the background-only data.
  // Save the parameter point snapshot in the Workspace
  //
  // POI value under the background hypothesis
  Double_t poiValueForBModel = 0.0;
  pNll = BModel.GetPdf()->createNLL(*data);
  //const RooArgSet * poi = BModel.GetParametersOfInterest();
  //pProfile = pNll->createProfile(*poi);
  RooArgSet poiAndGlobalObs("poiAndGlobalObs");
  poiAndGlobalObs.add( poi );
  poiAndGlobalObs.add( globalObs );
  pProfile = pNll->createProfile( poiAndGlobalObs );
  ((RooRealVar *)poi.first())->setVal(poiValueForBModel);
  pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
  pPoiAndNuisance = new RooArgSet("poiAndNuisance");
  if(BModel.GetNuisanceParameters())
    pPoiAndNuisance->add(*BModel.GetNuisanceParameters());
  pPoiAndNuisance->add(*BModel.GetParametersOfInterest());
  std::cout << "\nShould use these parameter points to generate pseudo data for bkg only" << std::endl;
  pPoiAndNuisance->Print("v");
  BModel.SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;

  // import the model configs, has to be after all snapshots are saved
  pWs->import(SbModel);
  pWs->import(BModel);

  // restore RooFit messaging level
  RooMsgService::instance().setGlobalKillBelow(msglevel);

  // We also need to set up parameter snapshots for the models
  // but we need data for that, so it is done in makeData()

  //pWs->Print();



  return pWs;
}


RooAbsData * CL95Calc::makeData( Int_t n ){
  //
  // make the dataset owned by the class
  // the current one is deleted
  //
  
  // create data
  pWs->var("n")         ->setVal(n);
  delete data;
  data = new RooDataSet("data","",*(SbModel.GetObservables()));
  data->add( *(SbModel.GetObservables()));

  return data;
}



RooStats::ModelConfig * CL95Calc::GetModelConfig( std::string mcName ){
  //
  // Return a pointer to the ModelConfig or 0 if not found
  // User does NOT take ownership.
  //

  if (pWs){
    RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)pWs->obj(mcName.c_str());
    _mc -> SetWorkspace(*pWs);
    //_mc->Print();
    //_mc->GetWorkspace()->Print();
    return _mc;
  }
  else return 0;
}



LikelihoodInterval * CL95Calc::GetPlrInterval( double conf_level ){
  //
  // Profile likelihood ratio interval calculation
  //

  delete pPlrInt;
  
  RooStats::ModelConfig * _mc = GetModelConfig();
  //_mc->Print();

  ProfileLikelihoodCalculator plc(*data, *_mc);
  plc.SetConfidenceLevel(conf_level);
  pPlrInt = plc.GetInterval();

  return pPlrInt;
}



MCMCInterval * CL95Calc::GetMcmcInterval(double conf_level,
					int n_iter,
					int n_burn,
					double left_side_tail_fraction,
					int n_bins){
  // use MCMCCalculator  (takes about 1 min)
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit
  
  RooFitResult * fit = pWs->pdf("model")->fitTo(*data,Save(),
					       Verbose(kFALSE),
					       PrintLevel(-1),
					       Warnings(0),
					       PrintEvalErrors(-1));
  // silence warnings about unused variables
  (void)fit;

  /*
  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  ProposalFunction* pf = ph.GetProposalFunction();
  */

  // this proposal function seems fairly robust
  SequentialProposal sp(10.0);
  
  MCMCCalculator mcmc( *data, SbModel );
  mcmc.SetConfidenceLevel(conf_level);
  mcmc.SetNumIters(n_iter);          // Metropolis-Hastings algorithm iterations
  //mcmc.SetProposalFunction(*pf);
  mcmc.SetProposalFunction(sp);
  mcmc.SetNumBurnInSteps(n_burn); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(left_side_tail_fraction);
  mcmc.SetNumBins(n_bins);
  
  delete mcInt;
  mcInt = mcmc.GetInterval();

  return mcInt;
}


void CL95Calc::makeMcmcPosteriorPlot( std::string filename ){
  
  TCanvas c1("c1");
  MCMCIntervalPlot plot(*mcInt);
  plot.Draw();
  c1.SaveAs(filename.c_str());
  
  // Markov chain scatter plots
  if (pWs->var("nsig_nuis")){
    TCanvas c2("c2");
    plot.DrawChainScatter(*pWs->var("xsec"),*pWs->var("nsig_nuis"));
    c2.SaveAs("scatter_mcmc_poi_vs_nsig_nuis.png");
  }
  if (pWs->var("nbkg")){
    TCanvas c3("c3");
    plot.DrawChainScatter(*pWs->var("xsec"),*pWs->var("nbkg"));
    c3.SaveAs("scatter_mcmc_poi_vs_nbkg.png");
  }
  
  return;
}


double CL95Calc::printMcmcUpperLimit( std::string filename ){
  //
  // print out the upper limit on the first Parameter of Interest
  //

  RooRealVar * firstPOI = (RooRealVar*) SbModel.GetParametersOfInterest()->first();
  double _limit = mcInt->UpperLimit(*firstPOI);
  cout << "\n95% upper limit on " <<firstPOI->GetName()<<" is : "<<
    _limit <<endl;

  if (filename.size()!=0){
    
    std::ofstream aFile;

    // append to file if exists
    aFile.open(filename.c_str(), std::ios_base::app);

    char buf[1024];
    sprintf(buf, "%7.6f", _limit);

    aFile << buf << std::endl;

    // close outfile here so it is safe even if subsequent iterations crash
    aFile.close();

  }

  return _limit;
}



Double_t CL95Calc::FC_calc(int Nbins, float conf_int, float ULprecision, bool UseAdaptiveSampling, bool CreateConfidenceBelt){


  Double_t upper_limit = 0;
  int cnt = 0;
  bool verbose = true; //Set to true to see the output of each FC step

  std::cout << "[roostats_cl95]: FC calculation is still experimental in this context!!!" << std::endl;
      
  std::cout << "[roostats_cl95]: Range of allowed cross section values: [" 
	    << pWs->var("xsec")->getMin() << ", " 
	    << pWs->var("xsec")->getMax() << "]" << std::endl;


  //prepare Feldman-Cousins Calulator

  delete fcCalc;
  fcCalc = new FeldmanCousins(*data,SbModel);
      
  fcCalc->SetConfidenceLevel(conf_int); // confidence interval
  //fcCalc->AdditionalNToysFactor(0.1); // to speed up the result 
  fcCalc->UseAdaptiveSampling(UseAdaptiveSampling); // speed it up a bit
  fcCalc->SetNBins(Nbins); // set how many points per parameter of interest to scan
  fcCalc->CreateConfBelt(CreateConfidenceBelt); // save the information in the belt for plotting

      
  if(!SbModel.GetPdf()->canBeExtended()){
    if(data->numEntries()==1)     
      fcCalc->FluctuateNumDataEntries(false);
    else
      cout <<"Not sure what to do about this model" <<endl;
  }

  RooRealVar* firstPOI = (RooRealVar*) SbModel.GetParametersOfInterest()->first();
  
  double max = firstPOI->getMax();
  double min = firstPOI->getMin();
  double med = (max + min)/2.0;
      
  double maxPerm = firstPOI->getMax();
  double minPerm = firstPOI->getMin();
    
  double UpperLimit = 0;
  
  PointSetInterval* interval = 0;

  while ( 1 ){
    
    ++cnt;
    firstPOI->setMax( max );
    firstPOI->setMin( min );
    
    if ( verbose ) std::cout << "[FeldmanCousins]: Setting max/min/med to = " << max << " / " << min << " / " << med <<  std::endl;
	
    interval = fcCalc->GetInterval();
    interval -> Delete();
	
    UpperLimit = interval -> UpperLimit(*firstPOI);
    if ( verbose ) std::cout <<"[FeldmanCousins]: Updating Upper Limt to = "<< UpperLimit << std::endl;

    if ( UpperLimit > 0.000001 ){

      min = med;
      med = (max + min)/2.0;
      
    }
    else{
      
      max = med;
      med = (max + min)/2.0;
      
    }
    
    if (  ( UpperLimit > 0.000001 ) && ( (max - min) < ULprecision)  ) {
      upper_limit = UpperLimit;
      std::cout <<"[FeldmanCousins]: In "<< cnt << " steps Upper Limt converged to " << upper_limit << std::endl;
      break; 
    }
    
    if ( cnt > 50 ) {
      upper_limit = -1;
      std::cout << std::endl;
      std::cout <<"[FeldmanCousins     WARNING!!!!!!!!!!!!       ]: Calculator could not converge in under 50 steps. Returning Upper Limit of -1." << std::endl;
      std::cout << std::endl;
      break;
    }

  }
      
  pWs->var("xsec")->setMax( maxPerm );
  pWs->var("xsec")->setMin( minPerm );

  return upper_limit;

}





Double_t CL95Calc::cl95( std::string method, LimitResult * result ){
  //
  // Compute the observed limit
  // For some methods - CLs - compute the expected limts too.
  // Extended results are returned via reference as LimitResul object
  //
  // this method assumes that the workspace,
  // data and model config are ready
  //

  std::string legend = "[CL95Calc::cl95]: ";

  Double_t upper_limit = -1.0;

  // make RooFit quiet
  // cash the current message level first
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  // get ugly RooFit print out of the way
  // FIXME: uncomment
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  Int_t _attempt = 0; // allow several attempts for limit calculation, stop after that
  while(1){

    ++_attempt;
    
    // too many attempts
    if (_attempt > 5){
      std::cout << "[roostats_cl95]: limit calculation did not converge, exiting..." << std::endl;
      return -1.0;
    }

    if (method.find("bayesian") != std::string::npos){
      
      std::cout << "[roostats_cl95]: Range of allowed cross section values: [" 
		<< pWs->var("xsec")->getMin() << ", " 
		<< pWs->var("xsec")->getMax() << "]" << std::endl;

      //prepare Bayesian Calulator
      //delete bcalc;
      bcalc = new BayesianCalculator(*data, SbModel);
      TString namestring = "mybc";
      bcalc->SetName(namestring);
      bcalc->SetConfidenceLevel(mConfidenceLevel);
      bcalc->SetLeftSideTailFraction(0.0);
      //bcalc->SetIntegrationType("ROOFIT");
      
      delete sInt;
      sInt = bcalc->GetInterval();
      upper_limit = sInt->UpperLimit();
      delete sInt;
      sInt = 0;

      if (result) result->_observed_limit = upper_limit;
    }
    else if (method.find("mcmc") != std::string::npos){
      
      std::cout << "[roostats_cl95]: Bayesian MCMC calculation is still experimental in this context!!!" << std::endl;
      
      std::cout << "[roostats_cl95]: Range of allowed cross section values: [" 
		<< pWs->var("xsec")->getMin() << ", " 
		<< pWs->var("xsec")->getMax() << "]" << std::endl;

      //prepare Bayesian Markov Chain MC Calulator
      mcInt = GetMcmcInterval(mConfidenceLevel, 1000000, 500, 0.0, 40);
      upper_limit = printMcmcUpperLimit();

      if (result) result->_observed_limit = upper_limit;
    }
    else if (method.find("plr") != std::string::npos){
      
      std::cout << "[roostats_cl95]: Range of allowed cross section values: [" 
		<< pWs->var("xsec")->getMin() << ", " 
		<< pWs->var("xsec")->getMax() << "]" << std::endl;

      //prepare profile likelihood ratio Calulator
      pPlrInt = GetPlrInterval(mConfidenceLevel);
      upper_limit = pPlrInt->UpperLimit(*pWs->var("xsec"));
    }
    else if (method.find("cls") != std::string::npos){
      //
      // testing CLs
      //
      
      std::cout << "[roostats_cl95]: CLs calculation is still experimental in this context!!!" << std::endl;
      
      std::cout << "[roostats_cl95]: Range of allowed cross section values: [" 
		<< pWs->var("xsec")->getMin() << ", " 
		<< pWs->var("xsec")->getMax() << "]" << std::endl;
      // timer
      TStopwatch t;
      t.Start();

      // load parameter point with the best fit to data
      SbModel.LoadSnapshot();

      // estimate upper range boundary using quick PLR limit
      GetPlrInterval(mConfidenceLevel);
      upper_limit = pPlrInt->UpperLimit( *pWs->var("xsec") );
      //Double_t upper_range = ((double)(int)(4.0 * upper_limit*100.0))/100.0; // round to ~1% precision
      Double_t upper_range = 2.0 * upper_limit;

      // debug output
      std::cout << legend
		<< "CLs scan range: [0, " << upper_range << "]" 
		<< std::endl;

      RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);


      //calc.SetParameter("PlotHypoTestResult", plotHypoTestResult);
      
      HypoTestInverterResult * _res = 
	RunHypoTestInverter( pWs,
			 "SbModel",
			 "BModel",
			 "observed_data",
			 mCalculatorType, // calculator type, 0-freq, 1-hybrid, 2-asymptotic CLs
			 mTestStatType, // test statistic, 0-lep, 1-tevatron, 2-PL, 3-PL 1-sided
			 true, // useCls
			 mNClsSteps, // npoints in the scan
			 0, // poimin: use default is poimin >= poimax
			 upper_range,
			 mNToys,// ntoys
			 true );


      t.Stop();
      t.Print();

      if (result){
	result->_observed_limit = _res->UpperLimit();
	result->_observed_limit_error = _res->UpperLimitEstimatedError();
	result->_expected_limit = _res->GetExpectedUpperLimit(0);
	result->_low68  = _res->GetExpectedUpperLimit(-1);
	result->_high68 = _res->GetExpectedUpperLimit(1);
	result->_low95  = _res->GetExpectedUpperLimit(-2);
	result->_high95 = _res->GetExpectedUpperLimit(2);
	result->_cover68 = -1.0;
	result->_cover95 = -1.0;
      }

      upper_limit = _res->UpperLimit();

    } // end of the CLs block
    else if (method.find("fc") != std::string::npos){

      int Nbins = 1;
      float conf_int = mConfidenceLevel;
      float ULprecision = 0.1;
      bool UseAdaptiveSampling = true;
      bool CreateConfidenceBelt = true;
      
      
      upper_limit = FC_calc(Nbins, conf_int, ULprecision, UseAdaptiveSampling, CreateConfidenceBelt);
	
    } // end of the FC block
    else{

      std::cout << "[roostats_cl95]: method " << method 
		<< " is not implemented, exiting" <<std::endl;
      return -1.0;

    } // end of choice of method block

    
    // adaptive range in case the POI range was not guessed properly
    Double_t _poi_max_range = pWs->var("xsec")->getMax();

    if (method.find("cls")!=std::string::npos) break;
    if (method.find("fc") != std::string::npos ) break;
    // range too wide
    else if (upper_limit < _poi_max_range/10.0){
      std::cout << "[roostats_cl95]: POI range is too wide, will narrow the range and rerun" << std::endl;
      pWs->var("xsec")->setMax(RoundUpperBound(_poi_max_range/2.0));
    }
    // range too narrow
    else if (upper_limit > _poi_max_range/2.0){
      std::cout << "[roostats_cl95]: upper limit is too narrow, will widen the range and rerun" << std::endl;
      pWs->var("xsec")->setMax(RoundUpperBound(2.0*_poi_max_range));
    }
    // all good, limit is ready
    else{
      break;
    }
    
  } // end of while(1) loop
  
  // restore RooFit messaging level
  RooMsgService::instance().setGlobalKillBelow(msglevel);

  return upper_limit;
  
}


Double_t CL95Calc::cla( Double_t ilum, Double_t slum,
			Double_t eff, Double_t seff,
			Double_t bck, Double_t sbck,
			Int_t nuisanceModel,
			std::string method ){

  MakeWorkspace( ilum, slum,
		 eff, seff,
		 bck, sbck,
		 (int)(bck+0.5), // supply expected rate as data
		 kFALSE,
		 nuisanceModel );
  
  Double_t CL95A = 0, precision = 1.e-4;

  Int_t i;
  for (i = int(bck+0.5); i >= 0; i--)
    {
      makeData( i );

      Double_t s95 = cl95( method );
      Double_t s95w =s95*TMath::Poisson( (Double_t)i, bck );
      CL95A += s95w;
      cout << "[roostats_cla]: n = " << i << "; 95% C.L. = " << s95 << " pb; weighted 95% C.L. = " << s95w << " pb; running <s95> = " << CL95A << " pb" << endl;

      if (s95w < CL95A*precision) break;
    }
  cout << "[roostats_cla]: Lower bound on n has been found at " << i+1 << endl;

  for (i = int(bck+0.5)+1; ; i++)
    {
      makeData( i );
      Double_t s95 = cl95( method );
      Double_t s95w =s95*TMath::Poisson( (Double_t)i, bck );
      CL95A += s95w;
      cout << "[roostats_cla]: n = " << i << "; 95% C.L. = " << s95 << " pb; weighted 95% C.L. = " << s95w << " pb; running <s95> = " << CL95A << " pb" << endl;

      if (s95w < CL95A*precision) break;
    }
  cout << "[roostats_cla]: Upper bound on n has been found at " << i << endl;
  cout << "[roostats_cla]: Average upper 95% C.L. limit = " << CL95A << " pb" << endl;

  return CL95A;
}



LimitResult CL95Calc::clm( Double_t ilum, Double_t slum,
			   Double_t eff, Double_t seff,
			   Double_t bck, Double_t sbck,
			   Int_t nit, Int_t nuisanceModel,
			   std::string method ){

  std::string _legend = "[CL95Calc::clm()]: ";
  
  MakeWorkspace( ilum, slum,
		 eff, seff,
		 bck, sbck,
		 (int)(bck+0.5),
		 kFALSE,
		 nuisanceModel );
  
  Double_t CLM = 0.0;
  LimitResult _result;

  Double_t b68[2] = {0.0, 0.0}; // 1-sigma expected band
  Double_t b95[2] = {0.0, 0.0}; // 2-sigma expected band

  std::vector<Double_t> pe;

  // timer
  TStopwatch t;
  t.Start(); // start timer
  Double_t _realtime = 0.0;
  Double_t _cputime = 0.0;
  Double_t _realtime_last = 0.0;
  Double_t _cputime_last = 0.0;
  Double_t _realtime_average = 0.0;
  Double_t _cputime_average = 0.0;

  // throw pseudoexperiments
  if (nit <= 0)return _result;
  std::map<Int_t,Double_t> cached_limit;
  for (Int_t i = 0; i < nit; i++)
    {
      std::cout << std::endl << _legend
		<< "Pseudoexperiment # " << i+1 
		<< " / " << nit << std::endl;
      // throw random nuisance parameter (bkg yield)
      //Double_t bmean = GetRandom("syst_nbkg", "nbkg");
      Double_t _beta_lumi = 0.0;
      Double_t _lumi_kappa = 1.0;
      if (pWs->pdf("constr_lumi")){
	_beta_lumi = GetRandom("constr_lumi", "beta_lumi");
	_lumi_kappa = pWs->var("lumi_kappa")->getVal();
      }
      Double_t _beta_nbkg = 0.0;
      Double_t _nbkg_kappa = 1.0;
      Double_t _nbkg_nom;
      if (pWs->pdf("constr_nbkg")){
	_beta_nbkg = GetRandom("constr_nbkg", "beta_nbkg");
	_nbkg_kappa = pWs->var("nbkg_kappa")->getVal();
	_nbkg_nom = pWs->var("nbkg_nom")->getVal();
      }
      else{
	_nbkg_nom = pWs->var("nbkg")->getVal();
      }
      Double_t _alpha_nbkg = pow(_nbkg_kappa,_beta_nbkg);
      Double_t _alpha_lumi = pow(_lumi_kappa,_beta_lumi);
      Double_t bmean;
      if (mbCorrelatedLumiSyst){
	bmean = _nbkg_nom*_alpha_lumi*_alpha_nbkg;
      }
      else{
	bmean = _nbkg_nom*_alpha_nbkg;
      }
      

      std::cout << "[roostats_clm]: generatin pseudo-data with bmean = " << bmean << std::endl;
      Int_t n = mRandom.Poisson(bmean);

      // check if the limit for this n is already cached
      Double_t _pe = -1.0;
      if (cached_limit.find(n)==cached_limit.end()){
	
	makeData( n );
	std::cout << "[roostats_clm]: invoking CL95 with n = " << n << std::endl;
	
	_pe = cl95( method );
	cached_limit[n] = _pe;
      }
      else{
	std::cout << "[roostats_clm]: returning previously cached limit for n = " << n << std::endl;
	_pe = cached_limit[n];
      }

      pe.push_back(_pe);
      CLM += pe[i];

      _realtime_last = t.RealTime() - _realtime;
      _cputime_last  = t.CpuTime() - _cputime;
      _realtime = t.RealTime();
      _cputime = t.CpuTime();
      t.Continue();
      _realtime_average = _realtime/((Double_t)(i+1));
      _cputime_average  = _cputime/((Double_t)(i+1));

      std::cout << "n = " << n << "; 95% C.L. = " << _pe << " pb; running <s95> = " << CLM/(i+1.) << std::endl;
      std::cout << "Real time (s), this iteration: " << _realtime_last << ", average per iteration: " << _realtime_average << ", total: " << _realtime << std::endl;
      std::cout << "CPU time (s),  this iteration: " << _cputime_last << ", average per iteration: " << _cputime_average << ", total: " << _cputime << std::endl << std::endl;
    }

  CLM /= nit;

  // sort the vector with limits
  std::sort(pe.begin(), pe.end());

  // median for the expected limit
  Double_t _median = TMath::Median(nit, &pe[0]);

  // quantiles for the expected limit bands
  Double_t _prob[4]; // array with quantile boundaries
  _prob[0] = 0.021;
  _prob[1] = 0.159;
  _prob[2] = 0.841;
  _prob[3] = 0.979;

  Double_t _quantiles[4]; // array for the results

  TMath::Quantiles(nit, 4, &pe[0], _quantiles, _prob); // evaluate quantiles

  b68[0] = _quantiles[1];
  b68[1] = _quantiles[2];
  b95[0] = _quantiles[0];
  b95[1] = _quantiles[3]; 

  // let's get actual coverages now

  Long64_t lc68 = LowBoundarySearch(&pe, _quantiles[1]);
  Long64_t uc68 = HighBoundarySearch(&pe, _quantiles[2]);
  Long64_t lc95 = LowBoundarySearch(&pe, _quantiles[0]);
  Long64_t uc95 = HighBoundarySearch(&pe, _quantiles[3]);

  Double_t _cover68 = (nit - lc68 - uc68)*100./nit;
  Double_t _cover95 = (nit - lc95 - uc95)*100./nit;

  std::cout << "[CL95Calc::clm()]: median limit: " << _median << std::endl;
  std::cout << "[CL95Calc::clm()]: 1 sigma band: [" << b68[0] << "," << b68[1] << 
    "]; actual coverage: " << _cover68 << 
    "%; lower/upper percentile: " << lc68*100./nit <<"/" << uc68*100./nit << std::endl;
  std::cout << "[CL95Calc::clm()]: 2 sigma band: [" << b95[0] << "," << b95[1] << 
    "]; actual coverage: " << _cover95 << 
    "%; lower/upper percentile: " << lc95*100./nit <<"/" << uc95*100./nit << std::endl;

  t.Print();

  _result._expected_limit = _median;
  _result._low68  = b68[0];
  _result._high68 = b68[1];
  _result._low95  = b95[0];
  _result._high95 = b95[1];
  _result._cover68 = _cover68;
  _result._cover95 = _cover95;

  return _result;
}



int CL95Calc::makePlot( std::string method,
			std::string plotFileName ){

  if (method.find("bayesian") != std::string::npos){

    std::cout << "[roostats_cl95]: making Bayesian posterior plot" << endl;
  
    TCanvas c1("posterior");
    bcalc->SetScanOfPosterior(100);
    RooPlot * plot = bcalc->GetPosteriorPlot();
    plot->Draw();
    c1.SaveAs(plotFileName.c_str());
  }
  else if (method.find("mcmc") != std::string::npos){

    std::cout << "[roostats_cl95]: making Bayesian MCMC posterior plot" << endl;

    makeMcmcPosteriorPlot(plotFileName);
  
  }
  else if (method.find("cls") != std::string::npos){

    // cls plot is made in AnalyzeResult(), do nothing here
  
  }
  else{
    std::cout << "[roostats_cl95]: plot for method " << method 
	      << " is not implemented" <<std::endl;
    return -1;
  }

  return 0;
}



Double_t CL95Calc::GetRandom( std::string pdf, std::string var ){
  //
  // generates a random number using a pdf in the workspace
  //
  
  // generate a dataset with one entry
  RooDataSet * _ds = pWs->pdf(pdf.c_str())->generate(*pWs->var(var.c_str()), 1);

  Double_t _result = ((RooRealVar *)(_ds->get(0)->first()))->getVal();
  delete _ds;

  return _result;
}


Long64_t CL95Calc::LowBoundarySearch(std::vector<Double_t> * cdf, Double_t value){
  //
  // return number of elements which are < value with precision 1e-10
  //

  Long64_t result = 0;
  std::vector<Double_t>::const_iterator i = cdf->begin();
  while( (*i<value) && fabs(*i-value)>1.0e-10 && (i!=cdf->end()) ){
    ++i;
    ++result;
  }
  return result;
}


Long64_t CL95Calc::HighBoundarySearch(std::vector<Double_t> * cdf, Double_t value){
  //
  // return number of elements which are > value with precision 1e-10
  //

  Long64_t result = 0;
  std::vector<Double_t>::const_iterator i = cdf->end();
  while(1){ // (*i<value) && (i!=cdf->begin()) ){
    --i;
    if (*i>value && fabs(*i-value)>1.0e-10 ){
      ++result;
    }
    else break;
    if (i==cdf->begin()) break;
  }
  return result;
}



Double_t CL95Calc::RoundUpperBound(Double_t bound){
  //
  // find a round upper bound for a floating point
  //
  Double_t power = log10(bound);
  Int_t int_power = power>0.0 ? (Int_t)power : (Int_t)(power-1.0);
  Int_t int_bound = (Int_t)(bound/pow(10,(Double_t)int_power) * 10.0 + 1.0);
  bound = (Double_t)(int_bound/10.0*pow(10,(Double_t)int_power));
  return bound;
}



void CL95Calc::PrintMethodInfo( std::string method ){
  //
  // Printout some info
  // 

  std::string legend = "[CL95Calc::ValidateInput]: ";

  std::cout << legend << "estimating 95% C.L. upper limit" << endl;
  if (method.find("bayesian") != std::string::npos){
    std::cout << legend << "using Bayesian calculation via numeric integration" << endl;
  }
  else if (method.find("plr") != std::string::npos){
    std::cout << legend << "using profile likelihood ratio calculation with Wilk's theorem" << endl;
  }
  else if (method.find("mcmc") != std::string::npos){
    std::cout << legend << "using Bayesian calculation via numeric integration" << endl;
  }
  else if (method.find("cls") != std::string::npos){
    std::cout << legend << "using CLs calculation" << endl;
  }
  else if (method.find("fc") != std::string::npos){
    std::cout << legend << "using Feldman-Cousins approach" << endl;
  }
  else if (method.find("workspace") != std::string::npos){
    std::cout << legend << "no interval calculation, only create and save workspace" << endl;
  }
  else{
    std::cout << legend << "method " << method 
	      << " is not implemented, exiting" <<std::endl;
    return;
  }

  return;
}



RooStats::HypoTestResult *
CL95Calc::GetHypoTest( std::string hypoName, std::string plotName ){
  //
  // Evaluate and return HypoTestResult for a hypothesis
  // defined as a model config object in the workspace
  // 
  // If non-empty plotName is specified, create a plot
  // of the sampling distribution
  //
  // User takes ownership of returned object
  //

  std::string legend = "[CL95Calc::GetHypoTest]: ";

  RooStats::HypoTestResult * pHtr = 0;

  // not used yet
  // name of prior PDF for hybrid calculator
  const char * nuisPriorName = 0;

  // silence warning
  (void)plotName;

  if (!pWs){
    std::cout << legend
	      << "workspace is not created, cannot do hypothesis test"
	      << std::endl;
  }

  // get data
  RooAbsData * pData = pWs->data("observed_data");

  // get s+b model config
  RooStats::ModelConfig * pSbModel = (RooStats::ModelConfig *)pWs->obj("SbModel");

  // get b-only model config
  RooStats::ModelConfig * pBModel = (RooStats::ModelConfig *)pWs->obj(hypoName.c_str());

  // parameter of interest
  const RooArgSet * pPoiSet = pSbModel->GetParametersOfInterest();
  RooRealVar *pPoi = (RooRealVar*)pPoiSet->first();

  // load parameter snapshot
  if ( !pBModel->GetSnapshot() ){
    Error(legend.c_str(),"B-only hypothesis has no snapshot saved, exiting" );
    std::exit(-1);
  }

  // test statistic
  SimpleLikelihoodRatioTestStat slrts(*pSbModel->GetPdf(),*pBModel->GetPdf());
  if (pSbModel->GetSnapshot()) slrts.SetNullParameters(*pSbModel->GetSnapshot());
  if (pBModel->GetSnapshot()) slrts.SetAltParameters(*pBModel->GetSnapshot());
  slrts.SetReuseNLL(mbOptimize);
  
  // ratio of profile likelihood - need to pass snapshot for the alt
  RatioOfProfiledLikelihoodsTestStat 
    ropl(*pSbModel->GetPdf(), *pBModel->GetPdf(), pBModel->GetSnapshot());
  ropl.SetSubtractMLE(false);
  //ropl.SetPrintLevel(mPrintLevel);
  //ropl.SetMinimizer(mMinimizerType.c_str());
  ropl.SetReuseNLL(mbOptimize);
  if (mbOptimize) ropl.SetStrategy(0);
  
  ProfileLikelihoodTestStat profll(*pSbModel->GetPdf());
  if (mTestStatType == 3) profll.SetOneSided(1);
  //profll.SetMinimizer(mMinimizerType.c_str());
  //profll.SetPrintLevel(mPrintLevel);
  profll.SetReuseNLL(mbOptimize);
  if (mbOptimize) profll.SetStrategy(0);

  MaxLikelihoodEstimateTestStat maxll(*pSbModel->GetPdf(),*pPoi); 

  /*
  ProfileLikelihoodTestStat profll(*pSbModel->GetPdf());
  profll.SetOneSided(0);
  profll.SetReuseNLL(1);
  profll.SetStrategy(0);
  //if (testStatType == 3) profll.SetOneSided(1);
  //profll.SetMinimizer(mMinimizerType.c_str());
  //profll.SetPrintLevel(mPrintLevel);
  */
  
  // create the HypoTest calculator object
  HypoTestCalculatorGeneric *  hc = 0;
  if (mCalculatorType == 0){
    hc = new FrequentistCalculator(*pData, *pSbModel, *pBModel);
    ((FrequentistCalculator*) hc)->SetToys(mNToys,0); 
  }
  else if (mCalculatorType == 1){
    hc = new HybridCalculator(*pData, *pSbModel, *pBModel);
    HybridCalculator *hhc = dynamic_cast<HybridCalculator*> (hc);
    assert(hhc);
    
    hhc->SetToys(mNToys,0); // can use less ntoys for b hypothesis 
    
    // remove global observables from ModelConfig (this is probably not needed anymore in 5.32)
    pBModel->SetGlobalObservables(RooArgSet() );
    pSbModel->SetGlobalObservables(RooArgSet() );
    
    // check for nuisance prior pdf in case of nuisance parameters 
    if (pBModel->GetNuisanceParameters() || pSbModel->GetNuisanceParameters() ) {
      RooAbsPdf * nuisPdf = 0; 
      if (nuisPriorName) nuisPdf = pWs->pdf(nuisPriorName);
      // use prior defined first in pBModel (then in SbModel)
      if (!nuisPdf)  { 
	Info("RunHypoTestInverter","No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the model");
	if (pBModel->GetPdf() && pBModel->GetObservables() ) 
	  nuisPdf = RooStats::MakeNuisancePdf(*pBModel,"nuisancePdf_bmodel");
	else 
	  nuisPdf = RooStats::MakeNuisancePdf(*pSbModel,"nuisancePdf_sbmodel");
      }   
      if (!nuisPdf ) {
	if (pBModel->GetPriorPdf())  { 
	  nuisPdf = pBModel->GetPriorPdf();
	  Info("RunHypoTestInverter","No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf->GetName());            
	}
	else { 
	  Error("RunHypoTestInverter","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived");
	  return 0;
	}
      }
      assert(nuisPdf);
      Info("RunHypoTestInverter","Using as nuisance Pdf ... " );
      nuisPdf->Print();
      
      const RooArgSet * nuisParams = (pBModel->GetNuisanceParameters() ) ? pBModel->GetNuisanceParameters() : pSbModel->GetNuisanceParameters();
      RooArgSet * np = nuisPdf->getObservables(*nuisParams);
      if (np->getSize() == 0) { 
	Warning("RunHypoTestInverter","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
      }
      delete np;
      
      hhc->ForcePriorNuisanceAlt(*nuisPdf);
      hhc->ForcePriorNuisanceNull(*nuisPdf);
      
      
    }
  }
  else if (mCalculatorType == 2){
    hc = new AsymptoticCalculator(*pData, *pSbModel, *pBModel);
    ((AsymptoticCalculator*) hc)->SetOneSided(true); 
    // ((AsymptoticCalculator*) hc)->SetQTilde(true); // not needed should be done automatically now
    //((AsymptoticCalculator*) hc)->SetPrintLevel(mPrintLevel+1); 
  }
  else {
    Error("RunHypoTestInverter","Invalid - calculator type = %d supported values are only :\n\t\t\t 0 (Frequentist) , 1 (Hybrid) , 2 (Asymptotic) ",mCalculatorType);
    return 0;
  }
   /*
   hc = new FrequentistCalculator(*pData, *pSbModel, *pBModel);
   */

  // set the test statistic 
  TestStatistic * testStat = 0;
  if (mTestStatType == 0) testStat = &slrts;
  if (mTestStatType == 1) testStat = &ropl;
  if (mTestStatType == 2 || mTestStatType == 3) testStat = &profll;
  if (mTestStatType == 4) testStat = &maxll;
  if (testStat == 0) { 
    Error("RunHypoTestInverter","Invalid - test statistic type = %d supported values are only :\n\t\t\t 0 (SLR) , 1 (Tevatron) , 2 (PLR), 3 (PLR1), 4(MLE)",mTestStatType);
    return 0;
  }
  /*
  testStat = &profll;
   */

  // toy MC sampler
  ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
  toymcs->SetNEventsPerToy(1);
  toymcs->SetTestStatistic(testStat);
  toymcs->SetUseMultiGen(1);

  // get hypo test
  pHtr = hc->GetHypoTest();
  
  // clean up
  delete hc;

  return pHtr;
}



Int_t banner(){
  //#define __ROOFIT_NOBANNER // banner temporary off
#ifndef __EXOST_NOBANNER
  std::cout << desc << std::endl;
#endif
  return 0 ;
}
static Int_t dummy_ = banner() ;






Double_t roostats_cl95(Double_t ilum, Double_t slum,
		       Double_t eff, Double_t seff,
		       Double_t bck, Double_t sbck,
		       Int_t n,
		       Bool_t gauss,
		       Int_t nuisanceModel,
		       std::string method,
		       std::string plotFileName,
		       UInt_t seed,
		       LimitResult * result){
  //
  // Global function to run the CL95 routine
  // 
  // If a non-null pointer to a LimitResult object is provided,
  // it will be filled, and the caller keeps the ownership of
  // the object. This is mainly an internal interface design solution,
  // users are not expected to use that (but they can of course)
  //

  // limit calculation
  //CL95Calc theCalc(seed);
  CL95Calc * theCalc = CL95Calc::GetInstance();
  theCalc->SetSeed(seed);

  //mbGaussianStatistics = gauss;
  //if (n == 0) gauss = kFALSE;

  theCalc->PrintMethodInfo( method );

  // container for computed limits
  LimitResult limitResult;

  RooWorkspace * pWs = theCalc->MakeWorkspace( ilum, slum,
					       eff, seff,
					       bck, sbck,
					       n,
					       gauss,
					       nuisanceModel );
  
  //RooDataSet * data = (RooDataSet *)( theCalc->makeData( n )->Clone() );
  //data->SetName("observed_data");
  //pWs->import(*data);

  pWs->Print();

  pWs->SaveAs("ws.root");

  // if only workspace requested, exit here
  if ( method.find("workspace") != std::string::npos ) return 0.0;

  Double_t limit = theCalc->cl95( method, &limitResult );
  std::cout << "[roostats_cl95]: 95% C.L. upper limit: " << limit << std::endl;

  // check if the plot is requested
  if (plotFileName.size() != 0){
    theCalc->makePlot(method, plotFileName);
  }

  if (result) *result = limitResult;

  return limit;
}



LimitResult
GetClsLimit( Double_t ilum, Double_t slum,
	     Double_t eff, Double_t seff,
	     Double_t bck, Double_t sbck,
	     Int_t n ){
  //
  // Compute observed and expected CLs limit
  //
  
  std::string legend = "[GetClsLimit()]: ";

  LimitResult limit_result;

  bool gaussian_statistics   = GetParameter("GaussianStatistics",bool());
  int nuisance_model         = 1;
  std::string method         = "cls";
  std::string plot_file_name = "dummy.pdf";
  int seed                   = GetParameter("RandomSeed",int());

  roostats_cl95(ilum, slum,
		eff, seff,
		bck, sbck,
		n,
		gaussian_statistics,
		nuisance_model,
		method,
		plot_file_name,
		seed,
		&limit_result);

  if ( GetParameter( "Verbosity", int() )>1 ){
    std::cout << legend << " expected limit (median) " << limit_result.GetExpectedLimit() << std::endl;
    std::cout << legend << " expected limit (-1 sig) " << limit_result.GetOneSigmaLowRange() << std::endl;
    std::cout << legend << " expected limit (+1 sig) " << limit_result.GetOneSigmaHighRange() << std::endl;
    std::cout << legend << " expected limit (-2 sig) " << limit_result.GetTwoSigmaLowRange() << std::endl;
    std::cout << legend << " expected limit (+2 sig) " << limit_result.GetTwoSigmaHighRange() << std::endl;
  }
  
  return limit_result;
}



double
GetBayesianLimit( Double_t ilum, Double_t slum,
		  Double_t eff, Double_t seff,
		  Double_t bck, Double_t sbck,
		  Int_t n,
		  std::string option ){
  //
  // Compute observed bayesian limit
  //
  // option: ""     - numerical integration (default)
  //         "mcmc" - Metropolis-Hastings sampling
  //
  
  std::string legend = "[GetBayesianLimit()]: ";

  LimitResult limit_result;

  bool gaussian_statistics   = GetParameter("GaussianStatistics",bool());
  int nuisance_model         = 1;
  std::string method         = "bayesian";
  if (option.find("mcmc")!=std::string::npos) method = "mcmc";
  std::string plot_file_name = "";
  if ( GetParameter("MakePlot", bool()) ) plot_file_name = "bayesian_posterior.pdf";
  int seed                   = GetParameter("RandomSeed",int());

  roostats_cl95(ilum, slum,
		eff, seff,
		bck, sbck,
		n,
		gaussian_statistics,
		nuisance_model,
		method,
		plot_file_name,
		seed,
		&limit_result);

  if ( GetParameter( "Verbosity", int() )>1 ){
    std::cout << legend << " limit = " << limit_result.GetObservedLimit() << std::endl;
  }
  
  return limit_result.GetObservedLimit();
}



LimitResult
GetExpectedLimit( Double_t ilum, Double_t slum,
		  Double_t eff, Double_t seff,
		  Double_t bck, Double_t sbck,
		  Int_t nit,
		  std::string method ){
  //
  // Compute expected limit quantiles for arbitrary method
  //
  
  std::string legend = "[GetExpectedLimit()]: ";

  LimitResult limit_result;

  //bool gaussian_statistics   = GetParameter("GaussianStatistics",bool());
  int nuisance_model         = 1;
  int seed                   = GetParameter("RandomSeed",int());

  limit_result = roostats_clm(ilum, slum,
			      eff, seff,
			      bck, sbck,
			      nit,
			      nuisance_model,
			      method,
			      seed);

  if ( GetParameter( "Verbosity", int() )>1 ){
    std::cout << legend << " expected limit (median) " << limit_result.GetExpectedLimit() << std::endl;
    std::cout << legend << " expected limit (-1 sig) " << limit_result.GetOneSigmaLowRange() << std::endl;
    std::cout << legend << " expected limit (+1 sig) " << limit_result.GetOneSigmaHighRange() << std::endl;
    std::cout << legend << " expected limit (-2 sig) " << limit_result.GetTwoSigmaLowRange() << std::endl;
    std::cout << legend << " expected limit (+2 sig) " << limit_result.GetTwoSigmaHighRange() << std::endl;
  }
  
  return limit_result;
}



LimitResult roostats_limit(Double_t ilum, Double_t slum,
			   Double_t eff, Double_t seff,
			   Double_t bck, Double_t sbck,
			   Int_t n,
			   Bool_t gauss,
			   Int_t nuisanceModel,
			   std::string method,
			   std::string plotFileName,
			   UInt_t seed){
  //
  // Global function to run the CL95 routine
  // 

  LimitResult limitResult;

  roostats_cl95(ilum, slum,
		eff, seff,
		bck, sbck,
		n,
		gauss,
		nuisanceModel,
		method,
		plotFileName,
		seed,
		&limitResult);

  std::cout << " expected limit (median) " << limitResult.GetExpectedLimit() << std::endl;
  std::cout << " expected limit (-1 sig) " << limitResult.GetOneSigmaLowRange() << std::endl;
  std::cout << " expected limit (+1 sig) " << limitResult.GetOneSigmaHighRange() << std::endl;
  std::cout << " expected limit (-2 sig) " << limitResult.GetTwoSigmaLowRange() << std::endl;
  std::cout << " expected limit (+2 sig) " << limitResult.GetTwoSigmaHighRange() << std::endl;
  
  return limitResult;
}



Double_t roostats_cla(Double_t ilum, Double_t slum,
		      Double_t eff, Double_t seff,
		      Double_t bck, Double_t sbck,
		      Int_t nuisanceModel,
		      std::string method,
		      UInt_t seed){
  //
  // Global function to run old-style average limit routine.
  // Please use roostats_clm() instead.
  //

  Double_t limit = -1.0;

  std::cout << "[roostats_cla]: estimating average 95% C.L. upper limit" << endl;
  if (method.find("bayesian") != std::string::npos){
    std::cout << "[roostats_cla]: using Bayesian calculation via numeric integration" << endl;
  }
  else if (method.find("mcmc") != std::string::npos){
    std::cout << "[roostats_cl95]: using Bayesian calculation via numeric integration" << endl;
  }
  else if (method.find("cls") != std::string::npos){
    std::cout << "[roostats_cl95]: using CLs calculation" << endl;
  }
  else if (method.find("fc") != std::string::npos){
    std::cout << "[roostats_cl95]: using Feldman-Cousins approach" << endl;
  }
  else{
    std::cout << "[roostats_cla]: method " << method 
	      << " is not implemented, exiting" <<std::endl;
    return -1.0;
  }

  std::cout << "[roostats_cla]: Poisson statistics used" << endl;
    
  //CL95Calc theCalc(seed);
  CL95Calc * theCalc = CL95Calc::GetInstance();
  theCalc->SetSeed(seed);
  limit = theCalc->cla( ilum, slum,
		       eff, seff,
		       bck, sbck,
		       nuisanceModel,
		       method );

  //std::cout << "[roostats_cla]: average 95% C.L. upper limit: " << limit << std::endl;

  return limit;
}



LimitResult roostats_clm(Double_t ilum, Double_t slum,
			 Double_t eff, Double_t seff,
			 Double_t bck, Double_t sbck,
			 Int_t nit, Int_t nuisanceModel,
			 std::string method,
			 UInt_t seed){
  //
  // Global function to evaluate median expected limit and 1/2 sigma bands.
  //
  
  LimitResult limit;

  std::cout << "[roostats_clm]: estimating expected 95% C.L. upper limit" << endl;
  if (method.find("bayesian") != std::string::npos){
    std::cout << "[roostats_clm]: using Bayesian calculation via numeric integration" << endl;
  }
  else if (method.find("mcmc") != std::string::npos){
    std::cout << "[roostats_cl95]: using Bayesian calculation via numeric integration" << endl;
  }
  else if (method.find("plr") != std::string::npos){
    std::cout << "[roostats_cl95]: using profile likelihood ratio calculation" << endl;
  }
  else if (method.find("cls") != std::string::npos){
    std::cout << "[roostats_cl95]: using CLs calculation" << endl;
  }
  else if (method.find("fc") != std::string::npos){
    std::cout << "[roostats_cl95]: using Feldman-Cousins approach" << endl;
  }
  else{
    std::cout << "[roostats_clm]: method " << method 
	      << "is not implemented, exiting" <<std::endl;
    return limit;
  }

  std::cout << "[roostats_clm]: Poisson statistics used" << endl;
    
  //CL95Calc theCalc(seed);
  CL95Calc * theCalc = CL95Calc::GetInstance();
  theCalc->SetSeed(seed);
  limit = theCalc->clm( ilum, slum,
		       eff, seff,
		       bck, sbck,
		       nit, nuisanceModel,
		       method );

  return limit;
}



double
roostats_zscore( Double_t ilum, Double_t slum,
		 Double_t eff, Double_t seff,
		 Double_t bck, Double_t sbck,
		 Int_t n,
		 Bool_t gauss,
		 Int_t nuisanceModel,
		 std::string method,
		 std::string plotFileName,
		 UInt_t seed ){
  //
  // Estimate z-score (a measure of signiicance)
  //

  std::string legend = "[roostats_zscore]: ";

  double zscore = std::numeric_limits<double>::min();

  CL95Calc * theCalc = CL95Calc::GetInstance();
  theCalc->SetSeed(seed);

  theCalc->PrintMethodInfo( method );

  theCalc->MakeWorkspace( ilum, slum,
			  eff, seff,
			  bck, sbck,
			  n,
			  gauss,
			  nuisanceModel );

  RooStats::HypoTestResult * pHtr = theCalc->GetHypoTest("BModel", "");

  double nullPValue = pHtr->NullPValue();
  zscore = pHtr->Significance();

  std::cout << legend
	    << "Null p-value = " << nullPValue << std::endl;

  std::cout << legend
	    << "Z-score = " << zscore << std::endl;

  (void)plotFileName;

  // clean up
  delete pHtr;

  return zscore;
}



void SetParameter(const char * name, const char * value){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  theCalc->SetParameter(name, value);
  
  return;
}



void SetParameter(const char * name, bool value){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  theCalc->SetParameter(name, value);
  
  return;
}



void SetParameter(const char * name, int value){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  theCalc->SetParameter(name, value);
  
  return;
}



void SetParameter(const char * name, double value){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  theCalc->SetParameter(name, value);
  
  return;
}



bool GetParameter(const char * name, bool){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  return theCalc->GetParameter(name, bool());
}



int GetParameter(const char * name, int){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  return theCalc->GetParameter(name, int());
}



double GetParameter(const char * name, double){
  CL95Calc * theCalc = CL95Calc::GetInstance();
  
  return theCalc->GetParameter(name, double());
}



/////////////////////////////////////////////////////////////////////////
//
// CLs helper methods from roostats macro StandardHypoTestInvDemo.C
// This is the core of the CLs calculation
//
// All this will be factored out in a ROOT class in near future
// and will disappear from here
//



// internal class to run the inverter and more

namespace RooStats { 

   class HypoTestInvTool{

   public:
      HypoTestInvTool();
      ~HypoTestInvTool(){};

      HypoTestInverterResult * 
      RunInverter(RooWorkspace * w, 
                  const char * modelSBName, const char * modelBName, 
                  const char * dataName,
                  int type,  int testStatType, 
                  bool useCLs, 
                  int npoints, double poimin, double poimax, int ntoys, 
                  bool useNumberCounting = false, 
                  const char * nuisPriorName = 0);



      void
      AnalyzeResult( HypoTestInverterResult * r,
                     int calculatorType,
                     int testStatType, 
                     bool useCLs,  
                     int npoints,
                     const char * fileNameBase = 0 );

      void SetParameter(const char * name, const char * value);
      void SetParameter(const char * name, bool value);
      void SetParameter(const char * name, int value);
      void SetParameter(const char * name, double value);

   private:

     // GENA: modified
      bool mNoSystematics;
      bool mPlotHypoTestResult;
      bool mWriteResult;
      bool mOptimize;
      bool mUseVectorStore;
      bool mGenerateBinned;
      bool mUseProof;
      bool mRebuild;
      int     mNWorkers;
      int     mNToyToRebuild;
      int     mPrintLevel;
      double  mNToysRatio;
      double  mMaxPoi;
      std::string mMassValue;
      std::string mMinimizerType;                  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()

   };

} // end namespace RooStats



// GENA: modified
RooStats::HypoTestInvTool::HypoTestInvTool() : mNoSystematics(false),
					       mPlotHypoTestResult(true),
                                               mWriteResult(false),
                                               mOptimize(true),
                                               mUseVectorStore(true),
                                               mGenerateBinned(false),
                                               mUseProof(false),
                                               mRebuild(false),
                                               mNWorkers(4),
                                               mNToyToRebuild(100),
                                               mPrintLevel(0),
                                               mNToysRatio(2),
                                               mMaxPoi(-1),
                                               mMassValue(""),
                                               mMinimizerType(""){
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, bool value){
   //
   // set boolean parameters
   //

   std::string s_name(name);

   if (s_name.find("PlotHypoTestResult") != std::string::npos) mPlotHypoTestResult = value;
   if (s_name.find("WriteResult") != std::string::npos) mWriteResult = value;
   if (s_name.find("Optimize") != std::string::npos) mOptimize = value;
   if (s_name.find("UseVectorStore") != std::string::npos) mUseVectorStore = value;
   if (s_name.find("GenerateBinned") != std::string::npos) mGenerateBinned = value;
   if (s_name.find("UseProof") != std::string::npos) mUseProof = value;
   if (s_name.find("Rebuild") != std::string::npos) mRebuild = value;

   return;
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, int value){
   //
   // set integer parameters
   //

   std::string s_name(name);

   if (s_name.find("NWorkers") != std::string::npos) mNWorkers = value;
   if (s_name.find("NToyToRebuild") != std::string::npos) mNToyToRebuild = value;
   if (s_name.find("PrintLevel") != std::string::npos) mPrintLevel = value;

   return;
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, double value){
   //
   // set double precision parameters
   //

   std::string s_name(name);

   if (s_name.find("NToysRatio") != std::string::npos) mNToysRatio = value;
   if (s_name.find("MaxPoi") != std::string::npos) mMaxPoi = value;

   return;
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, const char * value){
   //
   // set string parameters
   //

   std::string s_name(name);

   if (s_name.find("MassValue") != std::string::npos) mMassValue.assign(value);
   if (s_name.find("MinimizerType") != std::string::npos) mMinimizerType.assign(value);

   return;
}



HypoTestInverterResult * 
CL95Calc::RunHypoTestInverter(RooWorkspace * w,
			      const char * modelSBName,
			      const char * modelBName,
			      const char * dataName,
			      int calculatorType,
			      int testStatType,
			      bool useCLs,
			      int npoints,
			      double poimin,
			      double poimax,
			      int ntoys,
			      bool useNumberCounting,
			      const char * nuisPriorName){
  //
  // adopted from StandardHypoTestInvDemo()
  //
  
   HypoTestInvTool calc;

   // set parameters
   // GENA: modified
   calc.SetParameter("NoSystematics", false);
   calc.SetParameter("PlotHypoTestResult", GetParameter("MakePlot",bool()));
   calc.SetParameter("WriteResult", false);
   calc.SetParameter("Optimize", GetParameter("Optimize",bool()));
   calc.SetParameter("UseVectorStore", true);
   calc.SetParameter("GenerateBinned", false);
   calc.SetParameter("NToysRatio", GetParameter("NToysRatio",double()));
   calc.SetParameter("MaxPOI", -1.0);
   calc.SetParameter("UseProof", false);
   calc.SetParameter("Nworkers", 2);
   calc.SetParameter("Rebuild", false);
   calc.SetParameter("NToyToRebuild", 100);
   calc.SetParameter("MassValue", "");
   calc.SetParameter("MinimizerType", "");
   calc.SetParameter("PrintLevel", 0);


   //RooWorkspace * w = dynamic_cast<RooWorkspace*>( file->Get(wsName) );
   HypoTestInverterResult * r = 0;  
   //std::cout << w << "\t" << fileName << std::endl;
   if (w != NULL) {
      r = calc.RunInverter(w, modelSBName, modelBName,
                           dataName, calculatorType, testStatType, useCLs,
                           npoints, poimin, poimax,  
                           ntoys, useNumberCounting, nuisPriorName );    
      if (!r) { 
         std::cerr << "Error running the HypoTestInverter - Exit " << std::endl;
         return 0;          
      }
   }
  
   calc.AnalyzeResult( r, calculatorType, testStatType, useCLs, npoints, "cls.pdf" );
  
   return r;
}



void
RooStats::HypoTestInvTool::AnalyzeResult( HypoTestInverterResult * r,
                                          int calculatorType,
                                          int testStatType, 
                                          bool useCLs,  
                                          int npoints,
                                          const char * fileNameBase ){
  
  // analyize result produced by the inverter, optionally save it in a file 
  
  //double upperLimit = r->UpperLimit();
  //double ulError = r->UpperLimitEstimatedError();
  
   //std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
  
   // compute expected limit
   //std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
   //std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
   //std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
   //std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
   //std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;
  
  
   // write result in a file 
   if (r != NULL && mWriteResult) {
    
      // write to a file the results
      const char *  calcType = (calculatorType == 0) ? "Freq" : (calculatorType == 1) ? "Hybr" : "Asym";
      const char *  limitType = (useCLs) ? "CLs" : "Cls+b";
      const char * scanType = (npoints < 0) ? "auto" : "grid";
      TString resultFileName = TString::Format("%s_%s_%s_ts%d_",calcType,limitType,scanType,testStatType);      
      //strip the / from the filename
      if (mMassValue.size()>0) {
         resultFileName += mMassValue.c_str();
         resultFileName += "_";
      }
    
      TString name = fileNameBase; 
      name.Replace(0, name.Last('/')+1, "");
      resultFileName += name;
    
      TFile * fileOut = new TFile(resultFileName,"RECREATE");
      r->Write();
      fileOut->Close();                                                                     
   }   
  
  
   // plot the result ( p values vs scan points) 
   std::string typeName = "";
   if (calculatorType == 0 )
      typeName = "Frequentist";
   if (calculatorType == 1 )
      typeName = "Hybrid";   
   else if (calculatorType == 2 ) { 
      typeName = "Asymptotic";
      mPlotHypoTestResult = false; 
   }
  
   const char * resultName = r->GetName();
   TString plotTitle = TString::Format("%s CL Scan for workspace %s",typeName.c_str(),resultName);
   // GENA:
   if (mPlotHypoTestResult) { 
     TCanvas * c33 = new TCanvas();
     HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);
     plot->Draw("CLb 2CL");  // plot all and Clb
     c33->SaveAs("cls.pdf");
     delete c33;
     
     
     const int nEntries = r->ArraySize();
     
     // plot test statistics distributions for the two hypothesis 
     
     TCanvas * c2 = new TCanvas();
     if (nEntries > 1) { 
         int ny = TMath::CeilNint( sqrt(nEntries) );
         int nx = TMath::CeilNint(double(nEntries)/ny);
         c2->Divide( nx,ny);
     }
     for (int i=0; i<nEntries; i++) {
       if (nEntries > 1) c2->cd(i+1);
       SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
       pl->SetLogYaxis(true);
       pl->Draw();
     }
     c2->SaveAs("cls_sampling.pdf");
     delete c2;
   }
}



// internal routine to run the inverter
HypoTestInverterResult *
RooStats::HypoTestInvTool::RunInverter(RooWorkspace * w,
                                       const char * modelSBName, const char * modelBName, 
                                       const char * dataName, int type,  int testStatType, 
                                       bool useCLs, int npoints, double poimin, double poimax, 
                                       int ntoys,
                                       bool useNumberCounting,
                                       const char * nuisPriorName ){

   std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;
  
   w->Print();
  
  
   RooAbsData * data = w->data(dataName); 
   if (!data) { 
      Error("StandardHypoTestDemo","Not existing data %s",dataName);
      return 0;
   }
   else 
      std::cout << "Using data set " << dataName << std::endl;
  
   if (mUseVectorStore) { 
      RooAbsData::defaultStorageType = RooAbsData::Vector;
      data->convertToVectorStore() ;
   }
  
  
   // get models from WS
   // get the modelConfig out of the file
   ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
   ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);
  
   if (!sbModel) {
      Error("StandardHypoTestDemo","Not existing ModelConfig %s",modelSBName);
      return 0;
   }
   // check the model 
   if (!sbModel->GetPdf()) { 
      Error("StandardHypoTestDemo","Model %s has no pdf ",modelSBName);
      return 0;
   }
   if (!sbModel->GetParametersOfInterest()) {
      Error("StandardHypoTestDemo","Model %s has no poi ",modelSBName);
      return 0;
   }
   if (!sbModel->GetObservables()) {
      Error("RunHypoTestInverter","Model %s has no observables ",modelSBName);
      return 0;
   }
   if (!sbModel->GetSnapshot() ) { 
      Info("RunHypoTestInverter","Model %s has no snapshot  - make one using model poi",modelSBName);
      sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
   }
  
   // case of no systematics
   // remove nuisance parameters from model
   if (mNoSystematics) { 
      const RooArgSet * nuisPar = sbModel->GetNuisanceParameters();
      if (nuisPar && nuisPar->getSize() > 0) { 
         std::cout << "RunHypoTestInverter" << "  -  Switch off all systematics by setting them constant to their initial values" << std::endl;
         RooStats::SetAllConstant(*nuisPar);
      }
      if (bModel) { 
         const RooArgSet * bnuisPar = bModel->GetNuisanceParameters();
         if (bnuisPar) 
            RooStats::SetAllConstant(*bnuisPar);
      }
   }
  
   if (!bModel || bModel == sbModel) {
      Info("RunHypoTestInverter","The background model %s does not exist",modelBName);
      Info("RunHypoTestInverter","Copy it from ModelConfig %s and set POI to zero",modelSBName);
      bModel = (ModelConfig*) sbModel->Clone();
      bModel->SetName(TString(modelSBName)+TString("_with_poi_0"));      
      RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
      if (!var) return 0;
      double oldval = var->getVal();
      var->setVal(0);
      bModel->SetSnapshot( RooArgSet(*var)  );
      var->setVal(oldval);
   }
   else { 
      if (!bModel->GetSnapshot() ) { 
         Info("RunHypoTestInverter","Model %s has no snapshot  - make one using model poi and 0 values ",modelBName);
         RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
         if (var) { 
            double oldval = var->getVal();
            var->setVal(0);
            bModel->SetSnapshot( RooArgSet(*var)  );
            var->setVal(oldval);
         }
         else { 
            Error("RunHypoTestInverter","Model %s has no valid poi",modelBName);
            return 0;
         }         
      }
   }
  
   //-----> GENA:
   // estimate limit with PLR and set range
//   RooStats::ModelConfig * _mc = sbModel;
//   RooAbsData * _d = w->data("asimovData");
//   RooStats::ProfileLikelihoodCalculator plc(*_d, *_mc);
//   plc.SetConfidenceLevel(mConfidenceLevel);
//   RooStats::LikelihoodInterval * pPlrInt = plc.GetInterval();
//   Double_t _poimax = pPlrInt->UpperLimit( *w->var("xsec") );
//   std::cout << "POIMAX!!!!!!!!!!!!!*******************: "
//	     << _poimax << std::endl;
//   w->var("xsec")->setRange(0,_poimax);
//   w->var("xsec")->setRange(0,75);
   //std::exit(-1);

   //------------------------->

   // GENA:
   //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit","migradimproved");

   // run first a data fit 
  
   const RooArgSet * poiSet = sbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar*)poiSet->first();
  
   std::cout << "RunHypoTestInverter : POI initial value:   " << poi->GetName() << " = " << poi->getVal()   << std::endl;  
  
   // fit the data first (need to use constraint )
   Info( "RunHypoTestInverter"," Doing a first fit to the observed data ");
   // GENA: modified
   if (mMinimizerType.size()==0) mMinimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
   else 
      ROOT::Math::MinimizerOptions::SetDefaultMinimizer(mMinimizerType.c_str());
   Info("RunHypoTestInverter","Using %s as minimizer for computing the test statistic",
        ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str() );
   RooArgSet constrainParams;
   if (sbModel->GetNuisanceParameters() ) constrainParams.add(*sbModel->GetNuisanceParameters());
   RooStats::RemoveConstantParameters(&constrainParams);
   TStopwatch tw; 
   tw.Start(); 
   RooFitResult * fitres = sbModel->GetPdf()->fitTo(*data,InitialHesse(false),
						    Hesse(false),
						    Minimizer(mMinimizerType.c_str(),"Migrad"),
						    Strategy(0), PrintLevel(mPrintLevel+1),
						    Constrain(constrainParams), Save(true) );
   if (fitres->status() != 0) { 
      Warning("RunHypoTestInverter","Fit to the model failed - try with strategy 1 and perform first an Hesse computation");
      fitres = sbModel->GetPdf()->fitTo(*data,InitialHesse(true), Hesse(false),Minimizer(mMinimizerType.c_str(),"Migrad"), Strategy(1), PrintLevel(mPrintLevel+1), Constrain(constrainParams), Save(true) );
   }
   if (fitres->status() != 0) 
      Warning("RunHypoTestInverter"," Fit still failed - continue anyway.....");
  
  
   double poihat  = poi->getVal();
   std::cout << "RunHypoTestInverter - Best Fit value : " << poi->GetName() << " = "  
             << poihat << " +/- " << poi->getError() << std::endl;
   std::cout << "Time for fitting : "; tw.Print(); 
  
   //save best fit value in the poi snapshot 
   sbModel->SetSnapshot(*sbModel->GetParametersOfInterest());
   std::cout << "StandardHypoTestInvo: snapshot of S+B Model " << sbModel->GetName() 
             << " is set to the best fit value" << std::endl;
  
   // build test statistics and hypotest calculators for running the inverter 
  
   SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(),*bModel->GetPdf());
  
   if (sbModel->GetSnapshot()) slrts.SetNullParameters(*sbModel->GetSnapshot());
   if (bModel->GetSnapshot()) slrts.SetAltParameters(*bModel->GetSnapshot());
  
   // ratio of profile likelihood - need to pass snapshot for the alt
   RatioOfProfiledLikelihoodsTestStat 
      ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   ropl.SetPrintLevel(mPrintLevel);
   ropl.SetMinimizer(mMinimizerType.c_str());
  
   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   if (testStatType == 3) profll.SetOneSided(1);
   profll.SetMinimizer(mMinimizerType.c_str());
   profll.SetPrintLevel(mPrintLevel);

   profll.SetReuseNLL(mOptimize);
   slrts.SetReuseNLL(mOptimize);
   ropl.SetReuseNLL(mOptimize);

   if (mOptimize) { 
      profll.SetStrategy(0);
      ropl.SetStrategy(0);
   }
  
   if (mMaxPoi > 0) poi->setMax(mMaxPoi);  // increase limit
  
   MaxLikelihoodEstimateTestStat maxll(*sbModel->GetPdf(),*poi); 
  
   // create the HypoTest calculator class 
   HypoTestCalculatorGeneric *  hc = 0;
   if (type == 0) hc = new FrequentistCalculator(*data, *bModel, *sbModel);
   else if (type == 1) hc = new HybridCalculator(*data, *bModel, *sbModel);
   else if (type == 2) hc = new AsymptoticCalculator(*data, *bModel, *sbModel);
   else {
      Error("RunHypoTestInverter","Invalid - calculator type = %d supported values are only :\n\t\t\t 0 (Frequentist) , 1 (Hybrid) , 2 (Asymptotic) ",type);
      return 0;
   }
  
   // set the test statistic 
   TestStatistic * testStat = 0;
   if (testStatType == 0) testStat = &slrts;
   if (testStatType == 1) testStat = &ropl;
   if (testStatType == 2 || testStatType == 3) testStat = &profll;
   if (testStatType == 4) testStat = &maxll;
   if (testStat == 0) { 
      Error("RunHypoTestInverter","Invalid - test statistic type = %d supported values are only :\n\t\t\t 0 (SLR) , 1 (Tevatron) , 2 (PLR), 3 (PLR1), 4(MLE)",testStatType);
      return 0;
   }
  
  
   ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
   if (toymcs) { 
      if (useNumberCounting) toymcs->SetNEventsPerToy(1);
      toymcs->SetTestStatistic(testStat);
    
      if (data->isWeighted() && !mGenerateBinned) { 
         Info("RunHypoTestInverter","Data set is weighted, nentries = %d and sum of weights = %8.1f but toy generation is unbinned - it would be faster to set mGenerateBinned to true\n",data->numEntries(), data->sumEntries());
      }
      toymcs->SetGenerateBinned(mGenerateBinned);
    
      toymcs->SetUseMultiGen(mOptimize);
    
      if (mGenerateBinned &&  sbModel->GetObservables()->getSize() > 2) { 
         Warning("RunHypoTestInverter","generate binned is activated but the number of ovservable is %d. Too much memory could be needed for allocating all the bins",sbModel->GetObservables()->getSize() );
      }
    
   }
  
  
   if (type == 1) { 
      HybridCalculator *hhc = dynamic_cast<HybridCalculator*> (hc);
      assert(hhc);
    
      hhc->SetToys(ntoys,ntoys/mNToysRatio); // can use less ntoys for b hypothesis 
    
      // remove global observables from ModelConfig (this is probably not needed anymore in 5.32)
      bModel->SetGlobalObservables(RooArgSet() );
      sbModel->SetGlobalObservables(RooArgSet() );
    
    
      // check for nuisance prior pdf in case of nuisance parameters 
      if (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters() ) {
         RooAbsPdf * nuisPdf = 0; 
         if (nuisPriorName) nuisPdf = w->pdf(nuisPriorName);
         // use prior defined first in bModel (then in SbModel)
         if (!nuisPdf)  { 
            Info("RunHypoTestInverter","No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the model");
            if (bModel->GetPdf() && bModel->GetObservables() ) 
               nuisPdf = RooStats::MakeNuisancePdf(*bModel,"nuisancePdf_bmodel");
            else 
               nuisPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisancePdf_sbmodel");
         }   
         if (!nuisPdf ) {
            if (bModel->GetPriorPdf())  { 
               nuisPdf = bModel->GetPriorPdf();
               Info("RunHypoTestInverter","No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf->GetName());            
            }
            else { 
               Error("RunHypoTestInverter","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived");
               return 0;
            }
         }
         assert(nuisPdf);
         Info("RunHypoTestInverter","Using as nuisance Pdf ... " );
         nuisPdf->Print();
      
         const RooArgSet * nuisParams = (bModel->GetNuisanceParameters() ) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
         RooArgSet * np = nuisPdf->getObservables(*nuisParams);
         if (np->getSize() == 0) { 
            Warning("RunHypoTestInverter","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
         }
         delete np;
      
         hhc->ForcePriorNuisanceAlt(*nuisPdf);
         hhc->ForcePriorNuisanceNull(*nuisPdf);
      
      
      }
   } 
   else if (type == 2) { 
      ((AsymptoticCalculator*) hc)->SetOneSided(true); 
      // ((AsymptoticCalculator*) hc)->SetQTilde(true); // not needed should be done automatically now
      ((AsymptoticCalculator*) hc)->SetPrintLevel(mPrintLevel+1); 
   }
   else if (type != 2) 
      ((FrequentistCalculator*) hc)->SetToys(ntoys,ntoys/mNToysRatio); 
  
   // Get the result
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  
  
  
   HypoTestInverter calc(*hc);
   // GENA: conf level
   calc.SetConfidenceLevel(GetParameter("ConfidenceLevel",double()));
  
  
   calc.UseCLs(useCLs);
   calc.SetVerbose(true);
  
   // can speed up using proof-lite
   if (mUseProof && mNWorkers > 1) { 
      ProofConfig pc(*w, mNWorkers, "", kFALSE);
      toymcs->SetProofConfig(&pc);    // enable proof
   }
  
  
   if (npoints > 0) {
      if (poimin > poimax) { 
         // if no min/max given scan between MLE and +4 sigma 
         poimin = int(poihat);
         poimax = int(poihat +  4 * poi->getError());
      }
      std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
      calc.SetFixedScan(npoints,poimin,poimax);
   }
   else { 
      //poi->setMax(10*int( (poihat+ 10 *poi->getError() )/10 ) );
      std::cout << "Doing an  automatic scan  in interval : " << poi->getMin() << " , " << poi->getMax() << std::endl;
   }
  
   //------> GENA:
   //minimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
   //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit","migrad");

   tw.Start();
   HypoTestInverterResult * r = calc.GetInterval();
   std::cout << "Time to perform limit scan \n";
   tw.Print();
  
   if (mRebuild) {
      calc.SetCloseProof(1);
      tw.Start();
      SamplingDistribution * limDist = calc.GetUpperLimitDistribution(true,mNToyToRebuild);
      std::cout << "Time to rebuild distributions " << std::endl;
      tw.Print();
    
      if (limDist) { 
         std::cout << "expected up limit " << limDist->InverseCDF(0.5) << " +/- " 
                   << limDist->InverseCDF(0.16) << "  " 
                   << limDist->InverseCDF(0.84) << "\n"; 
      
         //update r to a new updated result object containing the rebuilt expected p-values distributions
         // (it will not recompute the expected limit)
         if (r) delete r;  // need to delete previous object since GetInterval will return a cloned copy
         r = calc.GetInterval();
      
      }
      else 
         std::cout << "ERROR : failed to re-build distributions " << std::endl; 
   }
  
   return r;
}



