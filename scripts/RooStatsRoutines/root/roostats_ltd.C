static const char* desc =
"=====================================================================\n"
"|                                                                    \n"
"|\033[1m        roostats_ltd.C  alpha                         \033[0m\n"
"|                                                                    \n"
"| Standard: RooStats-based c++ routine for 95% C.L. limit calculation\n"
"|                                                                    \n"
"| Model-independent: accepts any properly configured RooFit model    \n"
"|                                                                    \n"
"| Standalone: depends only on ROOT                                   \n"
"|                                                                    \n"
"|\033[1m Gena Kukartsev                                       \033[0m\n"
"|\033[1m Lorenzo Moneta (CLs core)                            \033[0m\n"
"|                                                                    \n"
"| September 2011: first version                                      \n"
"|                                                                    \n"
"=====================================================================\n"
"                                                                     \n"
"Prerequisites:                                                       \n"
"                ROOT version 5.30.02 or higher                       \n"
"                                                                     \n"
"                                                                     \n"
"                                                                    \n"
"The code should be compiled in ROOT:                                 \n"
"                                                                     \n"
"root -l                                                              \n"
"                                                                     \n"
".L roostats_ltd.C+                                                   \n"
"                                                                     \n"
"Usage:                                                               \n"
"                                                                     \n"
"Inputs:                                                              \n"
"                                                                     \n"
"For more details see                                                 \n"
"        https://twiki.cern.ch/twiki/bin/view/CMS/RooStatsLtd         \n"
"                                                                     \n"
"\033[1m       Note!                                           \033[0m\n"
"If you are running nonstandard ROOT environment, e.g. in CMSSW,      \n"
"you need to make sure that the RooFit and RooStats header files      \n"
"can be found since they might be in a nonstandard location.          \n"
"                                                                     \n"
"For CMSSW_4_2_0_pre8 and later, add the following line to your       \n"
"rootlogon.C:                                                         \n"
"      gSystem -> SetIncludePath( \"-I$ROOFITSYS/include\" );         \n";


#include <algorithm>

#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TBox.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/ProposalHelper.h"
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


class LimitResult; // container class for limit results

LimitResult limit( const char * method = "plr",
		   const char * inFileName = 0,
		   const char * workspaceName = 0,
		   const char * datasetName = 0 );

// ---> implementation below --------------------------------------------


class LimitResult{

  friend class LimitCalc;
  
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



class LimitCalc{

public:

  // no public constructors - this is a singleton class
  // use static method LimitCalc::GetInstance() to get the instance
  ~LimitCalc();

  static LimitCalc * GetInstance(void){
    if (!mspInstance) mspInstance = new LimitCalc();
    return mspInstance;
  }

  LimitResult       GetClsLimit( int nToys,
				 bool printResult = true );

  LimitResult       GetClsLimit( int nPoiScanPoints,
				 int nToys,
				 bool printResult = true );

  // run single point in the cls scan, with a given precision
  std::pair<double, double>
  GetClsSinglePoint( RooStats::HypoTestInverter & calc,
		     double poi,
		     double precision,
		     bool observed,
		     double bgSigma,
		     int maxNToys );

  RooWorkspace *    GetWorkspace(){ return mpWs;}
  RooAbsData *      GetData( void ){return mpData;};
  const RooArgSet * GetPoiSet( void );
  RooRealVar *      GetFirstPoi( void );
  Double_t          GetFirstPoiMin( void ){return GetFirstPoi()->getMin();};
  Double_t          GetFirstPoiMax( void ){return GetFirstPoi()->getMax();};

  RooAbsData *            LoadData( std::string datasetName );
  RooStats::ModelConfig * LoadModelConfig( std::string mcName );
  RooWorkspace *          LoadWorkspace( std::string wsFileName,
					 std::string wsName );

  Double_t GuessNextPoiStep( std::vector<double> & vPoi,
			     std::vector<double> & vCls,
			     std::vector<double> & vClsErr,
			     double cls );
    
  void SetSeed(UInt_t seed);
  void SetInverterCalcType(int type){mInverterCalcType=type;};
  void SetTestStatType(int type){mTestStatType=type;};
  void SetSbModelConfig( std::string modelConfigName );
  void SetBModelConfig( std::string modelConfigName );
  void SetSbModelConfig( RooStats::ModelConfig * pSbModelConfig );
  void SetBModelConfig( RooStats::ModelConfig * pBModelConfig );
  void SetData( RooAbsData * data );
  void SetWorkspace( RooWorkspace * data );
  void SetPlot(bool doPlot){mDoPlotHypoTestResult=doPlot;};
  void SetUseProof(bool useProof){mUseProof=useProof;};
  void SetOptimize(bool optimize){mOptimize=optimize;};
  void SetWriteResult(bool writeResult){mWriteResult=writeResult;};
  void SetNProofWorkers(int nProofWorkers){mNProofWorkers=nProofWorkers;};
  bool SetFirstPoiValue( double value );
  bool SetFirstPoiMax( double value );
  bool SetFirstPoiMin( double value );
  void SetTestMode( bool mode ){mTestMode = mode;};
  void SetRebuildSamplingDist( bool rebuild ){mRebuildSamplingDist = rebuild;};
  void SetNToysToRebuild( bool nToys ){mNToysToRebuild = nToys;};
  void SetNPoiScanPoints( bool nPoiScanPoints ){mNPoiScanPoints = nPoiScanPoints;};
  void SetMaxZeroToys( int ntoys ){mMaxZeroToys = ntoys;};

private:

  LimitCalc();
  LimitCalc(const LimitCalc &); // stop default

  void init( UInt_t seed ); //  to be called by constructor

  Double_t GetRandom( std::string pdfName, std::string varName );
  Long64_t FindLowBoundary(std::vector<Double_t> * pCdf, Double_t value);
  Long64_t FindHighBoundary(std::vector<Double_t> * pCdf, Double_t value);

  LimitResult ComputeInverterLimit( bool useCls,
				    int nPoints,
				    double poiMin,  // use default is poimin >= poimax
				    double poiMax,
				    int nToys );

  // internal routine to run the inverter
  // uses: workspace, mpData, sbmodel, bmodel, 
  //       calc type, test stat type
  RooStats::HypoTestInverterResult * RunInverter( int nPoints,
						  double poiMin,
						  double poiMax, 
						  int nSbToys,
						  int nBToys,
						  bool useCls );

  RooStats::HypoTestInverter * GetHypoTestInverter( void );
 
  // get the expected p-values for given quantiles
  // from hypotestinvresult
  Double_t GetExpectedPValue( RooStats::HypoTestInverterResult * pResult,
			      Int_t index,
			      Double_t nSigma );

  // get the uncertainty for a given expected p-value
  // from hypotestinvresult
  Double_t GetExpectedPValueErr( RooStats::HypoTestInverterResult * pResult,
				 Int_t index,
				 Double_t pValue );

  
  RooStats::MCMCInterval * GetMcmcInterval(double confLevel,
					   int nIter,
					   int nBurnIn,
					   double leftSideTailFraction,
					   int nBins);

  void PlotMcmcPosterior( std::string plotFileName );
  double PrintMcmcUpperLimit( std::string limitFileName = "" );

  Double_t RoundUpperBound(Double_t bound);


  // attributes
  bool mTestMode;
  int mInverterCalcType;
  int mTestStatType;

  // toy mc parameters
  int mNEventsPerToy; // ignored if S+B model PDF is extended

  // pointers whose objects we do not own
  RooStats::ModelConfig * mpSbModel;
  RooStats::ModelConfig * mpBModel;

  // data members (own those objects)
  RooWorkspace * mpWs;
  RooAbsData * mpData;
  RooStats::BayesianCalculator * mpBayesCalc;
  RooStats::SimpleInterval * mpSimpleInterval;

  // for Bayesian MCMC calculation
  RooStats::MCMCInterval * mpMcmcInterval;
  
  // for Feldman-Cousins Calculator
  RooStats::FeldmanCousins * mpFcCalc;

  // random numbers
  TRandom3 mRandom;

  // HypoTestInverter attr and objects (CLs, FC...)
  RooStats::HypoTestInverter * mpHypoTestInvCalc; // own it
  RooStats::HypoTestCalculatorGeneric *  mpHypoTestCalc;
  RooStats::FrequentistCalculator * mpFreqCalc;
  RooStats::HybridCalculator * mpHybridCalc;
  RooStats::SimpleLikelihoodRatioTestStat * slrts;
  RooStats::RatioOfProfiledLikelihoodsTestStat * ropl;
  RooStats::ProfileLikelihoodTestStat * profll;
  RooStats::MaxLikelihoodEstimateTestStat * maxll;
  bool mDoPlotHypoTestResult;
  bool mUseProof;
  bool mOptimize;
  bool mWriteResult;
  int  mNProofWorkers;
  std::string mNuisPriorName;
  bool mRebuildSamplingDist;
  int  mNToysToRebuild;
  int mNPoiScanPoints;
  int mMaxZeroToys; // max toys to throw when tail has zero

  // pointer to class instance
  static LimitCalc * mspInstance;
};



LimitCalc * LimitCalc::mspInstance = 0;



// default constructor
LimitCalc::LimitCalc(){
  init(12345);
}



void LimitCalc::init(UInt_t randomSeed){

  // set test mode
  mTestMode = false;

  // toy MC settings
  mNEventsPerToy = 0;

  // 0 - freq, 1 - hybrid
  mInverterCalcType = 0;

  // 0 - SimpleLikelihood (LEP), 
  // 1 - Profile likelihood ratio (Tevatron),
  // 2 - Two-sided profile likelihood,
  // 3 - One-sided profile likelihood
  mTestStatType = 3;

  // workspace-related
  mpWs = 0;
  mpSbModel = 0;
  mpBModel = 0;
  mpData = 0;

  // Bayesian atttributes
  mpBayesCalc = 0;
  mpSimpleInterval = 0;
  mpMcmcInterval = 0;

  // Feldman-Cousins attributes
  mpFcCalc = 0;

  // set random seed
  SetSeed(randomSeed);

  // HypoTestInverter attributes (CLs, FC...)
  mpHypoTestInvCalc = 0;
  mpHypoTestCalc = 0;
  mpFreqCalc = 0;
  mpHybridCalc = 0;
  slrts = 0;
  ropl = 0;
  profll = 0;
  maxll = 0;
  mDoPlotHypoTestResult =  false; 
  mUseProof =  false;
  mOptimize =  false;
  mWriteResult =  false;
  mNProofWorkers =  1;
  mNuisPriorName.clear();
  mRebuildSamplingDist = false;
  mNToysToRebuild = -1;
  mNPoiScanPoints = 0;
  mMaxZeroToys = 10000;
}


LimitCalc::~LimitCalc(){
  delete mpWs;
  delete mpData;
  delete mpBayesCalc;
  delete mpSimpleInterval;
  delete mpMcmcInterval;
  delete mpFcCalc;
  delete mpHypoTestInvCalc;
  delete mpHypoTestCalc;
  delete mpFreqCalc;
  delete mpHybridCalc;
  delete slrts;
  delete ropl;
  delete profll;
  delete maxll;
}



void LimitCalc::SetSeed( UInt_t seed ){
  //
  // Set random seed. If 0, set unique random.
  //
  std::string _legend = "[LimitCalc::SetSeed]: ";

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



const RooArgSet * LimitCalc::GetPoiSet( void ){
  //
  // Return a pointer to the set of parameters of interest
  // No ownership is transferred
  // Return null pinter if failed
  //
  if (!mpSbModel){
    return 0;
  }
  
  return mpSbModel->GetParametersOfInterest();
}



RooRealVar * LimitCalc::GetFirstPoi( void ){
  //
  // Return a pointer to the first parameters of interest
  // (often we only have one POI)
  // No ownership is transferred
  // Return null pinter if failed
  //

  const RooArgSet * p_poi_set = GetPoiSet();

  if (!p_poi_set){
    return 0;
  }

  RooAbsArg * first_poi = p_poi_set->first();

  return (RooRealVar *)first_poi;
}



bool LimitCalc::SetFirstPoiValue( double value ){
  //
  // Set value of the first POI
  // Return false if fail, true if success
  //
  RooRealVar * poi = GetFirstPoi();

  if (!poi) return false;

  poi->setVal(value);

  return true;
}



bool LimitCalc::SetFirstPoiMax( double value ){
  //
  // Set upper boundary of the range for the first POI
  // Return false if fail, true if success
  //
  RooRealVar * poi = GetFirstPoi();

  if (!poi) return false;

  poi->setMax(value);

  return true;
}



bool LimitCalc::SetFirstPoiMin( double value ){
  //
  // Set lower boundary of the range for the first POI
  // Return false if fail, true if success
  //
  RooRealVar * poi = GetFirstPoi();

  if (!poi) return false;

  poi->setMin(value);

  return true;
}



LimitResult LimitCalc::GetClsLimit( int nToys,
				    bool printResult ){
  //
  // Compute CLs limit
  //

  return GetClsLimit(mNPoiScanPoints, nToys, printResult);

}



LimitResult LimitCalc::GetClsLimit( int nPoiScanPoints,
				    int nToys,
				    bool printResult ){
  //
  // Compute CLs limit
  //

  std::string _legend = "[LimitCalc::GetClsLimit]: ";

  LimitResult result;

  // check necessary components

  // check workspace
  if (!mpWs){
    std::cout << _legend 
	      << "workspace not found, no model config loaded, can't continue"
	      << std::endl;
    return result;
  }
  else{
    std::cout << _legend 
	      << "workspace found..."
	      << std::endl;
  }

  // check data
  if (!mpData){
    std::cout << _legend 
	      << "dataset not loaded, can't continue"
	      << std::endl;
    return result;
  }
  else{
    std::cout << _legend 
	      << "data found..."
	      << std::endl;
  }

  // check S+B model config
  if (!mpSbModel){
    std::cout << _legend 
	      << "S+B ModelConfig not loaded, can't continue"
	      << std::endl;
    return result;
  }
  else{
    std::cout << _legend 
	      << "S+B ModelConfig found..."
	      << std::endl;
  }

  // check S+B snapshot
  const RooArgSet * _ss = mpSbModel->GetSnapshot();
  if (!_ss){
    std::cout << _legend 
	      << "found no snapshot for " 
	      << "S+B ModelConfig, can't continue" << std::endl;
    return result;
  }
  else{
    std::cout << _legend 
	      << "S+B snapshot found..."
	      << std::endl;
  }
  delete _ss;

  // check B-only model config
  if (!mpBModel){
    std::cout << _legend 
	      << "B-only ModelConfig not loaded, can't continue"
	      << std::endl;
    return result;
  }
  else{
    std::cout << _legend 
	      << "B-only ModelConfig found..."
	      << std::endl;
  }

  // check B-only snapshot
  _ss = mpBModel->GetSnapshot();
  if (!_ss){
    std::cout << _legend 
	      << "found no snapshot for " 
	      << "B-only ModelConfig, can't continue" << std::endl;
  }
  else{
    std::cout << _legend 
	      << "B-only snapshot found..."
	      << std::endl;
  }
  delete _ss;


  /*
    bool useCls,
    int nPoints,
    double poiMin,  // use default is poimin >= poimax
    double poiMax,
    int nToys );
  */

  // default values define automatic adaptive scan (max<min)
  double poi_min = 0.0;
  double poi_max = -1.0;

  RooRealVar * poi = GetFirstPoi();

  if (!poi){
    std::cout << _legend 
	      << "no parameter of interest found, can't continue " 
	      << std::endl;
  }

  if (nPoiScanPoints > 0){
    poi_min = poi->getMin();
    poi_max = poi->getMax();
  }
  else{
    std::cout << _legend 
	      << "Non-positive number of scan points requested - will attempt an adaptive scan instead" 
	      << std::endl;
  }

  // run CLs calculation
  result=ComputeInverterLimit( true,
			       nPoiScanPoints,
			       poi_min,
			       poi_max,
			       nToys );

  if (printResult){
    std::cout << " observed limit: " << result.GetObservedLimit() << std::endl;
    std::cout << " observed limit uncertainty: " << result.GetObservedLimitError() << std::endl;
    std::cout << " expected limit (median): " << result.GetExpectedLimit() << std::endl;
    std::cout << " expected limit (-1 sig): " << result.GetOneSigmaLowRange() << std::endl;
    std::cout << " expected limit (+1 sig): " << result.GetOneSigmaHighRange() << std::endl;
    std::cout << " expected limit (-2 sig): " << result.GetTwoSigmaLowRange() << std::endl;
    std::cout << " expected limit (+2 sig): " << result.GetTwoSigmaHighRange() << std::endl;
  }

  return result;
}



RooStats::ModelConfig * 
LimitCalc::LoadModelConfig( std::string mcName ){
  //
  // Return a pointer to the ModelConfig in the loaded
  // workspace, or 0 if no workspace is loaded.
  // User does NOT take ownership.
  //
  // Also load the snapshot if available
  //

  std::string _legend = "[LimitCalc::LoadModelConfig]: ";

  if (!mpWs){
    std::cout << _legend 
	      << "workspace not found, no model config loaded"
	      << std::endl;
    return 0;
  }

  RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)mpWs->obj(mcName.c_str());

  if (!_mc){
    std::cout << _legend 
	      << "ModelConfig "
	      << mcName
	      <<" not found"
	      << std::endl;
    return 0;
  }
  
  std::cout << _legend 
	    << "successfully loaded ModelConfig " 
	    << mcName << std::endl;
  
  // need to set workspace so ModelCOnfig knows where to find objects
  _mc -> SetWorkspace(*mpWs);
  
  // try to load parameter snapshot
  const RooArgSet * _ss = _mc->GetSnapshot();
  if (_ss){
    std::cout << _legend 
	      << "successfully loaded snapshot for " 
	      << "ModelConfig " << mcName << std::endl;
  }
  else{
    std::cout << _legend 
	      << "found no snapshot for " 
	      << "ModelConfig " << mcName << std::endl;
  }
  delete _ss;
  
  return _mc;
}



void
LimitCalc::SetSbModelConfig( std::string ModelConfigName ){
  // 
  // Set class member mpSbModelConfig
  // to point to a model config from the workspace
  //
  SetSbModelConfig( LoadModelConfig(ModelConfigName) );

  return;
}



void
LimitCalc::SetBModelConfig( std::string ModelConfigName ){
  // 
  // Set class member mpBModelConfig
  // to point to a model config from the workspace
  //
  SetBModelConfig( LoadModelConfig(ModelConfigName) );

  return;
}



void
LimitCalc::SetSbModelConfig( RooStats::ModelConfig * pSbModel ){
  // 
  // Set class member mpSbModelConfig - the signal plus background model
  //

  std::string _legend = "[LimitCalc::SetSbModelConfig]: ";

  mpSbModel = pSbModel;

  if (!mpSbModel){
    std::cout << _legend 
	      << "signal+background model config not set - null pointer specified"
	      << std::endl;
  }

  return;
}



void
LimitCalc::SetBModelConfig( RooStats::ModelConfig * pBModel ){
  // 
  // Set class member mpBModelConfig - the background-only model
  //

  std::string _legend = "[LimitCalc::SetBModelConfig]: ";

  mpBModel = pBModel;

  if (!mpBModel){
    std::cout << _legend 
	      << "background-only model config not set - null pointer specified"
	      << std::endl;
  }

  return;
}



RooWorkspace * 
LimitCalc::LoadWorkspace( std::string wsFileName,
			  std::string wsName ){
  //
  // load a workspace from a file
  // a copy of the workspace is kept as a class member,
  // the input file is immediately closed
  //
  
  std::string _legend = "[LimitCalc::LoadWorkspace]: ";

  TFile * p_infile = new TFile(wsFileName.c_str(), "read");
  if (p_infile->IsZombie()){
    std::cout << _legend 
	      << "no file " << wsFileName << " found!"
	      << std::endl;
    return 0;
  }
  else{
    std::cout << _legend 
	      << "file " << wsFileName << " open"
	      << std::endl;
  }

  RooWorkspace * _ws = (RooWorkspace *)p_infile->Get(wsName.c_str());
  if (!_ws){
    std::cout << _legend 
	      << "no workspace " << wsName << " found!"
	      << std::endl;
    return 0;
  }
  else{
    std::cout << _legend 
	      << "workspace " << wsName << " loaded"
	      << std::endl;
  }

  // delete current workspace or null
  delete mpWs;

  mpWs = (RooWorkspace *)_ws->Clone();

  delete p_infile;

  // try to load S+B and B-only model configs and their snapshots
  RooStats::ModelConfig * pMc = LoadModelConfig("SbModel");
  if (pMc) SetSbModelConfig( pMc );
  pMc = LoadModelConfig("BModel");
  if (pMc) SetBModelConfig( pMc );
  
  return mpWs;
}




void
LimitCalc::SetData( RooAbsData * data ){
  //
  // Copy data object and point mpData to it.
  // The class takes ownership of the new object.
  //

  std::string _legend = "[LimitCalc::SetData]: ";

  if (!data){
    std::cout << _legend 
	      << "dataset not found!"
	      << std::endl;
    return;
  }

  RooAbsData * _data = (RooAbsData *)data->Clone();

  if (!_data){
    std::cout << _legend 
	      << "failed to copy dataset!"
	      << std::endl;
    return;
  }

  // delete current dataset (or null)
  delete mpData;

  mpData = _data;

  return;
}



void
LimitCalc::SetWorkspace( RooWorkspace * ws ){
  //
  // Copy workspace object and point mpWs to it.
  // The class takes ownership of the new object.
  //

  std::string _legend = "[LimitCalc::SetWorkspace]: ";

  if (!ws){
    std::cout << _legend 
	      << "workspace not found!"
	      << std::endl;
    return;
  }

  RooWorkspace * _ws = (RooWorkspace *)ws->Clone();

  if (!_ws){
    std::cout << _legend 
	      << "failed to copy workspace!"
	      << std::endl;
    return;
  }

  // delete current workspace, if any
  delete mpWs;

  mpWs = _ws;

  return;
}



RooAbsData *
LimitCalc::LoadData( std::string datasetName ){
  //
  // Load dataset from workspace by name, copy and point mpData to it.
  // Leave mpData unchanged if new data object not found.
  // The class takes ownership of the new object.
  // Return the pointer to the new object.
  //

  std::string _legend = "[LimitCalc::LoadData]: ";

  if (!mpWs){
    std::cout << _legend 
	      << "workspace not found, no data loaded"
	      << std::endl;
    return 0;
  }

  RooAbsData * _data = (RooAbsData *)mpWs->data(datasetName.c_str());
  if (!_data){
    std::cout << _legend 
	      << "dataset " << datasetName
	      << " not found, no data loaded"
	      << std::endl;
    return 0;
  }

  // clean up old data and copy new one
  delete mpData;
  mpData = (RooAbsData *)_data->Clone();

  std::cout << _legend 
	    << "dataset " << datasetName
	    << " loaded"
	    << std::endl;

  return mpData;
}



RooStats::MCMCInterval * 
LimitCalc::GetMcmcInterval(double conf_level,
			   int n_iter,
			   int n_burn,
			   double left_side_tail_fraction,
			   int n_bins){
  // use MCMCCalculator  (takes about 1 min)
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit
  
  RooFitResult * fit = mpWs->pdf("model")->fitTo(*mpData,RooFit::Save(),
					       RooFit::Verbose(kFALSE),
					       RooFit::PrintLevel(-1),
					       RooFit::Warnings(0),
					       RooFit::PrintEvalErrors(-1));
  RooStats::ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  RooStats::ProposalFunction* pf = ph.GetProposalFunction();
  
  RooStats::MCMCCalculator mcmc( *mpData, *mpSbModel );
  mcmc.SetConfidenceLevel(conf_level);
  mcmc.SetNumIters(n_iter);          // Metropolis-Hastings algorithm iterations
  mcmc.SetProposalFunction(*pf);
  mcmc.SetNumBurnInSteps(n_burn); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(left_side_tail_fraction);
  mcmc.SetNumBins(n_bins);
  
  delete mpMcmcInterval;
  mpMcmcInterval = mcmc.GetInterval();

  return mpMcmcInterval;
}


void LimitCalc::PlotMcmcPosterior( std::string filename ){
  
  TCanvas c1("c1");
  RooStats::MCMCIntervalPlot plot(*mpMcmcInterval);
  plot.Draw();
  c1.SaveAs(filename.c_str());
  
  return;
}


double LimitCalc::PrintMcmcUpperLimit( std::string filename ){
  //
  // print out the upper limit on the first Parameter of Interest
  //

  RooRealVar * firstPOI = (RooRealVar*) mpSbModel->GetParametersOfInterest()->first();
  double _limit = mpMcmcInterval->UpperLimit(*firstPOI);
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



Double_t 
LimitCalc::GetRandom( std::string pdfName,
		      std::string varName ){
  //
  // generates a random number using a pdf in the workspace
  //
  
  // generate a dataset with one entry
  RooDataSet * _ds = mpWs->pdf(pdfName.c_str())->generate(*mpWs->var(varName.c_str()), 1);

  Double_t _result = ((RooRealVar *)(_ds->get(0)->first()))->getVal();
  delete _ds;

  return _result;
}



Long64_t
LimitCalc::FindLowBoundary(std::vector<Double_t> * pCdf,
			   Double_t value){
  //
  // return number of elements which are < value with precision 1e-10
  //

  Long64_t result = 0;
  std::vector<Double_t>::const_iterator i = pCdf->begin();
  while( (*i<value) && fabs(*i-value)>1.0e-10 && (i!=pCdf->end()) ){
    ++i;
    ++result;
  }
  return result;
}



Long64_t
LimitCalc::FindHighBoundary(std::vector<Double_t> * pCdf,
			    Double_t value){
  //
  // return number of elements which are > value with precision 1e-10
  //

  Long64_t result = 0;
  std::vector<Double_t>::const_iterator i = pCdf->end();
  while(1){ // (*i<value) && (i!=pCdf->begin()) ){
    --i;
    if (*i>value && fabs(*i-value)>1.0e-10 ){
      ++result;
    }
    else break;
    if (i==pCdf->begin()) break;
  }
  return result;
}



Double_t LimitCalc::RoundUpperBound(Double_t bound){
  //
  // find a round upper bound for a floating point
  //
  Double_t power = log10(bound);
  Int_t int_power = power>0.0 ? (Int_t)power : (Int_t)(power-1.0);
  Int_t int_bound = (Int_t)(bound/pow(10,(Double_t)int_power) * 10.0 + 1.0);
  bound = (Double_t)(int_bound/10.0*pow(10,(Double_t)int_power));
  return bound;
}



Int_t banner(){
  //#define __NOBANNER // banner temporary off
#ifndef __NOBANNER
  std::cout << desc << std::endl;
#endif
  return 0 ;
}
static Int_t dummy_ = banner() ;





/////////////////////////////////////////////////////////////////////////
//
// CLs helper methods from Lorenzo Moneta
// This is the core of the CLs calculation
//


LimitResult
LimitCalc::ComputeInverterLimit( bool useCls,
				 int nPoints,
				 double poiMin,  // use default is poimin >= poimax
				 double poiMax,
				 int nToys )
{

  std::string _legend = "[LimitCalc::ComputeInverterLimit]: ";

/*

   Other Parameter to pass in tutorial
   apart from standard for filename, ws, modelconfig and data

    type = 0 Freq calculator 
    type = 1 Hybrid 

    testStatType = 0 LEP
                 = 1 Tevatron 
                 = 2 Profile Likelihood
                 = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)

    useCLs          scan for CLs (otherwise for CLs+b)    

    npoints:        number of points to scan , for autoscan set npoints = -1 

    poimin,poimax:  min/max value to scan in case of fixed scans 
                    (if min >= max, try to find automatically)                           

    ntoys:         number of toys to use 

    extra options are available as global paramters of the macro. They are: 

    mPlotHypoTestResult   plot result of tests at each point (TS distributions) 
    mUseProof = true;
    mWriteResult = true;
    nworkers = 4;


   */

// FIXME:
  std::string suffix = "_test";

  // result
  //std::vector<Double_t> result;
  LimitResult result;

  // check that workspace is present
  if (!mpWs){
    std::cout << "No workspace found, null pointer" << std::endl;
    return result;
  }
  
  RooStats::HypoTestInverterResult * r = 0;
  RooStats::HypoTestInverterResult * r2 = 0;
  
  // FIXME: terrible hack to check appending results
  if (suffix.find("merge")!=std::string::npos){
    std::string resFile = "Freq_CLs_grid_ts2_test_1.root";
    std::string resFile2 = "Freq_CLs_grid_ts2_test_2.root";
    std::string resName = "result_xsec";
    //std::cout << "Reading an HypoTestInverterResult with name " << resName << " from file " << resFile << std::endl;
    TFile * file = new TFile(resFile.c_str(), "read");
    TFile * file2 = new TFile(resFile2.c_str(), "read");
    r = dynamic_cast<RooStats::HypoTestInverterResult*>( file->Get(resName.c_str()) ); 
    r2 = dynamic_cast<RooStats::HypoTestInverterResult*>( file2->Get(resName.c_str()) ); 
    r->Add(*r2);
  }
  else{
    r = RunInverter( nPoints, poiMin, poiMax,  nToys, nToys, useCls );    
    if (!r) { 
      std::cerr << "Error running the HypoTestInverter - Exit " << std::endl;
      return result;
    }
  }
      		

   double upperLimit = r->UpperLimit();
   double ulError = r->UpperLimitEstimatedError();
   result.SetObservedLimit(upperLimit);
   result.SetObservedLimitError(ulError);

   const char *  limitType = (useCls) ? "CLs" : "Cls+b";
   const char * scanType = (nPoints < 0) ? "auto" : "grid";
   const char *  typeName = (mInverterCalcType == 0) ? "Frequentist" : "Hybrid";
   const char * resultName = (mpWs) ? mpWs->GetName() : r->GetName();
   TString plotTitle = TString::Format("%s CL Scan for workspace %s",typeName,resultName);

   /*
     RooStats::HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);
   TCanvas c1;
   //plot->Draw("CLb 2CL");  // plot all and Clb
   plot->Draw("2CL");  // plot all and Clb
   TString resultFileName = TString::Format("%s_%s_ts%d_scan_",limitType,scanType,mTestStatType);
   resultFileName += suffix;
   resultFileName += ".pdf";
   c1.SaveAs(resultFileName);

   if (mPlotHypoTestResult) { 
      TCanvas * c2 = new TCanvas();
      c2->Divide( 2, TMath::Ceil(nEntries/2));
      for (int i=0; i<nEntries; i++) {
         c2->cd(i+1);
         SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
         pl->SetLogYaxis(true);
         pl->Draw();
      }
   }
   */

   Double_t q[5];
   q[0] = r->GetExpectedUpperLimit(0);
   q[1] = r->GetExpectedUpperLimit(-1);
   q[2] = r->GetExpectedUpperLimit(1);
   q[3] = r->GetExpectedUpperLimit(-2);
   q[4] = r->GetExpectedUpperLimit(2);
   result.SetExpectedLimit(q[0]);
   result.SetOneSigmaLowRange(q[1]);
   result.SetOneSigmaHighRange(q[2]);
   result.SetTwoSigmaLowRange(q[3]);
   result.SetTwoSigmaHighRange(q[4]);


   if (mpWs != NULL && mWriteResult) {

      // write to a file the results
      const char *  calcType = (mInverterCalcType == 0) ? "Freq" : "Hybr";
      //const char *  limitType = (useCls) ? "CLs" : "Cls+b";
      //const char * scanType = (nPoints < 0) ? "auto" : "grid";
      TString resultFileName = TString::Format("%s_%s_%s_ts%d_",calcType,limitType,scanType,mTestStatType);      
      //resultFileName += fileName;
      resultFileName += suffix;
      resultFileName += ".root";

      TFile * fileOut = new TFile(resultFileName,"RECREATE");
      r->Write();
      fileOut->Close();                                                                     
   }   

   // FIXME: get adaptive result out in some way
   if (nPoints <= 0 || poiMin >= poiMax){
     result.Clear();
     std::cout << _legend 
	       << "can't return the adaptive scan result yet, returning empty for now"
	       << std::endl;
   }

   return result;
}



Double_t
LimitCalc::GetExpectedPValue( RooStats::HypoTestInverterResult * pResult,
			      Int_t index,
			      Double_t nSigma ){
  //
  // Find a p-value in the expected p-value sampling distribution
  // that corresponds to a specified deviation from median,
  // for a POI value specified by index
  // Sampling distributions are taken from the provided result object
  //
  
  double p[1];
  double q[1];
  p[0] = ROOT::Math::normal_cdf(nSigma,1);
  RooStats::SamplingDistribution * s = pResult->GetExpectedPValueDist(index);
  const std::vector<double> & values = s->GetSamplingDistribution();
  double * x = const_cast<double *>(&values[0]); // cast for TMath::Quantiles
  TMath::Quantiles(values.size(), 1, x, q, p, false);
  delete s;
  
  return q[0];
}



Double_t
LimitCalc::GetExpectedPValueErr( RooStats::HypoTestInverterResult * pResult,
				 Int_t index, 
				 Double_t pValue ){
  //
  // Find the uncertainty for a specified p-value from the expected distribution.
  // index - index of a scan point in *pResult
  // Sampling distributions are taken from the provided result object
  //
  // try to estimate the corresponding p-value uncertainty

  RooStats::SamplingDistribution * bDistribution = pResult->GetBackgroundTestStatDist(index);
  RooStats::SamplingDistribution * sbDistribution = pResult->GetSignalAndBackgroundTestStatDist(index);
  if (!bDistribution || !sbDistribution) return -1.0;

  RooStats::HypoTestResult * result = pResult->GetResult(index);

  // create a new HypoTestResult
  RooStats::HypoTestResult tempResult; 
  tempResult.SetPValueIsRightTail( result->GetPValueIsRightTail() );
  tempResult.SetBackgroundAsAlt( true );
  tempResult.SetNullDistribution( sbDistribution );
  tempResult.SetAltDistribution( bDistribution );

  std::vector<double> p_values(bDistribution->GetSize()); 
  std::vector<double> p_value_errors(bDistribution->GetSize()); 
  for (int i = 0; i < bDistribution->GetSize(); ++i) { 
    tempResult.SetTestStatisticData( bDistribution->GetSamplingDistribution()[i] );
    p_values[i] = tempResult.CLs();
    p_value_errors[i] = tempResult.CLsError();
  }

  // get index array of sorted p_values vector
  std::vector<int> v_sorted_index( p_values.size() );
  TMath::Sort((int)(v_sorted_index.size()),
	      &((const std::vector<double>)p_values)[0],
	      &v_sorted_index[0]);

  double p_value_error = 1.0;
  bool first_entry = true;
  int previous_index = -1;
  for (std::vector<int>::const_iterator i = v_sorted_index.begin(); 
       i != v_sorted_index.end();
       ++i){

    // test
    //std::cout << *i << ": p = "
    //	      << p_values[*i] << ", error = "
    //	      << p_value_errors[*i] << std::endl;

    // find two p-values that contain pValue between them,
    // take as error their average errors
    if (p_values[*i]<pValue){

      if (first_entry){
	// special case when pValue > all p-values
	p_value_error = p_value_errors[0];
      }
      else{
	p_value_error = 0.5*(p_value_errors[*i]+p_value_errors[previous_index]);
      }

      // found error, no need to loop further
      break;
    }
    
    previous_index = *i;
    first_entry = false;
  }

  return p_value_error;
}



std::pair<double,double>
LimitCalc::GetClsSinglePoint( RooStats::HypoTestInverter & calc,
			      double poi,
			      double precision,
			      bool observed,
			      double bgSigma,
			      int maxNToys ){
  //
  // Run single HypoTestInverter calculationt 
  // with adaptive number of toys for precision
  // and speed
  //
  //   precision - realtive (fraction)

  std::pair<double,double> result;
  result.first = -1.0;
  result.second = -1.0;

  //int max_ntoys_zero = 1000;
  int max_ntoys_zero = mMaxZeroToys;

  // result and its current precision
  double _clsb = -1.0;
  double _clsb_err = -1.0;
  double _clb = -1.0;
  double _clb_err = -1.0;
  double _cls = -1.0;
  double _cls_err = -1.0;

  // initial number of toys
  double p_value_bg = ROOT::Math::normal_cdf(-bgSigma,1);
  int add_toys_sb = (int)(1.5/p_value_bg/p_value_bg + 1.0);
  int add_toys_b = add_toys_sb;

  // keep adding toys until precision is reached or something
  // is hopelessly zero, and max_toys_zero are reached, or
  // maxNToys is reached
  if (mTestMode) std::cout << "[DEBUG]: starting loop with adding toys until precision is reached" << std::endl;
  while(1){
    
    if (mTestMode){
      std::cout << "[DEBUG]: add S+B toys: " << add_toys_sb << std::endl;
      std::cout << "[DEBUG]: add B toys: " << add_toys_b << std::endl;
    }

    // check if too many toys requested
    if ( add_toys_sb > maxNToys ) add_toys_sb = maxNToys;
    if ( add_toys_b > maxNToys )  add_toys_b = maxNToys;

    if (mInverterCalcType == 0){
      ((RooStats::FrequentistCalculator *)(calc.GetHypoTestCalculator()))->SetToys(add_toys_sb,add_toys_b);
    }
    else{
      ((RooStats::HybridCalculator *)(calc.GetHypoTestCalculator()))->SetToys(add_toys_sb,add_toys_b);
    }

    // reset
    add_toys_sb = 0;
    add_toys_b = 0;

    // run hypo test inverter
    calc.RunOnePoint(poi);
    RooStats::HypoTestInverterResult * pResult = calc.GetInterval();
    
    // get the index of the current point
    int poi_index = pResult->FindIndex(poi);
    if (poi_index < 0) poi_index = pResult->ArraySize()-1;

    if (mTestMode) std::cout << "[DEBUG]: ArraySize(): " << pResult->ArraySize() << std::endl;
    
    // sampling distributions
    // (we do not own these)
    RooStats::SamplingDistribution * pBDist = pResult->GetBackgroundTestStatDist(poi_index);
    RooStats::SamplingDistribution * pSbDist = pResult->GetSignalAndBackgroundTestStatDist(poi_index);
    if (!pBDist || !pSbDist) return result;
    
    // 
    RooStats::HypoTestResult * result = pResult->GetResult(poi_index);
    
    // create a new HypoTestResult
    RooStats::HypoTestResult temp_result; 
    temp_result.SetPValueIsRightTail( result->GetPValueIsRightTail() );
    temp_result.SetBackgroundAsAlt( true );
    temp_result.SetNullDistribution( pSbDist );
    temp_result.SetAltDistribution( pBDist );
    
    // get test statistic value for CLs calculation
    double test_stat;
    double test_stat_rms = 0.0;
    double test_stat_err;
    if (observed){
      test_stat = result->GetTestStatisticData();
    }
    else{
      // expected limits requested
      // get test statistic from quantile of the b-only sampling dist
      double p[1];
      double q[1];
      p[0] = p_value_bg;
      const std::vector<double> & values = pBDist->GetSamplingDistribution();
      double * x = const_cast<double *>(&values[0]); // cast for TMath::Quantiles
      TMath::Quantiles(values.size(), 1, x, q, p, false);
      test_stat = q[0];
      test_stat_rms = TMath::RMS(values.size(), &values[0]);
    }
    
    std::vector<double> vBSorted(pBDist->GetSamplingDistribution());
    std::vector<double> vSbSorted(pSbDist->GetSamplingDistribution());
    int n_b_toys = vBSorted.size();
    int n_sb_toys = vSbSorted.size();
    //int n_b_toys = pBDist->GetSamplingDistribution().size();
    //int n_sb_toys = pSbDist->GetSamplingDistribution().size();

    // evaluate tails and their errors
    if (observed) test_stat_err = 0.0; // test stat known from data
    else test_stat_err = test_stat_rms / sqrt((double)n_b_toys);

    // sort sampling dist vectors
    std::sort(vBSorted.begin(), vBSorted.end());
    std::sort(vSbSorted.begin(), vSbSorted.end());

    // find index of sampling dist element <= test_stat
    // and of +/- uncertainty deviations from test_stat
    //int indexB = TMath::BinarySearch(n_b_toys, &vBSorted[0], test_stat);
    int indexSb = TMath::BinarySearch(n_sb_toys, &vSbSorted[0], test_stat);
    int indexSbDown = TMath::BinarySearch(n_sb_toys, &vSbSorted[0], test_stat-test_stat_err);
    int indexSbUp = TMath::BinarySearch(n_sb_toys, &vSbSorted[0], test_stat+test_stat_err);

    //double _alt_clb = ((double)indexB+1.0)/((double)n_b_toys);
    double _alt_clb = 1.0 - ROOT::Math::normal_cdf(-bgSigma,1);
    double _alt_clsb = 1.0 - ((double)indexSb+1.0)/((double)n_sb_toys);
    double _alt_clsb_up = 1.0 - ((double)indexSbDown+1.0)/((double)n_sb_toys);
    double _alt_clsb_down = 1.0 - ((double)indexSbUp+1.0)/((double)n_sb_toys);
    //double _alt_clb_err = std::max(_alt_clb, 1.0);
    double _alt_clb_err = 0.0;
    double _alt_clsb_err1 = std::max(_alt_clsb, 1.0); // two components: from Poisson
    double _alt_clsb_err2; // and from test stat value fluctuation
    double _alt_clsb_err = std::max(_alt_clsb, 1.0);
    //double _alt_clsb_err2 = std::max(_alt_clsb, 1.0);
    double _alt_cls = 0.0;
    double _alt_cls_err = 1.0;
    //if ( indexB+1 > 0) _alt_clb_err = _alt_clb/sqrt((double)(indexB+1));
    //if ( n_sb_toys-indexSb-1 > 0) _alt_clsb_err = _alt_clsb/sqrt((double)(n_sb_toys-indexSb-1));

    if ( n_sb_toys-indexSb-1 > 0){
      _alt_clsb_err1 = _alt_clsb/sqrt((double)(n_sb_toys-indexSb-1));
    }
    //_alt_clsb_err2 = _alt_clsb_up - _alt_clsb;
    _alt_clsb_err2 = std::max( fabs(_alt_clsb_up-_alt_clsb), fabs(_alt_clsb-_alt_clsb_down) );
    _alt_clsb_err = _alt_clsb_err1 + _alt_clsb_err2;

    if (_alt_clb > 0){
      _alt_cls = _alt_clsb/_alt_clb;
      _alt_cls_err = _alt_cls * sqrt(_alt_clsb_err*_alt_clsb_err/_alt_clsb/_alt_clsb + _alt_clb_err*_alt_clb_err/_alt_clb/_alt_clb);
    }

    temp_result.SetTestStatisticData( test_stat );
    _clsb = temp_result.CLsplusb();
    _clsb_err = temp_result.CLsplusbError();
    _clb = temp_result.CLb();
    _clb_err = temp_result.CLbError();

    if (mTestMode){
      std::cout << "[DEBUG]: N of B toys: " << n_b_toys  << std::endl;
      std::cout << "[DEBUG]: N of S+B toys: " << n_sb_toys  << std::endl;
      std::cout << "[DEBUG]: CLsb: " << _clsb << " +/- " << _clsb_err  << std::endl;
      std::cout << "[DEBUG]: CLb: " << _clb << " +/- " << _clb_err  << std::endl;
      //std::cout << "[DEBUG]: alt CLsb: " << _alt_clsb << " +/- " << _alt_clsb_err  << std::endl;
      std::cout << "[DEBUG]: test stat: " << test_stat << " +/- " << test_stat_err  << std::endl;
      std::cout << "[DEBUG]: alt CLsb: " << _alt_clsb << " +/- " << _alt_clsb_err1 << " +/- " << _alt_clsb_err2  << std::endl;
      std::cout << "[DEBUG]: alt CLsb: " << _alt_clsb << " +/- " << _alt_clsb_err  << std::endl;
      //std::cout << "[DEBUG]: alt CLsb_up: " << _alt_clsb_up << std::endl;
      std::cout << "[DEBUG]: alt CLb: " << _alt_clb << " +/- " << _alt_clb_err  << std::endl;
      std::cout << "[DEBUG]: alt CLs: " << _alt_cls << " +/- " << _alt_cls_err  << std::endl;
    }

    if (_alt_clb < 0.1/(double)n_b_toys){
      // BG p-value iz zero
      _alt_cls = 1.0;
      _alt_cls_err = -1.0;

      if (mTestMode){
	std::cout << "[DEBUG]: CLb is zero..." << std::endl;
      }

      if (n_b_toys >= max_ntoys_zero){
	// too many toys and still zero
	break;
      }

      add_toys_b = (int)(((double)n_b_toys)/precision/precision);
      if ( (add_toys_b + n_b_toys) > max_ntoys_zero ){
	// too many toys
	add_toys_b = max_ntoys_zero - n_b_toys;
      }

      add_toys_sb = 1;
      continue;
    }

    if (_alt_clsb < 0.1/(double)n_sb_toys){
      // S+B p-value is zero
      _alt_cls = 0.0;
      _alt_cls_err = -1.0;

      if (mTestMode){
	std::cout << "[DEBUG]: CLsb is zero..." << std::endl;
      }

      if (n_sb_toys >= max_ntoys_zero){
	// too many toys and still zero
	break;
      }

      add_toys_b = 1;
      add_toys_sb = (int)(((double)n_sb_toys)/precision/precision);

      if ( (add_toys_sb + n_sb_toys) > max_ntoys_zero ){
	// too many toys
	add_toys_sb = max_ntoys_zero - n_sb_toys;
      }
      
      continue;
    }

    if ( _alt_clsb_err1/_alt_clsb > precision ){
      add_toys_sb = (int)(1.0
			  + ((double)n_sb_toys)
			  * _alt_clsb_err1*_alt_clsb_err1/_alt_clsb/_alt_clsb/precision/precision
			  - (double)n_sb_toys);

      if ( (double)add_toys_sb < 0.1*((double)n_sb_toys)){
	// protection against too small increments
	add_toys_sb = (int)(0.1*((double)n_sb_toys)+1.0);
      }
    }

    // how well we know test statistic cut depends on b-only
    if ( _alt_clsb_err2/_alt_clsb > precision ){
      add_toys_b = (int)(1.0
			 + ((double)n_b_toys)
			 * _alt_clsb_err2*_alt_clsb_err2/_alt_clsb/_alt_clsb/precision/precision
			 - (double)n_b_toys);

      if ( (double)add_toys_b < 0.1*((double)n_b_toys)){
	// protection against too small increments
	add_toys_b = (int)(0.1*((double)n_b_toys)+1.0);
      }
    }

    // this is not needed for expected limits
    if ( observed && _alt_clb_err/_alt_clb > precision ){
      add_toys_b = (int)(1.0
			 + ((double)n_b_toys)
			 * _alt_clb_err*_alt_clb_err/_alt_clb/_alt_clb/precision/precision
			 - (double)n_b_toys);

      if ( (double)add_toys_b < 0.1*((double)n_b_toys)){
	// protection against too small increments
	add_toys_b = (int)(0.1*((double)n_b_toys)+1.0);
      }
    }

    // FIXME: alternative version of cls calculation
    //_cls = temp_result.CLs();
    //_cls_err = temp_result.CLsError();
    _cls = _alt_cls;
    _cls_err = _alt_cls_err;

    if (mTestMode){
      std::cout << "[DEBUG]: CLs: " << _cls << " +/- " << _cls_err  << std::endl;
      std::cout << "[DEBUG]: alt CLs: " << _alt_cls << " +/- " << _alt_cls_err  << std::endl;
    }

    // precision achieved or far away from goal
    if ( _alt_cls_err/_alt_cls < precision ) break;

    // how much over the desired precision?
    double extra = _alt_cls_err*_alt_cls_err/_alt_cls/_alt_cls/precision/precision - 1.0;
    
    // more toys only if both cls and clb reached the precision
    if (_alt_clb_err/_alt_clb < precision && _alt_clsb_err/_alt_clsb < precision){
      add_toys_b = (int)(extra*((double)n_b_toys) + 1.0);
      add_toys_sb = (int)(extra*((double)n_sb_toys) + 1.0);
    }

    // protection against too small increments
    if ( (double)add_toys_sb < 0.1*((double)n_sb_toys)){
      add_toys_sb = (int)(0.1*((double)n_sb_toys)+1.0);
    }
    if ( (double)add_toys_b < 0.1*((double)n_b_toys)){
      add_toys_b = (int)(0.1*((double)n_b_toys)+1.0);
    }

    // max toys reached, exit
    if (n_sb_toys >= maxNToys || n_b_toys >= maxNToys){
      break;
    }

    // check if too many toys requested
    if ( (add_toys_sb + n_sb_toys) > maxNToys ){
      add_toys_sb = maxNToys - n_sb_toys;
    }
    if ( (add_toys_b + n_b_toys) > maxNToys ){
      add_toys_b = maxNToys - n_b_toys;
    }

  }
  
  //std::pair<double,double> result;
  result.first = _cls;
  result.second = _cls_err;

  return result;
}



RooStats::HypoTestInverter * 
LimitCalc::GetHypoTestInverter( void ){
  // 
  // Configure, create and return a pointer to a HypoTestInverter object
  // The object is owned by the class
  // 

  std::string _legend = "[LimitCalc::GetHypoTestInverter]: ";

  // check data
   if (!mpData) { 
     Error(_legend.c_str(),"dataset not found");
     return 0;
   }
   
   // check model config
   if (!mpSbModel) {
     Error(_legend.c_str(),"S+B ModelConfig does not exist...");
     return 0;
   }

   // check model pdf
   if (!mpSbModel->GetPdf()) { 
     Error(_legend.c_str(),"S+B model has no pdf...");
     return 0;
   }

   // check POI
   if (!mpSbModel->GetParametersOfInterest()) {
     Error(_legend.c_str(),"S+B model has no POI...");
     return 0;
   }

   // check S+B snapshot
   if (!mpSbModel->GetSnapshot() ) { 
      Info(_legend.c_str(),"S+B model has no snapshot  - make one using model POI");
      mpSbModel->SetSnapshot( *mpSbModel->GetParametersOfInterest() );
   }

   // check B-only model config
   if (!mpBModel || mpBModel == mpSbModel) {
      Info(_legend.c_str(),"The background model does not exist...");
      Info(_legend.c_str(),"Copy it from S+B ModelConfig and set POI to zero");
      mpBModel = (RooStats::ModelConfig*) mpSbModel->Clone();
      mpBModel->SetName(TString("SbModel")+TString("_with_poi_0"));      
      RooRealVar * var = dynamic_cast<RooRealVar*>(mpBModel->GetParametersOfInterest()->first());
      if (!var) return 0;
      double oldval = var->getVal();
      var->setVal(0);
      mpBModel->SetSnapshot( RooArgSet(*var)  );
      var->setVal(oldval);
   }
   else { 
     // check B-only snapshot
     if (!mpBModel->GetSnapshot() ) { 
       Info(_legend.c_str(),"B-only model has no snapshot  - make one using model poi and 0 values ");
       RooRealVar * var = dynamic_cast<RooRealVar*>(mpBModel->GetParametersOfInterest()->first());
       if (var) { 
	 double oldval = var->getVal();
	 var->setVal(0);
	 mpBModel->SetSnapshot( RooArgSet(*var)  );
	 var->setVal(oldval);
       }
       else { 
	 Error(_legend.c_str(),"B-only model has no valid poi...");
	 return 0;
       }         
     }
   }

   delete slrts;
   slrts = new RooStats::SimpleLikelihoodRatioTestStat(*mpSbModel->GetPdf(),*mpBModel->GetPdf());
   if (mpSbModel->GetSnapshot()) slrts->SetNullParameters(*mpSbModel->GetSnapshot());
   if (mpBModel->GetSnapshot()) slrts->SetAltParameters(*mpBModel->GetSnapshot());

   // ratio of profile likelihood - need to pass snapshot for the alt
   delete ropl;
   ropl = new RooStats::RatioOfProfiledLikelihoodsTestStat(*mpSbModel->GetPdf(), *mpBModel->GetPdf(), mpBModel->GetSnapshot());
   ropl->SetSubtractMLE(false);
   
   delete profll;
   profll = new RooStats::ProfileLikelihoodTestStat(*mpSbModel->GetPdf());
   if (mTestStatType == 3) profll->SetOneSided(1);
   if (mOptimize){ 
      profll->SetReuseNLL(true);
      slrts->setReuseNLL(true);
   }

   RooRealVar * pMu = dynamic_cast<RooRealVar*>(mpSbModel->GetParametersOfInterest()->first());
   assert(pMu != 0);
   delete maxll;
   maxll = new RooStats::MaxLikelihoodEstimateTestStat(*mpSbModel->GetPdf(),*pMu); 

   RooStats::TestStatistic * testStat = slrts;
   if (mTestStatType == 1) testStat = ropl;
   if (mTestStatType == 2 || mTestStatType == 3) testStat = profll;
   if (mTestStatType == 4) testStat = maxll;
  
   
   delete mpHypoTestCalc;
   if (mInverterCalcType == 0) mpHypoTestCalc = new RooStats::FrequentistCalculator(*mpData, *mpBModel, *mpSbModel);
   else mpHypoTestCalc = new RooStats::HybridCalculator(*mpData, *mpBModel, *mpSbModel);

   RooStats::ToyMCSampler * toymcs = (RooStats::ToyMCSampler*)mpHypoTestCalc->GetTestStatSampler();
   if(mpSbModel->GetPdf()->canBeExtended()){
     // if the PDF is extended, number of events per toy
     // will be taken from the PDF normalization
     std::cout << _legend
	       << "will take number of entries per toy from S+B model PDF"
	       << std::endl;
   }
   else{
     // PDF is not extended

     if(mNEventsPerToy > 0){
       // number of event per toy was specified
       toymcs->SetNEventsPerToy(mNEventsPerToy);
       std::cout << _legend
		 << "number of entries per toy explicitely set to "
		 << mNEventsPerToy
		 << std::endl;
     }
     else if(mpData->numEntries()==1){
       toymcs->SetNEventsPerToy(1);
       std::cout << _legend
		 << "guessing that this must be a counting experiment, "
		 << std::endl
		 << _legend
		 << "number of entries per toy explicitely set to 1"
		 << std::endl
		 << _legend
		 << "if you wish otherwise, set it with SetNEventsPerToy()"
		 << std::endl;
     }
     else{
       std::cout << _legend
		 << "number of entries per toy is undefined"
		 << std::endl;
     }
   }

   toymcs->SetTestStatistic(testStat);
   if (mOptimize){
     // Lorenzo: works only of b pdf and s+b pdf are the same
     // (perhaps this is not the case anymore)
     if (mpBModel->GetPdf() == mpSbModel->GetPdf() ) 
       toymcs->SetUseMultiGen(true);
   }


   if (mInverterCalcType == 1) { 
     RooStats::HybridCalculator *hhc = (RooStats::HybridCalculator*) mpHypoTestCalc;

      // remove global observables from ModelConfig
      mpBModel->SetGlobalObservables( RooArgSet() );
      mpSbModel->SetGlobalObservables( RooArgSet() );

      // check for nuisance prior pdf in case of nuisance parameters 
      if (mpBModel->GetNuisanceParameters() || mpSbModel->GetNuisanceParameters() ){
	RooAbsPdf * pNuisPdf = 0; 
	if (mNuisPriorName.length()!=0) pNuisPdf = mpWs->pdf(mNuisPriorName.c_str());
	// use prior defined first in bModel (then in SbModel)
	if (!pNuisPdf)  { 
            Info(_legend.c_str(),"No nuisance pdf given for the HybridCalculator - try to use the prior pdf from the model");
            pNuisPdf = (mpBModel->GetPriorPdf() ) ?  mpBModel->GetPriorPdf() : mpSbModel->GetPriorPdf();
         }
         if (!pNuisPdf) { 
            Error(_legend.c_str(),"Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified");
            return 0;
         }
         const RooArgSet * cpNuisParams = (mpBModel->GetNuisanceParameters() ) ? mpBModel->GetNuisanceParameters() : mpSbModel->GetNuisanceParameters();
         RooArgSet * pNp = pNuisPdf->getObservables(*cpNuisParams);
         if (pNp->getSize() == 0) { 
            Warning(_legend.c_str(),"Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
         }
         delete pNp;
         hhc->ForcePriorNuisanceAlt(*pNuisPdf);
         hhc->ForcePriorNuisanceNull(*pNuisPdf);
      }
   } 
   else {
     // just small default ntoys numbers
     ((RooStats::FrequentistCalculator*) mpHypoTestCalc)->SetToys(100,100); 
   }

   // Get the result
   // FIXME: verbosity?
   //RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

   // fit the data first
   mpSbModel->GetPdf()->fitTo(*mpData, 
			      RooFit::Verbose(0),
			      RooFit::PrintLevel(-1),
			      RooFit::Warnings(0),
			      RooFit::PrintEvalErrors(-1));

   delete mpHypoTestInvCalc;
   mpHypoTestInvCalc = new RooStats::HypoTestInverter(*mpHypoTestCalc);

   mpHypoTestInvCalc->SetConfidenceLevel(0.95);
   mpHypoTestInvCalc->UseCLs(true);
   mpHypoTestInvCalc->SetVerbose(true);

   return mpHypoTestInvCalc;
}


RooStats::HypoTestInverterResult * 
LimitCalc::RunInverter( int    npoints,
			double poimin,
			double poimax, 
			int    nSbToys,
			int    nBToys,
			bool   useCls ){
  // internal routine to run the inverter
  // uses dataset from the class mpData

  std::string _legend = "[LimitCalc::RunInverter]: ";

  // check data
   if (!mpData) { 
     Error("LimitCalc::RunInverter","dataset not found");
     return 0;
   }
   
   // check model config
   if (!mpSbModel) {
     Error("LimitCalc::RunInverter","S+B ModelConfig does not exist...");
     return 0;
   }

   // check model pdf
   if (!mpSbModel->GetPdf()) { 
     Error("LimitCalc::RunInverter","S+B model has no pdf...");
     return 0;
   }

   // check POI
   if (!mpSbModel->GetParametersOfInterest()) {
     Error("LimitCalc::RunInverter","S+B model has no POI...");
     return 0;
   }

   // check S+B snapshot
   if (!mpSbModel->GetSnapshot() ) { 
      Info("LimitCalc::RunInverter","S+B model has no snapshot  - make one using model POI");
      mpSbModel->SetSnapshot( *mpSbModel->GetParametersOfInterest() );
   }

   // check B-only model config
   if (!mpBModel || mpBModel == mpSbModel) {
      Info("LimitCalc::RunInverter","The background model does not exist...");
      Info("LimitCalc::RunInverter","Copy it from S+B ModelConfig and set POI to zero");
      mpBModel = (RooStats::ModelConfig*) mpSbModel->Clone();
      mpBModel->SetName(TString("SbModel")+TString("_with_poi_0"));      
      RooRealVar * var = dynamic_cast<RooRealVar*>(mpBModel->GetParametersOfInterest()->first());
      if (!var) return 0;
      double oldval = var->getVal();
      var->setVal(0);
      mpBModel->SetSnapshot( RooArgSet(*var)  );
      var->setVal(oldval);
   }
   else { 
     // check B-only snapshot
     if (!mpBModel->GetSnapshot() ) { 
       Info("LimitCalc::RunInverter","B-only model has no snapshot  - make one using model poi and 0 values ");
       RooRealVar * var = dynamic_cast<RooRealVar*>(mpBModel->GetParametersOfInterest()->first());
       if (var) { 
	 double oldval = var->getVal();
	 var->setVal(0);
	 mpBModel->SetSnapshot( RooArgSet(*var)  );
	 var->setVal(oldval);
       }
       else { 
	 Error("LimitCalc::RunInverter","B-only model has no valid poi...");
	 return 0;
       }         
     }
   }


   RooStats::SimpleLikelihoodRatioTestStat slrts(*mpSbModel->GetPdf(),*mpBModel->GetPdf());
   if (mpSbModel->GetSnapshot()) slrts.SetNullParameters(*mpSbModel->GetSnapshot());
   if (mpBModel->GetSnapshot()) slrts.SetAltParameters(*mpBModel->GetSnapshot());

   // ratio of profile likelihood - need to pass snapshot for the alt
   RooStats::RatioOfProfiledLikelihoodsTestStat 
     ropl(*mpSbModel->GetPdf(), *mpBModel->GetPdf(), mpBModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   
   RooStats::ProfileLikelihoodTestStat profll(*mpSbModel->GetPdf());
   if (mTestStatType == 3) profll.SetOneSided(1);
   if (mOptimize){ 
      profll.SetReuseNLL(true);
      slrts.setReuseNLL(true);
   }

   RooRealVar * pMu = dynamic_cast<RooRealVar*>(mpSbModel->GetParametersOfInterest()->first());
   assert(pMu != 0);
   RooStats::MaxLikelihoodEstimateTestStat maxll(*mpSbModel->GetPdf(),*pMu); 

   RooStats::TestStatistic * testStat = &slrts;
   if (mTestStatType == 1) testStat = &ropl;
   if (mTestStatType == 2 || mTestStatType == 3) testStat = &profll;
   if (mTestStatType == 4) testStat = &maxll;
  
   
   RooStats::HypoTestCalculatorGeneric *  hc = 0;
   if (mInverterCalcType == 0) hc = new RooStats::FrequentistCalculator(*mpData, *mpBModel, *mpSbModel);
   else hc = new RooStats::HybridCalculator(*mpData, *mpBModel, *mpSbModel);

   RooStats::ToyMCSampler * toymcs = (RooStats::ToyMCSampler*)hc->GetTestStatSampler();
   if(mpSbModel->GetPdf()->canBeExtended()){
     // if the PDF is extended, number of events per toy
     // will be taken from the PDF normalization
     std::cout << _legend
	       << "will take number of entries per toy from S+B model PDF"
	       << std::endl;
   }
   else{
     // PDF is not extended

     if(mNEventsPerToy > 0){
       // number of event per toy was specified
       toymcs->SetNEventsPerToy(mNEventsPerToy);
       std::cout << _legend
		 << "number of entries per toy explicitely set to "
		 << mNEventsPerToy
		 << std::endl;
     }
     else if(mpData->numEntries()==1){
       toymcs->SetNEventsPerToy(1);
       std::cout << _legend
		 << "guessing that this must be a counting experiment, "
		 << std::endl
		 << _legend
		 << "number of entries per toy explicitely set to 1"
		 << std::endl
		 << _legend
		 << "if you wish otherwise, set it with SetNEventsPerToy()"
		 << std::endl;
     }
     else{
       std::cout << _legend
		 << "number of entries per toy is undefined"
		 << std::endl;
     }
   }

   toymcs->SetTestStatistic(testStat);
   if (mOptimize){
     // Lorenzo: works only of b pdf and s+b pdf are the same
     if (mpBModel->GetPdf() == mpSbModel->GetPdf() ) 
       toymcs->SetUseMultiGen(true);
   }


   if (mInverterCalcType == 1) { 
     RooStats::HybridCalculator *hhc = (RooStats::HybridCalculator*) hc;
      hhc->SetToys(nSbToys,nBToys); 
      // FIXME: trying adaptive sampling
      hhc->SetNToysInTails(0,0); 

      // remove global observables from ModelConfig
      mpBModel->SetGlobalObservables( RooArgSet() );
      mpSbModel->SetGlobalObservables( RooArgSet() );

      // check for nuisance prior pdf in case of nuisance parameters 
      if (mpBModel->GetNuisanceParameters() || mpSbModel->GetNuisanceParameters() ){
	RooAbsPdf * pNuisPdf = 0; 
	if (mNuisPriorName.length()!=0) pNuisPdf = mpWs->pdf(mNuisPriorName.c_str());
	// use prior defined first in bModel (then in SbModel)
	if (!pNuisPdf)  { 
            Info("StandardHypoTestInvDemo","No nuisance pdf given for the HybridCalculator - try to use the prior pdf from the model");
            pNuisPdf = (mpBModel->GetPriorPdf() ) ?  mpBModel->GetPriorPdf() : mpSbModel->GetPriorPdf();
         }
         if (!pNuisPdf) { 
            Error("StandardHypoTestInvDemo","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified");
            return 0;
         }
         const RooArgSet * cpNuisParams = (mpBModel->GetNuisanceParameters() ) ? mpBModel->GetNuisanceParameters() : mpSbModel->GetNuisanceParameters();
         RooArgSet * pNp = pNuisPdf->getObservables(*cpNuisParams);
         if (pNp->getSize() == 0) { 
            Warning("StandardHypoTestInvDemo","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
         }
         delete pNp;
         hhc->ForcePriorNuisanceAlt(*pNuisPdf);
         hhc->ForcePriorNuisanceNull(*pNuisPdf);
      }
   } 
   else {
     ((RooStats::FrequentistCalculator*) hc)->SetToys(nSbToys,nBToys); 
     
     // FIXME: trying out adaptive sampling
     ((RooStats::FrequentistCalculator*) hc)->SetNToysInTails(0,0); 
   }

   // Get the result
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);


   TStopwatch tw; tw.Start(); 
   const RooArgSet * poiSet = mpSbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar*)poiSet->first();

   // fit the data first

   mpSbModel->GetPdf()->fitTo(*mpData, 
			      RooFit::Verbose(0),
			      RooFit::PrintLevel(-1),
			      RooFit::Warnings(0),
			      RooFit::PrintEvalErrors(-1));

   double poihat  = poi->getVal();

   RooStats::HypoTestInverter calc(*hc);
   calc.SetConfidenceLevel(0.95);

   calc.UseCLs(useCls);
   calc.SetVerbose(true);


   // can speed up using proof-lite
   if (mUseProof && mNProofWorkers > 1) { 
     RooStats::ProofConfig pc(*mpWs, mNProofWorkers, "", kFALSE);
     toymcs->SetProofConfig(&pc);    // enable proof
   }

   
   if (npoints > 0) {
      if (poimin >= poimax) { 
         // if no min/max given scan between MLE and +4 sigma 
         poimin = int(poihat);
         poimax = int(poihat +  4 * poi->getError());
      }
      std::cout << "Doing a fixed scan in interval : " << poimin << " , " << poimax << std::endl;
      calc.SetFixedScan(npoints,poimin,poimax);
   }
   else { 
     //poi->setMax(10*int( (poihat+ 10 *poi->getError() )/10 ) );
     std::cout << "Doing an  automatic scan in interval : " << poi->getMin() << " , " << poi->getMax() << std::endl;
   }

   RooStats::HypoTestInverterResult * r = 0;

   if (mTestMode){
     // test new functionality here
     
     // binary search for expected limit
     double _p = 0.0;
     double _p_err = 0.0;
     double _cls = 0.05;
     double _precision = 0.01;
     double _sigma = 0.0;

     // establish starting bounds
     // FIXME: check that low and high cover the interval
     std::pair<double,double> _high;
     std::pair<double,double> _low;
     _high.first = GetFirstPoiMax();
     _low.first= GetFirstPoiMin();
     _low.second = 1.0;
     _high.second = 0.0;

     // begin binary search
     double _current_poi;
     bool _adaptive_sampling = false;
     std::vector<double> v_poi;
     std::vector<double> v_cls;
     std::vector<double> v_cls_err;
     while ( fabs(_p - _cls) > _precision || !_adaptive_sampling ){

       // predict the next point
       std::cout << _low.first << std::endl;
       std::cout << _low.second << std::endl;
       std::cout << _high.first << std::endl;
       std::cout << _high.second << std::endl;

       // turn on adaptive sampling if close
       // FIXME: hardcoded constants
       if ( fabs(_p - _cls) < 0.04 ){
	 //((RooStats::FrequentistCalculator*) calc.GetHypoTestCalculator())->SetNToysInTails(400,400);
	 if (!_adaptive_sampling){
	   _adaptive_sampling = true;
	 }
       }
       else{
	 ((RooStats::FrequentistCalculator*) calc.GetHypoTestCalculator())->SetNToysInTails(0,0);
       }

       // guess the current step
       //
       // FIXME:
       // 1. save all points and fit with exp
       // 2. track if we're getting closer, stop if not
       // 3. estimate uncertainty on the limit
       //
       if (_adaptive_sampling){
	 //_current_poi = _low.first + (_high.first-_low.first)*(_low.second-_cls)/(_low.second-_high.second);
	 //_current_poi = (_low.first*(_cls-_high.second)+_high.first*(_low.second-_cls))/(_low.second-_high.second);
	 //_current_poi = GuessNextPoiStep(v_poi, v_cls, v_cls_err, _cls);
	 _current_poi = (_low.first + _high.first)/2.0;
       }
       else{
	 _current_poi = (_low.first + _high.first)/2.0;
       }
       
       // FIXME: adaptinve
       //calc.RunOnePoint(_current_poi);
       std::pair<double,double> point = GetClsSinglePoint(calc, _current_poi, _precision/_cls,
							  false, _sigma, 10000);
       _p = point.first;
       _p_err = point.second;
       r = calc.GetInterval();
       //_p = GetExpectedPValue(r, r->ArraySize()-1, _sigma);
       //_p_err = GetExpectedPValueErr(r, r->ArraySize()-1, _p);
       
       // turn off adaptive sampling if far away
       // FIXME: hardcoded constants
       if ( fabs(_p - _cls) > 0.04 ) _adaptive_sampling = false;
       
       // cash the point
       // FIXME: which points to save - useful for prediction?
       //if (_adaptive_sampling){
       if ( _p > 0.0001 ){
	 v_poi.push_back(_current_poi);
	 v_cls.push_back(_p);
	 v_cls_err.push_back(_p_err);
       }
       
       // adjust the estimate of the POI interval that covers cls=0.05
       // FIXME: this is a pretty rough estimate
       //        if this is useful at all, should be done with uncertainty
       //        i.e. using some sort of a normalized residual
       if (_p > _cls){
	 _low.first = _current_poi;
	 _low.second = _p;
       }
       else{
	 if (_high.first > _current_poi){
	   _high.first = _current_poi;
	   _high.second = _p;
	 }
       }
       
       std::cout << "POI = " << _current_poi
		 << ", CLs = " << _p 
		 << " +/- " << _p_err
		 << std::endl;
     }
     
     // FIXME: debug: plot points that supposed to search for cls=0.05
     std::cout << "******************************" << std::endl << std::endl;
     std::vector<double> v_index;
     for (unsigned int i=0; i!=v_poi.size(); ++i){
       v_index.push_back((double)i);
       std::cout << i << ": cls = "
		 << v_cls[i] << " +/- "
		 << v_cls_err[i] << ", relative: "
		 << v_cls_err[i]/v_cls[i]
		 << std::endl;
     }
     std::cout << "******************************" << std::endl << std::endl;

     //std::cout << "fabs(_p - _cls) = " << fabs(_p - _cls) << std::endl;
     //std::cout << "_precision      = " << _precision << std::endl;

     //calc.RunOnePoint(0.15);
     //r = calc.GetInterval();
     std::cout << "ArraySize: " << r->ArraySize() << std::endl;
     std::cout << "CLsb: " << r->CLsplusb(r->ArraySize()-1) << std::endl;
     std::cout << "CLb: " << r->CLb(r->ArraySize()-1) << std::endl;
     std::cout << "CLs: " << r->CLs(r->ArraySize()-1) << std::endl;

     //double _p = GetExpectedPValue(r, 0, nSigma );
     std::cout << "POI = " << _current_poi << std::endl;
     std::cout << "exp median: " << GetExpectedPValue(r, r->ArraySize()-1, 2) << std::endl;
     std::cout << "exp median: " << GetExpectedPValue(r, r->ArraySize()-1, 1) << std::endl;
     std::cout << "exp median: " << GetExpectedPValue(r, r->ArraySize()-1, 0) << std::endl;
     std::cout << "exp median: " << GetExpectedPValue(r, r->ArraySize()-1, -1) << std::endl;
     std::cout << "exp median: " << GetExpectedPValue(r, r->ArraySize()-1, -2) << std::endl;

     // stop here for testing
     return r;
   }

   r = calc.GetInterval();

   if (mRebuildSamplingDist) {
     calc.SetCloseProof(1);
     RooStats::SamplingDistribution * pLimDist = calc.GetUpperLimitDistribution(true,mNToysToRebuild);
     if (pLimDist) { 
       std::cout << "expected up limit " << pLimDist->InverseCDF(0.5) << " +/- " 
		 << pLimDist->InverseCDF(0.16) << "  " 
		 << pLimDist->InverseCDF(0.84) << "\n"; 
     }
     else 
       std::cout << "ERROR : failed to re-build distributions " << std::endl; 
   }

   //update r to a new re-freshed copied
   r = calc.GetInterval();

   return r; 
}



Double_t
LimitCalc::GuessNextPoiStep( std::vector<double> & vPoi,
			     std::vector<double> & vCls,
			     std::vector<double> & vClsErr,
			     double cls ){
  //
  // Use points in a POI CLs scan to predict the most likely
  // POI value that would correspond to the given cls value
  //

  // default value
  Double_t next_poi = 0.5*( GetFirstPoiMax()+GetFirstPoiMin() );

  // only use the last 4 points at most
  int input_size = vPoi.size();
  int n_points = vPoi.size();
  if (n_points > 4) n_points = 4;

  std::vector<double> _poi;
  std::vector<double> _cls;
  std::vector<double> _poi_err;
  std::vector<double> _cls_err;
  for (int i = n_points; i != 0; --i){
    _poi.push_back(vPoi[input_size-i]);
    _cls.push_back(vCls[input_size-i]);
    _poi_err.push_back(0.0);
    _cls_err.push_back(vClsErr[input_size-i]);

    // test
    std::cout << vPoi[input_size-n_points] << std::endl;
    std::cout << _poi[n_points-i]
	      << ": poi = " << _poi[n_points-i]
	      << " /- "  << _poi_err[n_points-i]
	      << ", cls = " << _cls[n_points-i]
	      << " /- "  << _cls_err[n_points-i]
	      << std::endl;
  }
  
  double poi_max = TMath::MaxElement(n_points, &_poi[0]);
  double poi_min = TMath::MinElement(n_points, &_poi[0]);

  TGraphErrors graph(n_points,
		     &_poi[0], &_cls[0],
		     &_poi_err[0], &_cls_err[0]);

  TF1 f1("f1", "[0]+[1]*exp([2]*x+[3])", poi_min, poi_max);
  TF1 f2("f2", "[0]+[1]*x", poi_min, poi_max);

  
  TFitResultPtr fit_result;
  if (n_points < 2){
    // safety: step 10% down if not enough points to fit
    double _poi = 0.9 * vPoi[n_points-1];
    if (_poi > GetFirstPoiMin()) next_poi = _poi;
  }
  else if (n_points < 40){
    // linear fit
    fit_result = graph.Fit("f2", "MEWS", "", poi_min, poi_max);
    double p0 = fit_result->Value(0);
    double p1 = fit_result->Value(1);
    next_poi = (cls-p0)/p1;

    // DEBUG: draw a plot with the last few points
    /*
    TCanvas c("c","c");
    c.DrawFrame(0.0, 0.0, 0.1, 0.2);
    //TBox f(poi_min, 0.0, poi_max, 0.2);
    //f.Draw("l");
    graph.Draw("P");
    char buf[256];
    sprintf(buf, "plot_%i.png", input_size);
    c.SaveAs(buf);
    */
  }
  else{
    // exponential fit
    fit_result = graph.Fit("f1", "MEWS", "", poi_min, poi_max);
    double p0 = fit_result->Value(0);
    double p1 = fit_result->Value(1);
    double p2 = fit_result->Value(2);
    double p3 = fit_result->Value(3);
    next_poi = (TMath::Log((cls-p0)/p1)-p3)/p2;
  }

  return next_poi;
}



//--------> global functions --------------------------------------
//
std::pair<double,double> 
get_cls_value( bool isObserved,
	       double poi,
	       double precision,
	       double sigma,
	       int nToys){
  //
  // Computes the CLs VALUE (not limit!) with a requested
  // precision or until the max number of pseudoexperiments
  // is reached.
  // 
  // Returns the pair: value-uncertainty
  //
  
  // instantiate calculator
  LimitCalc * pCalc = LimitCalc::GetInstance();

  RooStats::HypoTestInverter * calc = pCalc->GetHypoTestInverter();
  
  std::pair<double,double> result;
  result = pCalc->GetClsSinglePoint(*calc, poi, precision,
				    isObserved, sigma, nToys);

  return result;
}



int load( const char * inFileName,
	  const char * workspaceName,
	  const char * datasetName ){
  //
  // Load workspace and data from a file
  // Returns a bitset: 
  //   - the youngest bit is workspace loaded
  //   - the second bit is dataset loaded
  //

  std::string _legend = "[load]: ";

  int result = 0;

  // instantiate calculator
  LimitCalc * pCalc = LimitCalc::GetInstance();

  // load workspace
  if (inFileName && workspaceName){
    pCalc->LoadWorkspace(inFileName, workspaceName);
    ++result;
  }

  // load dataset
  if (datasetName){
    pCalc->LoadData(datasetName);
    result += 2;
  }

  return result;
}



LimitResult limit( const char * method,
		   const char * inFileName,
		   const char * workspaceName,
		   const char * datasetName ){
  //
  // Do one of the prepackaged limit calculations
  //
  std::string _legend = "[limit]: ";

  // instantiate calculator
  LimitCalc * pCalc = LimitCalc::GetInstance();

  // load workspace
  if (inFileName && workspaceName){
    pCalc->LoadWorkspace(inFileName, workspaceName);
  }

  // load dataset
  if (datasetName){
    pCalc->LoadData(datasetName);
  }

  LimitResult limitResult;
  if (!method || std::string(method).find("no_limit") != std::string::npos){
    std::cout << _legend
	      << "no limit calculation requested, doing nothing"
	      << std::endl;
  }
  else if (std::string(method).find("cls") != std::string::npos){
    limitResult = pCalc->GetClsLimit(0, 1000, true);
    //limitResult = pCalc->GetClsLimit(1000, true);
  }
  else if (std::string(method).find("plr") != std::string::npos){
    //limitResult = pCalc->GetPlrLimit(0, 1000, true);
    std::cout << _legend
	      << "profile likelihood ratio"
	      << std::endl;
  }
  else{
    std::cout << _legend
	      << "method " << method << "is unknown, exiting"
	      << std::endl;
    std::exit(-1);
  }

  return limitResult;
}



