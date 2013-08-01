//===============================================================>
//
// roostats_ltd_test.C
//
// Test for using roostats_ltd.C for computing CLs limits
//
// Usage:
//        .L roostats_ltd.C+
//        .x roostats_ltd_test.C
//
// 2011 Gena Kukartsev
//
//===============================================================>
{

  // load the library
  gSystem->Load("roostats_ltd_C.so");

  load( "ws_cl95.root",
	"ws",
	"observed_data" );

  LimitCalc * calc = LimitCalc::GetInstance();

  // load S+B model config
  calc->SetSbModelConfig("SbModel");

  // load B model config
  calc->SetBModelConfig("BModel");

  // verbose output
  calc->SetTestMode(false);
  
  // freq/hybrid calc
  calc->SetInverterCalcType(0);
  calc->SetTestStatType(3);
  calc->SetMaxZeroToys(1000);

  std::vector<double> vCls;
  std::vector<double> vClsE;
  for (int i=0; i<1; ++i){
    std::pair<double,double> paCls = get_cls_value(false,
						   0.15,
						   0.1,
						   1.0,
						   1000);
    
    vCls.push_back( paCls.first );

    double _rms = TMath::RMS(vCls.size(), &vCls[0]);

    std::cout << "[roostats_ltd_test]: iteration " << i
	      << ": running RMS(CLs) = "
	      << _rms
	      << std::endl;
    //<< std::endl;
  }
    
  std::cout << "[roostats_ltd_test]: CLs = "
	    << paCls.first << " +/- "
	    << paCls.second
	    << std::endl;

  // all done
  std::cout << "[roostats_ltd_test]: done." << std::endl;
}
