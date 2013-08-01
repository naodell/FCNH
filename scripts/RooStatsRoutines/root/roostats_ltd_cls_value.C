//===============================================================>
//
// roostats_ltd_test.C
//
// Demo for roostats_ltd.C
//   - compute a CLs value (not limit!) for either observed data
//     or expected, with control over the precision and
//     number of pseudoexperiments.
//     This way, a user can either guarantee the precision
//     or control the time it takes to compute the value
//     and still know the achieved precision
//
// Usage:
//        .L roostats_ltd.C+
//        .x roostats_ltd_cls_value.C
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

  // first parameter: true if observed, false for expected
  // second parameter: value of POI, for which CLs is evaluated
  // third parameter: desired precision, 0.1 means 10%
  // fourth parameter: quantile for expected, 0.0 for median, 1.0 for +1 sigma
  // fifth parameter: max number of toys for each S+B and B
  std::pair<double,double> 
    paCls = get_cls_value(false, 0.15, 0.1, 1.0, 1000);
  
    
  std::cout << "[roostats_ltd_test]: CLs = "
	    << paCls.first << " +/- "
	    << paCls.second
	    << std::endl;

  // all done
  std::cout << "[roostats_ltd_cls_value]: done." << std::endl;
}
