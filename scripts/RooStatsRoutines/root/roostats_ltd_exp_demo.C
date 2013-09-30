//===============================================================>
//
// roostats_ltd_exp_demo.C
//
// Demo for using roostats_ltd.C for computing CLs expected values
// (Does not compute limits by itself)
// A user needs to scan a range of POI computing CLs at each scan
// point and obtain the CLs(poi). The poi value that corresponds
// to CLs=0.05 represents the limit.
//
// Usage:
//        .L roostats_ltd.C+
//        .x roostats_ltd_exp_demo.C
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

  // verbose output
  calc->SetTestMode(true);

  // freq/hybrid calc
  calc->SetInverterCalcType(0);
  calc->SetTestStatType(3);

  // compute expected CLs for given POI value
  // 1st parameter: POI value
  // 2nd parameter: precision, 0.2 means 20%
  // 3rd parameter: quantile: 0 for median, 1 for +1 sigma etc
  // 4th parameter: max number of toys for each S+B and B
  std::pair<double,double> paCls = get_cls_value(0.1,
						 0.2,
						 0.0,
						 10000);

  std::cout << "[roostats_ltd_test]: CLs = "
	    << paCls.first << " +/- "
	    << paCls.second
	    << std::endl;

  // all done
  std::cout << "[roostats_ltd_exp_demo]: done." << std::endl;
}
