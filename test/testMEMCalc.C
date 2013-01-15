#include <vector>
#include "TLorentzVector.h"
#include "../interface/MEMCalculators.h"
#include <iostream>

void testMEMCalc(){

  vector<TLorentzVector> p4;
  vector<int> id;

  TLorentzVector a( 1.4322060E+01,  1.2640476E+01,  1.4482903E+02,  1.4608337E+02);
  TLorentzVector b( -8.5089784E+00, -3.0103278E+01,  5.5726599E+00,  3.1775220E+01);
  TLorentzVector c( -1.1015321E+01,  3.1412915E+00, 9.9395327E+00,  1.5165729E+01);
  TLorentzVector d( 5.2022391E+00,  1.4321511E+01,  7.8841187E+00,  1.7155999E+01);

  p4.push_back(a);
  id.push_back(11);
  p4.push_back(b);
  id.push_back(-11);
  p4.push_back(c);
  id.push_back(13);
  p4.push_back(d);
  id.push_back(-13);

  MEMs test;

  double result;

  cout << "-------- SMHiggs --------" << endl;
  test.computeME(MEMs::kSMHiggs,MEMs::kAnalytical,p4,id,result);
  cout << "analytical: " << result << endl;
  test.computeME(MEMs::kSMHiggs,MEMs::kMEKD,p4,id,result);
  cout << "madgraph: " << result << endl;
  test.computeME(MEMs::kSMHiggs,MEMs::kJHUGen,p4,id,result);
  cout << "JHUGen: " << result << endl;
  test.computeME(MEMs::kSMHiggs,MEMs::kMCFM,p4,id,result);
  cout << "MCFM: " << result << endl;

  cout << "-------- pseudoscalar --------" << endl;
  test.computeME(MEMs::k0minus,MEMs::kAnalytical,p4,id,result);
  cout << "analytical: " << result << endl;
  test.computeME(MEMs::k0minus,MEMs::kMEKD,p4,id,result);
  cout << "madgraph: " << result << endl;
  test.computeME(MEMs::k0minus,MEMs::kJHUGen,p4,id,result);
  cout << "JHUGen: " << result << endl;
  test.computeME(MEMs::k0minus,MEMs::kMCFM,p4,id,result);
  cout << "MCFM: " << result << endl;

  cout << "-------- 0hplus --------" << endl;
  test.computeME(MEMs::k0hplus,MEMs::kAnalytical,p4,id,result);
  cout << "analytical: " << result << endl;
  test.computeME(MEMs::k0hplus,MEMs::kMEKD,p4,id,result);
  cout << "madgraph: " << result << endl;
  test.computeME(MEMs::k0hplus,MEMs::kJHUGen,p4,id,result);
  cout << "JHUGen: " << result << endl;
  test.computeME(MEMs::k0hplus,MEMs::kMCFM,p4,id,result);
  cout << "MCFM: " << result << endl;

}
