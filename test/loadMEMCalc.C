

{

  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/libmcfm.so");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMELA.so");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMEKD.so");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMEMCalculators.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h+");

}
