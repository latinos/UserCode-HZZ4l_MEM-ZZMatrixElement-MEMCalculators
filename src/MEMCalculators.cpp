/*************************************************************************
 *  Authors:   MEKD & MELA fans
 *  Contact:   ...
 *  Created:   11.01.2013.
 *************************************************************************/
#ifndef MEMCalc_MEMCalc_cpp
#define MEMCalc_MEMCalc_cpp

/// MEMs header
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "ZZMatrixElement/MEKD/interface/MEKD.h"
#include "../interface/MEMCalculators.h"
#include "ZZMatrixElement/MELA/src/computeAngles.h"

using namespace std;
using namespace MEMNames;

//////////////////////////////////////////////////////////////////////////
///  MEMs interface class to MELA & MEKD packages.
///
///  Provides interface to the MEKD & MELA packages to
///  computes MEs and KDs for the process specified by the user.
///
//////////////////////////////////////////////////////////////////////////


/// MEKD process name translation - initialisation (to be updated)
const TString MEMs::m_processNameMEKD[MEMNames::NUM_PROCESSES] = {"SMHiggs", "Spin0PH", "Higgs0M", "Spin1P", "", "Spin1M", "", "Graviton2PM", "qqGraviton2PM", "", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "ZZ", "", ""};

///----------------------------------------------------------------------------------------------
/// MEMs::MEMs - constructor
///----------------------------------------------------------------------------------------------
MEMs::MEMs(double collisionEnergy, string PDFName, bool debug_)
{

    /// Mapping between MEMs process enums and MELA process enums 
    /// - initialisation (to be updated)
    MELAprocMap[kSMHiggs]         =TVar::HZZ_4l;
    MELAprocMap[k0hplus]          =TVar::HDHZZ_4l;
    MELAprocMap[k0minus]          =TVar::PSHZZ_4l;
    MELAprocMap[k1plus]           =TVar::AVZZ_4l;
    MELAprocMap[k1plus_prodIndep] =TVar::AVZZ_4l;
    MELAprocMap[k1minus]          =TVar::VZZ_4l;
    MELAprocMap[k1minus_prodIndep]=TVar::VZZ_4l;
    MELAprocMap[k2mplus_gg]       =TVar::TZZ_4l;
    MELAprocMap[k2mplus_qqbar]    =TVar::TZZ_4l;
    MELAprocMap[k2mplus_prodIndep]=TVar::TZZ_4l;
    MELAprocMap[k2hplus]          =TVar::TZZ_2hplus_4l;
    MELAprocMap[k2hminus]         =TVar::PTZZ_2hminus_4l;
    MELAprocMap[k2bplus]          =TVar::TZZ_2bplus_4l;
    MELAprocMap[kqqZZ]            =TVar::ZZ_4e;
    MELAprocMap[kqqZZ_prodIndep]  =TVar::ZZ_4e;
    MELAprocMap[kggZZ]            =TVar::GGZZ_4l;

    /// Mapping between MEMs process enums and MELA production enums 
    /// - initialisation (to be updated)
    MELAprodMap[kSMHiggs]         =TVar::GG;
    MELAprodMap[k0hplus]          =TVar::GG;
    MELAprodMap[k0minus]          =TVar::GG;
    MELAprodMap[k1plus]           =TVar::QQB;
    MELAprodMap[k1plus_prodIndep] =TVar::INDEPENDENT;
    MELAprodMap[k1minus]          =TVar::QQB;
    MELAprodMap[k1minus_prodIndep]=TVar::INDEPENDENT;
    MELAprodMap[k2mplus_gg]       =TVar::GG;
    MELAprodMap[k2mplus_qqbar]    =TVar::QQB;
    MELAprodMap[k2mplus_prodIndep]=TVar::INDEPENDENT;
    MELAprodMap[k2hplus]          =TVar::GG;
    MELAprodMap[k2hminus]         =TVar::GG;
    MELAprodMap[k2bplus]          =TVar::GG;
    MELAprodMap[kqqZZ]            =TVar::QQB;
    MELAprodMap[kqqZZ_prodIndep]  =TVar::INDEPENDENT;
    MELAprodMap[kggZZ]            =TVar::GG;

    /// Mapping between MEMs calculator enums and MELA MatrixElement enums 
    /// - initialisation (to be updated)
    MELAcalcMap[kMCFM]      =TVar::MCFM;
    MELAcalcMap[kJHUGen]    =TVar::JHUGen;
    MELAcalcMap[kAnalytical]=TVar::ANALYTICAL;
    MELAcalcMap[kMELA_HCP]  =TVar::ANALYTICAL; 
    MELAcalcMap[kMEKD]      =TVar::MadGraph; 

    debug=debug_;

    /// Initialise MEKD
    m_MEKD = new MEKD(collisionEnergy, PDFName);
    /// Initialise MELA
    m_MELA = new Mela(collisionEnergy);
    
    /// Set some non-physical values for MEs initially
    for(int iMemCalc = 0; iMemCalc < NUM_MEMCALCS; iMemCalc++ )
        for(int iProcess = 0; iProcess < NUM_PROCESSES; iProcess++ )
            m_computedME[iProcess][iMemCalc] = -999.;

    m_weight = 0.0;
}


///----------------------------------------------------------------------------------------------
/// MEMs::computeME - Compute ME for the specified process.
///----------------------------------------------------------------------------------------------
int MEMs::computeME(Processes process, MEMCalcs calculator, vector<TLorentzVector> partP, vector<int> partId, double& me2process)
{

  if(debug)
    std::cout << "MEMs::computeME" << std::endl;

  //std::cout << m_processNameMEKD[process] << " " << calculator <<std::endl;
    /// check if process is supported
    if (!isProcSupported[process][calculator]) return ERR_PROCESS;
    
    if(debug)
      std::cout << "process is supported!" << endl;

    /// perform computation according to the specified process and MEM package
    switch ( calculator ) {
    case kMEKD:         /// compute ME with MEKD
      if( (m_MEKD->computeME(m_processNameMEKD[process], partP, partId, me2process)) != 0 ) return ERR_COMPUTE;
      //std::cout <<"MEKD: " << m_processNameMEKD[process] << " : " << me2process << std::endl;
      break;
    case kAnalytical:  	  /// compute ME with MELA
      if(debug)
	std::cout << "Analytical -> process: " << process << std::endl;
      return cacheMELAcalculation(process,calculator,partP,partId,me2process); 
      break;
    case kJHUGen:       /// compute ME with JHUGen
      if(debug)
	std::cout << "JHUGen -> process: " << process << std::endl;
      return cacheMELAcalculation(process,calculator,partP,partId,me2process);
      break;	  
    case kMCFM:         /// compute ME with MCFM
      return cacheMELAcalculation(process,calculator,partP,partId,me2process);  
      break;	  
    case kMELA_HCP:     /// compute ME with MELA_HCP
      if(debug)
	std::cout << "MELA_HCP -> process: " << process << std::endl;
      return ERR_PROCESS;
      break;
    default:
      if(debug)
	std::cout << " default case hit... don't recognize calculator" << std::endl;
      return ERR_PROCESS;
      break;
    }

    return NO_ERR;
}

///----------------------------------------------------------------------------------------------
/// MEMs::computeKD - Compute KD and MEs for the specified processes and MEM calculator.
///----------------------------------------------------------------------------------------------
int MEMs::computeKD(Processes processA, Processes processB, MEMCalcs calculator, vector<TLorentzVector> partP, vector<int> partId, double& kd, double& me2processA, double& me2processB )
{    
    /// check if processes are supported
    if (!isProcSupported[processA][calculator]) return ERR_PROCESS;
    if (!isProcSupported[processB][calculator]) return ERR_PROCESS;

    /// perform computation according to the specified process and MEM package
    switch ( calculator ) {
    case kMEKD:         /// compute KD with MEKD
      //std::cout << m_processNameMEKD[processA] << " " << m_processNameMEKD[processB] << std::endl;
            if( (m_MEKD->computeKD(m_processNameMEKD[processA], m_processNameMEKD[processB], partP, partId, kd, me2processA, me2processB)) != 0 ) return ERR_COMPUTE;
            break;
        case kAnalytical:   /// compute KD with MELA
          if(cacheMELAcalculation(processA,calculator,partP,partId,me2processA) || 
	     cacheMELAcalculation(processA,calculator,partP,partId,me2processA) )
	    return ERR_COMPUTE;
	  else{
	    kd=me2processA/(me2processA+me2processB);
	    return NO_ERR;
	  }
	  break;
        case kJHUGen:       /// compute KD with JHUGen
          if(cacheMELAcalculation(processA,calculator,partP,partId,me2processA) || 
	     cacheMELAcalculation(processA,calculator,partP,partId,me2processA) )
	    return ERR_COMPUTE;
	  else{
	    kd=me2processA/(me2processA+me2processB);
	    return NO_ERR;
	  }
          break;
        case kMCFM:         /// compute KD with MCFM
          if(cacheMELAcalculation(processA,calculator,partP,partId,me2processA) || 
	     cacheMELAcalculation(processA,calculator,partP,partId,me2processA) )
	    return ERR_COMPUTE;
	  else{
	    kd=me2processA/(me2processA+me2processB);
	    return NO_ERR;
	  }
          break;
        case kMELA_HCP:         /// compute KD with MCFM
	  return ERR_PROCESS;
          break;
        default:
            return ERR_PROCESS;
            break;
    }
    return NO_ERR;
}

///----------------------------------------------------------------------------------------------
/// MEMs::computeMEs - Compute MEs for the supported set of processes.
///----------------------------------------------------------------------------------------------
int MEMs::computeMEs(vector<TLorentzVector> partP, vector<int> partId)
{
    int returnERR = NO_ERR;
    //loop over MEMCalcs and loop over Processes
    for(int iMemCalc = 0; iMemCalc < NUM_MEMCALCS; iMemCalc++ ) {
        for(int iProcess = 0; iProcess < NUM_PROCESSES; iProcess++ ) {
            if (!isProcSupported[iProcess][iMemCalc]) continue;
            double me2process;
            if( (computeME(static_cast<Processes>(iProcess), static_cast<MEMCalcs>(iMemCalc), partP, partId, me2process)) != 0 ) {returnERR = ERR_COMPUTE; continue;}
            m_computedME[iProcess][iMemCalc] = me2process;
        }
    }
    //return NO_ERR only if all computations were successful
    return returnERR;
}

///----------------------------------------------------------------------------------------------
/// MEMs::retrieveME - Retrieve ME for specified process and specified calculator tool.
///----------------------------------------------------------------------------------------------
int MEMs::retrieveME(Processes process, MEMCalcs calculator, double& me2process)
{
    /// check if process is supported
    if (!isProcSupported[process][calculator]) return ERR_PROCESS;
    
    /// retrieve ME
    me2process = m_computedME[process][calculator];
    
    return NO_ERR;
}

///----------------------------------------------------------------------------------------------
/// MEMs::computeKD - Compute KD for process A and process B, for specified calculator.
///----------------------------------------------------------------------------------------------
int MEMs::computeKD(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB, double (MEMs::*funcKD)(double, double), double& kd, double& me2processA, double& me2processB )
{
    /// check if processes are supported
    if (!isProcSupported[processA][calculatorA]) return ERR_PROCESS;
    if (!isProcSupported[processB][calculatorB]) return ERR_PROCESS;

    //// compute KD
    me2processA = m_computedME[processA][calculatorA];
    me2processB = m_computedME[processB][calculatorB];
    kd = (*this.*funcKD)(me2processA, me2processB);
    
    return NO_ERR;
}

///----------------------------------------------------------------------------------------------
/// MEMs::computeKD - Compute KD for process A and process B, for specified calculator.
///----------------------------------------------------------------------------------------------
int MEMs::computeKD(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB, double (MEMs::*funcKD)(Processes, MEMCalcs, Processes, MEMCalcs), double& kd, double& me2processA, double& me2processB )
{
    /// check if processes are supported
    if (!isProcSupported[processA][calculatorA]) return ERR_PROCESS;
    if (!isProcSupported[processB][calculatorB]) return ERR_PROCESS;
    
    /// retrieve already computed MEs
    me2processA = m_computedME[processA][calculatorA];
    me2processB = m_computedME[processB][calculatorB];
    /// compute KD
    kd = (*this.*funcKD)(processA, calculatorA, processB, calculatorB);
    
    return NO_ERR;
}

///----------------------------------------------------------------------------------------------
/// MEMs::logRatio - KD function which returns ln( me2processA / me2processB )
///----------------------------------------------------------------------------------------------
double MEMs::logRatio(double me2processA, double me2processB){
    if (me2processB == 0) return -999.;
    return log( me2processA / me2processB );
}

///----------------------------------------------------------------------------------------------
/// MEMs::probRatio - KD function which returns me2processA / ( me2processA + c * me2processB )
///----------------------------------------------------------------------------------------------
double MEMs::probRatio(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB)
{
    /// check if processes are supported
    if (!isProcSupported[processA][calculatorA]) return ERR_PROCESS;
    if (!isProcSupported[processB][calculatorB]) return ERR_PROCESS;

    /// retrieve already computed MEs
    double me2processA = m_computedME[processA][calculatorA];
    double me2processB = m_computedME[processB][calculatorB];

    // compute KD per case basis (for the time being, find more elegant solution later)
    double c;
    // determine c
    // keep c = 1. for all MELA_HCP calculations
    // if case is not known, use c = 1. - this is not necessarily wrong though.
    if( (calculatorA==kJHUGen || calculatorA==kMEKD) && (calculatorB==kJHUGen || calculatorB==kMEKD) ){ // (JHUGen or MEKD)
        if( processA==kSMHiggs && processB==k0minus ){
            c = 1.; // for JHUGen or MEKD when 0+ vs 0-
        }else if( processA==kSMHiggs && processB==k2mplus_gg ){
            c = 1.; // for JHUGen or MEKD when 0+ vs 2m+
        }else if ( processB==kqqZZ ){
            c = 1.;
//            me2processB = qqZZ_MCFMNorm; // qqZZ_MCFMNorm/qqZZ_MCFM should be used for (JHUGen or MEKD) signal vs MEKD bkg
        }else{
            c = 1.; // default for all "non-known" cases
        }
    }else if( (calculatorA==kJHUGen || calculatorA==kMEKD) && (calculatorB==kMCFM) ){ // (JHUGen or MEKD) vs. MCFM
        if( processB==kqqZZ ){
            c = 1.;
//            me2processB = qqZZ_MCFMNorm; // qqZZ_MCFMNorm/qqZZ_MCFM should be used for (JHUGen or MEKD) signal vs MCFM bkg
        }else{
            c = 1.; // default for all "non-known" cases
        }
    }else{
        c = 1.; // default for all "non-known" cases
    }
    
    if (me2processA + c * me2processB == 0) return -999.;
    return me2processA/( me2processA + c * me2processB );
}

//////////////////////////////////////////////////////////////////////////
///----------------------------------------------------------------------------------------------
/// MEMCalculators::cacheMELAcalculation - method to interface with Mela::computeP and cache results
///----------------------------------------------------------------------------------------------
int  MEMs::cacheMELAcalculation(Processes process, MEMCalcs calculator,vector<TLorentzVector> partP, vector<int> partId,double& me2process){

  if(debug){
    std::cout << " MEMs::cacheMELAcalculation " << std::endl;
    std::cout << " process: " << process << std::endl;
    std::cout << " calculator: " << calculator << std::endl;
  }

  partPCache = partP;
  partIdCache = partId;

  // NEED TO INCLUDE SOME PROTECTION SO THAT USER CANT 
  // PASS FOUR-VECTORS IN WRONG ORDER.  FOR NOW ASSUMING
  // THEY ARE PASSED AS e-,e+,mu-,mu+
  // ------------------ channel ------------------------
  int flavor;

  if(abs(partId[0])==abs(partId[1])&&
     abs(partId[0])==abs(partId[2])&&
     abs(partId[0])==abs(partId[3])){

    if(abs(partId[0])==11) flavor=1;
    else flavor=2;

  }else flavor=3;

  if(debug)
    std::cout << "flavor: "  << flavor << std::endl;

  // ---------------------------------------------------
  // ---------- COMPUTE ANGLES and MASSES --------------

  if(debug)
    std::cout << "computing angles" << std::endl;

  float costheta1, costheta2, costhetastar;
  float phi, phi1;

  mela::computeAngles(partP[0], partId[0], partP[1], partId[1], 
		      partP[2], partId[2], partP[3], partId[3],
		      costhetastar,costheta1,costheta2,phi,phi1);

  // protections against NaNs, return ERR_COMPUTE and skip ME calc
  if(TMath::IsNaN(costhetastar)||
     TMath::IsNaN(costhetastar)||
     TMath::IsNaN(costhetastar)||
     TMath::IsNaN(costhetastar)||
     TMath::IsNaN(costhetastar) ) return ERR_COMPUTE; 

  float m1=(partP[0] + partP[1]).M();
  float m2=(partP[2] + partP[3]).M();

  TLorentzVector ZZ = (partP[0] + partP[1] + partP[2] + partP[3]);
  float mzz = ZZ.M();

  float pt4l  = ZZ.Pt();
  float Y4l   = ZZ.Rapidity(); // Fixme: should probably protect against NaN?

  if(debug){
    std::cout << "mzz: " << mzz << std::endl;
    std::cout << "m1: "  << m1 << std::endl;
    std::cout << "m2: " << m2 << std::endl;
    std::cout << "costheta1: " << costheta1 << std::endl;
    std::cout << "costheta2: " << costheta2 << std::endl;
    std::cout << "costhetastar: " << costhetastar << std::endl;
    std::cout << "phi: " << phi << std::endl;
    std::cout << "phi1: " << phi1 << std::endl;
    std::cout << "pT4l: " << pt4l << std::endl;
    std::cout << "Y4l: " << Y4l << std::endl;
  }

  // retrieve ME calculations
  // ---------------------------------------------------

  float me2process_float;
  //                 TVar::Process        TVar::MatrixElement     TVar::Production
  m_MELA->setProcess(MELAprocMap[process],MELAcalcMap[calculator],MELAprodMap[process]);

  // check if ZZ_4e is configured and event is 2e2mu event
  if(MELAprocMap[process]==TVar::ZZ_4e && MELAcalcMap[calculator]==TVar::MCFM && flavor==3)
    m_MELA->setProcess(TVar::ZZ_2e2m,MELAcalcMap[calculator],MELAprodMap[process]);

  m_MELA->computeP(mzz, m1, m2,
		   costhetastar,costheta1,costheta2,phi,phi1,
		   flavor,
		   me2process_float);

  me2process = (double) me2process_float;
  
  if(debug) cout << "MEMs::cacheMELACalculation - me2process: " << me2process << endl;

  /// Add weight calculation
  m_MELA->computeWeight(mzz, m1,  m2,
                        costhetastar,
                        costheta1,
                        costheta2,
                        phi,
                        phi1,
                        // return variables:
                        m_weight);

  if(debug)
    std::cout << "Done!" << std::endl;

  return NO_ERR;

}  

/// interface for calculating P(m4l) for superKD
void MEMs::computePm4l(vector<TLorentzVector> partP, 
		  vector<int> partId,
		  SuperKDsyst syst,
		  double& sigProb,
		  double& bkgProb){
  
  float prob_float;

  m_MELA->setProcess(TVar::SummedBackgrounds,TVar::JHUGen,TVar::GG);

  m_MELA->computePM4l(partP[0],partId[0],
		      partP[1],partId[1],
		      partP[2],partId[2],
		      partP[3],partId[3],
		      (TVar::SuperMelaSyst)syst,prob_float);

  bkgProb = (double) prob_float;

  m_MELA->setProcess(TVar::HZZ_4l,TVar::JHUGen,TVar::GG);

  m_MELA->computePM4l(partP[0],partId[0],
		      partP[1],partId[1],
		      partP[2],partId[2],
		      partP[3],partId[3],
		      (TVar::SuperMelaSyst)syst,prob_float);

  sigProb = (double) prob_float;

}

#endif
