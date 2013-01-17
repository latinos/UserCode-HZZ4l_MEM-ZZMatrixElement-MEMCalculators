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


///----------------------------------------------------------------------------------------------
/// MEMs::MEMs - constructor
///----------------------------------------------------------------------------------------------
MEMs::MEMs(double collisionEnergy, string PDFName, bool debug_)
{

    debug=debug_;

    /// Initialise MEKD
    m_MEKD = new MEKD(collisionEnergy, PDFName);
    /// Initialise MELA
    m_MELA = new Mela(false, collisionEnergy);
    
    /// Set some non-physical values for MEs initially
    for(int iMemCalc = 0; iMemCalc < NUM_MEMCALCS; iMemCalc++ )
        for(int iProcess = 0; iProcess < NUM_PROCESSES; iProcess++ )
            m_computedME[iProcess][iMemCalc] = -999.;

}


///----------------------------------------------------------------------------------------------
/// MEMs::computeME - Compute ME for the specified process.
///----------------------------------------------------------------------------------------------
int MEMs::computeME(Processes process, MEMCalcs calculator, vector<TLorentzVector> partP, vector<int> partId, double& me2process)
{
    /// check if process is supported
    if (!isProcSupported[process][calculator]) return ERR_PROCESS;
    
    /// perform computation according to the specified process and MEM package
    switch ( calculator ) {
    case kMEKD:         /// compute ME with MEKD
      if( (m_MEKD->computeME(m_processNameMEKD[process], partP, partId, me2process)) != 0 ) return ERR_COMPUTE;
      break;
    case kAnalytical:  	  /// compute ME with MELA
      cacheMELAcalculation(partP,partId); 
      me2process=m_computedME[process][calculator];
      break;
    case kJHUGen:       /// compute ME with JHUGen
      cacheMELAcalculation(partP,partId); 
      me2process=m_computedME[process][calculator];
      break;	  
    case kMCFM:         /// compute ME with MCFM
      cacheMELAcalculation(partP,partId); 
      me2process=m_computedME[process][calculator];
      break;	  
    case kMELA_HCP:
      cacheMELAcalculation(partP,partId); 
      me2process=m_computedME[process][calculator];
      break;
    default:
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
            if( (m_MEKD->computeKD(m_processNameMEKD[processA], m_processNameMEKD[processB], partP, partId, kd, me2processA, me2processB)) != 0 ) return ERR_COMPUTE;
            break;
        case kAnalytical:   /// compute KD with MELA
          cacheMELAcalculation(partP,partId); 
          me2processA=m_computedME[processA][calculator];
          me2processB=m_computedME[processB][calculator];
	  kd=me2processA/(me2processA+me2processB);
          break;
        case kJHUGen:       /// compute KD with JHUGen
          cacheMELAcalculation(partP,partId); 
          me2processA=m_computedME[processA][calculator];
          me2processB=m_computedME[processB][calculator];
	  kd=me2processA/(me2processA+me2processB);
	  break;
        case kMCFM:         /// compute KD with MCFM
	  cacheMELAcalculation(partP,partId); 
          me2processA=m_computedME[processA][calculator];
          me2processB=m_computedME[processB][calculator];
	  kd=me2processA/(me2processA+me2processB);
	  break;
        case kMELA_HCP:         /// compute KD with MCFM
	  cacheMELAcalculation(partP,partId); 
          me2processA=m_computedME[processA][calculator];
          me2processB=m_computedME[processB][calculator];
	  kd=me2processA/(me2processA+me2processB);
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
    //loop over MEMCalcs and loop over Processes
    for(int iMemCalc = 0; iMemCalc < NUM_MEMCALCS; iMemCalc++ ) {
        for(int iProcess = 0; iProcess < NUM_PROCESSES; iProcess++ ) {
            if (!isProcSupported[iProcess][iMemCalc]) continue;
            double me2process;
            if( (computeME(static_cast<Processes>(iProcess), static_cast<MEMCalcs>(iMemCalc), partP, partId, me2process)) != 0 ) return ERR_COMPUTE;
            m_computedME[iProcess][iMemCalc] = me2process;
        }
    }
    return NO_ERR;
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
/// MEMs::logRatio - KD function which returns ln( me2processA / me2processB )
///----------------------------------------------------------------------------------------------
double MEMs::logRatio(double me2processA, double me2processB){
    if (me2processB == 0) return -999.;
    return log( me2processA / me2processB );
}

///----------------------------------------------------------------------------------------------
/// MEMs::probRatio - KD function which returns me2processA / ( me2processA + c * me2processB )
///----------------------------------------------------------------------------------------------
double MEMs::probRatio(double me2processA, double me2processB){
    double c = 3000;
    if (me2processA + c * me2processB == 0) return -999.;
    return me2processA/( me2processA + c * me2processB );
}

//////////////////////////////////////////////////////////////////////////
///----------------------------------------------------------------------------------------------
/// MEMCalculators::cacheMELAcalculation - method to interface with Mela::computeP and cache results
///----------------------------------------------------------------------------------------------
void  MEMs::cacheMELAcalculation(vector<TLorentzVector> partP, vector<int> partId){

  if(debug)
    std::cout << " MEMs::cacheMELAcalculation " << std::endl;
  
  if( partPCache==partP &&  partIdCache==partId) // do nothing if we have already processed these
    return;

  partPCache = partP;
  partIdCache = partId;

  // NEED TO INCLUDE SOME PROTECTION SO THAT USER CANT 
  // PASS FOUR-VECTORS IN WRONG ORDER.  FOR NOW ASSUMING
  // THEY ARE PASSED AS e-,e+,mu-,mu+
  // ------------------ channel ------------------------
  int flavor;

  if(debug)
    std::cout << "initializing variables" << std::endl;

  float p0plus_melaNorm,p0plus_mela,p0minus_mela,p0hplus_mela; // analytical spin0 
  float p1_mela,p1plus_mela; // analytical spin1
  float p2_mela,p2qqb_mela; // analytical spin2
  float p0plus_VAJHU,p0minus_VAJHU, p0hplus_VAJHU; // JHUGen spin0 
  float p0plus_VAMCFM; // MCFM spin0 
  float p1_VAJHU, p1plus_VAJHU; // JHUGen spin1
  float p2_VAJHU,p2qqb_VAJHU; // JHUGen spin2
  float bkg_mela, bkg_VAMCFM,ggzz_VAMCFM,bkg_VAMCFMNorm; // background
  float p0_pt,p0_y,bkg_pt,bkg_y;                        // rapidity/pt
  float p0plus_m4l,bkg_m4l; //supermela
  

  if(abs(partId[0])==abs(partId[1])&&
     abs(partId[0])==abs(partId[2])&&
     abs(partId[0])==abs(partId[3])){

    if(abs(partId[0])==11) flavor=1;
    else flavor=2;

  }else flavor=3;

  // ---------------------------------------------------
  // ---------- COMPUTE ANGLES and MASSES --------------

  if(debug)
    std::cout << "computing angles" << std::endl;

  float costheta1, costheta2, costhetastar;
  float phi, phi1;

  mela::computeAngles(partP[0], partId[0], partP[1], partId[1], 
		      partP[2], partId[2], partP[3], partId[3],
		      costhetastar,costheta1,costheta2,phi,phi1);

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

  // ---------------------------------------------------
  m_MELA->computeP(mzz, m1, m2,
		   costhetastar,costheta1,costheta2,phi,phi1,
		   //signal probabilities
		   p0plus_melaNorm,   // higgs, analytic distribution, normalized
		   p0plus_mela,   // higgs, analytic distribution
		   p0minus_mela,  // pseudoscalar, analytic distribution
		   p0hplus_mela,  // pseudoscalar, analytic distribution
		   p0plus_VAJHU,  // higgs, vector algebra, JHUgen
		   p0minus_VAJHU, // pseudoscalar, vector algebra, JHUgen
		   p0plus_VAMCFM,// higgs, vector algebra, MCFM
		   p0hplus_VAJHU, // 0h+ (high dimensional operator), vector algebra, JHUgen
		   p1_mela,  // zprime, analytic distribution
		   p1plus_mela,  // zprime, analytic distribution
		   p1_VAJHU, // zprime, vector algebra, JHUgen,
		   p1plus_VAJHU, // 1+ (axial vector), vector algebra, JHUgen,
		   p2_mela , // graviton, analytic distribution,
		   p2qqb_mela , // graviton, analytic distribution,
		   p2_VAJHU, // graviton, vector algebra, JHUgen,
		   p2qqb_VAJHU, // graviton produced by qqbar, vector algebra, JHUgen,
		   //backgrounds
		   bkg_mela,  // background,  analytic distribution,
		   bkg_VAMCFM, // background, vector algebra, MCFM,
		   ggzz_VAMCFM, // background, vector algebra, MCFM for ggzz,
		   bkg_VAMCFMNorm, // background, vector algebra, MCFM,
		   //pt/rapidity
		   p0_pt, // multiplicative probability for signal pt,
		   p0_y, // multiplicative probability for signal y,
		   bkg_pt, // multiplicative probability for bkg pt,
		   bkg_y, // multiplicative probability for bkg y,
		   // supermela
		   p0plus_m4l,  // signal m4l probability as in datacards,
		   bkg_m4l,     // backgroun m4l probability as in datacards,
		   //optional input parameters
		   pt4l,Y4l,flavor // 1:4e, 2:4mu, 3:2e2mu (for interference effects)
		   );

  if(debug)
    std::cout << "got MEs" << std::endl;

  m_computedME[kSMHiggs][kAnalytical]      = p0plus_mela;
  m_computedME[k0minus][kAnalytical]       = p0minus_mela;
  m_computedME[k0hplus][kAnalytical]       = p0hplus_mela;
  m_computedME[k1minus][kAnalytical]       = p1_VAJHU;
  m_computedME[k1plus][kAnalytical]        = p1plus_VAJHU;
  m_computedME[k2mplus_gg][kAnalytical]    = p2_mela;
  m_computedME[k2mplus_qqbar][kAnalytical] = p2qqb_mela;
  m_computedME[kqqZZ][kAnalytical]         = bkg_mela;
  
  m_computedME[kSMHiggs][kJHUGen]       = p0plus_VAJHU;
  m_computedME[k0minus][kJHUGen]        = p0minus_VAJHU;
  m_computedME[k0hplus][kJHUGen]        = p0hplus_VAJHU;
  m_computedME[k1minus][kJHUGen]        = p1_VAJHU;
  m_computedME[k1plus][kJHUGen]         = p1plus_VAJHU;
  m_computedME[k2mplus_gg][kJHUGen]     = p2_VAJHU;
  m_computedME[k2mplus_qqbar][kJHUGen]  = p2qqb_VAJHU;

  m_computedME[kSMHiggs][kMCFM]         = p0plus_VAMCFM;
  m_computedME[kqqZZ][kMCFM]            = bkg_VAMCFM;
  m_computedME[kggZZ][kMCFM]            = ggzz_VAMCFM;

  m_computedME[kSMHiggs][kMELA_HCP]         = p0plus_melaNorm;
  m_computedME[kqqZZ][kMELA_HCP]            = bkg_mela;

  if(debug)
    std::cout << "Done!" << std::endl;

}  

#endif
