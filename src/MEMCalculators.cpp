/*************************************************************************
 *  Authors:   MEKD & MELA fans
 *  Contact:   ...
 *  Created:   11.01.2013.
 *************************************************************************/
#ifndef MEMCalc_MEMCalc_cpp
#define MEMCalc_MEMCalc_cpp

/// MEMs header
#include "../interface/MEMCalculators.h"

using namespace std;

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
MEMs::MEMs(double collisionEnergy, string PDFName)
{
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
        case kAnalytical:   /// compute ME with MELA
        case kJHUGen:       /// compute ME with JHUGen
        case kMCFM:         /// compute ME with MCFM
            // ...
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
        case kJHUGen:       /// compute KD with JHUGen
        case kMCFM:         /// compute KD with MCFM
            // ...
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

#endif
