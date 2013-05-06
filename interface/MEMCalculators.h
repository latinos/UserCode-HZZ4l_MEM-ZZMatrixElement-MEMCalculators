/*************************************************************************
*  Authors:   MEKD & MELA fans
*************************************************************************/
#ifndef MEMCalc_MEMCalc_h
#define MEMCalc_MEMCalc_h

// C++ includes
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

// ROOT includes
#include "TString.h"
#include "TLorentzVector.h"

// MELA
//#include "ZZMatrixElement/MELA/interface/Mela.h"
class Mela;
// MEKD
//#include "ZZMatrixElement/MEKD/interface/MEKD.h"
class MEKD;

using namespace std;

//////////////////////////////////////////////////////////////////////////
///
///  MEMNames namespace provides enum types for names of processes and
///  names of tools/calculators supported by MELA and MEKD packages.
///
///  More details can be found at the TWiki:
///    https://twiki.cern.ch/twiki/bin/view/CMS/HZZ4lME
///
//////////////////////////////////////////////////////////////////////////
namespace MEMNames {
    /// Enum type for supported processes in MELA and MEKD packages
  enum Processes    {kSMHiggs, k0hplus, k0minus, k1plus,k1plus_prodIndep, k1minus,k1minus_prodIndep, k2mplus_gg, k2mplus_qqbar,k2mplus_prodIndep, k2hplus, k2hminus, k2bplus, kqqZZ,kqqZZ_prodIndep, kggZZ, NUM_PROCESSES};
  
  /// Enum type for supported MEM calculators from MELA and MEKD packages
  enum MEMCalcs    {kAnalytical, kMEKD, kJHUGen, kMCFM, kMELA_HCP, NUM_MEMCALCS};

  enum SuperKDsyst {kNone, kScaleUp, kScaleDown, kResolUp, kResolDown, NUM_SuperKDsyst};
        
}

//////////////////////////////////////////////////////////////////////////
///
///  MEMs class provides an interface to the MEKD & MELA packages necessary
///  to compute///  for the processes and by tools specified by the user.
///
//////////////////////////////////////////////////////////////////////////
using namespace MEMNames;

class MEMs {
public:
    ///
    /// Constructor. Can specify the PDF to be use (ony CTEQ6L available at the moment).
    ///
    /// \param collisionEnergy              the sqrt(s) value in TeV (DEFAULT = 8).
	/// \param PDFName                      the name of the parton density functions to be used (DEFAULT = "", Optional: "CTEQ6L").
    ///
    MEMs(double collisionEnergy = 8, string PDFName = "", bool debug_=false);
    
    ///
    /// Compute individual ME for the specified process.
    ///
    /// \param[in]  process                 names of the process for which the ME should be retrieved.
    /// \param[in]  calculator              name of the calculator tool to be used.
    /// \param[in]  partP                   the input vector with TLorentzVectors for 4 leptons and 1 photon.
    /// \param[in]  partId                  the input vecor with IDs (PDG) for 4 leptons and 1 photon.
    /// \param[out] me2process              retrieved |ME|^2 for the specified process and calculator.
    /// \return                             error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS, 2 = ERR_COMPUTE
    ///
    int computeME(Processes process, MEMCalcs calculator, vector<TLorentzVector> partP, vector<int> partId, double& me2process);
    
    ///
    /// Compute individual KD and MEs for process A and process B, obtained with the specified calculator tool.
    ///
    /// \param[in]  processA, processB      names of the processes A and B for which the KDs and MEs are computed.
    /// \param[in]  calculator              name of the calculator tool to be used.
    /// \param[in]  partP                   the input vector with TLorentzVectors for 4 leptons and 1 photon.
    /// \param[in]  partId                  the input vecor with IDs (PDG) for 4 leptons and 1 photon.
    /// \param[out] kd                      computed KD value for discrimination of processes A and B.
    /// \param[out] me2processA             computed |ME|^2 for process A.
    /// \param[out] me2processB             computed |ME|^2 for process B.
    /// \return                             error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS, 2 = ERR_COMPUTE
    ///
    int computeKD(Processes processA, Processes processB, MEMCalcs calculator, vector<TLorentzVector> partP, vector<int> partId, double& kd, double& me2processA, double& me2processB );

    ///
    /// Compute MEs for all supported processes.
    ///
    /// Individual MEs and KDs can be retrieved using retrieveME(Processes,MEMCalcs,double&) and computeKD(Processes,MEMCalcs,Processes,MEMCalcs,double(*)(double,double),double&,double&,double&).
    ///
    /// \param[in]  partP                   the input vector with TLorentzVectors for 4 leptons and 1 photon.
    /// \param[in]  partId                  the input vecor with IDs (PDG) for 4 leptons and 1 photon.
    /// \return                             error code of the computation: 0 = NO_ERR, 2 = ERR_COMPUTE
    ///
    int computeMEs(vector<TLorentzVector> partP, vector<int> partId);
    
    ///
    /// Retrieve ME for specified process and specified calculator tool.
    ///
    /// Method should be called only after running computeMEs(vector<TLorentzVector> partP,vector<int> partId).
    ///
    /// \param[in]  process                 names of the process for which the ME should be retrieved.
    /// \param[in]  calculator              name of the calculator tool to be used.
    /// \param[out] me2process              retrieved |ME|^2 for the specified process and calculator.
    /// \return                             error codes: 0 = NO_ERR, 1 = ERR_PROCESS
    ///
    int retrieveME(Processes process, MEMCalcs calculator, double& me2process);

    ///
    /// Compute KD and retrieve MEs for process A and process B, obtained with the specified calculator tool.
    /// The KD is computed using KD function specified by the user as kd = funcKD(me2processA, me2processB).
    ///
    /// Method should be called only after running computeMEs(vector<TLorentzVector> partP,vector<int> partId).
    ///
    /// \param[in]  processA, processB          names of the processes for which the KD and MEs are computed.
    /// \param[in]  calculatorA, calculatorB    names of the calculator tools to be used.
    /// \param[in]  funcKD                      name of the method to be used for KD computation.
    /// \param[out] kd                          computed KD value for discrimination of processes A and B.
    /// \param[out] me2processA                 computed |ME|^2 for process A.
    /// \param[out] me2processB                 computed |ME|^2 for process B.
    /// \return                                 error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS
    ///
    int computeKD(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB, double (MEMs::*funcKD)(double, double), double& kd, double& me2processA, double& me2processB );

    ///
    /// Compute KD and retrieve MEs for process A and process B, obtained with the specified calculator tool.
    /// The KD is computed using KD function specified by the user which has defferent implementations for
    /// different combinations of processes and calculator tools. Functions is of the form
    /// kd = funcKD(processA, calculatorA, processB, calculatorB).
    ///
    /// Method should be called only after running computeMEs(vector<TLorentzVector> partP,vector<int> partId).
    ///
    /// \param[in]  processA, processB          names of the processes for which the KD and MEs are computed.
    /// \param[in]  calculatorA, calculatorB    names of the calculator tools to be used.
    /// \param[in]  funcKD                      name of the method to be used for KD computation.
    /// \param[out] kd                          computed KD value for discrimination of processes A and B.
    /// \param[out] me2processA                 computed |ME|^2 for process A.
    /// \param[out] me2processB                 computed |ME|^2 for process B.
    /// \return                                 error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS
    ///
    int computeKD(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB, double (MEMs::*funcKD)(Processes, MEMCalcs, Processes, MEMCalcs), double& kd, double& me2processA, double& me2processB );

    ///
    /// Retrieve the interference reweighting factor for the given event, computed using the JHUGen.
    ///
    /// Method should be called only after running computeMEs(vector<TLorentzVector> partP,vector<int> partId).
    ///
    /// \return                                 interference reweighting factor for the given event.
    ///
    double getMELAWeight() { return m_weight;}
    
    //////////////////////////////////////////////////////////////////////////
    ///----------------------------------------------------------------------------------------------
    /// interface for calculating P(m4l) for superKD 
    /// syst controls whether PDF mean or width is adjusted
    /// to gauge systematic effects
    /// process: should be either kSMHiggs or kqqZZ else ERR_PROCESS
    /// partP: 4-momentum of leptons (with photons added)
    /// partId: PDG Id of leptons
    /// kNone: nominal shape is used
    /// kScaleUp/kScaleDown: mean mass shifted up/down appropriate scale error
    /// kResolUp/kResolDown: width is varied by appropriate resolution error
    
    /// sigProb: variable return with P(m4l) for signal
    /// bkgProb: variable return with P(m4l) for background
    ///----------------------------------------------------------------------------------------------
    void computePm4l(vector<TLorentzVector> partP, 
		     vector<int> partId,
		     SuperKDsyst syst,
		     double& sigProb,
		     double& bkgProb);

    /// Simple KD function: kd = log(me2processA / me2processB).
    double logRatio(double me2processA, double me2processB);
    
    /// Case-dependent KD function of a general form: kd = me2processA / (me2processA + c * me2processB).
    double probRatio(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB);

    /// Matrix of supproted processes
    static const bool isProcSupported[NUM_PROCESSES][NUM_MEMCALCS];

    /// mapping of process enums between MEMNames and MELA (defined in TVar.hh
    map<Processes,TVar::Process> MELAprocMap;
    map<Processes,TVar::Production> MELAprodMap;
    map<MEMCalcs,TVar::MatrixElement> MELAcalcMap;

    /// enums for supported return values/errors
    enum ERRCodes    {NO_ERR, ERR_PROCESS, ERR_COMPUTE, NUM_ERRORS};

private:
    /// MEM calculators: MEKD (Madgraph) and MELA (Analytic, JHUGen, MCFM)
    MEKD* m_MEKD;
    Mela* m_MELA;

    // debug flag
    bool debug;

    /// MEKD process name translation
    static const TString m_processNameMEKD[NUM_PROCESSES];

    /// stored results of MEs computed with computeMEs(...)
    double m_computedME[NUM_PROCESSES][NUM_MEMCALCS];

    ///Mike : Also add weights for interference in caching
    float m_weight;

    
    /// MELA calculation
    int cacheMELAcalculation(Processes process, MEMCalcs calculator,vector<TLorentzVector> partP, vector<int> partId,double& me2process);

    /// for calculating JHUGen/MCFM signal vs background KD
    double qqZZ_MCFMNorm;
    
    // caches to avoid multiplemela computations
    std::vector<TLorentzVector> partPCache;
    std::vector<int> partIdCache;
};


/// Matrix of supproted processes - initialisation (to be updated)
const bool MEMs::isProcSupported[MEMNames::NUM_PROCESSES][MEMNames::NUM_MEMCALCS] = {
// kAnalytical   kMEKD       kJHUGen     kMCFM       kMELA_HCP
  {1,            1,          1,          1,          1},      // kSMHiggs
  {1,            1,          1,          0,          0},      // k0hplus
  {1,            1,          1,          0,          0},      // k0minus
  {1,            1,          1,          0,          0},      // k1plus
  {1,            0,          1,          0,          0},      // k1plus_prodIndep
  {1,            1,          1,          0,          0},      // k1minus
  {1,            0,          1,          0,          0},      // k1minus_prodIndep
  {1,            1,          1,          0,          0},      // k2mplus_gg
  {1,            1,          1,          0,          0},      // k2mplus_qqbar
  {1,            0,          1,          0,          0},      // k2mplus_prodIndep
  {1,            1,          1,          0,          0},      // k2hplus
  {1,            1,          1,          0,          0},      // k2hminus
  {1,            1,          1,          0,          0},      // k2bplus
  {1,            1,          0,          1,          1},      // kqqZZ
  {1,            0,          0,          1,          0},      // kqqZZ_prodIndep
  {0,            0,          0,          1,          0}};     // kggZZ

//////////////////////////////////////////////////////////////////////////

#endif
