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
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
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
//////////////////////////////////////////////////////////////////////////
namespace MEMNames {
    /// Enum type for supported processes in MELA and MEKD packages
    enum Processes    {kSMHiggs, k0hplus, k0minus, k1plus, k1minus, k2mplus_gg, k2mplus_qqbar, kqqZZ, kggZZ, NUM_PROCESSES};
    
    /// Enum type for supported MEM calculators from MELA and MEKD packages
    enum MEMCalcs    {kAnalytical, kMEKD, kJHUGen, kMCFM, kMELA_HCP, NUM_MEMCALCS};
}

//////////////////////////////////////////////////////////////////////////
///
///  MEMs class provides an interface to the MEKD & MELA packages necessary
///  to computes LO Matrix Elements (MEs) and Kinematic Discriminants (KDs)
///  for the processes and by tools specified by the user.
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


    /// compute KD as me2processA/(me2processA + c*me2processB)  
    ///    c will be determined on a case by case basis with 
    ///    the default as 1.  If case is not found 
    /// \param[in]  processA (B)    - name of process to be calculated or numerator (denominator) 
    ///                               (kSMHiggs, k0minus, etc.).
    /// \param[in]  calculatorA (B) - name of calculator to be used for processA (B)
    ///                            (kAnalytical, kMCFM, kJHUGen, kMEKD, kMELA_HCP)
    /// \param[out] kd              - kinematic discriminant
    /// \param[out] me2processA (B) - result of processA (B) calculation, |ME|^2
    ///
    /// \return                     - error code of the computation: 
    ///                               0 = NO_ERR, 1 = ERR_PROCESS, 2 = ERR_COMPUTE
    
    int computeKD(Processes processA, MEMCalcs calculatorA, 
		  Processes processB, MEMCalcs calculatorB,
		  vector<TLorentzVector> partP, vector<int> partId,
		  double& kd, double& me2processA, double& me2processB);

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
    /// \param[in]  funcKD                      name of the function to be used for KD computation.
    /// \param[out] kd                          computed KD value for discrimination of processes A and B.
    /// \param[out] me2processA                 computed |ME|^2 for process A.
    /// \param[out] me2processB                 computed |ME|^2 for process B.
    /// \return                                 error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS
    ///
    int computeKD(Processes processA, MEMCalcs calculatorA, Processes processB, MEMCalcs calculatorB, double (MEMs::*funcKD)(double, double), double& kd, double& me2processA, double& me2processB );
    
    /// Supproted KD functions, kd = f_KD(me2processA, me2processB).
    double logRatio(double me2processA, double me2processB);
    double probRatio(double me2processA, double me2processB);

    /// Matrix of supproted processes
    static const bool isProcSupported[NUM_PROCESSES][NUM_MEMCALCS];
    
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
    
    /// cache MELA calculiation from old interface
    void cacheMELAcalculation(vector<TLorentzVector> partP, vector<int> partId);

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
  {1,            0,          1,          0,          0},      // k1plus
  {1,            0,          1,          0,          0},      // k1minus
  {1,            1,          1,          0,          0},      // k2mplus_gg
  {1,            0,          1,          0,          0},      // k2mplus_qqbar
  {1,            1,          0,          1,          1},      // kqqZZ
  {0,            0,          0,          1,          0}};     // kggZZ


/// MEKD process name translation - initialisation (to be updated)
const TString MEMs::m_processNameMEKD[MEMNames::NUM_PROCESSES] = {"SMHiggs", "Spin0PH", "Higgs0M", "", "", "Graviton2PM", "", "ZZ", ""};

//////////////////////////////////////////////////////////////////////////

#endif
