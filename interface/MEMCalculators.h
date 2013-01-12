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
#include "ZZMatrixElement/MELA/interface/Mela.h"
// MEKD
#include "ZZMatrixElement/MEKD/interface/MEKD.h"


//////////////////////////////////////////////////////////////////////////
///  MEMs interface class to MELA & MEKD packages.
///
///  Provides interface to the MEKD & MELA packages to
///  computes MEs and KDs for the process specified by the user.
///
//////////////////////////////////////////////////////////////////////////

using namespace std;


class MEMs {
public:
    ///
    /// Constructor. Can specify the PDF to be use (ony CTEQ6L available at the moment).
    ///
    /// \param collisionEnergy              the sqrt(s) value in TeV (DEFAULT = 8).
	/// \param PDFName                      the name of the parton density functions to be used (DEFAULT = "", Optional: "CTEQ6L").
    ///
    MEMs(double collisionEnergy = 8, string PDFName = "");

    /// Enum type for supported processes in MELA and MEKD packages
    enum Processes    {kSMHiggs, k0hplus, k0minus, k1plus, k1minus, k2mplus_gg, k2mplus_qqbar, kqqZZ, kggZZ, NUM_PROCESSES};

    /// Enum type for supported MEM calculators from MELA and MEKD packages
    enum MEMCalcs    {kAnalytical, kMEKD, kJHUGen, kMCFM, kMELA_HCP, NUM_MEMCALCS};
    
    ///
    /// Compute individual ME for the specified process.
    ///
    /// \param[in]  process                 names of the process for which the ME should be retrieved (REQUIRED).
    /// \param[in]  calculator              name of the calculator tool to be used (REQUIRED).
    /// \param[in]  partP                   the input vector with TLorentzVectors for 4 leptons and 1 photon (REQUIRED).
    /// \param[in]  partId                  the input vecor with IDs (PDG) for 4 leptons and 1 photon (REQUIRED).
    /// \param[out] me2process              retrieved |ME|^2 for the specified process and calculator.
    /// \return                             error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS, 2 = ERR_COMPUTE
    ///
    int computeME(Processes process, MEMCalcs calculator, vector<TLorentzVector> partP, vector<int> partId, double& me2process);
    
    ///
    /// Compute individual KD and MEs for process A and process B, obtained with the specified calculator tool.
    ///
    /// \param[in]  processA, processB      names of the processes A and B for which the KDs and MEs are computed (REQUIRED).
    /// \param[in]  calculator              name of the calculator tool to be used (REQUIRED).
    /// \param[in]  partP                   the input vector with TLorentzVectors for 4 leptons and 1 photon (REQUIRED).
    /// \param[in]  partId                  the input vecor with IDs (PDG) for 4 leptons and 1 photon (REQUIRED).
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
    /// \param[in]  partP                   the input vector with TLorentzVectors for 4 leptons and 1 photon (REQUIRED).
    /// \param[in]  partId                  the input vecor with IDs (PDG) for 4 leptons and 1 photon (REQUIRED).
    /// \return                             error code of the computation: 0 = NO_ERR, 1 = ERR_PROCESS, 2 = ERR_COMPUTE
    ///
    int computeMEs(vector<TLorentzVector> partP, vector<int> partId);
    
    ///
    /// Retrieve ME for specified process and specified calculator tool.
    ///
    /// Method should be called only after running computeMEs(vector<TLorentzVector> partP,vector<int> partId).
    ///
    /// \param[in]  process                 names of the process for which the ME should be retrieved (REQUIRED).
    /// \param[in]  calculator              name of the calculator tool to be used (REQUIRED).
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
    /// \param[in]  processA, processB          names of the processes for which the KD and MEs are computed (REQUIRED).
    /// \param[in]  calculatorA, calculatorB    names of the calculator tools to be used (REQUIRED).
    /// \param[in]  funcKD                      name of the function to be used for KD computation (REQUIRED).
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

private:
    /// MEM calculators: MEKD (Madgraph) and MELA (Analytic, JHUGen, MCFM)
    MEKD* m_MEKD;
    Mela* m_MELA;

    /// MEKD process name translation
    static const TString m_processNameMEKD[NUM_PROCESSES];

    /// stored results of MEs computed with computeMEs(...)
    double m_computedME[NUM_PROCESSES][NUM_MEMCALCS];

    /// internal enums for supported return values/errors
    enum ERRCodes    {NO_ERR, ERR_PROCESS, ERR_COMPUTE, NUM_ERRORS};
    
};

/// Matrix of supproted processes - initialisation (to be updated)
const bool MEMs::isProcSupported[NUM_PROCESSES][NUM_MEMCALCS] = {// kAnalytical   kMEKD       kJHUGen     kMCFM       kMELA_HCP
                                                                             {0,            1,          0,          0,          0},      // kSMHiggs
                                                                             {0,            0,          0,          0,          0},      // k0hplus
                                                                             {0,            1,          0,          0,          0},      // k0minus
                                                                             {0,            0,          0,          0,          0},      // k1plus
                                                                             {0,            0,          0,          0,          0},      // k1minus
                                                                             {0,            1,          0,          0,          0},      // k2mplus_gg
                                                                             {0,            0,          0,          0,          0},      // k2mplus_qqbar
                                                                             {0,            1,          0,          0,          0},      // kqqZZ
                                                                             {0,            0,          0,          0,          0}};     // kggZZ

/// MEKD process name translation - initialisation (to be updated)
const TString MEMs::m_processNameMEKD[NUM_PROCESSES] = {"SMHiggs", "Higgs0P", "Higgs0M", "", "", "Graviton2PM", "", "ZZ", ""};

//////////////////////////////////////////////////////////////////////////

#endif
