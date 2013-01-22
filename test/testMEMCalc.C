#include <vector>
#include "TLorentzVector.h"
#include "../interface/MEMCalculators.h"
#include "TLorentzRotation.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

using namespace MEMNames;
using namespace std;

vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi,double pt4l, double Y4l);

void testMEMCalc(TString fileName,int channel){

  TChain* t = new TChain("SelectedTree");
  t->Add(fileName);
  
  fileName.ReplaceAll(".root","_test.root");
  TFile* file = new TFile(fileName,"RECREATE");

  TTree* tree = (TTree*) t->CloneTree(0,"fast");

  MEMs test(8,"",false);

  double p0plus_mela_test,p2_mela_test,p2qqb_mela_test;
  double bkg_mela_test,p0plus_melaNorm_test;
  double p0hplus_mela_test,p1plus_mela_test;
  double p0minus_mela_test,p1_mela_test;

  double p0plus_VAJHU_test,p2_VAJHU_test,p2qqb_VAJHU_test;
  double p0hplus_VAJHU_test,p1plus_VAJHU_test;
  double p0minus_VAJHU_test,p1_VAJHU_test;

  double qqZZ_VAMCFM_test,ggZZ_VAMCFM_test,p0plus_VAMCFM_test;
  double qqZZ_VAMCFMNorm_test;

  double p0plus_mad, p0minus_mad, p2_mad, qqZZ_mad;

  double KDsigbkg, KD0minus, KD0hplus, KD2mplus;

  double dummy1,dummy2;

  // new branches
  tree->Branch("p0plus_mela_test",&p0plus_mela_test,"p0plus_mela_test/D");
  tree->Branch("p0hplus_mela_test",&p0hplus_mela_test,"p0hplus_mela_test/D");
  tree->Branch("p0minus_mela_test",&p0minus_mela_test,"p0minus_mela_test/D");
  tree->Branch("p1_mela_test",&p1_mela_test,"p1_mela_test/D");
  tree->Branch("p1plus_mela_test",&p1plus_mela_test,"p1plus_mela_test/D");
  tree->Branch("p2_mela_test",&p2_mela_test,"p2_mela_test/D");
  tree->Branch("p2qqb_mela_test",&p2qqb_mela_test,"p2qqb_mela_test/D");

  tree->Branch("p0plus_melaNorm_test",&p0plus_melaNorm_test,"p0plus_melaNorm_test/D");
  tree->Branch("bkg_mela_test",&bkg_mela_test,"bkg_mela_test/D");

  tree->Branch("p0plus_VAMCFM_test",&p0plus_VAMCFM_test,"p0plus_VAMCFM_test/D");
  tree->Branch("qqZZ_VAMCFM_test",&qqZZ_VAMCFM_test,"qqZZ_VAMCFM_test/D");
  tree->Branch("qqZZ_VAMCFMNorm_test",&qqZZ_VAMCFMNorm_test,"qqZZ_VAMCFMNorm_test/D");
  tree->Branch("ggZZ_VAMCFM_test",&ggZZ_VAMCFM_test,"ggZZ_VAMCFM_test/D");

  tree->Branch("p0plus_VAJHU_test",&p0plus_VAJHU_test,"p0plus_VAJHU_test/D");
  tree->Branch("p0minus_VAJHU_test",&p0minus_VAJHU_test,"p0minus_VAJHU_test/D");
  tree->Branch("p0hplus_VAJHU_test",&p0hplus_VAJHU_test,"p0hplus_VAJHU_test/D");
  tree->Branch("p1_VAJHU_test",&p1_VAJHU_test,"p1_VAJHU_test/D");
  tree->Branch("p1plus_VAJHU_test",&p1plus_VAJHU_test,"p1plus_VAJHU_test/D");
  tree->Branch("p2_VAJHU_test",&p2_VAJHU_test,"p2_VAJHU_test/D");
  tree->Branch("p2qqb_VAJHU_test",&p2qqb_VAJHU_test,"p2qqb_VAJHU_test/D");

  tree->Branch("p2_mad",&p2_mad,"p2_mad/D");
  tree->Branch("p0plus_mad",&p0plus_mad,"p0plus_mad/D");
  tree->Branch("p0minus_mad",&p0minus_mad,"p0minus_mad/D");
  tree->Branch("qqZZ_mad",&qqZZ_mad,"qqZZ_mad/D");

  tree->Branch("KDsigbkg",&KDsigbkg,"KDsigbkg/D");
  tree->Branch("KD0minus",&KD0minus,"KD0minus/D");
  tree->Branch("KD0hplus",&KD0hplus,"KD0hplus/D");
  tree->Branch("KD2mplus",&KD2mplus,"KD2mplus/D");

  // old branches

  float mzz,m1,m2,h1,h2,hs,phi,phi1,pt4l,Y4l;

  t->SetBranchAddress("ZZMass",&mzz);
  t->SetBranchAddress("Z1Mass",&m1);
  t->SetBranchAddress("Z2Mass",&m2);
  t->SetBranchAddress("helcosthetaZ1",&h1);
  t->SetBranchAddress("helcosthetaZ2",&h2);
  t->SetBranchAddress("costhetastar",&hs);
  t->SetBranchAddress("helphi",&phi);
  t->SetBranchAddress("phistarZ1",&phi1);
  t->SetBranchAddress("ZZPt",&pt4l);
  t->SetBranchAddress("ZZRapidity",&Y4l);

  vector<TLorentzVector> p4;
  vector<int> id;

  // channel 1: 4e 2: 4mu 2: 2e2mu
  if(channel==1){
    id.push_back(11);
    id.push_back(-11);
    id.push_back(11);
    id.push_back(-11);  
  }
  if(channel==2){
    id.push_back(13);
    id.push_back(-13);
    id.push_back(13);
    id.push_back(-13);  
  }
  if(channel==3){
    id.push_back(11);
    id.push_back(-11);
    id.push_back(13);
    id.push_back(-13);  
  }

  for(int i=0 ; i<t->GetEntries(); i++){

    t->GetEntry(i);
    //t->Show(i);
    if ( i % 1000 == 0 ) cout << i << "/" << t->GetEntries() << endl;

    if(mzz>140) continue;

    p4 = Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi,pt4l,Y4l);
    
    test.computeME(kSMHiggs,kAnalytical,p4,id,p0plus_mela_test);
    test.computeME(kSMHiggs,kMELA_HCP,p4,id,p0plus_melaNorm_test);
    test.computeME(kqqZZ,kMELA_HCP,p4,id,bkg_mela_test);

    test.computeME(k0hplus,kAnalytical,p4,id,p0hplus_mela_test);
    test.computeME(k0minus,kAnalytical,p4,id,p0minus_mela_test);
    test.computeME(k1minus,kAnalytical,p4,id,p1_mela_test);
    test.computeME(k1plus,kAnalytical,p4,id,p1plus_mela_test);
    test.computeME(k2mplus_gg,kAnalytical,p4,id,p2_mela_test);
    test.computeME(k2mplus_qqbar,kAnalytical,p4,id,p2qqb_mela_test);
   
    test.computeME(kggZZ,kMCFM,p4,id,ggZZ_VAMCFM_test);
    test.computeME(kSMHiggs,kMCFM,p4,id,p0plus_VAMCFM_test);
    test.computeME(kqqZZ,kMCFM,p4,id,qqZZ_VAMCFM_test);

    test.computeME(kSMHiggs,kJHUGen,p4,id,p0plus_VAJHU_test);
    test.computeME(k0hplus,kJHUGen,p4,id,p0hplus_VAJHU_test);
    test.computeME(k0minus,kJHUGen,p4,id,p0minus_VAJHU_test);
    test.computeME(k1minus,kJHUGen,p4,id,p1_VAJHU_test);
    test.computeME(k1plus,kJHUGen,p4,id,p1plus_VAJHU_test);
    test.computeME(k2mplus_gg,kJHUGen,p4,id,p2_VAJHU_test);
    test.computeME(k2mplus_qqbar,kJHUGen,p4,id,p2qqb_VAJHU_test);

    test.computeME(kSMHiggs,kMEKD,p4,id,p0plus_mad);
    test.computeME(k2mplus_gg,kMEKD,p4,id,p2_mad);
    test.computeME(k0minus,kMEKD,p4,id,p0minus_mad);
    test.computeME(kqqZZ,kMEKD,p4,id,qqZZ_mad);

    qqZZ_VAMCFMNorm_test=test.qqZZ_MCFMNorm;

    test.computeMEs(p4,id);

    test.computeKD(kSMHiggs,kJHUGen,kqqZZ,kMCFM,&MEMs::probRatio,KDsigbkg,dummy1,dummy2);
    test.computeKD(kSMHiggs,kJHUGen,k0minus,kJHUGen,&MEMs::probRatio,KD0minus,dummy1,dummy2);
    test.computeKD(kSMHiggs,kJHUGen,k0hplus,kJHUGen,&MEMs::probRatio,KD0hplus,dummy1,dummy2);
    test.computeKD(kSMHiggs,kJHUGen,k2mplus_gg,kJHUGen,&MEMs::probRatio,KD2mplus,dummy1,dummy2);
    double KD;
    test.computeKD(kSMHiggs,kMEKD,kqqZZ,kMCFM,&MEMs::probRatio,KD,dummy1,dummy2);

    if(KD>=1.0 || KD<=0.0)
      cout << "KD: " << KD << endl;
    /*
    cout << "p0plusJHU: " << p0plus_VAJHU_test << endl;
    cout << "qqZZMCFM: " << qqZZ_VAMCFMNorm_test << endl;
    cout << "sigbkg: " << KDsigbkg << endl;
    */

    /*
    cout << "p0plus: " << p0plus_mad << endl;
    cout << "p0minus: " << p0minus_mad << endl;
    cout << "p2: " << p2_mad << endl;
    cout << "qqZZ: " << qqZZ_mad << endl;

    double KD,pA,pB;
    test.computeKD(kSMHiggs,k0minus,kAnalytical,p4,id,KD,pA,pB);
    cout << "pseudoMELA: " << KD << endl;
    test.computeKD(kSMHiggs,k2mplus_gg,kAnalytical,p4,id,KD,pA,pB);
    cout << "graviMELA: " << KD << endl;
    test.computeKD(kSMHiggs,kqqZZ,kMELA_HCP,p4,id,KD,pA,pB);
    cout << "MELA: " << KD << endl;
    */

    tree->Fill();

  } // end loop over events

  file->cd();
  tree->Write();
  file->Close();

}


vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi,double pt4l, double Y4l)
{
  double phi1,phi2;
  phi1=TMath::Pi()-Phi1;
  phi2=Phi1+Phi;
    
    
  double gamma1,gamma2,beta1,beta2;
    
  gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
  gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
  beta1=sqrt(1-1/(gamma1*gamma1));
  beta2=sqrt(1-1/(gamma2*gamma2));
    
    
  //gluon 4 vectors
  TLorentzVector p1CM(0,0,Mx/2,Mx/2);
  TLorentzVector p2CM(0,0,-Mx/2,Mx/2);
    
  //vector boson 4 vectors
  TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);   
  TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);
    
  //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.
  
  TLorentzRotation Z1ToZ,Z2ToZ;
    
  Z1ToZ.Boost(0,0,beta1);
  Z2ToZ.Boost(0,0,beta2);
  Z1ToZ.RotateY(theta);
  Z2ToZ.RotateY(TMath::Pi()+theta);
    
    
  //fermions 4 vectors in vector boson rest frame
    
  TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
       
  TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
      
  TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);
    
  TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);
      
    
  // fermions 4 vectors in CM frame
    
  TLorentzVector p3CM,p4CM,p5CM,p6CM;
    
  p3CM=Z1ToZ*p3Z1;
  p4CM=Z1ToZ*p4Z1;
  p5CM=Z2ToZ*p5Z2;
  p6CM=Z2ToZ*p6Z2;

  // boost to lab frame:
  //TLorentzRotation labFrame;
  //labFrame.Boost(pt4l/Mx,0,sqrt(1.-1./(cosh(Y4l)*cosh(Y4l)));
		 
  vector<TLorentzVector> p;
  /*
*labFrame
*labFrame
*labFrame
*labFrame
  */
  p.push_back(p3CM); 
  p.push_back(p4CM);
  p.push_back(p5CM);
  p.push_back(p6CM);

  return p;
}

