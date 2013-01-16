#include <vector>
#include "TLorentzVector.h"
#include "../interface/MEMCalculators.h"
#include "TLorentzRotation.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi,double pt4l, double Y4l);

void testMEMCalc(){

  TChain* t = new TChain("SelectedTree");
  t->Add("HZZ4lTree_jhuPseH125_withProbabilities_4e.root");
  
  TString fileName = "HZZ4lTree_jhuPseH125_withProbabilities_4e.root";
  fileName.ReplaceAll(".root","_test.root");
  TFile* file = new TFile(fileName,"RECREATE");

  TTree* tree = (TTree*) t->CloneTree(0,"fast");

  MEMs test(8,"",false);

  double p0plus_mela_test,p2_mela_test,p2qqb_mela_test;
  double p0hplus_mela_test,p1plus_mela_test;
  double p0minus_mela_test,p1_mela_test;

  double p0plus_VAJHU_test,p2_VAJHU_test,p2qqb_VAJHU_test;
  double p0hplus_VAJHU_test,p1plus_VAJHU_test;
  double p0minus_VAJHU_test,p1_VAJHU_test;

  double qqZZ_VAMCFM_test,ggZZ_VAMCFM_test,p0plus_VAMCFM_test;

  // new branches
  tree->Branch("p0plus_mela_test",&p0plus_mela_test,"p0plus_mela_test/D");
  tree->Branch("p0hplus_mela_test",&p0hplus_mela_test,"p0hplus_mela_test/D");
  tree->Branch("p0minus_mela_test",&p0minus_mela_test,"p0minus_mela_test/D");
  tree->Branch("p1_mela_test",&p1_mela_test,"p1_mela_test/D");
  tree->Branch("p1plus_mela_test",&p1plus_mela_test,"p1plus_mela_test/D");
  tree->Branch("p2_mela_test",&p2_mela_test,"p2_mela_test/D");
  tree->Branch("p2qqb_mela_test",&p2qqb_mela_test,"p2qqb_mela_test/D");

  tree->Branch("p0plus_VAMCFM_test",&p0plus_VAMCFM_test,"p0plus_VAMCFM_test/D");
  tree->Branch("qqZZ_VAMCFM_test",&qqZZ_VAMCFM_test,"qqZZ_VAMCFM_test/D");
  tree->Branch("ggZZ_VAMCFM_test",&ggZZ_VAMCFM_test,"ggZZ_VAMCFM_test/D");

  tree->Branch("p0plus_VAJHU_test",&p0plus_VAJHU_test,"p0plus_VAJHU_test/D");
  tree->Branch("p0minus_VAJHU_test",&p0minus_VAJHU_test,"p0minus_VAJHU_test/D");
  tree->Branch("p0hplus_VAJHU_test",&p0hplus_VAJHU_test,"p0hplus_VAJHU_test/D");
  tree->Branch("p1_VAJHU_test",&p1_VAJHU_test,"p1_VAJHU_test/D");
  tree->Branch("p1plus_VAJHU_test",&p1plus_VAJHU_test,"p1plus_VAJHU_test/D");
  tree->Branch("p2_VAJHU_test",&p2_VAJHU_test,"p2_VAJHU_test/D");
  tree->Branch("p2qqb_VAJHU_test",&p2qqb_VAJHU_test,"p2qqb_VAJHU_test/D");

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
  id.push_back(11);
  id.push_back(-11);
  id.push_back(11);
  id.push_back(-11);  

  for(int i=0 ; i<t->GetEntries(); i++){

    t->GetEntry(i);
    t->Show(i);

    p4 = Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi,pt4l,Y4l);
    
    test.computeME(test.kSMHiggs,test.kAnalytical,p4,id,p0plus_mela_test);
    test.computeME(test.k0hplus,test.kAnalytical,p4,id,p0hplus_mela_test);
    test.computeME(test.k0minus,test.kAnalytical,p4,id,p0minus_mela_test);
    test.computeME(test.k1minus,test.kAnalytical,p4,id,p1_mela_test);
    test.computeME(test.k1plus,test.kAnalytical,p4,id,p1plus_mela_test);
    test.computeME(test.k2mplus_gg,test.kAnalytical,p4,id,p2_mela_test);
    test.computeME(test.k2mplus_qqbar,test.kAnalytical,p4,id,p2qqb_mela_test);
   
    test.computeME(test.kggZZ,test.kMCFM,p4,id,ggZZ_VAMCFM_test);
    test.computeME(test.kSMHiggs,test.kMCFM,p4,id,p0plus_VAMCFM_test);
    test.computeME(test.kqqZZ,test.kMCFM,p4,id,qqZZ_VAMCFM_test);

    test.computeME(test.kSMHiggs,test.kJHUGen,p4,id,p0plus_VAJHU_test);
    test.computeME(test.k0hplus,test.kJHUGen,p4,id,p0hplus_VAJHU_test);
    test.computeME(test.k0minus,test.kJHUGen,p4,id,p0minus_VAJHU_test);
    test.computeME(test.k1minus,test.kJHUGen,p4,id,p1_VAJHU_test);
    test.computeME(test.k1plus,test.kJHUGen,p4,id,p1plus_VAJHU_test);
    test.computeME(test.k2mplus_gg,test.kJHUGen,p4,id,p2_VAJHU_test);
    test.computeME(test.k2mplus_qqbar,test.kJHUGen,p4,id,p2qqb_VAJHU_test);

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

