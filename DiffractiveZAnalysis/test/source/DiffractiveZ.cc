//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsDiffractiveZsAnalysis#Macro_Analysis
//

//#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TMath.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "statusbar.h"
#include "DiffractiveZ.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveZEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace diffractiveAnalysis;
using namespace diffractiveZAnalysis;
using namespace eventInfo;
using namespace reweight;

void DiffractiveZ::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventdiff = new DiffractiveEvent();
  eventdiffZ = new DiffractiveZEvent();
  eventinfo = new EventInfoEvent();
  diff = tr->GetBranch("DiffractiveAnalysis");
  diffZ = tr->GetBranch("DiffractiveZAnalysis");
  info = tr->GetBranch("EventInfo");
  diff->SetAddress(&eventdiff);
  diffZ->SetAddress(&eventdiffZ);
  info->SetAddress(&eventinfo);

}

void DiffractiveZ::CreateHistos(std::string type){

  double binarrayplus[19] = {0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,5.2,6.2};
  double binarrayminus[19] = {-6.2,-5.2,-4.,-3.75,-3.5,-3.25,-3.,-2.75,-2.5,-2.25,-2.,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.};
  double binarraydelta[15] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,9.};
  double xi_bin[18]={0.0003,0.002,0.0045,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};

  std::string step0 = "without_cuts"; 
  std::string step1 = "with_trigger"; 
  std::string step2 = "step2"; 
  std::string step3 = "step3"; 
  std::string step4 = "step4"; 
  std::string step5 = "step5"; 
  std::string step6 = "step6"; 
  std::string step7 = "step7"; 
  std::string step8 = "NGapCMS";
  std::string step9 = "PGapCMS";
  std::string step10 = "NGapCMSAndCASTOR";
  std::string step11 = "PGapCMSAndCastorActivity";
  std::string step12 = "NGapCMSAndZKinP";
  std::string step13 = "PGapCMSAndZKinN";
  std::string step14 = "NGapCMSAndCASTORAndZKinP";
  std::string step15 = "PGapCMSAndCastorActivityAndZKinN";
  std::string step16 = "NGapCASTOR";
  std::string step17 = "NGapCASTORAndZKinP";
  std::string step18 = "PGapCMSAndCASTOR";
  std::string step19 = "PGapCMSAndCASTORAndZKinP";

  Folders.push_back(step0);
  Folders.push_back(step1);
  Folders.push_back(step2);
  Folders.push_back(step3);
  Folders.push_back(step4);
  Folders.push_back(step5);
  Folders.push_back(step6);
  Folders.push_back(step7);
  Folders.push_back(step8);
  Folders.push_back(step9);
  Folders.push_back(step10);
  Folders.push_back(step11);
  Folders.push_back(step12);
  Folders.push_back(step13);
  Folders.push_back(step14);
  Folders.push_back(step15);
  Folders.push_back(step16);
  Folders.push_back(step17);
  Folders.push_back(step18);
  Folders.push_back(step19);

  int nloop;

  if (type=="multiple_pileup"){
    nloop=21;
  }else if (type=="no_multiple_pileup"){
    nloop=1;
  }else{
    std::cout << "Multiple pile-up error" << std::endl;
    exit(EXIT_FAILURE);
  }

  char tag[300];

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    // Kinematics 
    m_hVector_DiLeptonMass.push_back( std::vector<TH1F*>() );
    m_hVector_DiLeptonEta.push_back( std::vector<TH1F*>() );
    m_hVector_DiLeptonPt.push_back( std::vector<TH1F*>() );
    m_hVector_DiLeptonPhi.push_back( std::vector<TH1F*>() );

    m_hVector_BosonZPt.push_back( std::vector<TH1F*>() );
    m_hVector_BosonZEta.push_back( std::vector<TH1F*>() );
    m_hVector_BosonZPhi.push_back( std::vector<TH1F*>() );
    m_hVector_BosonZMass.push_back( std::vector<TH1F*>() );

    m_hVector_LeptonsPt.push_back( std::vector<TH1F*>() );
    m_hVector_LeptonsEta.push_back( std::vector<TH1F*>() );
    m_hVector_LeptonsPhi.push_back( std::vector<TH1F*>() );

    m_hVector_LeadingLeptonPt.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonEta.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonPhi.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonCharge.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonTkDr03.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonEcalDr03.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonHcalDr03.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonIsolation.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonInnerHits.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonDCot.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonDist.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonDeltaEtaTkClu.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonDeltaPhiTkClu.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonSigmaIeIe.push_back( std::vector<TH1F*>() );
    m_hVector_LeadingLeptonHE.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonPt.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonEta.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonPhi.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonCharge.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonTkDr03.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonEcalDr03.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonHcalDr03.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonIsolation.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonInnerHits.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonDCot.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonDist.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonDeltaEtaTkClu.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonDeltaPhiTkClu.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonSigmaIeIe.push_back( std::vector<TH1F*>() );
    m_hVector_SecondLeptonHE.push_back( std::vector<TH1F*>() );
    m_hVector_deltaphiLeptons.push_back( std::vector<TH1F*>() );
    m_hVector_deltapTLeptons.push_back( std::vector<TH1F*>() );
    m_hVector_deltaetaLeptons.push_back( std::vector<TH1F*>() );
    m_hVector_tracksOutLeptonsCone03.push_back( std::vector<TH1F*>() );
    m_hVector_tracksOutLeptonsCone04.push_back( std::vector<TH1F*>() );
    m_hVector_tracksOutLeptonsCone05.push_back( std::vector<TH1F*>() );    
    m_hVector_ElectronsN.push_back( std::vector<TH1F*>() );
    m_hVector_MuonsN.push_back( std::vector<TH1F*>() );

    // Event Info
    m_hVector_RunNumber.push_back( std::vector<TH1F*>() );
    m_hVector_RunNumberZeroCastor.push_back( std::vector<TH1F*>() );
    m_hVector_RunNumberHighCastor.push_back( std::vector<TH1F*>() );
    m_hVector_vertex.push_back( std::vector<TH1F*>() );
    m_hVector_lumi.push_back( std::vector<TH1F*>() );
    m_hVector_vertexvslumi.push_back( std::vector<TH2F*>() );
    m_hVector_tracks.push_back( std::vector<TH1F*>() );
    m_hVector_tracksLow.push_back( std::vector<TH1F*>() );

    // Detector
    m_hVector_ECaloVsEta.push_back( std::vector<TH2F*>() );
    m_hVector_ECaloVsEtaTProf.push_back( std::vector<TProfile*>() );
    m_hVector_EnergyVsEtaBin1D.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHFplus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHFminus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHEplus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHEminus.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFplus_S.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFminus_S.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFplus_L.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFminus_L.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFMax.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFMin.push_back( std::vector<TH1F*>() );
    m_hVector_EnergyHFPlusVsEnergyHFMinus.push_back( std::vector<TH2F*>() );
    m_hVector_EnergyEEPlusVsEnergyEEMinus.push_back( std::vector<TH2F*>() );
    m_hVector_sumEEEminus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEEEplus.push_back( std::vector<TH1F*>() );
    m_hVector_multhf.push_back( std::vector<TH2F*>() );
    m_hVector_etcalos_p.push_back( std::vector<TH2F*>() );
    m_hVector_etcalos_n.push_back( std::vector<TH2F*>() );
    m_hVector_ECastorSector.push_back( std::vector<TH2F*>() );
    m_hVector_ECastorSectorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_CastorMultiplicity.push_back( std::vector<TH1F*>() );
    m_hVector_CastorMultiplicityVsLumi.push_back( std::vector<TH2F*>() );
    m_hVector_SectorVsTotalCastorEnergy.push_back( std::vector<TH2F*>() );
    m_hVector_SectorVsTotalCastorEnergyTProf.push_back( std::vector<TProfile*>() );
    m_hVector_ECastorSectorBin1D.push_back( std::vector<TH1F*>() );
    m_hVector_sumECastorMinus.push_back( std::vector<TH1F*>() );
    m_hVector_sumECastorMinusLow.push_back( std::vector<TH1F*>() );
    m_hVector_sumECastorAndHFMinus.push_back( std::vector<TH1F*>() );

    // Diffraction
    m_hVector_asumE.push_back( std::vector<TH1F*>() );
    m_hVector_AEcastor.push_back( std::vector<TH1F*>() );
    m_hVector_etasignedHF.push_back( std::vector<TH1F*>() );
    m_hVector_etasignedCASTOR.push_back( std::vector<TH1F*>() );
    m_hVector_XiPlusPF.push_back( std::vector<TH1F*>() );
    m_hVector_XiMinusPF.push_back( std::vector<TH1F*>() );
    m_hVector_XiPF.push_back( std::vector<TH1F*>() );
    m_hVector_pfetamincastor.push_back( std::vector<TH1F*>() );
    m_hVector_absdeltaEtaPF.push_back( std::vector<TH1F*>() );
    m_hVector_deltaEtaPF.push_back( std::vector<TH1F*>() );
    m_hVector_maxetagap.push_back( std::vector<TH1F*>() );
    m_hVector_LimPlusgap.push_back( std::vector<TH1F*>() );
    m_hVector_LimMinusgap.push_back( std::vector<TH1F*>() );
    m_hVector_SumPTLimPlusgap.push_back( std::vector<TH1F*>() );
    m_hVector_SumPTLimMinusgap.push_back( std::vector<TH1F*>() );
    m_hVector_pfetamax.push_back( std::vector<TH1F*>() );
    m_hVector_pfetamin.push_back( std::vector<TH1F*>() );
    m_hVector_asumE.push_back( std::vector<TH1F*>() );

    for (int k=0;k<nloop;k++){

      if (type=="multiple_pileup"){
	sprintf(tag,"multiple_pileup_%i",k);
      }
      else{
	sprintf(tag,"single");
      }

      char name[300];

      // Kinematics
      sprintf(name,"ElectronsN_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_ElectronsN = new TH1F(name,"Electrons per Event Distribution; Number of Electrons; N events",100,0,100);
      m_hVector_ElectronsN[j].push_back(histo_ElectronsN);

      sprintf(name,"MuonsN_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_MuonsN = new TH1F(name,"Muons per Event Distribution; Number of Muons; N events",100,0,100);
      m_hVector_MuonsN[j].push_back(histo_MuonsN);

      sprintf(name,"DiLeptonMass_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_DiLeptonMass = new TH1F(name,"DiLepton Invariant Mass Distribution; M_{ll} [GeV]; N events",500,0,500);
      m_hVector_DiLeptonMass[j].push_back(histo_DiLeptonMass);

      sprintf(name,"DiLeptonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_DiLeptonPt = new TH1F(name,"DiLepton Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_DiLeptonPt[j].push_back(histo_DiLeptonPt);

      sprintf(name,"DiLeptonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_DiLeptonEta = new TH1F(name,"DiLepton #eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_DiLeptonEta[j].push_back(histo_DiLeptonEta);

      sprintf(name,"DiLeptonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_DiLeptonPhi = new TH1F(name,"DiLepton #phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_DiLeptonPhi[j].push_back(histo_DiLeptonPhi);

      sprintf(name,"BosonZPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_BosonZPt = new TH1F(name,"Boson Z Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_BosonZPt[j].push_back(histo_BosonZPt);

      sprintf(name,"BosonZEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_BosonZEta = new TH1F(name,"Boson Z #eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_BosonZEta[j].push_back(histo_BosonZEta);

      sprintf(name,"BosonZPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_BosonZPhi = new TH1F(name,"Boson Z #phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_BosonZPhi[j].push_back(histo_BosonZPhi);

      sprintf(name,"BosonZMass_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_BosonZMass = new TH1F(name,"Boson Z Mass Distribution; M_{Z} [GeV]; N events",500,0,500);
      m_hVector_BosonZMass[j].push_back(histo_BosonZMass);

      sprintf(name,"LeptonsPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeptonsPt = new TH1F(name,"Leptons - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",500,0,500);
      m_hVector_LeptonsPt[j].push_back(histo_LeptonsPt);

      sprintf(name,"LeptonsEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeptonsEta = new TH1F(name,"Leptons - #eta Distribution; #eta; N events",60,-6,6);
      m_hVector_LeptonsEta[j].push_back(histo_LeptonsEta);

      sprintf(name,"LeptonsPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeptonsPhi = new TH1F(name,"Leptons - #phi Distribution; #phi [rad]; N events",16,-3.2,3.2);
      m_hVector_LeptonsPhi[j].push_back(histo_LeptonsPhi);





      sprintf(name,"LeadingLeptonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonPt = new TH1F(name,"Leading Lepton - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",500,0,500);
      m_hVector_LeadingLeptonPt[j].push_back(histo_LeadingLeptonPt);

      sprintf(name,"LeadingLeptonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonEta = new TH1F(name,"Leading Lepton - #eta Distribution; #eta; N events",60,-6,6);
      m_hVector_LeadingLeptonEta[j].push_back(histo_LeadingLeptonEta);

      sprintf(name,"LeadingLeptonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonPhi = new TH1F(name,"Leading Lepton - #phi Distribution; #phi [rad]; N events",16,-3.2,3.2);
      m_hVector_LeadingLeptonPhi[j].push_back(histo_LeadingLeptonPhi);

      sprintf(name,"LeadingLeptonTkDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonTkDr03 = new TH1F(name,"Leading Lepton: Tracker Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonTkDr03[j].push_back(histo_LeadingLeptonTkDr03);

      sprintf(name,"LeadingLeptonEcalDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonEcalDr03 = new TH1F(name,"Leading Lepton: ECAL Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonEcalDr03[j].push_back(histo_LeadingLeptonEcalDr03);

      sprintf(name,"LeadingLeptonHcalDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonHcalDr03 = new TH1F(name,"Leading Lepton: HCAL Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonHcalDr03[j].push_back(histo_LeadingLeptonHcalDr03);

      sprintf(name,"LeadingLeptonIsolation_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonIsolation = new TH1F(name,"Leading Lepton: Isolation; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonIsolation[j].push_back(histo_LeadingLeptonIsolation);

      sprintf(name,"LeadingLeptonInnerHits_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonInnerHits = new TH1F(name,"Leading Lepton; Number Of Expected Inner Hits; N events",500,0,50);
      m_hVector_LeadingLeptonInnerHits[j].push_back(histo_LeadingLeptonInnerHits);

      sprintf(name,"LeadingLeptonDCot_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonDCot = new TH1F(name,"Leading Lepton; DCot [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDCot[j].push_back(histo_LeadingLeptonDCot);

      sprintf(name,"LeadingLeptonDist_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonDist = new TH1F(name,"Leading Lepton; Dist [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDist[j].push_back(histo_LeadingLeptonDist);

      sprintf(name,"LeadingLeptonDeltaEtaTkClu_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonDeltaEtaTkClu = new TH1F(name,"Leading Lepton; #Delta#eta_{TK, Cluster} [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDeltaEtaTkClu[j].push_back(histo_LeadingLeptonDeltaEtaTkClu);

      sprintf(name,"LeadingLeptonDeltaPhiTkClu_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonDeltaPhiTkClu = new TH1F(name,"Leading Lepton; #Delta#phi_{TK, Cluster} [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDeltaPhiTkClu[j].push_back(histo_LeadingLeptonDeltaPhiTkClu);

      sprintf(name,"LeadingLeptonSigmaIeIe_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonSigmaIeIe = new TH1F(name,"Leading Lepton; #sigma_{i#etai#eta} [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonSigmaIeIe[j].push_back(histo_LeadingLeptonSigmaIeIe);

      sprintf(name,"LeadingLeptonHE_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LeadingLeptonHE = new TH1F(name,"Leading Lepton; HE; N events",100,0,1);
      m_hVector_LeadingLeptonHE[j].push_back(histo_LeadingLeptonHE);

      sprintf(name,"SecondLeptonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonPt = new TH1F(name,"Second Lepton - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",500,0,500);
      m_hVector_SecondLeptonPt[j].push_back(histo_SecondLeptonPt);

      sprintf(name,"SecondLeptonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonEta = new TH1F(name,"Second Lepton - #eta Distribution; #eta; N events",60,-6,6);
      m_hVector_SecondLeptonEta[j].push_back(histo_SecondLeptonEta);

      sprintf(name,"SecondLeptonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonPhi = new TH1F(name,"Second Lepton - #phi Distribution; #phi [rad]; N events",16,-3.2,3.2);
      m_hVector_SecondLeptonPhi[j].push_back(histo_SecondLeptonPhi);

      sprintf(name,"SecondLeptonTkDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonTkDr03 = new TH1F(name,"Second Lepton: Tracker Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_SecondLeptonTkDr03[j].push_back(histo_SecondLeptonTkDr03);

      sprintf(name,"SecondLeptonEcalDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonEcalDr03 = new TH1F(name,"Second Lepton: ECAL Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_SecondLeptonEcalDr03[j].push_back(histo_SecondLeptonEcalDr03);

      sprintf(name,"SecondLeptonHcalDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonHcalDr03 = new TH1F(name,"Second Lepton: HCAL Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_SecondLeptonHcalDr03[j].push_back(histo_SecondLeptonHcalDr03);

      sprintf(name,"SecondLeptonIsolation_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonIsolation = new TH1F(name,"Second Lepton: Isolation; # Isolation; [u]", 100, 0., 1.);
      m_hVector_SecondLeptonIsolation[j].push_back(histo_SecondLeptonIsolation);

      sprintf(name,"SecondLeptonInnerHits_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonInnerHits = new TH1F(name,"Second Lepton; Number Of Expected Inner Hits; N events",500,0,50);
      m_hVector_SecondLeptonInnerHits[j].push_back(histo_SecondLeptonInnerHits);

      sprintf(name,"SecondLeptonDCot_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonDCot = new TH1F(name,"Second Lepton; DCot [a.u.]; N events",100,0,1);
      m_hVector_SecondLeptonDCot[j].push_back(histo_SecondLeptonDCot);

      sprintf(name,"SecondLeptonDist_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonDist = new TH1F(name,"Second Lepton; Dist [a.u.]; N events",100,0,1);
      m_hVector_SecondLeptonDist[j].push_back(histo_SecondLeptonDist);

      sprintf(name,"SecondLeptonDeltaEtaTkClu_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonDeltaEtaTkClu = new TH1F(name,"Second Lepton; #Delta#eta_{TK, Cluster} [a.u.]; N events",100,0,1);
      m_hVector_SecondLeptonDeltaEtaTkClu[j].push_back(histo_SecondLeptonDeltaEtaTkClu);

      sprintf(name,"SecondLeptonDeltaPhiTkClu_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonDeltaPhiTkClu = new TH1F(name,"Second Lepton; #Delta#phi_{TK, Cluster} [a.u.]; N events",100,0,1);
      m_hVector_SecondLeptonDeltaPhiTkClu[j].push_back(histo_SecondLeptonDeltaPhiTkClu);

      sprintf(name,"SecondLeptonSigmaIeIe_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonSigmaIeIe = new TH1F(name,"Second Lepton; #sigma_{i#etai#eta} [a.u.]; N events",100,0,1);
      m_hVector_SecondLeptonSigmaIeIe[j].push_back(histo_SecondLeptonSigmaIeIe);

      sprintf(name,"SecondLeptonHE_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SecondLeptonHE = new TH1F(name,"Second Lepton; HE; N events",100,0,1);
      m_hVector_SecondLeptonHE[j].push_back(histo_SecondLeptonHE);

      sprintf(name,"deltaphiLeptons_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_deltaphiLeptons = new TH1F(name,"#Delta#phi_{ee} Distribution; #Delta#phi_{ee}; N events",20,0.0,3.2);
      m_hVector_deltaphiLeptons[j].push_back(histo_deltaphiLeptons);

      sprintf(name,"deltapTLeptons_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_deltapTLeptons = new TH1F(name,"#Delta pT_{ee} Distribution; #Delta pT_{ee} [GeV.c^{-1}]; N events",50,0.0,150);
      m_hVector_deltapTLeptons[j].push_back(histo_deltapTLeptons);

      sprintf(name,"deltaetaLeptons_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_deltaetaLeptons = new TH1F(name,"#Delta#eta_{ee} Distribution; #Delta#eta_{ee}; N events",50,-11,11);
      m_hVector_deltaetaLeptons[j].push_back(histo_deltaetaLeptons);

      sprintf(name,"tracksOutLeptonsCone03_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_tracksOutLeptonsCone03 = new TH1F(name,"Tracks Outside Leading and Second Lepton Cone; Number of Tracks, #DeltaR>0.3; N events",500,0,500);
      m_hVector_tracksOutLeptonsCone03[j].push_back(histo_tracksOutLeptonsCone03);

      sprintf(name,"tracksOutLeptonsCone04_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_tracksOutLeptonsCone04 = new TH1F(name,"Tracks Outside Leading and Second Lepton Cone; Number of Tracks, #DeltaR>0.4; N events",500,0,500);
      m_hVector_tracksOutLeptonsCone04[j].push_back(histo_tracksOutLeptonsCone04);

      sprintf(name,"tracksOutLeptonsCone05_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_tracksOutLeptonsCone05 = new TH1F(name,"Tracks Outside Leading and Second Lepton Cone; Number of Tracks, #DeltaR>0.5; N events",500,0,500);
      m_hVector_tracksOutLeptonsCone05[j].push_back(histo_tracksOutLeptonsCone05);


      // Event Info
      sprintf(name,"RunNumber_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_RunNumber = new TH1F(name,"Run Number; Run Number; N Event",16000,134000,150000);
      m_hVector_RunNumber[j].push_back(histo_RunNumber);

      sprintf(name,"RunNumberZeroCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_RunNumberZeroCastor = new TH1F(name,"Run Number; Run Number; N Event",16000,134000,150000);
      m_hVector_RunNumberZeroCastor[j].push_back(histo_RunNumberZeroCastor);

      sprintf(name,"RunNumberHighCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_RunNumberHighCastor = new TH1F(name,"Run Number; Run Number; N Event",16000,134000,150000);
      m_hVector_RunNumberHighCastor[j].push_back(histo_RunNumberHighCastor);

      sprintf(name,"vertex_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_vertex = new TH1F(name,"Number of Vertex; # Vertex; N events",25,0,25);
      m_hVector_vertex[j].push_back(histo_vertex);

      sprintf(name,"lumi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_lumi = new TH1F(name,"Luminosity per Bunch; L_{Bunch} [#mub^{-1}s^{-1}]; N events",25,0,2);
      m_hVector_lumi[j].push_back(histo_lumi);

      sprintf(name,"VertexVsLuminosity_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_vertexvslumi = new TH2F(name,"Vertex vs Luminosity; # Vertex; Luminosity per Bunch [#mub^{-1}s^{-1}]", 25., 0., 25., 25, 0., 2.);
      m_hVector_vertexvslumi[j].push_back(histo_vertexvslumi);

      sprintf(name,"Tracks_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_Tracks = new TH1F(name,"Tracks Multiplicity; n Tracks; N events",150,0,150);
      m_hVector_tracks[j].push_back(histo_Tracks);

      sprintf(name,"TracksLow_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_TracksLow = new TH1F(name,"Tracks Multiplicity; n Tracks; N events",40,0,80);
      m_hVector_tracksLow[j].push_back(histo_TracksLow);


      // Detector
      sprintf(name,"ECaloVsEta_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ECaloVsEta = new TH2F(name,"Calorimeter Energy X #eta; #eta; Energy [GeV]", 500, -8, 8, 100, 0., 1000.);
      m_hVector_ECaloVsEta[j].push_back(histo_ECaloVsEta);

      sprintf(name,"ECaloVsEtaTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_ECaloVsEtaTProf = new TProfile(name,"Calorimeter Energy X #eta; #eta; <Energy> [GeV]", 100, -8, 8, 0., 1000.);
      m_hVector_ECaloVsEtaTProf[j].push_back(histo_ECaloVsEtaTProf);

      sprintf(name,"EnergyVsEtaBin1D_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_EnergyVsEtaBin1D = new TH1F(name,"Calorimeter Energy X #eta; #eta; #sum Energy_{calotower} [GeV]", 500, -8, 8);
      m_hVector_EnergyVsEtaBin1D[j].push_back(histo_EnergyVsEtaBin1D);

      sprintf(name,"sumEHFplus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEHFplus = new TH1F(name,"HF^{+} - Sum of Energy; #sum E_{HF^{+}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHFplus[j].push_back(histo_sumEHFplus);

      sprintf(name,"sumEHFminus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEHFminus = new TH1F(name,"HF^{-} - Sum of Energy; #sum E_{HF^{-}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHFminus[j].push_back(histo_sumEHFminus);

      sprintf(name,"sumEHEplus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEHEplus = new TH1F(name,"HE^{+} - Sum of Energy; #sum E_{HE^{+}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHEplus[j].push_back(histo_sumEHEplus);

      sprintf(name,"sumEHEminus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEHEminus = new TH1F(name,"HE^{-} - Sum of Energy; #sum E_{HE^{-}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHEminus[j].push_back(histo_sumEHEminus);

      sprintf(name,"sumEHFplus_S_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFplus_S = new TH1F(name,"HF^{+} - Sum of Energy, Short Fibers; #sum E_{HF^{+},Short} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFplus_S[j].push_back(histo_SumEHFplus_S);

      sprintf(name,"sumEHFminus_S_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFminus_S = new TH1F(name,"HF^{-} - Sum of Energy, Short Fibers; #sum E_{HF^{-},Short} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFminus_S[j].push_back(histo_SumEHFminus_S);

      sprintf(name,"sumEHFplus_L_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFplus_L = new TH1F(name,"HF^{+} - Sum of Energy, Long Fibers; #sum E_{HF^{+},Long} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFplus_L[j].push_back(histo_SumEHFplus_L);

      sprintf(name,"sumEHFminus_L_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFminus_L = new TH1F(name,"HF^{-} - Sum of Energy, Long Fibers; #sum E_{HF^{-},Long} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFminus_L[j].push_back(histo_SumEHFminus_L);

      sprintf(name,"sumEHFMax_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFMax = new TH1F(name,"HF - Sum of Energy; #sum E_{HF,Max} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFMax[j].push_back(histo_SumEHFMax);

      sprintf(name,"sumEHFMin_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFMin = new TH1F(name,"HF - Sum of Energy; #sum E_{HF,Min} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFMin[j].push_back(histo_SumEHFMin);

      sprintf(name,"EnergyHFPlusVsEnergyHFMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_EnergyHFPlusVsEnergyHFMinus = new TH2F(name,"HF^{+} and HF^{-}; #sum Energy HF^{+} [GeV]; #sum Energy HF^{-} [GeV]; N events",1000,0.,1000.,1000,0.,1000.);
      m_hVector_EnergyHFPlusVsEnergyHFMinus[j].push_back(histo_EnergyHFPlusVsEnergyHFMinus);

      sprintf(name,"EnergyEEPlusVsEnergyEEMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_EnergyEEPlusVsEnergyEEMinus = new TH2F(name,"EE^{+} and EE^{-}; #sum Energy EE^{+} [GeV]; #sum Energy EE^{-} [GeV]; N events",1000,0.,500.,1000,0.,500.);
      m_hVector_EnergyEEPlusVsEnergyEEMinus[j].push_back(histo_EnergyEEPlusVsEnergyEEMinus);

      sprintf(name,"sumEEEplus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEEEplus = new TH1F(name,"EE^{+} - Sum of Energy; #sum E_{EE^{+}} [GeV]; N events",1000,0.,500.);
      m_hVector_sumEEEplus[j].push_back(histo_sumEEEplus);

      sprintf(name,"sumEEEminus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEEEminus = new TH1F(name,"EE^{-} - Sum of Energy; #sum E_{EE^{-}} [GeV]; N events",1000,0.,500.);
      m_hVector_sumEEEminus[j].push_back(histo_sumEEEminus);

      sprintf(name,"mHF_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_MultHF = new TH2F(name,"HF^{+} and HF^{-} Multiplicity; n HF^{+}; n HF^{-}; N events", 100, 0., 100., 100, 0., 100. );
      m_hVector_multhf[j].push_back(histo_MultHF);

      sprintf(name,"EnergyHFMinusVsCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ET_Calos_n = new TH2F(name,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 6000, 0., 3000. );
      m_hVector_etcalos_n[j].push_back(histo_ET_Calos_n);

      sprintf(name,"EnergyHFPlusVsCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ET_Calos_p = new TH2F(name,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 6000, 0., 3000. );
      m_hVector_etcalos_p[j].push_back(histo_ET_Calos_p);

      sprintf(name,"sumECastorMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumECastorMinus = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",6000,0,3000);
      m_hVector_sumECastorMinus[j].push_back(histo_sumECastorMinus);

      sprintf(name,"ECastorSector_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ECastorSector = new TH2F(name,"Castor Energy X Sector; Sector; Energy [GeV]", 17, 0, 17, 44, 0., 220.);
      m_hVector_ECastorSector[j].push_back(histo_ECastorSector);

      sprintf(name,"ECastorSectorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_ECastorSectorTProf = new TProfile(name,"Castor Energy X Sector; Sector; <Energy> [GeV]", 100, 0, 17, 0., 7000.);
      m_hVector_ECastorSectorTProf[j].push_back(histo_ECastorSectorTProf);

      sprintf(name,"ECastorSectorBin1D_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_ECastorSectorBin1D = new TH1F(name,"Castor Energy X Sector; Sector; Energy [GeV]", 17, 0, 17);
      m_hVector_ECastorSectorBin1D[j].push_back(histo_ECastorSectorBin1D);

      sprintf(name,"CastorMultiplicity_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_CastorMultiplicity = new TH1F(name,"Castor: number of sectors with activity; #Sectors; N events",81,0,81);
      m_hVector_CastorMultiplicity[j].push_back(histo_CastorMultiplicity);

      sprintf(name,"CastorMultiplicityVsLumi_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_CastorMultiplicityVsLumi = new TH2F(name,"CastorMultiplicity Vs Luminosity; Luminosity per Bunch [#mub^{-1}s^{-1}]; Castor Multiplicity",5000,0,2,81,0,81);
      m_hVector_CastorMultiplicityVsLumi[j].push_back(histo_CastorMultiplicityVsLumi);

      sprintf(name,"EnergyHFPlusVsCastorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_EnergyHFPlusVsCastorTProf = new TProfile(name,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 3000. );
      m_hVector_EnergyHFPlusVsCastorTProf[j].push_back(histo_EnergyHFPlusVsCastorTProf);

      sprintf(name,"EnergyHFMinusVsCastorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_EnergyHFMinusVsCastorTProf = new TProfile(name,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 3000. );
      m_hVector_EnergyHFMinusVsCastorTProf[j].push_back(histo_EnergyHFMinusVsCastorTProf);

      sprintf(name,"sumECastorAndSumHFMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumECastorAndHFMinus = new TH1F(name,"HF^{-} and Castor Sum of Energy; Energy [GeV]; N events",6000,0,3000);
      m_hVector_sumECastorAndHFMinus[j].push_back(histo_sumECastorAndHFMinus);

      sprintf(name,"SectorVsTotalCastorEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_SectorVsTotalCastorEnergy = new TH2F(name,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",17,0,17,1500,0,1500);
      m_hVector_SectorVsTotalCastorEnergy[j].push_back(histo_SectorVsTotalCastorEnergy);

      sprintf(name,"SectorVsTotalCastorEnergyTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_SectorVsTotalCastorEnergyTProf = new TProfile(name,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",17,0,17,0,1500);
      m_hVector_SectorVsTotalCastorEnergyTProf[j].push_back(histo_SectorVsTotalCastorEnergyTProf);

      sprintf(name,"sumECastorMinusLowEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumECastorMinusLow = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",1000,0.,500.);
      m_hVector_sumECastorMinusLow[j].push_back(histo_sumECastorMinusLow);


      // Diffraction
      sprintf(name,"aEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_aSumE = new TH1F(name,"Forward Backward Asymmetry Distribution ; (#sum HF^{+} - #sum HF^{-})x(#sum HF^{+} + #sum HF^{-})^{-1}; N events",100,-2,2);
      m_hVector_asumE[j].push_back(histo_aSumE);

      sprintf(name,"pfetamax_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_PFEtamax = new TH1F(name,"Particle Flow #eta_{max} Distribution; #eta; N events",18,binarrayplus);
      m_hVector_pfetamax[j].push_back(histo_PFEtamax);

      sprintf(name,"pfetamin_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_PFEtamin = new TH1F(name,"Particle Flow #eta_{min} Distribution; #eta; N events",18,binarrayminus);
      m_hVector_pfetamin[j].push_back(histo_PFEtamin);

      sprintf(name,"maxetagap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_maxetagap = new TH1F(name,"Particle Flow #Delta#eta_{max} Distribution; #Delta#eta_{max}; N events",40,0.,6.);
      m_hVector_maxetagap[j].push_back(histo_maxetagap);

      sprintf(name,"deltaEtamaxminPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_deltaEtaPF = new TH1F(name,"#Delta#eta_{PF} Distribution; #eta_{max}-#eta_{min}; N events",14,binarraydelta);
      m_hVector_deltaEtaPF[j].push_back(histo_deltaEtaPF);

      sprintf(name,"absdeltaEtamaxminPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_absdeltaEtaPF = new TH1F(name,"|#Delta#eta_{PF}| Distribution; |#eta_{max}-#eta_{min}|; N events",14,binarraydelta);
      m_hVector_absdeltaEtaPF[j].push_back(histo_absdeltaEtaPF);

      sprintf(name,"xiPlusPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_XiPlusPF = new TH1F(name,"#xi_{plus} Particle Flow; #xi_{plus}; N Event",17,xi_bin);
      m_hVector_XiPlusPF[j].push_back(histo_XiPlusPF);

      sprintf(name,"xiMinusPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_XiMinusPF = new TH1F(name,"#xi_{minus} Particle Flow; #xi_{minus}; N Event",17,xi_bin);
      m_hVector_XiMinusPF[j].push_back(histo_XiMinusPF);

      sprintf(name,"xiPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_XiPF = new TH1F(name,"#xi Particle Flow; #xi; N Event",17,xi_bin);
      m_hVector_XiPF[j].push_back(histo_XiPF);

      sprintf(name,"pfetamincastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_PFEtamincastor = new TH1F(name,"Particle Flow #eta_{min} Distribution; #eta; N events",18,binarrayminus);
      m_hVector_pfetamincastor[j].push_back(histo_PFEtamincastor);

      sprintf(name,"LimPlusgap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LimPlusgap = new TH1F(name,"Particle Flow #eta_{max} Gap upper limit; #eta; N events",50,-5.2,5.2);
      m_hVector_LimPlusgap[j].push_back(histo_LimPlusgap);

      sprintf(name,"SumPTMaxgap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumPTLimPlusgap = new TH1F(name,"Particle Flow P_{T}; #sum pT_{max} [GeV]; N events",1200,0,600);
      m_hVector_SumPTLimPlusgap[j].push_back(histo_SumPTLimPlusgap);

      sprintf(name,"LimMinusgap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_LimMinusgap = new TH1F(name,"Particle Flow #eta_{max} Gap lower limit; #eta; N events",50,-5.2,5.2);
      m_hVector_LimMinusgap[j].push_back(histo_LimMinusgap);

      sprintf(name,"SumPTMingap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumPTLimMinusgap = new TH1F(name,"Particle Flow P_{T}; #sum pT_{min} [GeV]; N events",1200,0,600);
      m_hVector_SumPTLimMinusgap[j].push_back(histo_SumPTLimMinusgap);

      sprintf(name,"AEcastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_AEcastor = new TH1F(name,"A_{Energy}; A_{Energy}; N events",1200,-1.5,1.5);
      m_hVector_AEcastor[j].push_back(histo_AEcastor);

      sprintf(name,"etasignedHF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_etasignedHF = new TH1F(name,"Signed #eta_{l} Distribution; #tilde{#eta_{l}}; N events",100,-10.,10.);
      m_hVector_etasignedHF[j].push_back(histo_etasignedHF);

      sprintf(name,"etasignedCASTOR_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_etasignedCASTOR = new TH1F(name,"Signed #eta_{l} Distribution; #tilde{#eta_{l}}; N events",100,-10.,10.);
      m_hVector_etasignedCASTOR[j].push_back(histo_etasignedCASTOR);

    }
  }
}

void DiffractiveZ::FillHistos(int index, int pileup, double totalweight){

  m_hVector_DiElectronMass[index].at(pileup)->Fill(eventdiffZ->GetDiElectronMass(),totalweight);
  m_hVector_DiElectronEta[index].at(pileup)->Fill(eventdiffZ->GetDiElectronEta(),totalweight);
  m_hVector_DiElectronPhi[index].at(pileup)->Fill(eventdiffZ->GetDiElectronPhi(),totalweight);
  m_hVector_DiElectronPt[index].at(pileup)->Fill(eventdiffZ->GetDiElectronPt(),totalweight);
  m_hVector_LeadingElectronPt[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronPt(),totalweight);
  m_hVector_LeadingElectronEta[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronEta(),totalweight);
  m_hVector_LeadingElectronPhi[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronPhi(),totalweight);
  m_hVector_LeadingElectronCharge[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronCharge(),totalweight);
  m_hVector_SecondElectronPt[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronPt(),totalweight);
  m_hVector_SecondElectronEta[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronEta(),totalweight);
  m_hVector_SecondElectronPhi[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronPhi(),totalweight);
  m_hVector_SecondElectronCharge[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronCharge(),totalweight);
  m_hVector_ElectronsN[index].at(pileup)->Fill(eventdiffZ->GetElectronsN(),totalweight);
  m_hVector_DiMuonMass[index].at(pileup)->Fill(eventdiffZ->GetDiMuonMass(),totalweight);
  m_hVector_DiMuonEta[index].at(pileup)->Fill(eventdiffZ->GetDiMuonEta(),totalweight);
  m_hVector_DiMuonPhi[index].at(pileup)->Fill(eventdiffZ->GetDiMuonPhi(),totalweight);
  m_hVector_DiMuonPt[index].at(pileup)->Fill(eventdiffZ->GetDiMuonPt(),totalweight);
  m_hVector_LeadingMuonPt[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonPt(),totalweight);
  m_hVector_LeadingMuonEta[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonEta(),totalweight);
  m_hVector_LeadingMuonPhi[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonPhi(),totalweight);
  m_hVector_LeadingMuonCharge[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonCharge(),totalweight);
  m_hVector_SecondMuonPt[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonPt(),totalweight);
  m_hVector_SecondMuonEta[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonEta(),totalweight);
  m_hVector_SecondMuonPhi[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonPhi(),totalweight);
  m_hVector_SecondMuonCharge[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonCharge(),totalweight);
  m_hVector_MuonsN[index].at(pileup)->Fill(eventdiffZ->GetMuonsN(),totalweight);
  m_hVector_sumEHFplus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),totalweight);
  m_hVector_sumEHFminus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),totalweight);
  m_hVector_sumEHEplus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHEPlus(),totalweight);
  m_hVector_sumEHEminus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHEMinus(),totalweight);
  m_hVector_lumi[index].at(pileup)->Fill(eventinfo->GetInstLumiBunch(),totalweight);
  m_hVector_asumE[index].at(pileup)->Fill(aSumE,totalweight);
  m_hVector_multhf[index].at(pileup)->Fill(eventdiff->GetMultiplicityHFPlus(),eventdiff->GetMultiplicityHFMinus(),totalweight);
  m_hVector_tracks[index].at(pileup)->Fill(eventdiff->GetMultiplicityTracks(),totalweight);
  m_hVector_tracksLow[index].at(pileup)->Fill(eventdiff->GetMultiplicityTracks(),totalweight);
  m_hVector_vertex[index].at(pileup)->Fill(eventdiff->GetNVertex(),totalweight);
  m_hVector_deltaphielectrons[index].at(pileup)->Fill(deltaphielectrons,totalweight);
  m_hVector_deltaphimuons[index].at(pileup)->Fill(deltaphimuons,totalweight);
  m_hVector_deltaetaelectrons[index].at(pileup)->Fill(deltaetaelectrons,totalweight);
  m_hVector_deltaetamuons[index].at(pileup)->Fill(deltaetamuons,totalweight);
  m_hVector_deltapTelectrons[index].at(pileup)->Fill(deltapTelectrons,totalweight);
  m_hVector_deltapTmuons[index].at(pileup)->Fill(deltapTmuons,totalweight);
  m_hVector_vertexvslumi[index].at(pileup)->Fill(eventdiff->GetNVertex(),eventinfo->GetInstLumiBunch(),totalweight);
  m_hVector_patNElectron[index].at(pileup)->Fill(eventdiffZ->GetPatNElectron(),totalweight);
  m_hVector_patElectron1Pt[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Pt(),totalweight);
  m_hVector_patElectron1Eta[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Eta(),totalweight);
  m_hVector_patElectron1Phi[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Phi(),totalweight);
  m_hVector_patElectron1Charge[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Charge(),totalweight);
  m_hVector_patElectron1Et[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Et(),totalweight);
  m_hVector_patNMuon[index].at(pileup)->Fill(eventdiffZ->GetPatNMuon(),totalweight);
  m_hVector_patMuon1Pt[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Pt(),totalweight);
  m_hVector_patMuon1Eta[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Eta(),totalweight);
  m_hVector_patMuon1Phi[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Phi(),totalweight);
  m_hVector_patMuon1Charge[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Charge(),totalweight);
  m_hVector_patMuon1Et[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Et(),totalweight);
  m_hVector_patElectron2Pt[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Pt(),totalweight);
  m_hVector_patElectron2Eta[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Eta(),totalweight);
  m_hVector_patElectron2Phi[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Phi(),totalweight);
  m_hVector_patElectron2Charge[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Charge(),totalweight);
  m_hVector_patElectron2Et[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Et(),totalweight);
  m_hVector_patMuon2Pt[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Pt(),totalweight);
  m_hVector_patMuon2Eta[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Eta(),totalweight);
  m_hVector_patMuon2Phi[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Phi(),totalweight);
  m_hVector_patMuon2Charge[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Charge(),totalweight);
  m_hVector_patMuon2Et[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Et(),totalweight);
  m_hVector_patDiElectronMass[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronMass(),totalweight);
  m_hVector_patDiMuonMass[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonMass(),totalweight);
  m_hVector_patDiElectronEta[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronEta(),totalweight);
  m_hVector_patDiMuonEta[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonEta(),totalweight);
  m_hVector_patDiElectronPhi[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronPhi(),totalweight);
  m_hVector_patDiMuonPhi[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonPhi(),totalweight);
  m_hVector_patDiElectronPt[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronPt(),totalweight);
  m_hVector_patDiMuonPt[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonPt(),totalweight);
  m_hVector_patElectron1TkDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1TkDr03(),totalweight);
  m_hVector_patElectron1TkDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1TkDr04(),totalweight);
  m_hVector_patElectron1EcalDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1EcalDr03(),totalweight);
  m_hVector_patElectron1EcalDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1EcalDr04(),totalweight);
  m_hVector_patElectron1HcalDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1HcalDr03(),totalweight);
  m_hVector_patElectron1HcalDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1HcalDr04(),totalweight);
  m_hVector_patElectron2TkDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2TkDr03(),totalweight);
  m_hVector_patElectron2TkDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2TkDr04(),totalweight);
  m_hVector_patElectron2EcalDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2EcalDr03(),totalweight);
  m_hVector_patElectron2EcalDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2EcalDr04(),totalweight);
  m_hVector_patElectron2HcalDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2HcalDr03(),totalweight);
  m_hVector_patElectron2HcalDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1HcalDr04(),totalweight);
  m_hVector_patElectron1relIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1relIsoDr03(),totalweight);
  m_hVector_patElectron1relIsoDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1relIsoDr04(),totalweight);
  m_hVector_patElectron2relIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2relIsoDr03(),totalweight);
  m_hVector_patElectron2relIsoDr04[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2relIsoDr04(),totalweight);
  m_hVector_patMuon1SumPtR03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1SumPtR03(),totalweight);
  m_hVector_patMuon1SumPtR05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1SumPtR05(),totalweight);
  m_hVector_patMuon1EmEtR03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1EmEtR03(),totalweight);
  m_hVector_patMuon1EmEtR05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1EmEtR05(),totalweight);
  m_hVector_patMuon1HadEtR03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1HadEtR03(),totalweight);
  m_hVector_patMuon1HadEtR05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1HadEtR05(),totalweight);
  m_hVector_patMuon2SumPtR03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2SumPtR03(),totalweight);
  m_hVector_patMuon2SumPtR05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2SumPtR05(),totalweight);
  m_hVector_patMuon2EmEtR03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2EmEtR03(),totalweight);
  m_hVector_patMuon2EmEtR05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2EmEtR05(),totalweight);
  m_hVector_patMuon2HadEtR03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2HadEtR03(),totalweight);
  m_hVector_patMuon2HadEtR05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2HadEtR05(),totalweight);
  m_hVector_patMuon1relIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1relIsoDr03(),totalweight);
  m_hVector_patMuon1relIsoDr05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1relIsoDr05(),totalweight);
  m_hVector_patMuon2relIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2relIsoDr03(),totalweight);
  m_hVector_patMuon2relIsoDr05[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2relIsoDr05(),totalweight);
  m_hVector_patMuon1relIso[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1relIso(),totalweight);
  m_hVector_patMuon2relIso[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2relIso(),totalweight);
  m_hVector_LeadingElectronTkDr03[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronTkDr03(),totalweight);
  m_hVector_LeadingElectronTkDr04[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronTkDr04(),totalweight);
  m_hVector_LeadingElectronEcalDr03[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronEcalDr03(),totalweight);
  m_hVector_LeadingElectronEcalDr04[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronEcalDr04(),totalweight);
  m_hVector_LeadingElectronHcalDr03[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronHcalDr03(),totalweight);
  m_hVector_LeadingElectronHcalDr04[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronHcalDr04(),totalweight);
  m_hVector_SecondElectronTkDr03[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronTkDr03(),totalweight);
  m_hVector_SecondElectronTkDr04[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronTkDr04(),totalweight);
  m_hVector_SecondElectronEcalDr03[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronEcalDr03(),totalweight);
  m_hVector_SecondElectronEcalDr04[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronEcalDr04(),totalweight);
  m_hVector_SecondElectronHcalDr03[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronHcalDr03(),totalweight);
  m_hVector_SecondElectronHcalDr04[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronHcalDr04(),totalweight);
  m_hVector_LeadingElectronrelIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronrelIsoDr03(),totalweight);
  m_hVector_LeadingElectronrelIsoDr04[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronrelIsoDr04(),totalweight);
  m_hVector_SecondElectronrelIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronrelIsoDr03(),totalweight);
  m_hVector_SecondElectronrelIsoDr04[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronrelIsoDr04(),totalweight);
  m_hVector_LeadingMuonSumPtR03[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonSumPtR03(),totalweight);
  m_hVector_LeadingMuonSumPtR05[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonSumPtR05(),totalweight);
  m_hVector_LeadingMuonEmEtR03[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonEmEtR03(),totalweight);
  m_hVector_LeadingMuonEmEtR05[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonEmEtR05(),totalweight);
  m_hVector_LeadingMuonHadEtR03[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonHadEtR03(),totalweight);
  m_hVector_LeadingMuonHadEtR05[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonHadEtR05(),totalweight);
  m_hVector_SecondMuonSumPtR03[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonSumPtR03(),totalweight);
  m_hVector_SecondMuonSumPtR05[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonSumPtR05(),totalweight);
  m_hVector_SecondMuonEmEtR03[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonEmEtR03(),totalweight);
  m_hVector_SecondMuonEmEtR05[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonEmEtR05(),totalweight);
  m_hVector_SecondMuonHadEtR03[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonHadEtR03(),totalweight);
  m_hVector_SecondMuonHadEtR05[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonHadEtR05(),totalweight);
  m_hVector_LeadingMuonrelIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonrelIsoDr03(),totalweight);
  m_hVector_LeadingMuonrelIsoDr05[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonrelIsoDr05(),totalweight);
  m_hVector_SecondMuonrelIsoDr03[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonrelIsoDr03(),totalweight);
  m_hVector_SecondMuonrelIsoDr05[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonrelIsoDr05(),totalweight);

  m_hVector_correlationElectron1Pt[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronPt(),eventdiffZ->GetPatElectron1Pt());
  m_hVector_correlationElectron1Eta[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronEta(),eventdiffZ->GetPatElectron1Eta());
  m_hVector_correlationElectron1Phi[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronPhi(),eventdiffZ->GetPatElectron1Phi());
  m_hVector_correlationElectron1Iso[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronrelIsoDr03(),eventdiffZ->GetPatElectron1relIsoDr03());
  m_hVector_correlationMuon1Pt[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonPt(),eventdiffZ->GetPatMuon1Pt());
  m_hVector_correlationMuon1Eta[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonEta(),eventdiffZ->GetPatMuon1Eta());
  m_hVector_correlationMuon1Phi[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonPhi(),eventdiffZ->GetPatMuon1Phi());
  m_hVector_correlationMuon1Iso[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonrelIsoDr03(),eventdiffZ->GetPatMuon1relIsoDr03());
  m_hVector_correlationElectron2Pt[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronPt(),eventdiffZ->GetPatElectron2Pt());
  m_hVector_correlationElectron2Eta[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronEta(),eventdiffZ->GetPatElectron2Eta());
  m_hVector_correlationElectron2Phi[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronPhi(),eventdiffZ->GetPatElectron2Phi());
  m_hVector_correlationElectron2Iso[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronrelIsoDr03(),eventdiffZ->GetPatElectron2relIsoDr03());
  m_hVector_correlationMuon2Pt[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonPt(),eventdiffZ->GetPatMuon2Pt());
  m_hVector_correlationMuon2Eta[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonEta(),eventdiffZ->GetPatMuon2Eta());
  m_hVector_correlationMuon2Phi[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonPhi(),eventdiffZ->GetPatMuon2Phi());
  m_hVector_correlationMuon2Iso[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonrelIsoDr03(),eventdiffZ->GetPatMuon2relIsoDr03());

  for (l=0; l<16;l++){
    if (CastorEnergySector[l] >=1.){
      m_hVector_ECastorSector[index].at(pileup)->Fill(l+1,CastorEnergySector[l]);
      m_hVector_ECastorSectorTProf[index].at(pileup)->Fill(l+1,CastorEnergySector[l]);
      m_hVector_ECastorSectorBin1D[index].at(pileup)->Fill(l+1,CastorEnergySector[l]);
    }
    else{
      m_hVector_ECastorSector[index].at(pileup)->Fill(l+1,0);
    }
  }

  m_hVector_pfetamin[index].at(pileup)->Fill(eventdiff->GetEtaMinFromPFCands(),totalweight);
  m_hVector_pfetamincastor[index].at(pileup)->Fill(etamin_,totalweight);
  m_hVector_pfetamax[index].at(pileup)->Fill(eventdiff->GetEtaMaxFromPFCands(),totalweight);

  m_hVector_SectorVsTotalCastorEnergy[index].at(pileup)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);
  m_hVector_SectorVsTotalCastorEnergyTProf[index].at(pileup)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);
  m_hVector_sumECastorMinus[index].at(pileup)->Fill(sumCastorEnergy,totalweight);
  m_hVector_sumECastorMinusLow[index].at(pileup)->Fill(sumCastorEnergy,totalweight);
  sumCastorAndHFMinusEnergy = sumCastorEnergy+eventdiff->GetSumEnergyHFMinus();
  m_hVector_sumECastorAndHFMinus[index].at(pileup)->Fill(sumCastorAndHFMinusEnergy,totalweight);
  m_hVector_etcalos_p[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_etcalos_n[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_EnergyHFMinusVsCastorTProf[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_EnergyHFPlusVsCastorTProf[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_CastorMultiplicity[index].at(pileup)->Fill(counterHit,totalweight);
  m_hVector_CastorMultiplicityVsLumi[index].at(pileup)->Fill(eventinfo->GetInstLumiBunch(),counterHit,totalweight);

  if (SectorCastorHit < 1) m_hVector_RunNumberZeroCastor[index].at(pileup)->Fill(eventdiff->GetRunNumber());
  if (SectorCastorHit > 15) m_hVector_RunNumberHighCastor[index].at(pileup)->Fill(eventdiff->GetRunNumber());

  m_hVector_RunNumber[index].at(pileup)->Fill(eventdiff->GetRunNumber());

  for (k=0; k<eventdiffZ->GetEachTowerCounter();k++){
    m_hVector_ECaloVsEta[index].at(pileup)->Fill(eventdiffZ->GetEachTowerEta(k),eventdiffZ->GetEachTowerEnergy(k),totalweight);
    m_hVector_ECaloVsEtaTProf[index].at(pileup)->Fill(eventdiffZ->GetEachTowerEta(k),eventdiffZ->GetEachTowerEnergy(k),totalweight);
    m_hVector_EnergyVsEtaBin1D[index].at(pileup)->Fill(eventdiffZ->GetEachTowerEta(k),eventdiffZ->GetEachTowerEnergy(k)*totalweight);
  }

  m_hVector_ECaloVsEta[index].at(pileup)->Fill(-6.,sumCastorEnergy,totalweight);
  m_hVector_ECaloVsEtaTProf[index].at(pileup)->Fill(-6.,sumCastorEnergy,totalweight);
  m_hVector_EnergyVsEtaBin1D[index].at(pileup)->Fill(-6.,sumCastorEnergy*totalweight);

  if (eventdiff->GetSumEnergyHFPlus() > eventdiff->GetSumEnergyHFMinus()){
    m_hVector_SumEHFMax[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),totalweight);
    m_hVector_SumEHFMin[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),totalweight);
  }

  else{
    m_hVector_SumEHFMax[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),totalweight);
    m_hVector_SumEHFMin[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),totalweight);
  }

  m_hVector_tracksOutMuonsCone03[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConeMuon03(),totalweight);  
  m_hVector_tracksOutElectronsCone03[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConeElectron03(),totalweight);   
  m_hVector_tracksOutpatMuonsCone03[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConepatMuon03(),totalweight);
  m_hVector_tracksOutpatElectronsCone03[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConepatElectron03(),totalweight);

  m_hVector_tracksOutMuonsCone04[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConeMuon04(),totalweight);
  m_hVector_tracksOutElectronsCone04[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConeElectron04(),totalweight);
  m_hVector_tracksOutpatMuonsCone04[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConepatMuon04(),totalweight);
  m_hVector_tracksOutpatElectronsCone04[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConepatElectron04(),totalweight);

  m_hVector_tracksOutMuonsCone05[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConeMuon05(),totalweight);
  m_hVector_tracksOutElectronsCone05[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConeElectron05(),totalweight);
  m_hVector_tracksOutpatMuonsCone05[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConepatMuon05(),totalweight);
  m_hVector_tracksOutpatElectronsCone05[index].at(pileup)->Fill(eventdiffZ->GetTracksNonConepatElectron05(),totalweight);

  m_hVector_LeadingElectronInnerHits[index].at(pileup)->Fill(innerHits1,totalweight);
  m_hVector_LeadingElectronDCot[index].at(pileup)->Fill(Dcot1,totalweight);
  m_hVector_LeadingElectronDist[index].at(pileup)->Fill(Dist1,totalweight);
  m_hVector_LeadingElectronDeltaEtaTkClu[index].at(pileup)->Fill(DeltaEtaTkClu1,totalweight);
  m_hVector_LeadingElectronDeltaPhiTkClu[index].at(pileup)->Fill(DeltaPhiTkClu1,totalweight);
  m_hVector_LeadingElectronSigmaIeIe[index].at(pileup)->Fill(sigmaIeIe1,totalweight);
  m_hVector_LeadingElectronHE[index].at(pileup)->Fill(HE1,totalweight);
  m_hVector_SecondElectronInnerHits[index].at(pileup)->Fill(innerHits2,totalweight);
  m_hVector_SecondElectronDCot[index].at(pileup)->Fill(Dcot2,totalweight);
  m_hVector_SecondElectronDist[index].at(pileup)->Fill(Dist2,totalweight);
  m_hVector_SecondElectronDeltaEtaTkClu[index].at(pileup)->Fill(DeltaEtaTkClu2,totalweight);
  m_hVector_SecondElectronDeltaPhiTkClu[index].at(pileup)->Fill(DeltaPhiTkClu2,totalweight);
  m_hVector_SecondElectronSigmaIeIe[index].at(pileup)->Fill(sigmaIeIe2,totalweight);
  m_hVector_SecondElectronHE[index].at(pileup)->Fill(HE2,totalweight);
  m_hVector_SumEHFplus_S[index].at(pileup)->Fill(eventdiffZ->GetSumEHF_SPlus(),totalweight);
  m_hVector_SumEHFminus_S[index].at(pileup)->Fill(eventdiffZ->GetSumEHF_SMinus(),totalweight);
  m_hVector_SumEHFplus_L[index].at(pileup)->Fill(eventdiffZ->GetSumEHF_LPlus(),totalweight);
  m_hVector_SumEHFminus_L[index].at(pileup)->Fill(eventdiffZ->GetSumEHF_LMinus(),totalweight);
  m_hVector_maxetagap[index].at(pileup)->Fill(fabs(eventdiffZ->GetMaxGapPF()),totalweight);
  m_hVector_LimPlusgap[index].at(pileup)->Fill(eventdiffZ->GetLimPlusGapPF(),totalweight);
  m_hVector_LimMinusgap[index].at(pileup)->Fill(eventdiffZ->GetLimMinusGapPF(),totalweight);
  m_hVector_SumPTLimPlusgap[index].at(pileup)->Fill(eventdiffZ->GetPTMaxGapMaxPF(),totalweight);
  m_hVector_SumPTLimMinusgap[index].at(pileup)->Fill(eventdiffZ->GetPTMinGapMaxPF(),totalweight);

  m_hVector_ElectronsPt[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronPt(),totalweight);
  m_hVector_ElectronsEta[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronEta(),totalweight);
  m_hVector_ElectronsPhi[index].at(pileup)->Fill(eventdiffZ->GetLeadingElectronPhi(),totalweight);
  m_hVector_MuonsPt[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonPt(),totalweight);
  m_hVector_MuonsEta[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonEta(),totalweight);
  m_hVector_MuonsPhi[index].at(pileup)->Fill(eventdiffZ->GetLeadingMuonPhi(),totalweight);
  m_hVector_patElectronsPt[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Pt(),totalweight);
  m_hVector_patElectronsEta[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Eta(),totalweight);
  m_hVector_patElectronsPhi[index].at(pileup)->Fill(eventdiffZ->GetPatElectron1Phi(),totalweight);
  m_hVector_patMuonsPt[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Pt(),totalweight);
  m_hVector_patMuonsEta[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Eta(),totalweight);
  m_hVector_patMuonsPhi[index].at(pileup)->Fill(eventdiffZ->GetPatMuon1Phi(),totalweight);
  m_hVector_ElectronsPt[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronPt(),totalweight);
  m_hVector_ElectronsEta[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronEta(),totalweight);
  m_hVector_ElectronsPhi[index].at(pileup)->Fill(eventdiffZ->GetSecondElectronPhi(),totalweight);
  m_hVector_MuonsPt[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonPt(),totalweight);
  m_hVector_MuonsEta[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonEta(),totalweight);
  m_hVector_MuonsPhi[index].at(pileup)->Fill(eventdiffZ->GetSecondMuonPhi(),totalweight);
  m_hVector_patElectronsPt[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Pt(),totalweight);
  m_hVector_patElectronsEta[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Eta(),totalweight);
  m_hVector_patElectronsPhi[index].at(pileup)->Fill(eventdiffZ->GetPatElectron2Phi(),totalweight);
  m_hVector_patMuonsPt[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Pt(),totalweight);
  m_hVector_patMuonsEta[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Eta(),totalweight);
  m_hVector_patMuonsPhi[index].at(pileup)->Fill(eventdiffZ->GetPatMuon2Phi(),totalweight);

  m_hVector_absdeltaEtaPF[index].at(pileup)->Fill(absdeltaetapf,totalweight);
  m_hVector_deltaEtaPF[index].at(pileup)->Fill(deltaetapf,totalweight);

  m_hVector_XiPlusPF[index].at(pileup)->Fill(eventdiff->GetXiPlusFromPFCands(),totalweight);
  m_hVector_XiMinusPF[index].at(pileup)->Fill(eventdiff->GetXiMinusFromPFCands(),totalweight);
  m_hVector_XiPF[index].at(pileup)->Fill(eventdiff->GetXiPlusFromPFCands(),totalweight);
  m_hVector_XiPF[index].at(pileup)->Fill(eventdiff->GetXiMinusFromPFCands(),totalweight);

  m_hVector_AEcastor[index].at(pileup)->Fill(AEcastor,totalweight);
  m_hVector_etasignedHF[index].at(pileup)->Fill(etasignedHF,totalweight);
  m_hVector_etasignedCASTOR[index].at(pileup)->Fill(etasignedCASTOR,totalweight);

  m_hVector_sumEEEminus[index].at(pileup)->Fill(eventdiffZ->GetSumEEEMinus(),totalweight);
  m_hVector_sumEEEplus[index].at(pileup)->Fill(eventdiffZ->GetSumEEEPlus(),totalweight);

  m_hVector_EnergyHFPlusVsEnergyHFMinus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),eventdiff->GetSumEnergyHFMinus(),totalweight);
  m_hVector_EnergyEEPlusVsEnergyEEMinus[index].at(pileup)->Fill(eventdiffZ->GetSumEEEPlus(),eventdiffZ->GetSumEEEMinus(),totalweight);

  if (typesel == "RecoMuon"){
    m_hVector_BosonZPt[index].at(pileup)->Fill(eventdiffZ->GetDiMuonPt(),totalweight);
    m_hVector_BosonZEta[index].at(pileup)->Fill(eventdiffZ->GetDiMuonEta(),totalweight);
    m_hVector_BosonZPhi[index].at(pileup)->Fill(eventdiffZ->GetDiMuonPhi(),totalweight);
    m_hVector_BosonZMass[index].at(pileup)->Fill(eventdiffZ->GetDiMuonMass(),totalweight);
  }
  if (typesel == "PatMuon"){
    m_hVector_patBosonZPt[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonPt(),totalweight);
    m_hVector_patBosonZEta[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonEta(),totalweight);
    m_hVector_patBosonZPhi[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonPhi(),totalweight);
    m_hVector_patBosonZMass[index].at(pileup)->Fill(eventdiffZ->GetPatDiMuonMass(),totalweight);
  }
  if (typesel == "RecoElectron"){
    m_hVector_BosonZPt[index].at(pileup)->Fill(eventdiffZ->GetDiElectronPt(),totalweight);
    m_hVector_BosonZEta[index].at(pileup)->Fill(eventdiffZ->GetDiElectronEta(),totalweight);
    m_hVector_BosonZPhi[index].at(pileup)->Fill(eventdiffZ->GetDiElectronPhi(),totalweight);
    m_hVector_BosonZMass[index].at(pileup)->Fill(eventdiffZ->GetDiElectronMass(),totalweight);
  }
  if (typesel == "PatElectron"){
    m_hVector_patBosonZPt[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronPt(),totalweight);
    m_hVector_patBosonZEta[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronEta(),totalweight);
    m_hVector_patBosonZPhi[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronPhi(),totalweight);
    m_hVector_patBosonZMass[index].at(pileup)->Fill(eventdiffZ->GetPatDiElectronMass(),totalweight);
  }

}

void DiffractiveZ::SaveHistos(std::string type,std::string typesel){

  // Creating Correlation Histograms
  bool correlation = false;

  int ipileup;

  if (type=="multiple_pileup") ipileup=21;
  else ipileup=1;

  for (int i = 0; i < ipileup; i++){
    for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

      if(typesel=="RecoMuon" || typesel=="RecoElectron"){
	// Lepton Kinematics Folder
	foldersFile[0]->cd();
	m_hVector_DiElectronMass[j].at(i)->Write();
	m_hVector_LeadingElectronPt[j].at(i)->Write();
	m_hVector_LeadingElectronEta[j].at(i)->Write();
	m_hVector_LeadingElectronPhi[j].at(i)->Write();
	m_hVector_LeadingElectronCharge[j].at(i)->Write();
	m_hVector_SecondElectronPt[j].at(i)->Write();
	m_hVector_SecondElectronEta[j].at(i)->Write();
	m_hVector_SecondElectronPhi[j].at(i)->Write();
	m_hVector_SecondElectronCharge[j].at(i)->Write();
	m_hVector_ElectronsN[j].at(i)->Write();
	m_hVector_DiMuonMass[j].at(i)->Write();
	m_hVector_LeadingMuonPt[j].at(i)->Write();
	m_hVector_LeadingMuonEta[j].at(i)->Write();
	m_hVector_LeadingMuonPhi[j].at(i)->Write();
	m_hVector_LeadingMuonCharge[j].at(i)->Write();
	m_hVector_SecondMuonPt[j].at(i)->Write();
	m_hVector_SecondMuonEta[j].at(i)->Write();
	m_hVector_SecondMuonPhi[j].at(i)->Write();
	m_hVector_SecondMuonCharge[j].at(i)->Write();
	m_hVector_MuonsN[j].at(i)->Write();
	m_hVector_LeadingElectronTkDr03[j].at(i)->Write();
	m_hVector_LeadingElectronTkDr04[j].at(i)->Write();
	m_hVector_LeadingElectronEcalDr03[j].at(i)->Write();
	m_hVector_LeadingElectronEcalDr04[j].at(i)->Write();
	m_hVector_LeadingElectronHcalDr03[j].at(i)->Write();
	m_hVector_LeadingElectronHcalDr04[j].at(i)->Write();
	m_hVector_SecondElectronTkDr03[j].at(i)->Write();
	m_hVector_SecondElectronTkDr04[j].at(i)->Write();
	m_hVector_SecondElectronEcalDr03[j].at(i)->Write();
	m_hVector_SecondElectronEcalDr04[j].at(i)->Write();
	m_hVector_SecondElectronHcalDr03[j].at(i)->Write();
	m_hVector_SecondElectronHcalDr04[j].at(i)->Write();
	m_hVector_LeadingElectronrelIsoDr03[j].at(i)->Write();
	m_hVector_LeadingElectronrelIsoDr04[j].at(i)->Write();
	m_hVector_SecondElectronrelIsoDr03[j].at(i)->Write();
	m_hVector_SecondElectronrelIsoDr04[j].at(i)->Write();
	m_hVector_LeadingMuonSumPtR03[j].at(i)->Write();
	m_hVector_LeadingMuonSumPtR05[j].at(i)->Write();
	m_hVector_LeadingMuonEmEtR03[j].at(i)->Write();
	m_hVector_LeadingMuonEmEtR05[j].at(i)->Write();
	m_hVector_LeadingMuonHadEtR03[j].at(i)->Write();
	m_hVector_LeadingMuonHadEtR05[j].at(i)->Write();
	m_hVector_SecondMuonSumPtR03[j].at(i)->Write();
	m_hVector_SecondMuonSumPtR05[j].at(i)->Write();
	m_hVector_SecondMuonEmEtR03[j].at(i)->Write();
	m_hVector_SecondMuonEmEtR05[j].at(i)->Write();
	m_hVector_SecondMuonHadEtR03[j].at(i)->Write();
	m_hVector_SecondMuonHadEtR05[j].at(i)->Write();
	m_hVector_LeadingMuonrelIsoDr03[j].at(i)->Write();
	m_hVector_LeadingMuonrelIsoDr05[j].at(i)->Write();
	m_hVector_SecondMuonrelIsoDr03[j].at(i)->Write();
	m_hVector_SecondMuonrelIsoDr05[j].at(i)->Write();
	m_hVector_LeadingElectronInnerHits[j].at(i)->Write();
	m_hVector_LeadingElectronDCot[j].at(i)->Write();
	m_hVector_LeadingElectronDist[j].at(i)->Write();
	m_hVector_LeadingElectronDeltaEtaTkClu[j].at(i)->Write();
	m_hVector_LeadingElectronDeltaPhiTkClu[j].at(i)->Write();
	m_hVector_LeadingElectronSigmaIeIe[j].at(i)->Write();
	m_hVector_LeadingElectronHE[j].at(i)->Write();
	m_hVector_SecondElectronInnerHits[j].at(i)->Write();
	m_hVector_SecondElectronDCot[j].at(i)->Write();
	m_hVector_SecondElectronDist[j].at(i)->Write();
	m_hVector_SecondElectronDeltaEtaTkClu[j].at(i)->Write();
	m_hVector_SecondElectronDeltaPhiTkClu[j].at(i)->Write();
	m_hVector_SecondElectronSigmaIeIe[j].at(i)->Write();
	m_hVector_SecondElectronHE[j].at(i)->Write();
	m_hVector_DiElectronEta[j].at(i)->Write();
	m_hVector_DiElectronPhi[j].at(i)->Write();
	m_hVector_DiElectronPt[j].at(i)->Write();
	m_hVector_DiMuonEta[j].at(i)->Write();
	m_hVector_DiMuonPhi[j].at(i)->Write();
	m_hVector_DiMuonPt[j].at(i)->Write();
	m_hVector_ElectronsPt[j].at(i)->Write();
	m_hVector_ElectronsEta[j].at(i)->Write();
	m_hVector_ElectronsPhi[j].at(i)->Write();
	m_hVector_MuonsPt[j].at(i)->Write();
	m_hVector_MuonsEta[j].at(i)->Write();
	m_hVector_MuonsPhi[j].at(i)->Write();
	m_hVector_BosonZPt[j].at(i)->Write();
	m_hVector_BosonZEta[j].at(i)->Write();
	m_hVector_BosonZPhi[j].at(i)->Write();
	m_hVector_BosonZMass[j].at(i)->Write();
	m_hVector_tracksOutMuonsCone03[j].at(i)->Write();
	m_hVector_tracksOutElectronsCone03[j].at(i)->Write();
	m_hVector_tracksOutMuonsCone04[j].at(i)->Write();
	m_hVector_tracksOutElectronsCone04[j].at(i)->Write();
	m_hVector_tracksOutMuonsCone05[j].at(i)->Write();
	m_hVector_tracksOutElectronsCone05[j].at(i)->Write();
	m_hVector_deltaphielectrons[j].at(i)->Write();
	m_hVector_deltaphimuons[j].at(i)->Write();
	m_hVector_deltaetaelectrons[j].at(i)->Write();
	m_hVector_deltaetamuons[j].at(i)->Write();
	m_hVector_deltapTelectrons[j].at(i)->Write();
	m_hVector_deltapTmuons[j].at(i)->Write();
      }

      if(typesel=="PatElectron" || typesel=="PatMuon"){
	// Lepton Kinematics Folder
	foldersFile[0]->cd();
	m_hVector_patNElectron[j].at(i)->Write();
	m_hVector_patElectron1Pt[j].at(i)->Write();
	m_hVector_patElectron1Eta[j].at(i)->Write();
	m_hVector_patElectron1Phi[j].at(i)->Write();
	m_hVector_patElectron1Charge[j].at(i)->Write();
	m_hVector_patElectron1Et[j].at(i)->Write();
	m_hVector_patNMuon[j].at(i)->Write();
	m_hVector_patMuon1Pt[j].at(i)->Write();
	m_hVector_patMuon1Eta[j].at(i)->Write();
	m_hVector_patMuon1Phi[j].at(i)->Write();
	m_hVector_patMuon1Charge[j].at(i)->Write();
	m_hVector_patMuon1Et[j].at(i)->Write();
	m_hVector_patElectron2Pt[j].at(i)->Write();
	m_hVector_patElectron2Eta[j].at(i)->Write();
	m_hVector_patElectron2Phi[j].at(i)->Write();
	m_hVector_patElectron2Charge[j].at(i)->Write();
	m_hVector_patElectron2Et[j].at(i)->Write();
	m_hVector_patMuon2Pt[j].at(i)->Write();
	m_hVector_patMuon2Eta[j].at(i)->Write();
	m_hVector_patMuon2Phi[j].at(i)->Write();
	m_hVector_patMuon2Charge[j].at(i)->Write();
	m_hVector_patMuon2Et[j].at(i)->Write();
	m_hVector_patDiElectronMass[j].at(i)->Write();
	m_hVector_patDiMuonMass[j].at(i)->Write();
	m_hVector_patElectron1TkDr03[j].at(i)->Write();
	m_hVector_patElectron1TkDr04[j].at(i)->Write();
	m_hVector_patElectron1EcalDr03[j].at(i)->Write();
	m_hVector_patElectron1EcalDr04[j].at(i)->Write();
	m_hVector_patElectron1HcalDr03[j].at(i)->Write();
	m_hVector_patElectron1HcalDr04[j].at(i)->Write();
	m_hVector_patElectron2TkDr03[j].at(i)->Write();
	m_hVector_patElectron2TkDr04[j].at(i)->Write();
	m_hVector_patElectron2EcalDr03[j].at(i)->Write();
	m_hVector_patElectron2EcalDr04[j].at(i)->Write();
	m_hVector_patElectron2HcalDr03[j].at(i)->Write();
	m_hVector_patElectron2HcalDr04[j].at(i)->Write();
	m_hVector_patElectron1relIsoDr03[j].at(i)->Write();
	m_hVector_patElectron1relIsoDr04[j].at(i)->Write();
	m_hVector_patElectron2relIsoDr03[j].at(i)->Write();
	m_hVector_patElectron2relIsoDr04[j].at(i)->Write();
	m_hVector_patMuon1SumPtR03[j].at(i)->Write();
	m_hVector_patMuon1SumPtR05[j].at(i)->Write();
	m_hVector_patMuon1EmEtR03[j].at(i)->Write();
	m_hVector_patMuon1EmEtR05[j].at(i)->Write();
	m_hVector_patMuon1HadEtR03[j].at(i)->Write();
	m_hVector_patMuon1HadEtR05[j].at(i)->Write();
	m_hVector_patMuon2SumPtR03[j].at(i)->Write();
	m_hVector_patMuon2SumPtR05[j].at(i)->Write();
	m_hVector_patMuon2EmEtR03[j].at(i)->Write();
	m_hVector_patMuon2EmEtR05[j].at(i)->Write();
	m_hVector_patMuon2HadEtR03[j].at(i)->Write();
	m_hVector_patMuon2HadEtR05[j].at(i)->Write();
	m_hVector_patMuon1relIsoDr03[j].at(i)->Write();
	m_hVector_patMuon1relIsoDr05[j].at(i)->Write();
	m_hVector_patMuon2relIsoDr03[j].at(i)->Write();
	m_hVector_patMuon2relIsoDr05[j].at(i)->Write();
	m_hVector_patMuon1relIso[j].at(i)->Write();
	m_hVector_patMuon2relIso[j].at(i)->Write();
	m_hVector_patDiElectronEta[j].at(i)->Write();
	m_hVector_patDiMuonEta[j].at(i)->Write();
	m_hVector_patDiElectronPhi[j].at(i)->Write();
	m_hVector_patDiMuonPhi[j].at(i)->Write();
	m_hVector_patDiElectronPt[j].at(i)->Write();
	m_hVector_patDiMuonPt[j].at(i)->Write();
	m_hVector_patBosonZPt[j].at(i)->Write();
	m_hVector_patBosonZEta[j].at(i)->Write();
	m_hVector_patBosonZPhi[j].at(i)->Write();
	m_hVector_patBosonZMass[j].at(i)->Write();
	m_hVector_patElectronsPt[j].at(i)->Write();
	m_hVector_patElectronsEta[j].at(i)->Write();
	m_hVector_patElectronsPhi[j].at(i)->Write();
	m_hVector_patMuonsPt[j].at(i)->Write();
	m_hVector_patMuonsEta[j].at(i)->Write();
	m_hVector_patMuonsPhi[j].at(i)->Write();
	m_hVector_tracksOutpatMuonsCone03[j].at(i)->Write();
	m_hVector_tracksOutpatElectronsCone03[j].at(i)->Write();
	m_hVector_tracksOutpatMuonsCone04[j].at(i)->Write();
	m_hVector_tracksOutpatElectronsCone04[j].at(i)->Write();
	m_hVector_tracksOutpatMuonsCone05[j].at(i)->Write();
	m_hVector_tracksOutpatElectronsCone05[j].at(i)->Write();
      }

      if(correlation){
	// Detector Folder
	foldersFile[1]->cd();
	m_hVector_correlationElectron1Pt[j].at(i)->Write();
	m_hVector_correlationElectron1Eta[j].at(i)->Write();
	m_hVector_correlationElectron1Phi[j].at(i)->Write();
	m_hVector_correlationElectron1Iso[j].at(i)->Write();
	m_hVector_correlationMuon1Pt[j].at(i)->Write();
	m_hVector_correlationMuon1Eta[j].at(i)->Write();
	m_hVector_correlationMuon1Phi[j].at(i)->Write();
	m_hVector_correlationMuon1Iso[j].at(i)->Write();
	m_hVector_correlationElectron2Pt[j].at(i)->Write();
	m_hVector_correlationElectron2Eta[j].at(i)->Write();
	m_hVector_correlationElectron2Phi[j].at(i)->Write();
	m_hVector_correlationElectron2Iso[j].at(i)->Write();
	m_hVector_correlationMuon2Pt[j].at(i)->Write();
	m_hVector_correlationMuon2Eta[j].at(i)->Write();
	m_hVector_correlationMuon2Phi[j].at(i)->Write();
	m_hVector_correlationMuon2Iso[j].at(i)->Write();
      }

      // Detector Folder
      foldersFile[1]->cd();
      m_hVector_sumEHFplus[j].at(i)->Write();
      m_hVector_sumEHFminus[j].at(i)->Write();
      m_hVector_sumEHEplus[j].at(i)->Write();
      m_hVector_sumEHEminus[j].at(i)->Write();
      m_hVector_multhf[j].at(i)->Write();
      m_hVector_etcalos_p[j].at(i)->Write();
      m_hVector_etcalos_n[j].at(i)->Write();
      m_hVector_tracks[j].at(i)->Write();
      m_hVector_tracksLow[j].at(i)->Write();
      m_hVector_vertex[j].at(i)->Write();
      m_hVector_ECaloVsEta[j].at(i)->Write();
      m_hVector_ECaloVsEtaTProf[j].at(i)->Write();
      m_hVector_sumECastorMinus[j].at(i)->Write();
      m_hVector_sumECastorMinusLow[j].at(i)->Write();
      m_hVector_ECastorSector[j].at(i)->Write();
      m_hVector_ECastorSectorTProf[j].at(i)->Write();
      m_hVector_SumEHFplus_S[j].at(i)->Write();
      m_hVector_SumEHFminus_S[j].at(i)->Write();
      m_hVector_SumEHFplus_L[j].at(i)->Write();
      m_hVector_SumEHFminus_L[j].at(i)->Write();
      m_hVector_SumEHFMax[j].at(i)->Write();
      m_hVector_SumEHFMin[j].at(i)->Write();
      m_hVector_EnergyVsEtaBin1D[j].at(i)->Write();
      m_hVector_ECastorSectorBin1D[j].at(i)->Write();
      m_hVector_sumECastorAndHFMinus[j].at(i)->Write();
      m_hVector_CastorMultiplicity[j].at(i)->Write();
      m_hVector_CastorMultiplicityVsLumi[j].at(i)->Write();
      m_hVector_SectorVsTotalCastorEnergy[j].at(i)->Write();
      m_hVector_SectorVsTotalCastorEnergyTProf[j].at(i)->Write();
      m_hVector_sumEEEminus[j].at(i)->Write();
      m_hVector_sumEEEplus[j].at(i)->Write();
      m_hVector_EnergyHFPlusVsEnergyHFMinus[j].at(i)->Write();
      m_hVector_EnergyEEPlusVsEnergyEEMinus[j].at(i)->Write();

      // EventInfo Folder
      foldersFile[2]->cd();
      m_hVector_lumi[j].at(i)->Write();
      m_hVector_vertexvslumi[j].at(i)->Write();
      m_hVector_RunNumberHighCastor[j].at(i)->Write();
      m_hVector_RunNumberZeroCastor[j].at(i)->Write();
      m_hVector_RunNumber[j].at(i)->Write();

      // Diffractive Variables Folder
      foldersFile[3]->cd();
      m_hVector_pfetamax[j].at(i)->Write();
      m_hVector_pfetamin[j].at(i)->Write();
      m_hVector_asumE[j].at(i)->Write();
      m_hVector_maxetagap[j].at(i)->Write();
      m_hVector_LimPlusgap[j].at(i)->Write();
      m_hVector_LimMinusgap[j].at(i)->Write();
      m_hVector_SumPTLimPlusgap[j].at(i)->Write();
      m_hVector_SumPTLimMinusgap[j].at(i)->Write();
      m_hVector_absdeltaEtaPF[j].at(i)->Write();
      m_hVector_deltaEtaPF[j].at(i)->Write();
      m_hVector_XiMinusPF[j].at(i)->Write();
      m_hVector_XiPlusPF[j].at(i)->Write();
      m_hVector_XiPF[j].at(i)->Write();
      m_hVector_pfetamincastor[j].at(i)->Write();
      m_hVector_AEcastor[j].at(i)->Write();
      m_hVector_etasignedHF[j].at(i)->Write();
      m_hVector_etasignedCASTOR[j].at(i)->Write();
    }
  }

}

void DiffractiveZ::Run(std::string filein_, std::string processname_, std::string savehistofile_, std::string switchtrigger_, int optTrigger_, double lepton1pt_, double lepton2pt_, int nVertex_, std::string type_, std::string switchlumiweight_, double mcweight_, std::string typesel_, double castorthreshold_, double channelsthreshold_, std::string castorcorrfile_, std::string gapseltype_, std::string pumfile_, std::string pudfile_){

  bool debug = false;

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  type = type_;
  switchtrigger = switchtrigger_;
  switchlumiweight = switchlumiweight_;
  mcweight = mcweight_;
  nVertex = nVertex_;
  optTrigger = optTrigger_;
  lepton1pt = lepton1pt_;
  lepton2pt = lepton2pt_;
  typesel = typesel_;
  castorthreshold = castorthreshold_;
  channelsthreshold = channelsthreshold_;
  castorcorrfile = castorcorrfile_;
  gapseltype = gapseltype_;
  pumfile = pumfile_;
  pudfile = pudfile_;

  std::string selStatus;
  std::string TriggerStatus;
  char selCastor[300];

  // Adding TTree Golden Events
  TString TTreeoutput, TTreeAllZ, TTreeCASTOR;
  TTreeoutput = "TTreeGoldenDiffZ_" + savehistofile;
  fOut = new TFile(TTreeoutput, "RECREATE");
  fOut->cd();

  trout = new TTree("Events", "Events");
  trout->Branch("RunNumber",&bRunNumber,"bRunNumber/I");
  trout->Branch("LumiSection",&bLumiSection,"bLumiSection/I");
  trout->Branch("EventNumber",&bEventNumber,"bEventNumber/I");
  trout->Branch("DiBosonPt",&bDiBosonPt,"bDiBosonPt/D");
  trout->Branch("DiBosonEta",&bDiBosonEta,"bDiBosonEta/D");
  trout->Branch("DiBosonPhi",&bDiBosonPhi,"bDiBosonPhi/D");
  trout->Branch("DiBosonMass",&bDiBosonMass,"bDiBosonMass/D");
  trout->Branch("MultiplicityTracks",&bMultiplicityTracks,"bMultiplicityTracks/I");
  trout->Branch("SumEEEMinus",&bSumEEEMinus,"bSumEEEMinus/D");
  trout->Branch("SumEEEPlus",&bSumEEEPlus,"bSumEEEPlus/D");
  trout->Branch("SumEnergyHFMinus",&bSumEnergyHFMinus,"bSumEnergyHFMinus/D");
  trout->Branch("SumEnergyHFPlus",&bSumEnergyHFPlus,"bSumEnergyHFPlus/D");
  trout->Branch("sumCastorEnergy",&bsumCastorEnergy,"bsumCastorEnergy/D");
  trout->Branch("SectorCastorHit",&bSectorCastorHit,"bSectorCastorHit/D");
  trout->Branch("deltaetapf",&bdeltaetapf,"bdeltaetapf/D");
  trout->Branch("AEcastor",&bAEcastor,"AEcastor/D");
  trout->Branch("etasignedHF",&betasignedHF,"betasignedHF/D");
  trout->Branch("etasignedCASTOR",&betasignedCASTOR,"betasignedCASTOR/D");
  trout->Branch("MaxGapPF",&bMaxGapPF,"bMaxGapPF/D");
  trout->Branch("PTMaxGapMaxPF",&bPTMaxGapMaxPF,"bPTMaxGapMaxPF/D");
  trout->Branch("PTMinGapMaxPF",&bPTMinGapMaxPF,"bPTMinGapMaxPF/D");
  trout->Branch("XiPlusFromPFCands",&bXiPlusFromPFCands,"bXiPlusFromPFCands/D");
  trout->Branch("XiMinusFromPFCands",&bXiMinusFromPFCands,"bXiMinusFromPFCands/D");
  trout->Branch("EtaMaxPF",&betamax,"betamax/D");
  trout->Branch("EtaMinPF",&betamin,"betamin/D");
  trout->Branch("EtaLimMinus",&betalimmin,"betalimmin/D");
  trout->Branch("EtaLimPlus",&betalimmax,"betalimmax/D");

  TTreeCASTOR = "TTreeCASTOR_" + savehistofile;
  fOutCASTOR = new TFile(TTreeCASTOR, "RECREATE");
  fOutCASTOR->cd();
  troutCASTOR = trout->CloneTree(0);

  TTreeAllZ = "TTreeAllZ_" + savehistofile;
  fOutZ = new TFile(TTreeAllZ, "RECREATE");
  fOutZ->cd();
  troutZ = trout->CloneTree(0);

  TFile check1(filein.c_str());

  if (check1.IsZombie()){

    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " There is no the file " << filein << " or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;

  }

  if (check1.GetDirectory(processname.c_str())){
    LoadFile(filein,processname);
  }

  else{
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << " There is no directory/path " << processname << std::endl;
    std::cout << " in the file." << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    return;
  }

  TFile check2(castorcorrfile.c_str());

  TFile *outf = new TFile(savehistofile.c_str(),"RECREATE");
  TString outtxt = savehistofile;
  outtxt.ReplaceAll("root","txt");
  std::ofstream outstring(outtxt);
  outf->cd();

  outstring << "" << std::endl;
  outstring << "<< Gold Events >>" << std::endl;
  outstring << "" << std::endl;
  outstring << "Please, insert this events in another text file to be used by PickEvent Tool. " << std::endl;
  outstring << "Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents " << std::endl;
  outstring << ">>---------------------------------------------------------------------- " << std::endl;

  int NEVENTS = tr->GetEntries();
  int pileup = -999;
  int triggercounter[20]={0};
  int totalT=0;

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  unsigned NEntries = tr->GetEntries();
  std::cout << "" << std::endl;
  std::cout<< "Reading Tree: "<< NEntries << " events"<<std::endl;
  std::cout << "" << std::endl;

  std::string status;  

  double energycorr[5][16];
  h_castor_channel = (TH2F*)check2.Get("channelcorrector");

  // Get CASTOR Corrections
  for(int i=1; i<6;i++){
    for(int j=1; j<17; j++){
      energycorr[i-1][j-1]=0;
      if (switchtrigger == "no_trigger_correction"){
	energycorr[i-1][j-1]=h_castor_channel->GetBinContent(i,j);
      }else{
	energycorr[i-1][j-1]=1.;
      }
    }
  }

  // Lumiweight PU
  edm::LumiReWeighting LumiWeights_(pumfile.c_str(),pudfile.c_str(),"pileUpBx0_complete_without_cuts","pileup");

  for(int i=0;i<NEVENTS;i++){

    tr->GetEntry(i);

    if ( type=="multiple_pileup" && (eventinfo->GetNPileUpBx0()==-1 && eventinfo->GetNPileUpBxm1()==-1 && eventinfo->GetNPileUpBxp1()==-1 )){
      std::cout << " " << std::endl; 
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << " There is no pile-up TTree information in your PATTuplefile."   << std::endl;
      std::cout << " Please, use another PATTuple with PU information to run mul- " << std::endl;
      std::cout << " tiple PU option." << std::endl;
      std::cout << "--------------------------------------------------------------" << std::endl;
      return;
    }

    if (type=="multiple_pileup"){
      pileup = eventinfo->GetNPileUpBx0();
    }
    else{
      pileup = 0;
    }

    for (int nt=0;nt<20;nt++){
      if(eventdiffZ->GetHLTPath(nt) > 0){
	triggercounter[nt]++;
      }
    }

    etasignedHF = -999.;
    etasignedCASTOR = -999;
    aSumE = -999.;
    AEcastor = -999.;

    double totalASum = eventdiff->GetSumEnergyHFPlus() + eventdiff->GetSumEnergyHFMinus();
    if (totalASum > 0.){
      aSumE = (eventdiff->GetSumEnergyHFPlus() - eventdiff->GetSumEnergyHFMinus())/(eventdiff->GetSumEnergyHFPlus() + eventdiff->GetSumEnergyHFMinus());
    }else{
      aSumE = 999.;
    }

    deltaphielectrons = fabs(eventdiffZ->GetLeadingElectronPhi() - eventdiffZ->GetSecondElectronPhi());
    deltaphimuons = fabs(eventdiffZ->GetLeadingMuonPhi() - eventdiffZ->GetSecondMuonPhi());
    deltaetaelectrons = eventdiffZ->GetLeadingElectronEta() - eventdiffZ->GetSecondElectronEta();
    deltaetamuons = eventdiffZ->GetLeadingMuonEta() - eventdiffZ->GetSecondMuonEta();
    deltapTelectrons = eventdiffZ->GetLeadingElectronPt() - eventdiffZ->GetSecondElectronPt();
    deltapTmuons = eventdiffZ->GetLeadingMuonPt() - eventdiffZ->GetSecondMuonPt();

    double puweight;
    double totalcommon;

    if (switchlumiweight == "mc_lumi_weight"){
      totalcommon = mcweight;
    }else if (switchlumiweight == "mc_lumi_pu_weight"){
      puweight = LumiWeights_.weight(eventinfo->GetNPileUpBx0());
      totalcommon = mcweight*puweight;
    }else if (switchlumiweight == "no_weight"){
      totalcommon = 1.;
    }
    else{
      std::cout << "\nPlease, re-run with a correct event weight option.\n" << std::endl;
      return;
    }
    if (!debug){
      if (i==0) {
	std::cout << "" << std::endl;
	std::cout<< "Status Bar" << std::endl;
	std::cout << "" << std::endl;
      }
      loadBar(i,NEVENTS,100,100);
    }

    if (debug){
      std::cout << " " << std::endl;
      std::cout << "# of active towers: " << eventdiffZ->GetEachTowerCounter() << std::endl;
      std::cout << " " <<std::endl;
      for (int ic=0; ic<eventdiffZ->GetEachTowerCounter(); ic++){
	std::cout << "Each Tower Energy: " << eventdiffZ->GetEachTowerEnergy(ic) << std::endl;
      }
      std::cout << " " <<std::endl;
    }

    if (debug){
      if( i % 1000 == 0 ){
	std::cout << "\nEvent " << i << std::endl;
	std::cout << "Nr. events Bx 0: " << eventinfo->GetNPileUpBx0() << std::endl;
	std::cout << "Lumi per Bunch: " << eventinfo->GetInstLumiBunch() << std::endl;
	std::cout << "Event Weight: " << totalcommon <<std::endl;
      }
    }

    bool triggerE_a = false;
    bool triggerE_b = false;
    bool triggerE_c = false;
    bool triggerE_d = false;
    bool triggerE_e = false;
    bool triggerE_f = false;
    bool triggerE_g = false;
    bool trigger = false;
    bool vertex = false;
    bool diffseln = false;
    bool diffselp = false;

    bool presel = false;
    bool charge = false;
    bool dimass = false;
    bool nSel = false;
    bool eleBarrel1 = false;
    bool eleEndCap1 = false;
    bool eleBarrel2 = false;
    bool eleEndCap2 = false; 
    bool candSel = false;
    bool isoBarrel1 = false;
    bool isoEndCap1 = false;
    bool isoBarrel2 = false;
    bool isoEndCap2 = false;
    bool isolation = false;
    bool ZKinN = false;
    bool ZKinP = false;

    if (switchtrigger == "trigger_all_electron"){
      if ( (eventdiff->GetRunNumber() >= 132440 && eventdiff->GetRunNumber() <= 137028) && eventdiffZ->GetHLTPath(0) > 0 ) triggerE_a = true;
      if ( (eventdiff->GetRunNumber() >= 138564 && eventdiff->GetRunNumber() <= 140401) && eventdiffZ->GetHLTPath(1) > 0) triggerE_b = true;
      if ( (eventdiff->GetRunNumber() >= 141956 && eventdiff->GetRunNumber() <= 144114) && eventdiffZ->GetHLTPath(2) > 0) triggerE_c = true;
      if ( (eventdiff->GetRunNumber() >= 144115 && eventdiff->GetRunNumber() <= 147145) && eventdiffZ->GetHLTPath(3) > 0) triggerE_d = true;
      if ( (eventdiff->GetRunNumber() >= 147146 && eventdiff->GetRunNumber() <= 148058) && eventdiffZ->GetHLTPath(4) > 0) triggerE_e = true;
      if ( (eventdiff->GetRunNumber() >= 148103 && eventdiff->GetRunNumber() <= 149065) && eventdiffZ->GetHLTPath(5) > 0) triggerE_f = true;
      if ( (eventdiff->GetRunNumber() >= 149180 && eventdiff->GetRunNumber() <= 149442) && eventdiffZ->GetHLTPath(6) > 0) triggerE_g = true;
      if (triggerE_a || triggerE_b || triggerE_c || triggerE_d || triggerE_e || triggerE_f || triggerE_g) trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Electron accepted." << std::endl;
      TriggerStatus = "trigger_all_electron";
    }
    else if (switchtrigger == "trigger_all_muon"){
      if (eventdiffZ->GetHLTPath(0) > 0|| eventdiffZ->GetHLTPath(1) > 0) trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Muon accepted." << std::endl;
      TriggerStatus = "trigger_all_muon";
    }
    else if (switchtrigger == "trigger"){
      if (eventdiffZ->GetHLTPath(optTrigger) > 0) trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger <<", trigger accepted." << std::endl;
      TriggerStatus = "trigger";
    }
    else if (switchtrigger == "no_trigger_nocorrection") {
      trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger << ", no trigger." << std::endl; 
      TriggerStatus = "no_trigger_nocorrection";
    }
    else if (switchtrigger == "no_trigger_correction") {
      trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger << ", no trigger." << std::endl; 
      TriggerStatus = "no_trigger_correction";
    }
    else{
      exit(EXIT_FAILURE);
    }

    sumCastorEnergy = 0.;
    sumCastorAndHFMinusEnergy = 0.;
    SectorCastorHit = 0;
    counterHit = 0;
    SectorZeroCastorCounter = 0;
    castorgap = false;
    castoractivity = false;

    if (channelsthreshold > 0. && castorthreshold < 0.){
      sprintf(selCastor,">> Castor Channel Threshold: %0.2f GeV",channelsthreshold);
      for (int i=0; i < 16; i++){
	CastorEnergySector[i]=0.;
	if (i==4 || i==5){
	  if (eventdiffZ->GetCastorModule2Energy(i)*energycorr[1][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule2Energy(i)*energycorr[1][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule3Energy(i)*energycorr[2][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule3Energy(i)*energycorr[2][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule4Energy(i)*energycorr[3][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule4Energy(i)*energycorr[3][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule5Energy(i)*energycorr[4][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule5Energy(i)*energycorr[4][i]; ++counterHit;}
	}else{
	  if (eventdiffZ->GetCastorModule1Energy(i)*energycorr[0][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule1Energy(i)*energycorr[0][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule2Energy(i)*energycorr[1][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule2Energy(i)*energycorr[1][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule3Energy(i)*energycorr[2][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule3Energy(i)*energycorr[2][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule4Energy(i)*energycorr[3][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule4Energy(i)*energycorr[3][i]; ++counterHit;}
	  if (eventdiffZ->GetCastorModule5Energy(i)*energycorr[4][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffZ->GetCastorModule5Energy(i)*energycorr[4][i]; ++counterHit;}
	}
      }

      for (l=0; l<16;l++){
	sumCastorEnergy+=CastorEnergySector[l];
      }

      if (sumCastorEnergy>0.){
	++SectorCastorHit;
      }
      else{
	++SectorZeroCastorCounter;
      }

    }

    if (castorthreshold > 0. && channelsthreshold < 0. ) {

      sprintf(selCastor,">> Castor Sector Threshold: %0.2f GeV",castorthreshold);
      for (int i=0; i < 16; i++){
	CastorEnergySector[i]=0.;
	if (i==4 || i==5){
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule2Energy(i)*energycorr[1][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule3Energy(i)*energycorr[2][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule4Energy(i)*energycorr[3][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule5Energy(i)*energycorr[4][i];
	}else{
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule1Energy(i)*energycorr[0][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule2Energy(i)*energycorr[1][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule3Energy(i)*energycorr[2][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule4Energy(i)*energycorr[3][i];
	  CastorEnergySector[i]+=eventdiffZ->GetCastorModule5Energy(i)*energycorr[4][i];
	}
      }
      for (l=0; l<16;l++){
	if (CastorEnergySector[l] >= castorthreshold){
	  ++SectorCastorHit;
	  ++counterHit;
	  sumCastorEnergy+=CastorEnergySector[l];
	}
	else{
	  ++SectorZeroCastorCounter;
	}
      }
    } 

    if (SectorCastorHit > 0){
      castoractivity = true;
    }else{
      castorgap = true;
    }

    // Redefinition of etamin_
    //------------------------

    // It is possible to include CASTOR as the min. eta with activity.

    /*
    // CMS and CASTOR acceptance
    etamin_=0.;
    if (castoractivity) {
    etamin_ = -6.;
    }else{
    etamin_ = eventdiff->GetEtaMinFromPFCands();
    }
     */

    // Only CMS
    etamin_= eventdiff->GetEtaMinFromPFCands();

    //----->

    if (eventdiff->GetNVertex() == nVertex) vertex = true;

    if (gapseltype == "PF" || gapseltype == "pf"){
      if ((eventdiff->GetEtaMinFromPFCands() > -3.)) diffseln = true;
      if ((eventdiff->GetEtaMaxFromPFCands() <  3.)) diffselp = true;
    }
    else if(gapseltype == "HF" || gapseltype == "hf"){
      if ((eventdiff->GetSumEnergyHFMinus()<1.)) diffseln = true;
      if ((eventdiff->GetSumEnergyHFPlus()<1.)) diffselp = true;   
    }else{
      std::cout << "" << std::endl;
      std::cout << "Please insert the correct gap type selection: " << std::endl;
      std::cout << ">> HF or hf: for hadron forward calorimeter gap selection." << std::endl;
      std::cout << ">> PF or pf: for etamax, etamin gap selection." << std::endl;
      std::cout << "" << std::endl;
      exit(EXIT_FAILURE);
    }

    deltaetapf = eventdiff->GetEtaMaxFromPFCands() - etamin_;
    absdeltaetapf = fabs(eventdiff->GetEtaMaxFromPFCands() - etamin_);

    if (typesel == "RecoElectron"){
      selStatus = "Reco::Electron";
      isoTk1 = eventdiffZ->GetLeadingElectronTkDr03()/eventdiffZ->GetLeadingElectronPt();
      isoEcal1 = eventdiffZ->GetLeadingElectronEcalDr03()/eventdiffZ->GetLeadingElectronPt();
      isoHcal1 = eventdiffZ->GetLeadingElectronHcalDr03()/eventdiffZ->GetLeadingElectronPt();
      isoTk2 = eventdiffZ->GetSecondElectronTkDr03()/eventdiffZ->GetSecondElectronPt();
      isoEcal2 = eventdiffZ->GetSecondElectronEcalDr03()/eventdiffZ->GetSecondElectronPt();
      isoHcal2 = eventdiffZ->GetSecondElectronHcalDr03()/eventdiffZ->GetSecondElectronPt();
      innerHits1 = eventdiffZ->GetLeadingElectronInnerHits();
      Dcot1 = eventdiffZ->GetLeadingElectronDCot();
      Dist1 = eventdiffZ->GetLeadingElectronDist();
      DeltaEtaTkClu1 = eventdiffZ->GetLeadingElectronDeltaEtaTkClu();
      DeltaPhiTkClu1 = eventdiffZ->GetLeadingElectronDeltaPhiTkClu();
      sigmaIeIe1 = eventdiffZ->GetLeadingElectronSigmaIeIe();
      HE1 = eventdiffZ->GetLeadingElectronHE();
      innerHits2 = eventdiffZ->GetSecondElectronInnerHits();
      Dcot2 = eventdiffZ->GetSecondElectronDCot();
      Dist2 = eventdiffZ->GetSecondElectronDist();
      DeltaEtaTkClu2 = eventdiffZ->GetSecondElectronDeltaEtaTkClu();
      DeltaPhiTkClu2 = eventdiffZ->GetSecondElectronDeltaPhiTkClu();
      sigmaIeIe2 = eventdiffZ->GetSecondElectronSigmaIeIe();
      HE2 = eventdiffZ->GetSecondElectronHE();

      double totalASumCastor = eventdiff->GetSumEnergyHFMinus() + sumCastorEnergy;
      if(totalASumCastor > 0.){
	AEcastor = (eventdiff->GetSumEnergyHFMinus() - sumCastorEnergy)/(eventdiff->GetSumEnergyHFMinus() + sumCastorEnergy);
      }
      else{
	AEcastor = 999.;
      }

      if (eventdiffZ->GetLeadingElectronPt() > lepton1pt && eventdiffZ->GetSecondElectronPt() > lepton2pt) presel = true;
      if (eventdiffZ->GetLeadingElectronCharge()*eventdiffZ->GetSecondElectronCharge()==-1) charge = true;
      if (eventdiffZ->GetDiElectronMass() > 60. && eventdiffZ->GetDiElectronMass() < 110.) dimass = true;
      if (eventdiffZ->GetElectronsN() > 1) nSel = true;

      //Isolation Electron
      if ((fabs (eventdiffZ->GetLeadingElectronEta()) <= 1.4442) ){
	if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1 = true;
      }

      if ((fabs (eventdiffZ->GetLeadingElectronEta()) >= 1.5660) && (fabs (eventdiffZ->GetLeadingElectronEta()) <= 2.5)){
	if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1 = true;
      }

      if ((fabs (eventdiffZ->GetSecondElectronEta()) <= 1.4442) ){
	if (isoTk2<0.09 && isoEcal2<0.07 && isoHcal2<0.10) isoBarrel2 = true;
      }

      if ((fabs (eventdiffZ->GetSecondElectronEta()) >= 1.5660) && (fabs (eventdiffZ->GetSecondElectronEta()) <= 2.5)){
	if (isoTk2<0.04 && isoEcal2<0.05 && isoHcal2<0.025) isoEndCap2 = true;
      }

      if ((isoEndCap1 || isoBarrel1) && (isoEndCap2 || isoBarrel2)) isolation = true;

      //Electron Candidate Selection
      if ((fabs (eventdiffZ->GetLeadingElectronEta()) <= 1.4442) ){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.004 && fabs (DeltaPhiTkClu1) < 0.06 && sigmaIeIe1 < 0.01 && HE1 < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (eventdiffZ->GetLeadingElectronEta()) >= 1.5660) && (fabs (eventdiffZ->GetLeadingElectronEta()) <= 2.5)){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.007 && fabs (DeltaPhiTkClu1) < 0.03 && sigmaIeIe1 < 0.03 && HE1 < 0.025) eleEndCap1 = true;
      }

      if ((fabs (eventdiffZ->GetSecondElectronEta()) <= 1.4442) ){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.004 && fabs (DeltaPhiTkClu2) < 0.06 && sigmaIeIe2 < 0.01 && HE2 < 0.04 ) eleBarrel2 = true;
      }

      if ((fabs (eventdiffZ->GetSecondElectronEta()) >= 1.5660) && (fabs (eventdiffZ->GetSecondElectronEta()) <= 2.5)){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.007 && fabs (DeltaPhiTkClu2) < 0.03 && sigmaIeIe2 < 0.03 && HE2 < 0.025) eleEndCap2 = true;
      }

      if ((eleEndCap1 || eleBarrel1) && (eleEndCap2 || eleBarrel2)) candSel = true;

      if (eventdiffZ->GetDiElectronEta()>0.) ZKinP = true;
      if (eventdiffZ->GetDiElectronEta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffZ->GetDiElectronEta()); }
      else{ etasignedHF = fabs(eventdiffZ->GetDiElectronEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffZ->GetDiElectronEta());}
      else {etasignedCASTOR = fabs(eventdiffZ->GetDiElectronEta());}

    }

    else if (typesel == "RecoMuon"){
      selStatus = "Reco::Muon";
      if (eventdiffZ->GetLeadingMuonPt() > lepton1pt && eventdiffZ->GetSecondMuonPt() > lepton2pt) presel = true;
      if (eventdiffZ->GetLeadingMuonCharge()*eventdiffZ->GetSecondMuonCharge()==-1) charge = true;
      if (eventdiffZ->GetDiMuonMass() > 60. && eventdiffZ->GetDiMuonMass() < 110.) dimass = true;
      if (eventdiffZ->GetMuonsN() > 1) nSel = true;
      if (eventdiffZ->GetLeadingMuonSumPtR03() < 3 && eventdiffZ->GetSecondMuonSumPtR03() < 3 ) { 
	isolation = true;
	candSel = true;
      }

      if (eventdiffZ->GetDiMuonEta()>0.) ZKinP = true;
      if (eventdiffZ->GetDiMuonEta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffZ->GetDiMuonEta()); }
      else{ etasignedHF = fabs(eventdiffZ->GetDiMuonEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffZ->GetDiMuonEta());}
      else {etasignedCASTOR = fabs(eventdiffZ->GetDiMuonEta());}

    }

    else if (typesel == "PatElectron"){
      selStatus = "Pat::Electron";
      isoTk1 = eventdiffZ->GetPatElectron1TkDr03()/eventdiffZ->GetPatElectron1Pt();
      isoEcal1 = eventdiffZ->GetPatElectron1EcalDr03()/eventdiffZ->GetPatElectron1Pt();
      isoHcal1 = eventdiffZ->GetPatElectron1HcalDr03()/eventdiffZ->GetPatElectron1Pt();
      isoTk2 = eventdiffZ->GetPatElectron2TkDr03()/eventdiffZ->GetPatElectron2Pt();
      isoEcal2 = eventdiffZ->GetPatElectron2EcalDr03()/eventdiffZ->GetPatElectron2Pt();
      isoHcal2 = eventdiffZ->GetPatElectron2HcalDr03()/eventdiffZ->GetPatElectron2Pt();
      innerHits1 = eventdiffZ->GetPatElectron1InnerHits();
      Dcot1 = eventdiffZ->GetPatElectron1DCot();
      Dist1 = eventdiffZ->GetPatElectron1Dist();
      DeltaEtaTkClu1 = eventdiffZ->GetPatElectron1DeltaEtaTkClu();
      DeltaPhiTkClu1 = eventdiffZ->GetPatElectron1DeltaPhiTkClu();
      sigmaIeIe1 = eventdiffZ->GetPatElectron1SigmaIeIe();
      HE1 = eventdiffZ->GetPatElectron1HE();
      innerHits2 = eventdiffZ->GetPatElectron2InnerHits();
      Dcot2 = eventdiffZ->GetPatElectron2DCot();
      Dist2 = eventdiffZ->GetPatElectron2Dist();
      DeltaEtaTkClu2 = eventdiffZ->GetPatElectron2DeltaEtaTkClu();
      DeltaPhiTkClu2 = eventdiffZ->GetPatElectron2DeltaPhiTkClu();
      sigmaIeIe2 = eventdiffZ->GetPatElectron2SigmaIeIe();
      HE2 = eventdiffZ->GetPatElectron2HE();

      if (eventdiffZ->GetPatElectron1Pt() > lepton1pt && eventdiffZ->GetPatElectron2Pt() > lepton2pt) presel = true;
      if (eventdiffZ->GetPatElectron1Charge()*eventdiffZ->GetPatElectron2Charge()==-1) charge = true;
      if (eventdiffZ->GetPatDiElectronMass() > 60. && eventdiffZ->GetPatDiElectronMass() < 110.) dimass = true;
      if (eventdiffZ->GetPatNElectron() > 1) nSel = true;


      //Isolation Electron
      if ((fabs (eventdiffZ->GetPatElectron1Eta()) <= 1.4442) ){
	if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1 = true;
      }

      if ((fabs (eventdiffZ->GetPatElectron1Eta()) >= 1.5660) && (fabs (eventdiffZ->GetPatElectron1Eta()) <= 2.5)){
	if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1 = true;
      }

      if ((fabs (eventdiffZ->GetPatElectron2Eta()) <= 1.4442) ){
	if (isoTk2<0.09 && isoEcal2<0.07 && isoHcal2<0.10) isoBarrel2 = true;
      }

      if ((fabs (eventdiffZ->GetPatElectron2Eta()) >= 1.5660) && (fabs (eventdiffZ->GetPatElectron2Eta()) <= 2.5)){
	if (isoTk2<0.04 && isoEcal2<0.05 && isoHcal2<0.025) isoEndCap2 = true;
      }

      if ((isoEndCap1 || isoBarrel1) && (isoEndCap2 || isoBarrel2)) isolation = true;

      //Electron Candidate Selection
      if ((fabs (eventdiffZ->GetPatElectron1Eta()) <= 1.4442) ){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.004 && fabs (DeltaPhiTkClu1) < 0.06 && sigmaIeIe1 < 0.01 && HE1 < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (eventdiffZ->GetPatElectron1Eta()) >= 1.5660) && (fabs (eventdiffZ->GetPatElectron1Eta()) <= 2.5)){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.007 && fabs (DeltaPhiTkClu1) < 0.03 && sigmaIeIe1 < 0.03 && HE1 < 0.025) eleEndCap1 = true;
      }

      if ((fabs (eventdiffZ->GetPatElectron2Eta()) <= 1.4442) ){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.004 && fabs (DeltaPhiTkClu2) < 0.06 && sigmaIeIe2 < 0.01 && HE2 < 0.04 ) eleBarrel2 = true;
      }

      if ((fabs (eventdiffZ->GetPatElectron2Eta()) >= 1.5660) && (fabs (eventdiffZ->GetPatElectron2Eta()) <= 2.5)){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.007 && fabs (DeltaPhiTkClu2) < 0.03 && sigmaIeIe2 < 0.03 && HE2 < 0.025) eleEndCap2 = true;
      }

      if ((eleEndCap1 || eleBarrel1) && (eleEndCap2 || eleBarrel2)) candSel = true;

      if (eventdiffZ->GetPatDiElectronEta()>0.) ZKinP = true;
      if (eventdiffZ->GetPatDiElectronEta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffZ->GetPatDiElectronEta()); }
      else{ etasignedHF = fabs(eventdiffZ->GetPatDiElectronEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffZ->GetPatDiElectronEta());}
      else {etasignedCASTOR = fabs(eventdiffZ->GetPatDiElectronEta());}

    }

    else if(typesel == "PatMuon"){
      selStatus = "Pat::Muon";
      if (eventdiffZ->GetPatMuon1Pt() > lepton1pt && eventdiffZ->GetPatMuon2Pt() > lepton2pt) presel = true;
      if (eventdiffZ->GetPatMuon1Charge()*eventdiffZ->GetPatMuon2Charge()==-1) charge = true;
      if (eventdiffZ->GetPatDiMuonMass() > 60. && eventdiffZ->GetPatDiMuonMass() < 110.) dimass = true;
      if (eventdiffZ->GetPatNMuon() > 1) nSel = true; 
      if (eventdiffZ->GetPatMuon1SumPtR03() < 3 && eventdiffZ->GetPatMuon2SumPtR03() < 3 ) {
	candSel = true;
	isolation = true;
      }

      if (eventdiffZ->GetPatDiMuonEta()>0.) ZKinP = true;
      if (eventdiffZ->GetPatDiMuonEta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffZ->GetPatDiMuonEta()); }
      else{ etasignedHF = fabs(eventdiffZ->GetPatDiMuonEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffZ->GetPatDiMuonEta());}
      else {etasignedCASTOR = fabs(eventdiffZ->GetPatDiMuonEta());}

    }

    else{
      exit(EXIT_FAILURE);
    }

    //Branch Defining
    bRunNumber = eventdiff->GetRunNumber();
    bLumiSection = eventdiff->GetLumiSection();
    bEventNumber = eventdiff->GetEventNumber();
    if (typesel == "RecoMuon"){ 
      bDiBosonPt = eventdiffZ->GetDiMuonPt();
      bDiBosonEta = eventdiffZ->GetDiMuonEta();
      bDiBosonPhi = eventdiffZ->GetDiMuonPhi();
      bDiBosonMass = eventdiffZ->GetDiMuonMass();
    }
    if (typesel == "PatMuon"){
      bDiBosonPt = eventdiffZ->GetPatDiMuonPt();
      bDiBosonEta = eventdiffZ->GetPatDiMuonEta();
      bDiBosonPhi = eventdiffZ->GetPatDiMuonPhi();
      bDiBosonMass = eventdiffZ->GetPatDiMuonMass();
    }
    if (typesel == "RecoElectron"){
      bDiBosonPt = eventdiffZ->GetDiElectronPt();
      bDiBosonEta = eventdiffZ->GetDiElectronEta();
      bDiBosonPhi = eventdiffZ->GetDiElectronPhi();
      bDiBosonMass = eventdiffZ->GetDiElectronMass();
    }
    if (typesel == "PatElectron"){
      bDiBosonPt = eventdiffZ->GetPatDiElectronPt();
      bDiBosonEta = eventdiffZ->GetPatDiElectronEta();
      bDiBosonPhi = eventdiffZ->GetPatDiElectronPhi();
      bDiBosonMass = eventdiffZ->GetPatDiElectronMass();
    }
    bMultiplicityTracks = eventdiff->GetMultiplicityTracks();
    bSumEEEMinus = eventdiffZ->GetSumEEEMinus();
    bSumEEEPlus = eventdiffZ->GetSumEEEPlus();
    bSumEnergyHFMinus = eventdiff->GetSumEnergyHFMinus();
    bSumEnergyHFPlus = eventdiff->GetSumEnergyHFPlus();
    bsumCastorEnergy = sumCastorEnergy;
    bSectorCastorHit = counterHit;
    bMaxGapPF = eventdiffZ->GetMaxGapPF();
    bPTMaxGapMaxPF = eventdiffZ->GetPTMaxGapMaxPF();
    bPTMinGapMaxPF = eventdiffZ->GetPTMinGapMaxPF();
    bXiPlusFromPFCands = eventdiff->GetXiPlusFromPFCands();
    bXiMinusFromPFCands = eventdiff->GetXiMinusFromPFCands();
    betasignedHF = etasignedHF;
    betasignedCASTOR = etasignedCASTOR;
    bAEcastor = AEcastor;
    bdeltaetapf = deltaetapf;
    betamax = eventdiff->GetEtaMaxFromPFCands();
    betamin = etamin_;
    betalimmin = eventdiffZ->GetLimMinusGapPF();
    betalimmax = eventdiffZ->GetLimPlusGapPF();

    if(pileup < 21){ // Never comment this line. It is the program defense.

      if(switchtrigger == "trigger" || switchtrigger == "trigger_all_electron" || switchtrigger == "trigger_all_muon"){ 
	FillHistos(0,pileup,totalcommon); 
	if(trigger) {
	  ++totalT;
	  FillHistos(1,pileup,totalcommon);
	} 
	if(trigger && vertex && presel) FillHistos(2,pileup,totalcommon);
	if(trigger && vertex && presel && nSel) FillHistos(3,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge) FillHistos(4,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass) FillHistos(5,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation) FillHistos(6,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel) {
	  FillHistos(7,pileup,totalcommon);
	  fOutZ->cd();
	  troutZ->Fill();
	}
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln) FillHistos(8,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp) FillHistos(9,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln && castorgap) FillHistos(10,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castoractivity) FillHistos(11,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln && ZKinP){
	  outstring << "HF- Gap, Z Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(12,pileup,totalcommon);
	}
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && ZKinN){
	  outstring << "HF+ Gap, Z Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(13,pileup,totalcommon);
	}
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln && castorgap && ZKinP) FillHistos(14,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castoractivity && ZKinN) FillHistos(15,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && castorgap){
	  FillHistos(16,pileup,totalcommon);
	  fOutCASTOR->cd();
	  troutCASTOR->Fill();
	}
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && castorgap && ZKinP){
	  outstring << "CASTOR Gap, Z Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(17,pileup,totalcommon);
	  fOut->cd();
	  trout->Fill();
	}

	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castorgap) FillHistos(18,pileup,totalcommon);
	if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castorgap && ZKinP) FillHistos(19,pileup,totalcommon);

      }

      else if (switchtrigger =="no_trigger_nocorrection" || switchtrigger == "no_trigger_correction" ){
	--totalT;
	FillHistos(0,pileup,totalcommon);
	if(vertex && presel) FillHistos(2,pileup,totalcommon);
	if(vertex && presel && nSel) FillHistos(3,pileup,totalcommon);
	if(vertex && presel && nSel && charge) FillHistos(4,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass) FillHistos(5,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation) FillHistos(6,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel) {
	  FillHistos(7,pileup,totalcommon);
	  fOutZ->cd();
	  troutZ->Fill();
	}
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln) FillHistos(8,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp) FillHistos(9,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln && castorgap) FillHistos(10,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castoractivity) FillHistos(11,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln && ZKinP) FillHistos(12,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && ZKinN) FillHistos(13,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffseln && castorgap && ZKinP) FillHistos(14,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castoractivity && ZKinN) FillHistos(15,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && castorgap){ 
	  FillHistos(16,pileup,totalcommon);
	  fOutCASTOR->cd();
	  troutCASTOR->Fill();
	}
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && castorgap && ZKinP){
	  FillHistos(17,pileup,totalcommon);
	  fOut->cd();
	  trout->Fill();
	}
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castorgap) FillHistos(18,pileup,totalcommon);
	if(vertex && presel && nSel && charge && dimass && isolation && candSel && diffselp && castorgap && ZKinP) FillHistos(19,pileup,totalcommon);
      }

      else{
	exit(EXIT_FAILURE);
      }
    }  

  }

  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << ">> Input file: " << filein << std::endl;
  outstring << ">> Output file: " << savehistofile << std::endl;
  outstring << ">> TTree Name: " << processname << std::endl;
  outstring << " " << std::endl;
  outstring << "<< OPTIONS >>" << std::endl; 
  outstring << " " << std::endl;
  outstring << ">> Trigger: " << TriggerStatus << std::endl;
  outstring << ">> # Vertex: " << nVertex << std::endl;
  outstring << ">> Event Weight: " << switchlumiweight << " | MC Weight:" << mcweight << std::endl;
  outstring << ">> Lepton1(pT) > " << lepton1pt <<std::endl;
  outstring << ">> Lepton2(pT) > " << lepton2pt <<std::endl;
  outstring << ">> Type of Selection: " << selStatus << std::endl;
  outstring << ">> Gap Type of Selection: " << gapseltype << std::endl;
  outstring << selCastor << std::endl;
  outstring << " " << std::endl;
  outstring << "<< TRIGGER >> " << std::endl;
  outstring << " " << std::endl;
  outstring << "Total Number of Events: " <<  NEntries << std::endl;
  outstring << "Number of Triggered Events: " << totalT << " | Option: " << optTrigger << std::endl;
  outstring << "(if negative, trigger option is off)" << std::endl;
  outstring << "" << std::endl; 
  outstring << "Each Trigger Fired: " <<  std::endl;
  outstring << "Trigger 0: " << triggercounter[0] << std::endl;
  outstring << "Trigger 1: " << triggercounter[1] << std::endl;
  outstring << "Trigger 2: " << triggercounter[2] << std::endl;
  outstring << "Trigger 3: " << triggercounter[3] << std::endl;
  outstring << "Trigger 4: " << triggercounter[4] << std::endl;
  outstring << "Trigger 5: " << triggercounter[5] << std::endl;
  outstring << "Trigger 6: " << triggercounter[6] << std::endl;
  outstring << "Trigger 7: " << triggercounter[7] << std::endl;
  outstring << "Trigger 8: " << triggercounter[8] << std::endl;
  outstring << "Trigger 9: " << triggercounter[9] << std::endl;
  outstring << "Trigger 10: " << triggercounter[10] << std::endl;
  outstring << "Trigger 11: " << triggercounter[11] << std::endl;
  outstring << "Trigger 12: " << triggercounter[12] << std::endl;
  outstring << "Trigger 13: " << triggercounter[13] << std::endl;
  outstring << "Trigger 14: " << triggercounter[14] << std::endl;
  outstring << "Trigger 15: " << triggercounter[15] << std::endl;
  outstring << "Trigger 16: " << triggercounter[16] << std::endl;
  outstring << "Trigger 17: " << triggercounter[17] << std::endl;
  outstring << "Trigger 18: " << triggercounter[18] << std::endl;
  outstring << "Trigger 19: " << triggercounter[19] << std::endl;
  outstring << "" << std::endl;

  outf->cd();
  foldersFile[0] = outf->mkdir("LeptonsKinematics");
  foldersFile[1] = outf->mkdir("Detector");
  foldersFile[2] = outf->mkdir("Run");
  foldersFile[3] = outf->mkdir("Diffraction");
  SaveHistos(type,typesel);
  outf->Close();
  std::cout << "\n" << std::endl;

  fOut->cd();
  trout->Write();
  fOut->Close();

  fOutZ->cd();
  troutZ->Write();
  fOutZ->Close();

  fOutCASTOR->cd();
  troutCASTOR->Write();
  fOutCASTOR->Close();

}

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

int main(int argc, char **argv)
{

  AutoLibraryLoader::enable();

  const char *s1="*";
  std::string filein_;
  std::string processname_;
  std::string savehistofile_;
  double lepton1pt_;
  double lepton2pt_;
  int nVertex_;
  int optTrigger_;
  std::string switchtrigger_;
  double mcweight_;
  std::string switchlumiweight_;
  std::string type_;
  std::string typesel_;
  std::string pumfile_;
  std::string pudfile_;
  double castorthreshold_;
  double channelsthreshold_;
  std::string castorcorrfile_;
  std::string gapseltype_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0) filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0) processname_ = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0) savehistofile_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0) switchtrigger_ = argv[4];
  if (argc > 5 && strcmp(s1,argv[5]) != 0) optTrigger_   = atoi(argv[5]);
  if (argc > 6 && strcmp(s1,argv[6]) != 0) lepton1pt_ = atof(argv[6]);
  if (argc > 7 && strcmp(s1,argv[7]) != 0) lepton2pt_ = atof(argv[7]);
  if (argc > 8 && strcmp(s1,argv[8]) != 0) nVertex_ = atoi(argv[8]);
  if (argc > 9 && strcmp(s1,argv[9]) != 0) type_ = argv[9];
  if (argc > 10 && strcmp(s1,argv[10]) != 0) switchlumiweight_ = argv[10];
  if (argc > 11 && strcmp(s1,argv[11]) != 0) mcweight_ = atof(argv[11]);
  if (argc > 12 && strcmp(s1,argv[12]) != 0) typesel_ = argv[12];
  if (argc > 13 && strcmp(s1,argv[13]) != 0) castorthreshold_ = atof(argv[13]);
  if (argc > 14 && strcmp(s1,argv[14]) != 0) channelsthreshold_ = atof(argv[14]);
  if (argc > 15 && strcmp(s1,argv[15]) != 0) castorcorrfile_ = argv[15];
  if (argc > 16 && strcmp(s1,argv[16]) != 0) gapseltype_ = argv[16];
  if (argc > 17 && strcmp(s1,argv[17]) != 0) pumfile_ = argv[17];
  if (argc > 18 && strcmp(s1,argv[18]) != 0) pudfile_ = argv[18];

  std::cout << " " << std::endl;
  std::cout << ">>>>> Run with the Options <<<<< " << std::endl;
  std::cout << "Filename, PATTuple: " << filein_ << std::endl;
  std::cout << "ProcessName, PATTuple: " << processname_ << std::endl;
  std::cout << "Output name, histograms: " << savehistofile_ << std::endl;
  std::cout << "Trigger Switch: " << switchtrigger_ << std::endl;
  std::cout << "Trigger Path Option: " << optTrigger_ << std::endl;
  std::cout << "Lepton 1, pT [GeV]: " << lepton1pt_ << std::endl;
  std::cout << "Lepton 2, pT [GeV]: " << lepton2pt_ << std::endl;
  std::cout << "# Vertex: " << nVertex_ << std::endl;
  std::cout << "Multiple Pile-Up Option: " << type_ << std::endl;
  std::cout << "Luminosity Weight Option: " << switchlumiweight_ << std::endl;
  std::cout << "MC Weight: " << mcweight_ << std::endl;
  std::cout << "Type of Selection: " << typesel_ << std::endl;
  std::cout << "Castor Threshold: " << castorthreshold_ << std::endl;
  std::cout << "Channels Threshold: " << channelsthreshold_ << std::endl;
  std::cout << "Gap Type of Selection: " << gapseltype_ <<std::endl;
  std::cout << "" << std::endl;

  if (type_=="multiple_pileup" || type_=="no_multiple_pileup") {

    if (switchtrigger_=="trigger" || switchtrigger_=="no_trigger_nocorrection" || switchtrigger_=="no_trigger_correction" || switchtrigger_=="trigger_all_electron" || switchtrigger_=="trigger_all_muon") {}
    else{
      std::cout << "Please Insert type of Swithtrigger: " << std::endl;
      std::cout << "1) trigger: run with trigger. Need optTrigger >=0;" << std::endl;
      std::cout << "2) trigger_all_electron: all trigger electron path will be accepted. Do not require optTrigger." << std::endl;
      std::cout << "3) no_trigger: run without trigger." << std::endl;
      return 0;
    }

    if (switchlumiweight_=="mc_lumi_weight" || switchlumiweight_ == "mc_lumi_pu_weight" || switchlumiweight_ == "no_weight") {}
    else{
      std::cout << "Please Insert type of Event Weight (swithlumiweight): " << std::endl;
      std::cout << "1) no_weight: for data. Event weight = 1." << std::endl;
      std::cout << "2) mc_lumi_weight: Luminosity weight calculated by cross section. For Monte Carlo." << std::endl;
      std::cout << "3) mc_lumi_pu_weight: Luminosity weight and PU Reweight by Monte Carlo." << std::endl;
      return 0;
    }

    if(switchlumiweight_=="mc_lumi_pu_weight"){
      TFile pudatafile(pudfile_.c_str());
      TFile pumcfile(pumfile_.c_str());
      if (pudatafile.IsZombie()){
	std::cout << "-------------------------------------------" << std::endl;
	std::cout << " There is no file " << pudfile_ << " to be used for the" << std::endl;
	std::cout << " MC pile up correction." << std::endl;
	std::cout << "-------------------------------------------" << std::endl;
	return 0;
      }
      if (pumcfile.IsZombie()){
	std::cout << "-------------------------------------------" << std::endl;
	std::cout << " There is no file " << pumfile_ <<" to be used for the" << std::endl;
	std::cout << " MC pile up correction." << std::endl;
	std::cout << "-------------------------------------------" << std::endl;
	return 0;
      }
    }

    if (typesel_=="RecoMuon" || typesel_=="RecoElectron" || typesel_=="PatElectron" || typesel_=="PatMuon" ) {}
    else{
      std::cout << "Please Insert type of Selections: " << std::endl;
      std::cout << "1) RecoMuon: selections with Reco::Muon." << std::endl;
      std::cout << "2) RecoElectron: selections with Reco::Electron." << std::endl;
      std::cout << "3) PatMuon: selections with Pat::Muon." << std::endl;
      std::cout << "4) PatElectron: selections with Pat::Electron." << std::endl;
      return 0;
    }

    if (nVertex_ <= 0 || optTrigger_ < 0 || mcweight_ <= 0 || lepton1pt_ < 0 || lepton2pt_ < 0){
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << " Pay attention on the input numbers parameters" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << ">> Requirements:                             " << std::endl;
      std::cout << "I)   nVertex_ > 0 " << std::endl;
      std::cout << "II)  optTrigger >= 0" << std::endl;
      std::cout << "III) mcweight_ > 0" << std::endl;
      std::cout << "IV)  Lepton1pt_ and Lepton2pt_ >= 0" << std::endl;  
      std::cout << "----------------------------------------------" << std::endl; 
      return 0;
    }

    if ( (castorthreshold_ < 0 && channelsthreshold_ < 0) || (castorthreshold_ > 0 && channelsthreshold_ > 0) ){
      std::cout << "---------------------------------------------------" << std::endl;
      std::cout << "   Pay attention on the input numbers parameters" << std::endl;
      std::cout << "---------------------------------------------------" << std::endl;
      std::cout << ">> Requirements:                             " << std::endl;
      std::cout << "If castorthreshold_ and channelsthreshold are both" << std::endl;
      std::cout << "positive or negative, you can not run program. "           << std::endl; 
      std::cout << "" << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout << ">> Instructions: " << std::endl;
      std::cout << "I) castorthreshold_ > 0 and channelsthreshold_ < 0, allow Sector threshold." << std::endl;
      std::cout << "I) castorthreshold_ < 0 and channelsthreshold_ > 0, allow Channel threshold." << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      return 0;
    }

    if(switchtrigger_=="no_trigger_correction"){
      TFile castorfile(castorcorrfile_.c_str());
      if (castorfile.IsZombie()){
	std::cout << "-------------------------------------------" << std::endl;
	std::cout << " There is no the file " << castorcorrfile_ << " or the" << std::endl;
	std::cout << " path is not correct." << std::endl;
	std::cout << "-------------------------------------------" << std::endl;
	return 0;
      }
    }

    DiffractiveZ* diffZRun = new DiffractiveZ();
    diffZRun->CreateHistos(type_);
    diffZRun->Run(filein_, processname_, savehistofile_, switchtrigger_, optTrigger_, lepton1pt_, lepton2pt_, nVertex_, type_, switchlumiweight_, mcweight_, typesel_, castorthreshold_, channelsthreshold_, castorcorrfile_, gapseltype_, pumfile_, pudfile_);
    return 0;
  }

  else{
    std::cout << "Please Insert type of histograms: " << std::endl;
    std::cout << "1) multiple_pileup: create histograms for each pile-up event. It works only for MC with pile-up." << std::endl;
    std::cout << "2) no_multiple_pileup: create histograms without each pile-up event. It works for data and MC with/without pile-up." << std::endl;
    return 0;
  }
}

#endif
