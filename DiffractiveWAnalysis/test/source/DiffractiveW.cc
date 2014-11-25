//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsDiffractiveWsAnalysis#Macro_Analysis
//
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
#include "DiffractiveW.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveWEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#define isfinite(x) !std::isinf(x)

using namespace diffractiveAnalysis;
using namespace diffractiveWAnalysis;
using namespace eventInfo;
using namespace reweight;

void DiffractiveW::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventdiff = new DiffractiveEvent();
  eventdiffW = new DiffractiveWEvent();
  eventinfo = new EventInfoEvent();
  diff = tr->GetBranch("DiffractiveAnalysis");
  diffW = tr->GetBranch("DiffractiveWAnalysis");
  info = tr->GetBranch("EventInfo");
  diff->SetAddress(&eventdiff);
  diffW->SetAddress(&eventdiffW);
  info->SetAddress(&eventinfo);

}

void DiffractiveW::CleanVariables(){

  bosonWMass = -999.;
  bosonWEta = -999.;
  bosonWPt = -999.;
  bosonWPhi = -999.;
  bosonWCharge = -999;
  NElectrons = -999;
  NMuons = -999;
  isoTk1 = -999.;
  isoEcal1 = -999.;
  isoHcal1 = -999.;
  isoRec = -999.;
  innerHits1 = -999;
  Dcot1 = -999.;
  Dist1 = -999.;
  DeltaEtaTkClu1 = -999;
  DeltaPhiTkClu1 = -999.;
  sigmaIeIe1 = -999.;
  HE1 = -999.;
  etasignedHF = -999.;
  etasignedCASTOR = -999;
  aSumE = -999.;
  AEcastor = -999.;
  etamin_= -999.;
  deltaetapf = -999.;
  deltaetapfcastor = -999.;

}

void DiffractiveW::CreateHistos(std::string type){

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
  std::string step6 = "NGapCMS";
  std::string step7 = "PGapCMS";
  std::string step8 = "NGapCMSAndCASTOR";
  std::string step9 = "PGapCMSAndCastorActivity";
  std::string step10 = "NGapCMSAndWKinP";
  std::string step11 = "PGapCMSAndWKinN";
  std::string step12 = "NGapCMSAndCASTORAndWKinP";
  std::string step13 = "PGapCMSAndCastorActivityAndWKinN";
  std::string step14 = "NGapCASTOR";
  std::string step15 = "NGapCASTORAndWKinP";
  std::string step16 = "PGapCMSAndCASTOR";
  std::string step17 = "PGapCMSAndCASTORAndWKinP";

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
    m_hVector_NElectrons.push_back( std::vector<TH1F*>() );
    m_hVector_NMuons.push_back( std::vector<TH1F*>() );
    m_hVector_WMass.push_back( std::vector<TH1F*>() );
    m_hVector_WEta.push_back( std::vector<TH1F*>() );
    m_hVector_WPt.push_back( std::vector<TH1F*>() );
    m_hVector_WPhi.push_back( std::vector<TH1F*>() );
    m_hVector_WCharge.push_back( std::vector<TH1F*>() );
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

    // Detector
    m_hVector_sumEHFplus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHFminus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHEplus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEHEminus.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFplus_S.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFminus_S.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFplus_L.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFminus_L.push_back( std::vector<TH1F*>() );
    m_hVector_sumEEEplus.push_back( std::vector<TH1F*>() );
    m_hVector_sumEEEminus.push_back( std::vector<TH1F*>() );
    m_hVector_EnergyHFPlusVsEnergyHFMinus.push_back( std::vector<TH2F*>() );
    m_hVector_EnergyEEPlusVsEnergyEEMinus.push_back( std::vector<TH2F*>() );
    m_hVector_SumEHFMax.push_back( std::vector<TH1F*>() );
    m_hVector_SumEHFMin.push_back( std::vector<TH1F*>() );
    m_hVector_etcalos_p.push_back( std::vector<TH2F*>() );
    m_hVector_etcalos_n.push_back( std::vector<TH2F*>() );
    m_hVector_ECaloVsEta.push_back( std::vector<TH2F*>() );
    m_hVector_ECaloVsEtaTProf.push_back( std::vector<TProfile*>() );
    m_hVector_EnergyVsEtaBin1D.push_back( std::vector<TH1F*>() );
    m_hVector_sumECastorMinus.push_back( std::vector<TH1F*>() );
    m_hVector_ECastorSector.push_back( std::vector<TH2F*>() );
    m_hVector_ECastorSectorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_ECastorSectorBin1D.push_back( std::vector<TH1F*>() );
    m_hVector_EnergyHFMinusVsCastorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_EnergyHFPlusVsCastorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_sumECastorAndHFMinus.push_back( std::vector<TH1F*>() );
    m_hVector_CastorMultiplicity.push_back( std::vector<TH1F*>() );
    m_hVector_CastorMultiplicityVsLumi.push_back( std::vector<TH2F*>() );
    m_hVector_SectorVsTotalCastorEnergy.push_back( std::vector<TH2F*>() );
    m_hVector_SectorVsTotalCastorEnergyTProf.push_back( std::vector<TProfile*>() );


    // Event Information
    m_hVector_lumi.push_back( std::vector<TH1F*>() );
    m_hVector_tracks.push_back( std::vector<TH1F*>() );
    m_hVector_vertex.push_back( std::vector<TH1F*>() );


    // Diffraction
    m_hVector_asumE.push_back( std::vector<TH1F*>() );
    m_hVector_multhf.push_back( std::vector<TH2F*>() );
    m_hVector_pfetamax.push_back( std::vector<TH1F*>() );
    m_hVector_pfetamin.push_back( std::vector<TH1F*>() );
    m_hVector_pfetamincastor.push_back( std::vector<TH1F*>() );
    m_hVector_maxetagap.push_back( std::vector<TH1F*>() );
    m_hVector_SumPTMax.push_back( std::vector<TH1F*>() );
    m_hVector_SumPTMin.push_back( std::vector<TH1F*>() );
    m_hVector_absdeltaEtaPF.push_back( std::vector<TH1F*>() );
    m_hVector_deltaEtaPF.push_back( std::vector<TH1F*>() );
    m_hVector_absdeltaEtaPFCastor.push_back( std::vector<TH1F*>() );
    m_hVector_deltaEtaPFCastor.push_back( std::vector<TH1F*>() );
    m_hVector_XiPlusPF.push_back( std::vector<TH1F*>() );
    m_hVector_XiMinusPF.push_back( std::vector<TH1F*>() );
    m_hVector_XiPF.push_back( std::vector<TH1F*>() );
    m_hVector_AEcastor.push_back( std::vector<TH1F*>() );
    m_hVector_etasignedHF.push_back( std::vector<TH1F*>() );
    m_hVector_etasignedCASTOR.push_back( std::vector<TH1F*>() );

    for (int k=0;k<nloop;k++){

      if (type=="multiple_pileup"){
	sprintf(tag,"multiple_pileup_%i",k);
      }
      else{
	sprintf(tag,"single");
      }

      char name[300];

      // Kinematics
      sprintf(name,"NMuons_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_NMuons = new TH1F(name,"Number of Muons per Event; # Muons; Multiplicity", 100, 0., 100.);
      m_hVector_NMuons[j].push_back(histo_NMuons);

      sprintf(name,"NElectrons_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_NElectrons = new TH1F(name,"Number of Electrons per Event; # Electrons; Multiplicity", 100, 0., 100.);
      m_hVector_NElectrons[j].push_back(histo_NElectrons);

      sprintf(name,"BosonWMass_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WMass = new TH1F(name,"Boson W Transverse Mass Distribution; M_{T}(l#nu) [GeV]; N events",500,0,500);
      m_hVector_WMass[j].push_back(histo_WMass);

      sprintf(name,"BosonWPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WPt = new TH1F(name,"Boson W Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_WPt[j].push_back(histo_WPt);

      sprintf(name,"BosonWEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WEta = new TH1F(name,"Boson W #eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_WEta[j].push_back(histo_WEta);

      sprintf(name,"BosonWPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WPhi = new TH1F(name,"Boson W #phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_WPhi[j].push_back(histo_WPhi);

      sprintf(name,"BosonWCharge_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WCharge = new TH1F(name,"Boson W Charge Distribution; Charge; N events",6,-3,3);
      m_hVector_WCharge[j].push_back(histo_WCharge);

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


      // Detector

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

      sprintf(name,"sumEEEplus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEEEplus = new TH1F(name,"EE^{+} - Sum of Energy; #sum E_{EE^{+}} [GeV]; N events",1000,0.,500.);
      m_hVector_sumEEEplus[j].push_back(histo_sumEEEplus);

      sprintf(name,"sumEEEminus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumEEEminus = new TH1F(name,"EE^{-} - Sum of Energy; #sum E_{EE^{-}} [GeV]; N events",1000,0.,500.);
      m_hVector_sumEEEminus[j].push_back(histo_sumEEEminus);

      sprintf(name,"EnergyHFPlusVsEnergyHFMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_EnergyHFPlusVsEnergyHFMinus = new TH2F(name,"HF^{+} and HF^{-}; #sum Energy HF^{+} [GeV]; #sum Energy HF^{-} [GeV]; N events",1000,0.,1000.,1000,0.,1000.);
      m_hVector_EnergyHFPlusVsEnergyHFMinus[j].push_back(histo_EnergyHFPlusVsEnergyHFMinus);

      sprintf(name,"EnergyEEPlusVsEnergyEEMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_EnergyEEPlusVsEnergyEEMinus = new TH2F(name,"EE^{+} and EE^{-}; #sum Energy EE^{+} [GeV]; #sum Energy EE^{-} [GeV]; N events",1000,0.,500.,1000,0.,500.);
      m_hVector_EnergyEEPlusVsEnergyEEMinus[j].push_back(histo_EnergyEEPlusVsEnergyEEMinus);    

      sprintf(name,"sumEHFMax_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFMax = new TH1F(name,"HF - Sum of Energy; #sum E_{HF,Max} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFMax[j].push_back(histo_SumEHFMax);

      sprintf(name,"sumEHFMin_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumEHFMin = new TH1F(name,"HF - Sum of Energy; #sum E_{HF,Min} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFMin[j].push_back(histo_SumEHFMin);

      sprintf(name,"EnergyHFPlusVsCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ET_Calos_p = new TH2F(name,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 6000, 0., 3000. );
      m_hVector_etcalos_p[j].push_back(histo_ET_Calos_p);

      sprintf(name,"EnergyHFMinusVsCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ET_Calos_n = new TH2F(name,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 6000, 0., 3000. );
      m_hVector_etcalos_n[j].push_back(histo_ET_Calos_n);

      sprintf(name,"ECaloVsEta_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_ECaloVsEta = new TH2F(name,"Calorimeter Energy X #eta; #eta; Energy [GeV]", 500, -8, 8, 100, 0., 1000.);
      m_hVector_ECaloVsEta[j].push_back(histo_ECaloVsEta);

      sprintf(name,"ECaloVsEtaTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_ECaloVsEtaTProf = new TProfile(name,"Calorimeter Energy X #eta; #eta; <Energy> [GeV]", 100, -8, 8, 0., 1000.);
      m_hVector_ECaloVsEtaTProf[j].push_back(histo_ECaloVsEtaTProf);

      sprintf(name,"EnergyVsEtaBin1D_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_EnergyVsEtaBin1D = new TH1F(name,"Calorimeter Energy X #eta; #eta; #sum Energy_{calotower} [GeV]", 500, -8, 8);
      m_hVector_EnergyVsEtaBin1D[j].push_back(histo_EnergyVsEtaBin1D);

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

      sprintf(name,"EnergyHFPlusVsCastorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_EnergyHFPlusVsCastorTProf = new TProfile(name,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 3000. );
      m_hVector_EnergyHFPlusVsCastorTProf[j].push_back(histo_EnergyHFPlusVsCastorTProf);

      sprintf(name,"EnergyHFMinusVsCastorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_EnergyHFMinusVsCastorTProf = new TProfile(name,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 3000. );
      m_hVector_EnergyHFMinusVsCastorTProf[j].push_back(histo_EnergyHFMinusVsCastorTProf);

      sprintf(name,"sumECastorAndSumHFMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_sumECastorAndHFMinus = new TH1F(name,"HF^{-} and Castor Sum of Energy; Energy [GeV]; N events",6000,0,3000);
      m_hVector_sumECastorAndHFMinus[j].push_back(histo_sumECastorAndHFMinus);

      sprintf(name,"CastorMultiplicity_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_CastorMultiplicity = new TH1F(name,"Castor: number of sectors with activity; #Sectors; N events",81,0,81);
      m_hVector_CastorMultiplicity[j].push_back(histo_CastorMultiplicity);

      sprintf(name,"CastorMultiplicityVsLumi_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_CastorMultiplicityVsLumi = new TH2F(name,"CastorMultiplicity Vs Luminosity; Luminosity per Bunch [#mub^{-1}s^{-1}]; Castor Multiplicity",5000,0,2,81,0,81);
      m_hVector_CastorMultiplicityVsLumi[j].push_back(histo_CastorMultiplicityVsLumi);

      sprintf(name,"SectorVsTotalCastorEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_SectorVsTotalCastorEnergy = new TH2F(name,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",17,0,17,1500,0,1500);
      m_hVector_SectorVsTotalCastorEnergy[j].push_back(histo_SectorVsTotalCastorEnergy);

      sprintf(name,"SectorVsTotalCastorEnergyTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_SectorVsTotalCastorEnergyTProf = new TProfile(name,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",17,0,17,0,1500);
      m_hVector_SectorVsTotalCastorEnergyTProf[j].push_back(histo_SectorVsTotalCastorEnergyTProf);


      // Event Information
      sprintf(name,"lumi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_lumi = new TH1F(name,"Luminosity per Bunch; L_{Bunch} [#mub^{-1}s^{-1}]; N events",25,0,2);
      m_hVector_lumi[j].push_back(histo_lumi); 

      sprintf(name,"TracksMultiplicity_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_Tracks = new TH1F(name,"Tracks Multiplicity; n Tracks; N events",150,0,150);
      m_hVector_tracks[j].push_back(histo_Tracks);

      sprintf(name,"vertex_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_vertex = new TH1F(name,"Number of Vertex; # Vertex; N events",25,0,25);
      m_hVector_vertex[j].push_back(histo_vertex);


      // Diffraction
      sprintf(name,"aEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_aSumE = new TH1F(name,"Forward Backward Asymmetry Distribution ; (#sum HF^{+} - #sum HF^{-})x(#sum HF^{+} + #sum HF^{-})^{-1}; N events",100,-2,2);
      m_hVector_asumE[j].push_back(histo_aSumE);

      sprintf(name,"mHF_%s_%s",tag,Folders.at(j).c_str());
      TH2F *histo_MultHF = new TH2F(name,"HF^{+} and HF^{-} Multiplicity; n HF^{+}; n HF^{-}; N events", 100, 0., 100., 100, 0., 100. );
      m_hVector_multhf[j].push_back(histo_MultHF);

      sprintf(name,"pfetamax_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_PFEtamax = new TH1F(name,"Particle Flow #eta_{max} Distribution; #eta; N events",18,binarrayplus);
      m_hVector_pfetamax[j].push_back(histo_PFEtamax);

      sprintf(name,"pfetamin_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_PFEtamin = new TH1F(name,"Particle Flow #eta_{min} Distribution; #eta; N events",18,binarrayminus);
      m_hVector_pfetamin[j].push_back(histo_PFEtamin);

      sprintf(name,"pfetamincastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_PFEtamincastor = new TH1F(name,"Particle Flow #eta_{min} Distribution; #eta; N events",18,binarrayminus);
      m_hVector_pfetamincastor[j].push_back(histo_PFEtamincastor);

      sprintf(name,"maxetagap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_maxetagap = new TH1F(name,"Particle Flow #Delta#eta_{max} Distribution; #Delta#eta_{max}; N events",40,0.,6.);
      m_hVector_maxetagap[j].push_back(histo_maxetagap);

      sprintf(name,"SumPTMaxgap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumPTMax = new TH1F(name,"Particle Flow P_{T}; #sum pT_{max} [GeV]; N events",1200,0,600);
      m_hVector_SumPTMax[j].push_back(histo_SumPTMax);

      sprintf(name,"SumPTMingap_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_SumPTMin = new TH1F(name,"Particle Flow P_{T}; #sum pT_{min} [GeV]; N events",1200,0,600);
      m_hVector_SumPTMin[j].push_back(histo_SumPTMin);

      sprintf(name,"deltaEtamaxminPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_deltaEtaPF = new TH1F(name,"#Delta#eta_{PF} Distribution; #eta_{max}-#eta_{min}; N events",14,binarraydelta);
      m_hVector_deltaEtaPF[j].push_back(histo_deltaEtaPF);

      sprintf(name,"absdeltaEtamaxminPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_absdeltaEtaPF = new TH1F(name,"|#Delta#eta_{PF}| Distribution; |#eta_{max}-#eta_{min}|; N events",14,binarraydelta);
      m_hVector_absdeltaEtaPF[j].push_back(histo_absdeltaEtaPF);

      sprintf(name,"deltaEtamaxminPFCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_deltaEtaPFCastor = new TH1F(name,"#Delta#eta_{PF} Distribution; #eta_{max}-#eta_{min}; N events",14,binarraydelta);
      m_hVector_deltaEtaPFCastor[j].push_back(histo_deltaEtaPFCastor);

      sprintf(name,"absdeltaEtamaxminPFCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_absdeltaEtaPFCastor = new TH1F(name,"|#Delta#eta_{PF}| Distribution; |#eta_{max}-#eta_{min}|; N events",14,binarraydelta);
      m_hVector_absdeltaEtaPFCastor[j].push_back(histo_absdeltaEtaPFCastor);

      sprintf(name,"xiPlusPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_XiPlusPF = new TH1F(name,"#xi_{plus} Particle Flow; #xi_{plus}; N Event",17,xi_bin);
      m_hVector_XiPlusPF[j].push_back(histo_XiPlusPF);

      sprintf(name,"xiMinusPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_XiMinusPF = new TH1F(name,"#xi_{minus} Particle Flow; #xi_{minus}; N Event",17,xi_bin);
      m_hVector_XiMinusPF[j].push_back(histo_XiMinusPF);

      sprintf(name,"xiPF_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_XiPF = new TH1F(name,"#xi Particle Flow; #xi; N Event",17,xi_bin);
      m_hVector_XiPF[j].push_back(histo_XiPF);

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

void DiffractiveW::FillHistos(int index, int pileup, double totalweight){

  // Defense
  if (!isfinite(bosonWMass)) {
    bosonWMass = -999.;
  }

  // Kinematics
  m_hVector_WMass[index].at(pileup)->Fill(bosonWMass,totalweight);
  m_hVector_WEta[index].at(pileup)->Fill(bosonWEta,totalweight);
  m_hVector_WPt[index].at(pileup)->Fill(bosonWPt,totalweight);
  m_hVector_WPhi[index].at(pileup)->Fill(bosonWPhi,totalweight);
  m_hVector_WCharge[index].at(pileup)->Fill(bosonWCharge,totalweight);
  m_hVector_NElectrons[index].at(pileup)->Fill(NElectrons,totalweight);
  m_hVector_NMuons[index].at(pileup)->Fill(NMuons,totalweight);
  m_hVector_LeadingLeptonTkDr03[index].at(pileup)->Fill(isoTk1,totalweight);
  m_hVector_LeadingLeptonEcalDr03[index].at(pileup)->Fill(isoEcal1,totalweight);
  m_hVector_LeadingLeptonHcalDr03[index].at(pileup)->Fill(isoHcal1,totalweight);
  m_hVector_LeadingLeptonIsolation[index].at(pileup)->Fill(isoRec,totalweight);
  m_hVector_LeadingLeptonInnerHits[index].at(pileup)->Fill(innerHits1,totalweight);
  m_hVector_LeadingLeptonDCot[index].at(pileup)->Fill(Dcot1,totalweight);
  m_hVector_LeadingLeptonDist[index].at(pileup)->Fill(Dist1,totalweight);
  m_hVector_LeadingLeptonDeltaEtaTkClu[index].at(pileup)->Fill(DeltaEtaTkClu1,totalweight);
  m_hVector_LeadingLeptonDeltaPhiTkClu[index].at(pileup)->Fill(DeltaPhiTkClu1,totalweight);
  m_hVector_LeadingLeptonSigmaIeIe[index].at(pileup)->Fill(sigmaIeIe1,totalweight);
  m_hVector_LeadingLeptonHE[index].at(pileup)->Fill(HE1,totalweight);

  // Detector
  m_hVector_sumEHFplus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),totalweight);
  m_hVector_sumEHFminus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),totalweight);
  m_hVector_sumEHEplus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHEPlus(),totalweight);
  m_hVector_sumEHEminus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHEMinus(),totalweight);
  m_hVector_SumEHFplus_S[index].at(pileup)->Fill(eventdiffW->GetSumEHF_SPlus(),totalweight);
  m_hVector_SumEHFminus_S[index].at(pileup)->Fill(eventdiffW->GetSumEHF_SMinus(),totalweight);
  m_hVector_SumEHFplus_L[index].at(pileup)->Fill(eventdiffW->GetSumEHF_LPlus(),totalweight);
  m_hVector_SumEHFminus_L[index].at(pileup)->Fill(eventdiffW->GetSumEHF_LMinus(),totalweight);
  m_hVector_sumEEEminus[index].at(pileup)->Fill(eventdiffW->GetSumEEEMinus(),totalweight);
  m_hVector_sumEEEplus[index].at(pileup)->Fill(eventdiffW->GetSumEEEPlus(),totalweight);
  m_hVector_EnergyHFPlusVsEnergyHFMinus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),eventdiff->GetSumEnergyHFMinus(),totalweight);
  m_hVector_EnergyEEPlusVsEnergyEEMinus[index].at(pileup)->Fill(eventdiffW->GetSumEEEPlus(),eventdiffW->GetSumEEEMinus(),totalweight);

  for (int k=0; k<eventdiffW->GetEachTowerCounter();k++){
    m_hVector_ECaloVsEta[index].at(pileup)->Fill(eventdiffW->GetEachTowerEta(k),eventdiffW->GetEachTowerEnergy(k),totalweight);
    m_hVector_ECaloVsEtaTProf[index].at(pileup)->Fill(eventdiffW->GetEachTowerEta(k),eventdiffW->GetEachTowerEnergy(k),totalweight);
    m_hVector_EnergyVsEtaBin1D[index].at(pileup)->Fill(eventdiffW->GetEachTowerEta(k),eventdiffW->GetEachTowerEnergy(k)*totalweight);
  }

  if (eventdiff->GetSumEnergyHFPlus() > eventdiff->GetSumEnergyHFMinus()){
    m_hVector_SumEHFMax[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),totalweight);
    m_hVector_SumEHFMin[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),totalweight);
  }else{
    m_hVector_SumEHFMax[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),totalweight);
    m_hVector_SumEHFMin[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),totalweight);
  }

  m_hVector_etcalos_p[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_etcalos_n[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_sumECastorMinus[index].at(pileup)->Fill(sumCastorEnergy,totalweight);

  for (int k=0; k<16;k++){
    if (CastorEnergySector[k] >=1.){
      m_hVector_ECastorSector[index].at(pileup)->Fill(k+1,CastorEnergySector[k]);
      m_hVector_ECastorSectorTProf[index].at(pileup)->Fill(k+1,CastorEnergySector[k]);
      m_hVector_ECastorSectorBin1D[index].at(pileup)->Fill(k+1,CastorEnergySector[k]);
    }else{
      m_hVector_ECastorSector[index].at(pileup)->Fill(k+1,0);
    }
  }

  m_hVector_EnergyHFMinusVsCastorTProf[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_EnergyHFPlusVsCastorTProf[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_sumECastorAndHFMinus[index].at(pileup)->Fill(sumCastorAndHFMinusEnergy,totalweight);
  m_hVector_CastorMultiplicity[index].at(pileup)->Fill(counterHit,totalweight);
  m_hVector_CastorMultiplicityVsLumi[index].at(pileup)->Fill(eventinfo->GetInstLumiBunch(),counterHit,totalweight);
  m_hVector_SectorVsTotalCastorEnergy[index].at(pileup)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);
  m_hVector_SectorVsTotalCastorEnergyTProf[index].at(pileup)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);

  // Event Information
  m_hVector_lumi[index].at(pileup)->Fill(eventinfo->GetInstLumiBunch(),totalweight);
  m_hVector_tracks[index].at(pileup)->Fill(eventdiff->GetMultiplicityTracks(),totalweight);
  m_hVector_vertex[index].at(pileup)->Fill(eventdiff->GetNVertex(),totalweight);

  // Diffractive Variables
  m_hVector_asumE[index].at(pileup)->Fill(aSumE,totalweight);
  m_hVector_multhf[index].at(pileup)->Fill(eventdiff->GetMultiplicityHFPlus(),eventdiff->GetMultiplicityHFMinus(),totalweight);
  m_hVector_pfetamax[index].at(pileup)->Fill(eventdiff->GetEtaMaxFromPFCands(),totalweight);
  m_hVector_pfetamin[index].at(pileup)->Fill(eventdiff->GetEtaMinFromPFCands(),totalweight);
  m_hVector_pfetamincastor[index].at(pileup)->Fill(etamin_,totalweight);
  m_hVector_maxetagap[index].at(pileup)->Fill(fabs(eventdiffW->GetLrgPF()),totalweight);
  if(fabs(eventdiffW->GetSumptPFLeft())>fabs(eventdiffW->GetSumptPFRight())){
    m_hVector_SumPTMax[index].at(pileup)->Fill(eventdiffW->GetSumptPFLeft(),totalweight);
    m_hVector_SumPTMin[index].at(pileup)->Fill(eventdiffW->GetSumptPFRight(),totalweight);
  }else{
    m_hVector_SumPTMax[index].at(pileup)->Fill(eventdiffW->GetSumptPFRight(),totalweight);
    m_hVector_SumPTMin[index].at(pileup)->Fill(eventdiffW->GetSumptPFLeft(),totalweight);
  }
  m_hVector_absdeltaEtaPF[index].at(pileup)->Fill(fabs(deltaetapf),totalweight);
  m_hVector_deltaEtaPF[index].at(pileup)->Fill(deltaetapf,totalweight);
  m_hVector_absdeltaEtaPFCastor[index].at(pileup)->Fill(fabs(deltaetapfcastor),totalweight);
  m_hVector_deltaEtaPFCastor[index].at(pileup)->Fill(deltaetapfcastor,totalweight);
  m_hVector_XiPlusPF[index].at(pileup)->Fill(eventdiff->GetXiPlusFromPFCands(),totalweight);
  m_hVector_XiMinusPF[index].at(pileup)->Fill(eventdiff->GetXiMinusFromPFCands(),totalweight);
  m_hVector_XiPF[index].at(pileup)->Fill(eventdiff->GetXiPlusFromPFCands(),totalweight);
  m_hVector_XiPF[index].at(pileup)->Fill(eventdiff->GetXiMinusFromPFCands(),totalweight);
  m_hVector_AEcastor[index].at(pileup)->Fill(AEcastor,totalweight);
  m_hVector_etasignedHF[index].at(pileup)->Fill(etasignedHF,totalweight);
  m_hVector_etasignedCASTOR[index].at(pileup)->Fill(etasignedCASTOR,totalweight);

}

void DiffractiveW::SaveHistos(std::string type,std::string typesel){

  // Creating Correlation Histograms

  int ipileup;

  if (type=="multiple_pileup") ipileup=21;
  else ipileup=1;

  for (int i = 0; i < ipileup; i++){
    for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

      // Kinematics Folder
      foldersFile[0]->cd();
      m_hVector_WMass[j].at(i)->Write();
      m_hVector_WEta[j].at(i)->Write();
      m_hVector_WPt[j].at(i)->Write();
      m_hVector_WPhi[j].at(i)->Write();
      m_hVector_WCharge[j].at(i)->Write();
      m_hVector_NElectrons[j].at(i)->Write();
      m_hVector_NMuons[j].at(i)->Write();
      m_hVector_LeadingLeptonTkDr03[j].at(i)->Write();
      m_hVector_LeadingLeptonEcalDr03[j].at(i)->Write();
      m_hVector_LeadingLeptonHcalDr03[j].at(i)->Write();
      m_hVector_LeadingLeptonIsolation[j].at(i)->Write();
      m_hVector_LeadingLeptonInnerHits[j].at(i)->Write();
      m_hVector_LeadingLeptonDCot[j].at(i)->Write();
      m_hVector_LeadingLeptonDist[j].at(i)->Write();
      m_hVector_LeadingLeptonDeltaEtaTkClu[j].at(i)->Write();
      m_hVector_LeadingLeptonDeltaPhiTkClu[j].at(i)->Write();
      m_hVector_LeadingLeptonSigmaIeIe[j].at(i)->Write();
      m_hVector_LeadingLeptonHE[j].at(i)->Write();

      // Detector Folder
      foldersFile[1]->cd();
      m_hVector_sumEHFplus[j].at(i)->Write();
      m_hVector_sumEHFminus[j].at(i)->Write();
      m_hVector_sumEHEplus[j].at(i)->Write();
      m_hVector_sumEHEminus[j].at(i)->Write();
      m_hVector_SumEHFplus_S[j].at(i)->Write();
      m_hVector_SumEHFminus_S[j].at(i)->Write();
      m_hVector_SumEHFplus_L[j].at(i)->Write();
      m_hVector_SumEHFminus_L[j].at(i)->Write();
      m_hVector_sumEEEplus[j].at(i)->Write();
      m_hVector_sumEEEminus[j].at(i)->Write();
      m_hVector_EnergyHFPlusVsEnergyHFMinus[j].at(i)->Write();
      m_hVector_EnergyEEPlusVsEnergyEEMinus[j].at(i)->Write();
      m_hVector_SumEHFMax[j].at(i)->Write();
      m_hVector_SumEHFMin[j].at(i)->Write();
      m_hVector_etcalos_p[j].at(i)->Write();
      m_hVector_etcalos_n[j].at(i)->Write();
      m_hVector_ECaloVsEta[j].at(i)->Write();
      m_hVector_ECaloVsEtaTProf[j].at(i)->Write();  
      m_hVector_EnergyVsEtaBin1D[j].at(i)->Write();
      m_hVector_sumECastorMinus[j].at(i)->Write();
      m_hVector_ECastorSector[j].at(i)->Write();
      m_hVector_ECastorSectorTProf[j].at(i)->Write();
      m_hVector_ECastorSectorBin1D[j].at(i)->Write();
      m_hVector_EnergyHFMinusVsCastorTProf[j].at(i)->Write();
      m_hVector_EnergyHFPlusVsCastorTProf[j].at(i)->Write();
      m_hVector_sumECastorAndHFMinus[j].at(i)->Write();
      m_hVector_CastorMultiplicity[j].at(i)->Write();
      m_hVector_CastorMultiplicityVsLumi[j].at(i)->Write();
      m_hVector_SectorVsTotalCastorEnergy[j].at(i)->Write();
      m_hVector_SectorVsTotalCastorEnergyTProf[j].at(i)->Write(); 

      // Event Information
      foldersFile[2]->cd();
      m_hVector_lumi[j].at(i)->Write();
      m_hVector_tracks[j].at(i)->Write();
      m_hVector_vertex[j].at(i)->Write();

      // Diffraction
      foldersFile[3]->cd();
      m_hVector_asumE[j].at(i)->Write();
      m_hVector_multhf[j].at(i)->Write();
      m_hVector_pfetamax[j].at(i)->Write();
      m_hVector_pfetamin[j].at(i)->Write();
      m_hVector_pfetamincastor[j].at(i)->Write();
      m_hVector_maxetagap[j].at(i)->Write();
      m_hVector_SumPTMax[j].at(i)->Write();
      m_hVector_SumPTMin[j].at(i)->Write();
      m_hVector_absdeltaEtaPF[j].at(i)->Write();
      m_hVector_deltaEtaPF[j].at(i)->Write();
      m_hVector_absdeltaEtaPFCastor[j].at(i)->Write();
      m_hVector_deltaEtaPFCastor[j].at(i)->Write();
      m_hVector_XiPlusPF[j].at(i)->Write();
      m_hVector_XiMinusPF[j].at(i)->Write();
      m_hVector_XiPF[j].at(i)->Write();
      m_hVector_AEcastor[j].at(i)->Write();
      m_hVector_etasignedHF[j].at(i)->Write();
      m_hVector_etasignedCASTOR[j].at(i)->Write();

    }
  }

}

void DiffractiveW::Run(std::string filein_, std::string processname_, std::string savehistofile_, std::string switchtrigger_, int optTrigger_, double lepton1pt_, double lepton2pt_, int nVertex_, std::string type_, std::string switchlumiweight_, double mcweight_, std::string typesel_, double castorthreshold_, double channelsthreshold_, std::string castorcorrfile_, std::string gapseltype_, std::string pumfile_, std::string pudfile_){

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
  TString TTreeoutput, TTreeAllW, TTreeCASTOR;
  TTreeoutput = "TTreeGoldenDiffW_" + savehistofile;
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
  trout->Branch("AEcastor",&bAEcastor,"AEcastor/D");
  trout->Branch("etasignedHF",&betasignedHF,"betasignedHF/D");
  trout->Branch("etasignedCASTOR",&betasignedCASTOR,"betasignedCASTOR/D");
  trout->Branch("MaxGapPF",&bMaxGapPF,"bMaxGapPF/D");
  trout->Branch("SumPTMaxLrgPF",&SumPTMaxLrgPF,"SumPTMaxLrgPF/D");
  trout->Branch("SumPTMinLrgPF",&SumPTMinLrgPF,"SumPTMinLrgPF/D");
  trout->Branch("XiPlusFromPFCands",&bXiPlusFromPFCands,"bXiPlusFromPFCands/D");
  trout->Branch("XiMinusFromPFCands",&bXiMinusFromPFCands,"bXiMinusFromPFCands/D");
  trout->Branch("EtaMaxPF",&betamax,"betamax/D");
  trout->Branch("EtaMinPF",&betamin,"betamin/D");
  trout->Branch("EtaMinPFCastor",&betamincastor,"betamincastor/D");

  TTreeCASTOR = "TTreeCASTORW_" + savehistofile;
  fOutCASTOR = new TFile(TTreeCASTOR, "RECREATE");
  fOutCASTOR->cd();
  troutCASTOR = trout->CloneTree(0);

  TTreeAllW = "TTreeAllW_" + savehistofile;
  fOutW = new TFile(TTreeAllW, "RECREATE");
  fOutW->cd();
  troutW = trout->CloneTree(0);

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

    // Clean Variables
    CleanVariables();

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
      if(eventdiffW->GetHLTPath(nt) > 0){
	triggercounter[nt]++;
      }
    }

    double totalASum = eventdiff->GetSumEnergyHFPlus() + eventdiff->GetSumEnergyHFMinus();
    if (totalASum > 0.){
      aSumE = (eventdiff->GetSumEnergyHFPlus() - eventdiff->GetSumEnergyHFMinus())/(eventdiff->GetSumEnergyHFPlus() + eventdiff->GetSumEnergyHFMinus());
    }

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
      loadBar(i,NEVENTS);
    }

    if (debug){
      std::cout << " " << std::endl;
      std::cout << "# of active towers: " << eventdiffW->GetEachTowerCounter() << std::endl;
      std::cout << " " <<std::endl;
      for (int ic=0; ic<eventdiffW->GetEachTowerCounter(); ic++){
	std::cout << "Each Tower Energy: " << eventdiffW->GetEachTowerEnergy(ic) << std::endl;
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
    bool dimass = false;
    bool eleBarrel1 = false;
    bool eleEndCap1 = false;
    bool candSel = false;
    bool isoBarrel1 = false;
    bool isoEndCap1 = false;
    bool isolation = false;
    bool WKinN = false;
    bool WKinP = false;
    bool acceptEvt = false;

    if (switchtrigger == "trigger_all_electron"){
      if ( (eventdiff->GetRunNumber() >= 132440 && eventdiff->GetRunNumber() <= 137028) && eventdiffW->GetHLTPath(0) > 0 ) triggerE_a = true;
      if ( (eventdiff->GetRunNumber() >= 138564 && eventdiff->GetRunNumber() <= 140401) && eventdiffW->GetHLTPath(1) > 0) triggerE_b = true;
      if ( (eventdiff->GetRunNumber() >= 141956 && eventdiff->GetRunNumber() <= 144114) && eventdiffW->GetHLTPath(2) > 0) triggerE_c = true;
      if ( (eventdiff->GetRunNumber() >= 144115 && eventdiff->GetRunNumber() <= 147145) && eventdiffW->GetHLTPath(3) > 0) triggerE_d = true;
      if ( (eventdiff->GetRunNumber() >= 147146 && eventdiff->GetRunNumber() <= 148058) && eventdiffW->GetHLTPath(4) > 0) triggerE_e = true;
      if ( (eventdiff->GetRunNumber() >= 148103 && eventdiff->GetRunNumber() <= 149065) && eventdiffW->GetHLTPath(5) > 0) triggerE_f = true;
      if ( (eventdiff->GetRunNumber() >= 149180 && eventdiff->GetRunNumber() <= 149442) && eventdiffW->GetHLTPath(6) > 0) triggerE_g = true;
      if (triggerE_a || triggerE_b || triggerE_c || triggerE_d || triggerE_e || triggerE_f || triggerE_g) trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Electron accepted." << std::endl;
      TriggerStatus = "trigger_all_electron";
    }
    else if (switchtrigger == "trigger_all_muon"){
      if (eventdiffW->GetHLTPath(0) > 0|| eventdiffW->GetHLTPath(1) > 0) trigger = true;
      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Muon accepted." << std::endl;
      TriggerStatus = "trigger_all_muon";
    }
    else if (switchtrigger == "trigger"){
      if (eventdiffW->GetHLTPath(optTrigger) > 0) trigger = true;
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
	  if (eventdiffW->GetCastorModule2Energy(i)*energycorr[1][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule2Energy(i)*energycorr[1][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule3Energy(i)*energycorr[2][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule3Energy(i)*energycorr[2][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule4Energy(i)*energycorr[3][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule4Energy(i)*energycorr[3][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule5Energy(i)*energycorr[4][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule5Energy(i)*energycorr[4][i]; ++counterHit;}
	}else{
	  if (eventdiffW->GetCastorModule1Energy(i)*energycorr[0][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule1Energy(i)*energycorr[0][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule2Energy(i)*energycorr[1][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule2Energy(i)*energycorr[1][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule3Energy(i)*energycorr[2][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule3Energy(i)*energycorr[2][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule4Energy(i)*energycorr[3][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule4Energy(i)*energycorr[3][i]; ++counterHit;}
	  if (eventdiffW->GetCastorModule5Energy(i)*energycorr[4][i] > channelsthreshold){ CastorEnergySector[i]+=eventdiffW->GetCastorModule5Energy(i)*energycorr[4][i]; ++counterHit;}
	}
      }

      for (int l=0; l<16;l++){
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
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule2Energy(i)*energycorr[1][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule3Energy(i)*energycorr[2][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule4Energy(i)*energycorr[3][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule5Energy(i)*energycorr[4][i];
	}else{
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule1Energy(i)*energycorr[0][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule2Energy(i)*energycorr[1][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule3Energy(i)*energycorr[2][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule4Energy(i)*energycorr[3][i];
	  CastorEnergySector[i]+=eventdiffW->GetCastorModule5Energy(i)*energycorr[4][i];
	}
      }
      for (int l=0; l<16;l++){
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

    sumCastorAndHFMinusEnergy = sumCastorEnergy+eventdiff->GetSumEnergyHFMinus();

    // Redefinition of etamin_
    //------------------------

    // CMS and CASTOR acceptance
    if (castoractivity) {
      etamin_ = -6.;
    }else{
      etamin_ = eventdiff->GetEtaMinFromPFCands();
    }
    //----->

    deltaetapf = eventdiff->GetEtaMaxFromPFCands() - eventdiff->GetEtaMinFromPFCands();
    deltaetapfcastor = eventdiff->GetEtaMaxFromPFCands() - etamin_;

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

    if (typesel == "RecoElectron"){
      selStatus = "Reco::Electron";

      bosonWMass = eventdiffW->GetBosonElectronMass();
      bosonWEta = eventdiffW->GetLeadingElectronEta();
      bosonWPhi = eventdiffW->GetLeadingElectronPhi();
      bosonWPt = eventdiffW->GetLeadingElectronPt();
      bosonWCharge = eventdiffW->GetLeadingElectronCharge();
      NMuons = eventdiffW->GetMuonsN();
      NElectrons = eventdiffW->GetElectronsN();

      isoTk1 = eventdiffW->GetLeadingElectronTkDr03()/eventdiffW->GetLeadingElectronPt();
      isoEcal1 = eventdiffW->GetLeadingElectronEcalDr03()/eventdiffW->GetLeadingElectronPt();
      isoHcal1 = eventdiffW->GetLeadingElectronHcalDr03()/eventdiffW->GetLeadingElectronPt();
      innerHits1 = eventdiffW->GetLeadingElectronInnerHits();
      Dcot1 = eventdiffW->GetLeadingElectronDCot();
      Dist1 = eventdiffW->GetLeadingElectronDist();
      DeltaEtaTkClu1 = eventdiffW->GetLeadingElectronDeltaEtaTkClu();
      DeltaPhiTkClu1 = eventdiffW->GetLeadingElectronDeltaPhiTkClu();
      sigmaIeIe1 = eventdiffW->GetLeadingElectronSigmaIeIe();
      HE1 = eventdiffW->GetLeadingElectronHE();

      double totalASumCastor = eventdiff->GetSumEnergyHFMinus() + sumCastorEnergy;

      if(totalASumCastor > 0.){
	AEcastor = (eventdiff->GetSumEnergyHFMinus() - sumCastorEnergy)/(eventdiff->GetSumEnergyHFMinus() + sumCastorEnergy);
      }

      if (eventdiffW->GetLeadingElectronPt() > lepton1pt && eventdiffW->GetMETPt() > lepton2pt && eventdiffW->GetLeadingMuonPt()<10) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;

      if(NElectrons==1) acceptEvt = true;

      //Isolation Electron
      if ((fabs (eventdiffW->GetLeadingElectronEta()) <= 1.4442) ){
	if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1 = true;
      }

      if ((fabs (eventdiffW->GetLeadingElectronEta()) >= 1.5660) && (fabs (eventdiffW->GetLeadingElectronEta()) <= 2.5)){
	if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1 = true;
      }

      if ((isoEndCap1 || isoBarrel1)) isolation = true;

      //Electron Candidate Selection
      if ((fabs (eventdiffW->GetLeadingElectronEta()) <= 1.4442) ){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.004 && fabs (DeltaPhiTkClu1) < 0.06 && sigmaIeIe1 < 0.01 && HE1 < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (eventdiffW->GetLeadingElectronEta()) >= 1.5660) && (fabs (eventdiffW->GetLeadingElectronEta()) <= 2.5)){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.007 && fabs (DeltaPhiTkClu1) < 0.03 && sigmaIeIe1 < 0.03 && HE1 < 0.025) eleEndCap1 = true;
      }

      if ((eleEndCap1 || eleBarrel1)) candSel = true;

      if (eventdiffW->GetLeadingElectronEta()>0.) WKinP = true;
      if (eventdiffW->GetLeadingElectronEta()<0.) WKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetLeadingElectronEta()); }
      else{ etasignedHF = fabs(eventdiffW->GetLeadingElectronEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetLeadingElectronEta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetLeadingElectronEta());}

    }

    else if (typesel == "RecoMuon"){
      selStatus = "Reco::Muon";

      bosonWMass = eventdiffW->GetBosonMuonMass();
      bosonWEta = eventdiffW->GetLeadingMuonEta();
      bosonWPhi = eventdiffW->GetLeadingMuonPhi();
      bosonWPt = eventdiffW->GetLeadingMuonPt();
      bosonWCharge = eventdiffW->GetLeadingMuonCharge();
      NMuons = eventdiffW->GetMuonsN();
      NElectrons = eventdiffW->GetElectronsN();
      isoRec = eventdiffW->GetLeadingMuonSumPtR03();

      if(NMuons ==1) acceptEvt = true;

      if (eventdiffW->GetLeadingMuonPt() > lepton1pt && eventdiffW->GetMETPt() > lepton2pt && eventdiffW->GetLeadingElectronPt()<10) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;
      if (isoRec < 3) { 
	isolation = true;
	candSel = true;
      }

      if (eventdiffW->GetLeadingMuonEta()>0.) WKinP = true;
      if (eventdiffW->GetLeadingMuonEta()<0.) WKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetLeadingMuonEta()); }
      else{ etasignedHF = fabs(eventdiffW->GetLeadingMuonEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetLeadingMuonEta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetLeadingMuonEta());}

    }

    else if (typesel == "PatElectron"){
      selStatus = "Pat::Electron";

      bosonWMass = eventdiffW->GetPatBosonElectronMass();
      bosonWEta = eventdiffW->GetPatElectron1Eta();
      bosonWPhi = eventdiffW->GetPatElectron1Phi();
      bosonWPt = eventdiffW->GetPatElectron1Pt();
      bosonWCharge = eventdiffW->GetPatMuon1Charge();
      NMuons = eventdiffW->GetPatNMuon();
      NElectrons = eventdiffW->GetPatNElectron();

      isoTk1 = eventdiffW->GetPatElectron1TkDr03()/eventdiffW->GetPatElectron1Pt();
      isoEcal1 = eventdiffW->GetPatElectron1EcalDr03()/eventdiffW->GetPatElectron1Pt();
      isoHcal1 = eventdiffW->GetPatElectron1HcalDr03()/eventdiffW->GetPatElectron1Pt();
      innerHits1 = eventdiffW->GetPatElectron1InnerHits();
      Dcot1 = eventdiffW->GetPatElectron1DCot();
      Dist1 = eventdiffW->GetPatElectron1Dist();
      DeltaEtaTkClu1 = eventdiffW->GetPatElectron1DeltaEtaTkClu();
      DeltaPhiTkClu1 = eventdiffW->GetPatElectron1DeltaPhiTkClu();
      sigmaIeIe1 = eventdiffW->GetPatElectron1SigmaIeIe();
      HE1 = eventdiffW->GetPatElectron1HE();

      if (eventdiffW->GetPatElectron1Pt() > lepton1pt && eventdiffW->GetPatMETPt() > lepton2pt && eventdiffW->GetPatMuon1Pt()<10) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;

      if(NElectrons==1) acceptEvt = true;

      //Isolation Electron
      if ((fabs (eventdiffW->GetPatElectron1Eta()) <= 1.4442) ){
	if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1 = true;
      }

      if ((fabs (eventdiffW->GetPatElectron1Eta()) >= 1.5660) && (fabs (eventdiffW->GetPatElectron1Eta()) <= 2.5)){
	if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1 = true;
      }

      if ((isoEndCap1 || isoBarrel1)) isolation = true;

      //Electron Candidate Selection
      if ((fabs (eventdiffW->GetPatElectron1Eta()) <= 1.4442) ){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.004 && fabs (DeltaPhiTkClu1) < 0.06 && sigmaIeIe1 < 0.01 && HE1 < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (eventdiffW->GetPatElectron1Eta()) >= 1.5660) && (fabs (eventdiffW->GetPatElectron1Eta()) <= 2.5)){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.007 && fabs (DeltaPhiTkClu1) < 0.03 && sigmaIeIe1 < 0.03 && HE1 < 0.025) eleEndCap1 = true;
      }

      if ((eleEndCap1 || eleBarrel1)) candSel = true;

      if (eventdiffW->GetPatElectron1Eta()>0.) WKinP = true;
      if (eventdiffW->GetPatElectron1Eta()<0.) WKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetPatElectron1Eta()); }
      else{ etasignedHF = fabs(eventdiffW->GetPatElectron1Eta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetPatElectron1Eta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetPatElectron1Eta());}

    }

    else if(typesel == "PatMuon"){
      selStatus = "Pat::Muon";

      bosonWMass = eventdiffW->GetPatBosonMuonMass();
      bosonWEta = eventdiffW->GetPatMuon1Eta();
      bosonWPhi = eventdiffW->GetPatMuon1Phi();
      bosonWPt = eventdiffW->GetPatMuon1Pt();
      bosonWCharge = eventdiffW->GetPatMuon1Charge();
      NMuons = eventdiffW->GetPatNMuon();
      NElectrons = eventdiffW->GetPatNElectron();
      isoRec = eventdiffW->GetPatMuon1SumPtR03();

      if(NMuons ==1) acceptEvt = true;

      if (eventdiffW->GetPatMuon1Pt() > lepton1pt && eventdiffW->GetPatMETPt() > lepton2pt && eventdiffW->GetPatElectron1Pt()<10) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;
      if (isoRec < 3) {
	candSel = true;
	isolation = true;
      }

      if (eventdiffW->GetPatMuon1Eta()>0.) WKinP = true;
      if (eventdiffW->GetPatMuon1Eta()<0.) WKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetPatMuon1Eta()); }
      else{ etasignedHF = fabs(eventdiffW->GetPatMuon1Eta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetPatMuon1Eta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetPatMuon1Eta());}
    }

    else{
      exit(EXIT_FAILURE);
    }

    //Branch Defining
    bRunNumber = eventdiff->GetRunNumber();
    bLumiSection = eventdiff->GetLumiSection();
    bEventNumber = eventdiff->GetEventNumber();
    if (typesel == "RecoMuon"){
      bDiBosonPt = eventdiffW->GetLeadingMuonPt();
      bDiBosonEta = eventdiffW->GetLeadingMuonEta();
      bDiBosonPhi = eventdiffW->GetLeadingMuonPhi();
      bDiBosonMass = bosonWMass;
    }
    if (typesel == "PatMuon"){
      bDiBosonPt = eventdiffW->GetPatMuon1Pt();
      bDiBosonEta = eventdiffW->GetPatMuon1Eta();
      bDiBosonPhi = eventdiffW->GetPatMuon1Phi();
      bDiBosonMass = bosonWMass;
    }
    if (typesel == "RecoElectron"){
      bDiBosonPt = eventdiffW->GetLeadingElectronPt();
      bDiBosonEta = eventdiffW->GetLeadingElectronEta();
      bDiBosonPhi = eventdiffW->GetLeadingElectronPhi();
      bDiBosonMass = bosonWMass;
    }
    if (typesel == "PatElectron"){
      bDiBosonPt = eventdiffW->GetPatElectron1Pt();
      bDiBosonEta = eventdiffW->GetPatElectron1Eta();
      bDiBosonPhi = eventdiffW->GetPatElectron1Phi();
      bDiBosonMass = bosonWMass;
    }
    bMultiplicityTracks = eventdiff->GetMultiplicityTracks();
    bSumEEEMinus = eventdiffW->GetSumEEEMinus();
    bSumEEEPlus = eventdiffW->GetSumEEEPlus();
    bSumEnergyHFMinus = eventdiff->GetSumEnergyHFMinus();
    bSumEnergyHFPlus = eventdiff->GetSumEnergyHFPlus();
    bsumCastorEnergy = sumCastorEnergy;
    bSectorCastorHit = counterHit;
    bMaxGapPF = eventdiffW->GetLrgPF();
    if(fabs(eventdiffW->GetSumptPFLeft())>fabs(eventdiffW->GetSumptPFRight())){
      SumPTMaxLrgPF = eventdiffW->GetSumptPFLeft();
      SumPTMinLrgPF = eventdiffW->GetSumptPFRight();
    }else{
      SumPTMaxLrgPF = eventdiffW->GetSumptPFRight();
      SumPTMinLrgPF = eventdiffW->GetSumptPFLeft();
    }
    bXiPlusFromPFCands = eventdiff->GetXiPlusFromPFCands();
    bXiMinusFromPFCands = eventdiff->GetXiMinusFromPFCands();
    betasignedHF = etasignedHF;
    betasignedCASTOR = etasignedCASTOR;
    bAEcastor = AEcastor;
    betamax = eventdiff->GetEtaMaxFromPFCands();
    betamin = eventdiff->GetEtaMinFromPFCands();
    betamincastor = etamin_;

    if(pileup < 21){ // Never comment this line. It is the program defense.

      if(switchtrigger == "trigger" || switchtrigger == "trigger_all_electron" || switchtrigger == "trigger_all_muon"){ 
	FillHistos(0,pileup,totalcommon); 
	if(trigger) {
	  ++totalT;
	  FillHistos(1,pileup,totalcommon);
	} 
	if(trigger && vertex && presel) FillHistos(2,pileup,totalcommon);
	if(trigger && vertex && presel && isolation) FillHistos(3,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel) FillHistos(4,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt) {
	  FillHistos(5,pileup,totalcommon);
	  fOutW->cd();
	  troutW->Fill();
	}
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln) FillHistos(6,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp) FillHistos(7,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln && castorgap) FillHistos(8,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity) FillHistos(9,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln && WKinP){
	  outstring << "HF- Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(10,pileup,totalcommon);
	}
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && WKinN){
	  outstring << "HF+ Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(11,pileup,totalcommon);
	}
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln && castorgap && WKinP) FillHistos(12,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity && WKinN) FillHistos(13,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && castorgap){
	  FillHistos(14,pileup,totalcommon);
	  fOutCASTOR->cd();
	  troutCASTOR->Fill();
	}
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && castorgap && WKinP){
	  outstring << "CASTOR Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(15,pileup,totalcommon);
	  fOut->cd();
	  trout->Fill();
	}

	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castorgap) FillHistos(16,pileup,totalcommon);
	if(trigger && vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castorgap && WKinP) FillHistos(17,pileup,totalcommon);

      }

      else if (switchtrigger =="no_trigger_nocorrection" || switchtrigger == "no_trigger_correction" ){
	--totalT;
	FillHistos(0,pileup,totalcommon);
	if(vertex && presel) FillHistos(2,pileup,totalcommon);
	if(vertex && presel && isolation) FillHistos(3,pileup,totalcommon);
	if(vertex && presel && isolation && candSel) FillHistos(4,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt) {
	  FillHistos(5,pileup,totalcommon);
	  fOutW->cd();
	  troutW->Fill();
	}
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln) FillHistos(6,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp) FillHistos(7,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln && castorgap) FillHistos(8,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity) FillHistos(9,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln && WKinP) FillHistos(10,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && WKinN) FillHistos(11,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffseln && castorgap && WKinP) FillHistos(12,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity && WKinN) FillHistos(13,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && castorgap){ 
	  FillHistos(14,pileup,totalcommon);
	  fOutCASTOR->cd();
	  troutCASTOR->Fill();
	}
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && castorgap && WKinP){
	  FillHistos(15,pileup,totalcommon);
	  fOut->cd();
	  trout->Fill();
	}
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castorgap) FillHistos(16,pileup,totalcommon);
	if(vertex && presel && isolation && candSel && dimass && acceptEvt && diffselp && castorgap && WKinP) FillHistos(17,pileup,totalcommon);
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
  foldersFile[0] = outf->mkdir("BosonKinematics");
  foldersFile[1] = outf->mkdir("Detector");
  foldersFile[2] = outf->mkdir("Run");
  foldersFile[3] = outf->mkdir("Diffraction");
  SaveHistos(type,typesel);
  outf->Close();
  std::cout << "\n" << std::endl;

  fOut->cd();
  trout->SetWeight(totalweight);
  trout->Write();
  fOut->Close();

  fOutW->cd();
  troutW->SetWeight(totalweight);
  troutW->Write();
  fOutW->Close();

  fOutCASTOR->cd();
  troutCASTOR->SetWeight(totalweight);
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

    DiffractiveW* diffWRun = new DiffractiveW();
    clock_t tStart = clock();
    diffWRun->CreateHistos(type_);
    diffWRun->Run(filein_, processname_, savehistofile_, switchtrigger_, optTrigger_, lepton1pt_, lepton2pt_, nVertex_, type_, switchlumiweight_, mcweight_, typesel_, castorthreshold_, channelsthreshold_, castorcorrfile_, gapseltype_, pumfile_, pudfile_);
    std::cout<< "Time taken: " << (double)(clock() - tStart)/(60*CLOCKS_PER_SEC) << " min" << std::endl;
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
