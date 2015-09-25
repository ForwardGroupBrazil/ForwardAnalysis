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
#include <TF1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TMath.h>
#include <math.h>

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

  leptoneta = -999.;
  leptonpt = -999.;
  leptonphi = -999.;
  deltaphi = -999.;
  bosonWMass = -999.;
  bosonWEta = -999.;
  bosonWEtaGen = -999.;
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
  deltaeta = -999.;
  etamax = -999.;
  etamin = -999.;
  xi = -999.;
  xiplus = -999.;
  ximinus = -999.;
  xigen = -999.;
  xigenplus = -999.;
  xigenminus = -999.;
  ximinusunique = -999.;
  maxLRG = -999.;
  resLeadingPt = -999.;
  resLeadingEta = -999.;
  resLeadingPhi = -999.;
  resMETPt = -999.;
  resMETPhi = -999.;

}

void DiffractiveW::CreateHistos(std::string type){

  //double binarrayplus[19] = {0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,5.2,6.2};
  //double binarrayminus[19] = {-6.2,-5.2,-4.,-3.75,-3.5,-3.25,-3.,-2.75,-2.5,-2.25,-2.,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.};
  double binarraydelta[15] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,9.};
  double xi_bin[18]={0.0003,0.002,0.0045,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  int nbinxi = 17;

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
  std::string step12 = "NGapCMSAndWKinP";
  std::string step13 = "PGapCMSAndWKinN";
  std::string step14 = "NGapCMSAndCASTORAndWKinP";
  std::string step15 = "PGapCMSAndCastorActivityAndWKinN";
  std::string step16 = "NGapCASTOR";
  std::string step17 = "NGapCASTORAndWKinP";
  std::string step18 = "PGapCMSAndCASTOR";
  std::string step19 = "PGapCMSAndCASTORAndWKinP";
  std::string step20 = "generator";
  std::string step21 = "gencastorlow";
  std::string step22 = "gencastorhigh";
  std::string step23 = "zeropileup";
  std::string step24 = "zeropileup_NGapCASTOR";
  std::string step25 = "xi";
  std::string step26 = "generatorW";
  std::string step27 = "generatorWDiff";

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
  Folders.push_back(step20);
  Folders.push_back(step21);
  Folders.push_back(step22);
  Folders.push_back(step23);
  Folders.push_back(step24);
  Folders.push_back(step25);
  Folders.push_back(step26);
  Folders.push_back(step27);

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
    m_hVector_NElectrons.push_back( std::vector<TH1D*>() );
    m_hVector_NMuons.push_back( std::vector<TH1D*>() );
    m_hVector_WMass.push_back( std::vector<TH1D*>() );
    m_hVector_WEta.push_back( std::vector<TH1D*>() );
    m_hVector_WPt.push_back( std::vector<TH1D*>() );
    m_hVector_WPhi.push_back( std::vector<TH1D*>() );
    m_hVector_WCharge.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonTkDr03.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonEcalDr03.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonHcalDr03.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonIsolation.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonInnerHits.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonDCot.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonDist.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonDeltaEtaTkClu.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonDeltaPhiTkClu.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonSigmaIeIe.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonHE.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonPt.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonEta.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonPhi.push_back( std::vector<TH1D*>() );
    m_hVector_SecondLeptonPt.push_back( std::vector<TH1D*>() );
    m_hVector_SecondLeptonEta.push_back( std::vector<TH1D*>() );
    m_hVector_SecondLeptonPhi.push_back( std::vector<TH1D*>() );
    m_hVector_DeltaPhi.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonDxy.push_back( std::vector<TH1D*>() );
    m_hVector_LeadingLeptonDz.push_back( std::vector<TH1D*>() );
    m_hVector_METPt.push_back( std::vector<TH1D*>() );
    m_hVector_METSignificance.push_back( std::vector<TH1D*>() );
    m_hVector_METPhi.push_back( std::vector<TH1D*>() );
    m_hVector_DVtxMuon.push_back( std::vector<TH1D*>() );
    m_hVector_DVtxMuonZ.push_back( std::vector<TH1D*>() );
    m_hVector_DVtxElectron.push_back( std::vector<TH1D*>() );
    m_hVector_DVtxElectronZ.push_back( std::vector<TH1D*>() );
    m_hVector_DMuonElectron.push_back( std::vector<TH1D*>() );
    m_hVector_DMuonElectronZ.push_back( std::vector<TH1D*>() );
    m_hVector_DMuons.push_back( std::vector<TH1D*>() );
    m_hVector_DElectrons.push_back( std::vector<TH1D*>() );


    // Detector
    m_hVector_sumEHFplus.push_back( std::vector<TH1D*>() );
    m_hVector_sumEHFminus.push_back( std::vector<TH1D*>() );
    m_hVector_sumEHEplus.push_back( std::vector<TH1D*>() );
    m_hVector_sumEHEminus.push_back( std::vector<TH1D*>() );
    m_hVector_SumEHFplus_S.push_back( std::vector<TH1D*>() );
    m_hVector_SumEHFminus_S.push_back( std::vector<TH1D*>() );
    m_hVector_SumEHFplus_L.push_back( std::vector<TH1D*>() );
    m_hVector_SumEHFminus_L.push_back( std::vector<TH1D*>() );
    m_hVector_sumEEEplus.push_back( std::vector<TH1D*>() );
    m_hVector_sumEEEminus.push_back( std::vector<TH1D*>() );
    m_hVector_EnergyHFPlusVsEnergyHFMinus.push_back( std::vector<TH2D*>() );
    m_hVector_EnergyEEPlusVsEnergyEEMinus.push_back( std::vector<TH2D*>() );
    m_hVector_SumEHFMax.push_back( std::vector<TH1D*>() );
    m_hVector_SumEHFMin.push_back( std::vector<TH1D*>() );
    m_hVector_etcalos_p.push_back( std::vector<TH2D*>() );
    m_hVector_etcalos_n.push_back( std::vector<TH2D*>() );
    m_hVector_ECaloVsEta.push_back( std::vector<TH2D*>() );
    m_hVector_ECaloVsEtaTProf.push_back( std::vector<TProfile*>() );
    m_hVector_EnergyVsEtaBin1D.push_back( std::vector<TH1D*>() );
    m_hVector_sumECastorMinus.push_back( std::vector<TH1D*>() );
    m_hVector_ECastorSector.push_back( std::vector<TH2D*>() );
    m_hVector_ECastorSectorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_ECastorSectorBin1D.push_back( std::vector<TH1D*>() );
    m_hVector_EnergyHFMinusVsCastorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_EnergyHFPlusVsCastorTProf.push_back( std::vector<TProfile*>() );
    m_hVector_sumECastorAndHFMinus.push_back( std::vector<TH1D*>() );
    m_hVector_CastorMultiplicity.push_back( std::vector<TH1D*>() );
    m_hVector_CastorMultiplicityVsLumi.push_back( std::vector<TH2D*>() );
    m_hVector_SectorVsTotalCastorEnergy.push_back( std::vector<TH2D*>() );
    m_hVector_SectorVsTotalCastorEnergyTProf.push_back( std::vector<TProfile*>() );
    m_hVector_sumEHFplusPF.push_back( std::vector<TH1D*>() );
    m_hVector_sumEHFminusPF.push_back( std::vector<TH1D*>() );
    m_hVector_EPFVsTowerMinus.push_back( std::vector<TH2D*>() );
    m_hVector_EPFVsTowerPlus.push_back( std::vector<TH2D*>() );

    // Event Information
    m_hVector_lumi.push_back( std::vector<TH1D*>() );
    m_hVector_tracks.push_back( std::vector<TH1D*>() );
    m_hVector_vertex.push_back( std::vector<TH1D*>() );
    m_hVector_pu.push_back( std::vector<TH1D*>() );


    // Diffraction
    m_hVector_asumE.push_back( std::vector<TH1D*>() );
    m_hVector_multhf.push_back( std::vector<TH2D*>() );
    m_hVector_minusnhf.push_back( std::vector<TH1D*>() );
    m_hVector_plusnhf.push_back( std::vector<TH1D*>() );
    m_hVector_etaminimum.push_back( std::vector<TH1D*>() );
    m_hVector_etamax.push_back( std::vector<TH1D*>() );
    m_hVector_etamin.push_back( std::vector<TH1D*>() );
    m_hVector_maxetagap.push_back( std::vector<TH1D*>() );
    m_hVector_SumPTMax.push_back( std::vector<TH1D*>() );
    m_hVector_SumPTMin.push_back( std::vector<TH1D*>() );
    m_hVector_absdeltaEta.push_back( std::vector<TH1D*>() );
    m_hVector_deltaEta.push_back( std::vector<TH1D*>() );
    m_hVector_XiPlus.push_back( std::vector<TH1D*>() );
    m_hVector_XiMinus.push_back( std::vector<TH1D*>() );
    m_hVector_XiPlusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_XiMinusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_XiPlusAll.push_back( std::vector<TH1D*>() );
    m_hVector_XiMinusAll.push_back( std::vector<TH1D*>() );
    m_hVector_XiMinusUnique.push_back( std::vector<TH1D*>() );
    m_hVector_Xi.push_back( std::vector<TH1D*>() );
    m_hVector_AEcastor.push_back( std::vector<TH1D*>() );
    m_hVector_etasignedHF.push_back( std::vector<TH1D*>() );
    m_hVector_etasignedCASTOR.push_back( std::vector<TH1D*>() );

    // Generator Studies
    m_hVector_genProtonMinusXi.push_back( std::vector<TH1D*>() );
    m_hVector_genProtonPlusXi.push_back( std::vector<TH1D*>() );
    m_hVector_genXiPlus.push_back( std::vector<TH1D*>() );                                                                                                                
    m_hVector_genXiMinus.push_back( std::vector<TH1D*>() );
    m_hVector_genXiMinusMore.push_back( std::vector<TH1D*>() );
    m_hVector_genXiPlusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_genXiMinusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_genXiPlusAll.push_back( std::vector<TH1D*>() );
    m_hVector_genXiMinusAll.push_back( std::vector<TH1D*>() );
    m_hVector_genXi.push_back( std::vector<TH1D*>() );
    m_hVector_resLeadingLeptonPt.push_back( std::vector<TH1D*>() );
    m_hVector_resLeadingLeptonEta.push_back( std::vector<TH1D*>() );
    m_hVector_resLeadingLeptonPhi.push_back( std::vector<TH1D*>() );
    m_hVector_resMETPt.push_back( std::vector<TH1D*>() );
    m_hVector_resMETPhi.push_back( std::vector<TH1D*>() );
    m_hVector_resXiPlus.push_back( std::vector<TH1D*>() );
    m_hVector_resXiMinus.push_back( std::vector<TH1D*>() );
    m_hVector_resXiMinusMore.push_back( std::vector<TH1D*>() );
    m_hVector_resXiPlusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_resXiMinusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_resXiPlusAll.push_back( std::vector<TH1D*>() );
    m_hVector_resXiMinusAll.push_back( std::vector<TH1D*>() );
    m_hVector_resAbsXiPlus.push_back( std::vector<TH1D*>() );
    m_hVector_resAbsXiMinus.push_back( std::vector<TH1D*>() );
    m_hVector_resAbsXiPlusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_resAbsXiMinusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_resAbsXiPlusAll.push_back( std::vector<TH1D*>() );
    m_hVector_resAbsXiMinusAll.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiPlus.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiMinus.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiMinusMore.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiPlusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiMinusCastor.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiPlusAll.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiMinusAll.push_back( std::vector<TH1D*>() );
    m_hVector_ratioXiMinusUnique.push_back( std::vector<TH1D*>() );
    m_hVector_resHFEnergy.push_back( std::vector<TH1D*>() );
    m_hVector_resCASTOREnergy.push_back( std::vector<TH1D*>() );
    m_hVector_ratioCASTOREnergy.push_back( std::vector<TH1D*>() );
    m_hVector_ratioCASTORPionsEnergy.push_back( std::vector<TH1D*>() );
    m_hVector_ratioCASTOREPhotonEnergy.push_back( std::vector<TH1D*>() );
    m_hVector_ratioCASTOROthersEnergy.push_back( std::vector<TH1D*>() );
    m_hVector_correlHFEnergy.push_back( std::vector<TH2D*>() );
    m_hVector_correlCASTOREnergy.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiPlus.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinus.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusMore.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiPlusCastor.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusCastor.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusCastorDef.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiPlusAll.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusAll.push_back( std::vector<TH2D*>() );
    m_hVector_correlRatioCastor.push_back( std::vector<TH2D*>() );
    m_hVector_correlRatioCastorPion.push_back( std::vector<TH2D*>() );
    m_hVector_correlRatioCastorEPhoton.push_back( std::vector<TH2D*>() );
    m_hVector_correlRatioCastorOthers.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusRecoGen.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusCastorRecoGen.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusAllRecoGen.push_back( std::vector<TH2D*>() );
    m_hVector_correlXiMinusUniqueRecoGen.push_back( std::vector<TH2D*>() );
    m_hVector_gensumECastorMinus.push_back( std::vector<TH1D*>() );
    m_hVector_gensumECastorPionsMinus.push_back( std::vector<TH1D*>() );
    m_hVector_gensumECastorEPhotonMinus.push_back( std::vector<TH1D*>() );
    m_hVector_gensumECastorOthersMinus.push_back( std::vector<TH1D*>() );
    m_hVector_gensumEHFplus.push_back( std::vector<TH1D*>() );
    m_hVector_gensumEHFminus.push_back( std::vector<TH1D*>() );
    m_hVector_WEtaGen.push_back( std::vector<TH1D*>() );
    m_hVector_NumberOfProtons.push_back( std::vector<TH1D*>() );
    m_hVector_ProtonPz.push_back( std::vector<TH1D*>() );

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
      TH1D *histo_NMuons = new TH1D(name,"Number of Muons per Event; # Muons; Multiplicity", 100, 0., 100.);
      m_hVector_NMuons[j].push_back(histo_NMuons);

      sprintf(name,"NElectrons_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_NElectrons = new TH1D(name,"Number of Electrons per Event; # Electrons; Multiplicity", 100, 0., 100.);
      m_hVector_NElectrons[j].push_back(histo_NElectrons);

      sprintf(name,"BosonWMass_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_WMass = new TH1D(name,"Boson W Transverse Mass Distribution; M_{T}(l#nu) [GeV]; N events",500,0,500);
      m_hVector_WMass[j].push_back(histo_WMass);

      sprintf(name,"BosonWPt_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_WPt = new TH1D(name,"Boson W Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_WPt[j].push_back(histo_WPt);

      sprintf(name,"BosonWEta_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_WEta = new TH1D(name,"Boson W #eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_WEta[j].push_back(histo_WEta);

      sprintf(name,"BosonWPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_WPhi = new TH1D(name,"Boson W #phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_WPhi[j].push_back(histo_WPhi);

      sprintf(name,"BosonWCharge_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_WCharge = new TH1D(name,"Boson W Charge Distribution; Charge; N events",6,-3,3);
      m_hVector_WCharge[j].push_back(histo_WCharge);

      sprintf(name,"LeadingLeptonTkDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonTkDr03 = new TH1D(name,"Leading Lepton: Tracker Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonTkDr03[j].push_back(histo_LeadingLeptonTkDr03);

      sprintf(name,"LeadingLeptonEcalDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonEcalDr03 = new TH1D(name,"Leading Lepton: ECAL Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonEcalDr03[j].push_back(histo_LeadingLeptonEcalDr03);

      sprintf(name,"LeadingLeptonHcalDr03_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonHcalDr03 = new TH1D(name,"Leading Lepton: HCAL Isolation DR03; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonHcalDr03[j].push_back(histo_LeadingLeptonHcalDr03);

      sprintf(name,"LeadingLeptonIsolation_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonIsolation = new TH1D(name,"Leading Lepton: Isolation; # Isolation; [u]", 100, 0., 1.);
      m_hVector_LeadingLeptonIsolation[j].push_back(histo_LeadingLeptonIsolation);

      sprintf(name,"LeadingLeptonInnerHits_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonInnerHits = new TH1D(name,"Leading Lepton; Number Of Expected Inner Hits; N events",500,0,50);
      m_hVector_LeadingLeptonInnerHits[j].push_back(histo_LeadingLeptonInnerHits);

      sprintf(name,"LeadingLeptonDCot_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonDCot = new TH1D(name,"Leading Lepton; DCot [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDCot[j].push_back(histo_LeadingLeptonDCot);

      sprintf(name,"LeadingLeptonDist_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonDist = new TH1D(name,"Leading Lepton; Dist [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDist[j].push_back(histo_LeadingLeptonDist);

      sprintf(name,"LeadingLeptonDeltaEtaTkClu_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonDeltaEtaTkClu = new TH1D(name,"Leading Lepton; #Delta#eta_{TK, Cluster} [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDeltaEtaTkClu[j].push_back(histo_LeadingLeptonDeltaEtaTkClu);

      sprintf(name,"LeadingLeptonDeltaPhiTkClu_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonDeltaPhiTkClu = new TH1D(name,"Leading Lepton; #Delta#phi_{TK, Cluster} [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonDeltaPhiTkClu[j].push_back(histo_LeadingLeptonDeltaPhiTkClu);

      sprintf(name,"LeadingLeptonSigmaIeIe_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonSigmaIeIe = new TH1D(name,"Leading Lepton; #sigma_{i#etai#eta} [a.u.]; N events",100,0,1);
      m_hVector_LeadingLeptonSigmaIeIe[j].push_back(histo_LeadingLeptonSigmaIeIe);

      sprintf(name,"LeadingLeptonHE_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonHE = new TH1D(name,"Leading Lepton; HE; N events",100,0,1);
      m_hVector_LeadingLeptonHE[j].push_back(histo_LeadingLeptonHE);

      sprintf(name,"LeadingLeptonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonPt = new TH1D(name,"Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_LeadingLeptonPt[j].push_back(histo_LeadingLeptonPt);

      sprintf(name,"LeadingLeptonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonEta = new TH1D(name,"#eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_LeadingLeptonEta[j].push_back(histo_LeadingLeptonEta);

      sprintf(name,"LeadingLeptonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonPhi = new TH1D(name,"#phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_LeadingLeptonPhi[j].push_back(histo_LeadingLeptonPhi);

      sprintf(name,"DeltaPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DeltaPhi = new TH1D(name,"#Delta#phi Distribution; #Delta#phi [rad]; N events",60,0.,3.3);
      m_hVector_DeltaPhi[j].push_back(histo_DeltaPhi);

      sprintf(name,"LeadingLeptonDz_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonDz = new TH1D(name,"dz Distribution; dz [mm]; N events",200,0.,0.05);
      m_hVector_LeadingLeptonDz[j].push_back(histo_LeadingLeptonDz);

      sprintf(name,"LeadingLeptonDxy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_LeadingLeptonDxy = new TH1D(name,"dxy Distribution; dxy [mm]; N events",200,0.,0.05);
      m_hVector_LeadingLeptonDxy[j].push_back(histo_LeadingLeptonDxy);

      sprintf(name,"SecondLeptonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SecondLeptonPt = new TH1D(name,"Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_SecondLeptonPt[j].push_back(histo_SecondLeptonPt);

      sprintf(name,"SecondLeptonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SecondLeptonEta = new TH1D(name,"#eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_SecondLeptonEta[j].push_back(histo_SecondLeptonEta);

      sprintf(name,"SecondLeptonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SecondLeptonPhi = new TH1D(name,"#phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_SecondLeptonPhi[j].push_back(histo_SecondLeptonPhi);

      sprintf(name,"METPt_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_METPt = new TH1D(name,"Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0.,1000.);
      m_hVector_METPt[j].push_back(histo_METPt);

      sprintf(name,"METSignificance_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_METSignificance = new TH1D(name,"Significance Distribution; #sigma [GeV.c^{-1}]; N events",200,-500,500.);
      m_hVector_METSignificance[j].push_back(histo_METSignificance);

      sprintf(name,"METPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_METPhi = new TH1D(name,"#phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_METPhi[j].push_back(histo_METPhi);

      sprintf(name,"DVtxMuon_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DVtxMuon = new TH1D(name,"Distance Vertex - Leading Muon; d [cm]; N events",100,-1.,1.);
      m_hVector_DVtxMuon[j].push_back(histo_DVtxMuon);

      sprintf(name,"DVtxMuonZ_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DVtxMuonZ = new TH1D(name,"Distance Z, Vertex - Leading Muon; d [cm]; N events",100,-1.,1.);
      m_hVector_DVtxMuonZ[j].push_back(histo_DVtxMuonZ);

      sprintf(name,"DVtxElectron_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DVtxElectron = new TH1D(name,"Distance Vertex - Leading Electron; d [cm]; N events",100,-1.,1.);
      m_hVector_DVtxElectron[j].push_back(histo_DVtxElectron);

      sprintf(name,"DVtxElectronZ_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DVtxElectronZ = new TH1D(name,"Distance Z, Vertex - Leading Electron; d [cm]; N events",100,-1.,1.);
      m_hVector_DVtxElectronZ[j].push_back(histo_DVtxElectronZ);

      sprintf(name,"DMuonElectron_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DMuonElectron = new TH1D(name,"Distance Leading Muon - Leading Electron; d [cm]; N events",100,-1.,1.);
      m_hVector_DMuonElectron[j].push_back(histo_DMuonElectron);

      sprintf(name,"DMuonElectronZ_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DMuonElectronZ = new TH1D(name,"Distance Z, Leading Muon - Leading Electron; d [cm]; N events",100,-1.,1.);
      m_hVector_DMuonElectronZ[j].push_back(histo_DMuonElectronZ);

      sprintf(name,"DMuons_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DMuons = new TH1D(name,"Distance Leading and Second Leading Muon; d [cm]; N events",100,-1,1.);
      m_hVector_DMuons[j].push_back(histo_DMuons);

      sprintf(name,"DElectrons_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_DElectrons = new TH1D(name,"Distance Leading and Second Leading Electron; d [cm]; N events",100,-1.,1.);
      m_hVector_DElectrons[j].push_back(histo_DElectrons);


      // Detector

      sprintf(name,"sumEHFplus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEHFplus = new TH1D(name,"HF^{+} - Sum of Energy; #sum E_{HF^{+}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHFplus[j].push_back(histo_sumEHFplus);

      sprintf(name,"sumEHFminus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEHFminus = new TH1D(name,"HF^{-} - Sum of Energy; #sum E_{HF^{-}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHFminus[j].push_back(histo_sumEHFminus);

      sprintf(name,"sumEHEplus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEHEplus = new TH1D(name,"HE^{+} - Sum of Energy; #sum E_{HE^{+}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHEplus[j].push_back(histo_sumEHEplus);

      sprintf(name,"sumEHEminus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEHEminus = new TH1D(name,"HE^{-} - Sum of Energy; #sum E_{HE^{-}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHEminus[j].push_back(histo_sumEHEminus);

      sprintf(name,"sumEHFplus_S_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumEHFplus_S = new TH1D(name,"HF^{+} - Sum of Energy, Short Fibers; #sum E_{HF^{+},Short} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFplus_S[j].push_back(histo_SumEHFplus_S);

      sprintf(name,"sumEHFminus_S_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumEHFminus_S = new TH1D(name,"HF^{-} - Sum of Energy, Short Fibers; #sum E_{HF^{-},Short} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFminus_S[j].push_back(histo_SumEHFminus_S);

      sprintf(name,"sumEHFplus_L_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumEHFplus_L = new TH1D(name,"HF^{+} - Sum of Energy, Long Fibers; #sum E_{HF^{+},Long} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFplus_L[j].push_back(histo_SumEHFplus_L);

      sprintf(name,"sumEHFminus_L_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumEHFminus_L = new TH1D(name,"HF^{-} - Sum of Energy, Long Fibers; #sum E_{HF^{-},Long} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFminus_L[j].push_back(histo_SumEHFminus_L);

      sprintf(name,"sumEEEplus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEEEplus = new TH1D(name,"EE^{+} - Sum of Energy; #sum E_{EE^{+}} [GeV]; N events",1000,0.,500.);
      m_hVector_sumEEEplus[j].push_back(histo_sumEEEplus);

      sprintf(name,"sumEEEminus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEEEminus = new TH1D(name,"EE^{-} - Sum of Energy; #sum E_{EE^{-}} [GeV]; N events",1000,0.,500.);
      m_hVector_sumEEEminus[j].push_back(histo_sumEEEminus);

      sprintf(name,"mHF_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_MultHF = new TH2D(name,"HF^{+} and HF^{-} Multiplicity; n HF^{+}; n HF^{-}; N events", 100, 0., 100., 100, 0., 100. );
      m_hVector_multhf[j].push_back(histo_MultHF);

      sprintf(name,"minusNHF_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_minusNHF = new TH1D(name,"HF^{-} Multiplicity; n HF^{-}; N events", 100, 0., 100. );
      m_hVector_minusnhf[j].push_back(histo_minusNHF);

      sprintf(name,"plusNHF_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_plusNHF = new TH1D(name,"HF^{+} Multiplicity; n HF^{+}; N events", 100, 0., 100. );
      m_hVector_plusnhf[j].push_back(histo_plusNHF);

      sprintf(name,"EnergyHFPlusVsEnergyHFMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_EnergyHFPlusVsEnergyHFMinus = new TH2D(name,"HF^{+} and HF^{-}; #sum Energy HF^{+} [GeV]; #sum Energy HF^{-} [GeV]; N events",1000,0.,1000.,1000,0.,1000.);
      m_hVector_EnergyHFPlusVsEnergyHFMinus[j].push_back(histo_EnergyHFPlusVsEnergyHFMinus);

      sprintf(name,"EnergyEEPlusVsEnergyEEMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_EnergyEEPlusVsEnergyEEMinus = new TH2D(name,"EE^{+} and EE^{-}; #sum Energy EE^{+} [GeV]; #sum Energy EE^{-} [GeV]; N events",1000,0.,500.,1000,0.,500.);
      m_hVector_EnergyEEPlusVsEnergyEEMinus[j].push_back(histo_EnergyEEPlusVsEnergyEEMinus);    

      sprintf(name,"sumEHFMax_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumEHFMax = new TH1D(name,"HF - Sum of Energy; #sum E_{HF,Max} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFMax[j].push_back(histo_SumEHFMax);

      sprintf(name,"sumEHFMin_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumEHFMin = new TH1D(name,"HF - Sum of Energy; #sum E_{HF,Min} [GeV]; N events",2000,0,2000);
      m_hVector_SumEHFMin[j].push_back(histo_SumEHFMin);

      sprintf(name,"EnergyHFPlusVsCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_ET_Calos_p = new TH2D(name,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 2000, 0., 1000. );
      m_hVector_etcalos_p[j].push_back(histo_ET_Calos_p);

      sprintf(name,"EnergyHFMinusVsCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_ET_Calos_n = new TH2D(name,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 2000, 0., 1000. );
      m_hVector_etcalos_n[j].push_back(histo_ET_Calos_n);

      sprintf(name,"ECaloVsEta_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_ECaloVsEta = new TH2D(name,"Calorimeter Energy X #eta; #eta; Energy [GeV]", 500, -8, 8, 100, 0., 1000.);
      m_hVector_ECaloVsEta[j].push_back(histo_ECaloVsEta);

      sprintf(name,"ECaloVsEtaTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_ECaloVsEtaTProf = new TProfile(name,"Calorimeter Energy X #eta; #eta; <Energy> [GeV]", 100, -8, 8, 0., 1000.);
      m_hVector_ECaloVsEtaTProf[j].push_back(histo_ECaloVsEtaTProf);

      sprintf(name,"EnergyVsEtaBin1D_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_EnergyVsEtaBin1D = new TH1D(name,"Calorimeter Energy X #eta; #eta; #sum Energy_{calotower} [GeV]", 500, -8, 8);
      m_hVector_EnergyVsEtaBin1D[j].push_back(histo_EnergyVsEtaBin1D);

      sprintf(name,"sumECastorMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumECastorMinus = new TH1D(name,"Castor Sum of Energy; Energy [GeV]; N events",3000,0,1500); // 6000,0,3000
      m_hVector_sumECastorMinus[j].push_back(histo_sumECastorMinus);

      sprintf(name,"ECastorSector_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_ECastorSector = new TH2D(name,"Castor Energy X Sector; Sector; Energy [GeV]", 17, 0, 17, 44, 0., 220.);
      m_hVector_ECastorSector[j].push_back(histo_ECastorSector);

      sprintf(name,"ECastorSectorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_ECastorSectorTProf = new TProfile(name,"Castor Energy X Sector; Sector; <Energy> [GeV]", 100, 0, 17, 0., 7000.);
      m_hVector_ECastorSectorTProf[j].push_back(histo_ECastorSectorTProf); 

      sprintf(name,"ECastorSectorBin1D_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ECastorSectorBin1D = new TH1D(name,"Castor Energy X Sector; Sector; Energy [GeV]", 17, 0, 17);
      m_hVector_ECastorSectorBin1D[j].push_back(histo_ECastorSectorBin1D);

      sprintf(name,"EnergyHFPlusVsCastorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_EnergyHFPlusVsCastorTProf = new TProfile(name,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 1000. );
      m_hVector_EnergyHFPlusVsCastorTProf[j].push_back(histo_EnergyHFPlusVsCastorTProf);

      sprintf(name,"EnergyHFMinusVsCastorTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_EnergyHFMinusVsCastorTProf = new TProfile(name,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 1000. );
      m_hVector_EnergyHFMinusVsCastorTProf[j].push_back(histo_EnergyHFMinusVsCastorTProf);

      sprintf(name,"sumECastorAndSumHFMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumECastorAndHFMinus = new TH1D(name,"HF^{-} and Castor Sum of Energy; Energy [GeV]; N events",2000,0,1000);
      m_hVector_sumECastorAndHFMinus[j].push_back(histo_sumECastorAndHFMinus);

      sprintf(name,"CastorMultiplicity_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_CastorMultiplicity = new TH1D(name,"Castor: number of active cells; # Active Cells; N events",81,0,81);
      m_hVector_CastorMultiplicity[j].push_back(histo_CastorMultiplicity);

      sprintf(name,"CastorMultiplicityVsLumi_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_CastorMultiplicityVsLumi = new TH2D(name,"CastorMultiplicity Vs Luminosity; Luminosity per Bunch [#mub^{-1}s^{-1}]; Castor Multiplicity",5000,0,2,81,0,81);
      m_hVector_CastorMultiplicityVsLumi[j].push_back(histo_CastorMultiplicityVsLumi);

      sprintf(name,"SectorVsTotalCastorEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_SectorVsTotalCastorEnergy = new TH2D(name,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",17,0,17,1500,0,1500);
      m_hVector_SectorVsTotalCastorEnergy[j].push_back(histo_SectorVsTotalCastorEnergy);

      sprintf(name,"SectorVsTotalCastorEnergyTProf_%s_%s",tag,Folders.at(j).c_str());
      TProfile *histo_SectorVsTotalCastorEnergyTProf = new TProfile(name,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",17,0,17,0,1500);
      m_hVector_SectorVsTotalCastorEnergyTProf[j].push_back(histo_SectorVsTotalCastorEnergyTProf);

      sprintf(name,"sumEHFplusPF_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEHFplusPF = new TH1D(name,"PF HF^{+} - Sum of Energy; #sum E_{HF^{+}} [GeV]; N events",1000,0,1000);
      m_hVector_sumEHFplusPF[j].push_back(histo_sumEHFplusPF);

      sprintf(name,"sumEHFminusPF_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_sumEHFminusPF = new TH1D(name,"PF HF^{-} - Sum of Energy; #sum E_{HF^{-}} [GeV]; N events",1000,0,1000);
      m_hVector_sumEHFminusPF[j].push_back(histo_sumEHFminusPF);

      sprintf(name,"EPFVsTowerMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_EPFVsTowerMinus = new TH2D(name,"PF Energy Vs Tower Energy; #sum E_{PF^{-}} [GeV]; #sum E_{HF^{-}} [GeV]",1000,0,1000,1000,0,1000);
      m_hVector_EPFVsTowerMinus[j].push_back(histo_EPFVsTowerMinus);

      sprintf(name,"EPFVsTowerPlus_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_EPFVsTowerPlus = new TH2D(name,"PF Energy Vs Tower Energy; #sum E_{PF^{-}} [GeV]; #sum E_{HF^{-}} [GeV]",1000,0,1000,1000,0,1000);
      m_hVector_EPFVsTowerPlus[j].push_back(histo_EPFVsTowerPlus);


      // Event Information

      sprintf(name,"lumi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_lumi = new TH1D(name,"Luminosity per Bunch; L_{Bunch} [#mub^{-1}s^{-1}]; N events",25,0,2);
      m_hVector_lumi[j].push_back(histo_lumi); 

      sprintf(name,"TracksMultiplicity_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_Tracks = new TH1D(name,"Tracks Multiplicity; n Tracks; N events",150,0,150);
      m_hVector_tracks[j].push_back(histo_Tracks);

      sprintf(name,"vertex_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_vertex = new TH1D(name,"Number of Vertex; # Vertex; N events",25,0,25);
      m_hVector_vertex[j].push_back(histo_vertex);

      sprintf(name,"pileup_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_pu = new TH1D(name,"Pile Up; # pile up; N events",25,0,25);
      m_hVector_pu[j].push_back(histo_pu);


      // Diffraction
      sprintf(name,"aEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_aSumE = new TH1D(name,"Forward Backward Asymmetry Distribution ; (#sum HF^{+} - #sum HF^{-})x(#sum HF^{+} + #sum HF^{-})^{-1}; N events",100,-2,2);
      m_hVector_asumE[j].push_back(histo_aSumE);

      sprintf(name,"etaminimum_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_EtaMinimum = new TH1D(name,"; min(#eta_{max},-#eta_{min}); N events",20,-5,5);
      m_hVector_etaminimum[j].push_back(histo_EtaMinimum);

      sprintf(name,"etamax_%s_%s",tag,Folders.at(j).c_str());
      //TH1D *histo_Etamax = new TH1D(name,"#eta_{max} Distribution; #eta; N events",18,binarrayplus);
      TH1D *histo_Etamax = new TH1D(name,"#eta_{max} Distribution; #eta; N events",20,-5.,5.);
      m_hVector_etamax[j].push_back(histo_Etamax);

      sprintf(name,"etamin_%s_%s",tag,Folders.at(j).c_str());
      //TH1D *histo_Etamin = new TH1D(name,"#eta_{min} Distribution; #eta; N events",18,binarrayminus);
      TH1D *histo_Etamin = new TH1D(name,"#eta_{min} Distribution; #eta; N events",20,-5.,5.);
      m_hVector_etamin[j].push_back(histo_Etamin);

      sprintf(name,"maxetagap_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_maxetagap = new TH1D(name,"#Delta#eta_{max} Distribution; #Delta#eta_{max}; N events",40,0.,6.);
      m_hVector_maxetagap[j].push_back(histo_maxetagap);

      sprintf(name,"SumPTMaxgap_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumPTMax = new TH1D(name,"P_{T}; #sum pT_{max} [GeV]; N events",1200,0,600);
      m_hVector_SumPTMax[j].push_back(histo_SumPTMax);

      sprintf(name,"SumPTMingap_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_SumPTMin = new TH1D(name,"P_{T}; #sum pT_{min} [GeV]; N events",1200,0,600);
      m_hVector_SumPTMin[j].push_back(histo_SumPTMin);

      sprintf(name,"deltaEtamaxmin_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_deltaEta = new TH1D(name,"#Delta#eta Distribution; #eta_{max}-#eta_{min}; N events",14,binarraydelta);
      m_hVector_deltaEta[j].push_back(histo_deltaEta);

      sprintf(name,"absdeltaEtamaxmin_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_absdeltaEta = new TH1D(name,"|#Delta#eta| Distribution; |#eta_{max}-#eta_{min}|; N events",14,binarraydelta);
      m_hVector_absdeltaEta[j].push_back(histo_absdeltaEta);

      sprintf(name,"xiPlus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiPlus = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_XiPlus[j].push_back(histo_XiPlus);

      sprintf(name,"xiMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiMinus = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_XiMinus[j].push_back(histo_XiMinus);

      sprintf(name,"xiPlusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiPlusCastor = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_XiPlusCastor[j].push_back(histo_XiPlusCastor);

      sprintf(name,"xiMinusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiMinusCastor = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_XiMinusCastor[j].push_back(histo_XiMinusCastor);

      sprintf(name,"xiPlusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiPlusAll = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_XiPlusAll[j].push_back(histo_XiPlusAll);

      sprintf(name,"xiMinusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiMinusAll = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_XiMinusAll[j].push_back(histo_XiMinusAll);

      sprintf(name,"xiMinusUnique_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_XiMinusUnique = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_XiMinusUnique[j].push_back(histo_XiMinusUnique);

      sprintf(name,"xi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_Xi = new TH1D(name,"#xi; #xi; N Event",nbinxi,xi_bin);
      m_hVector_Xi[j].push_back(histo_Xi);

      sprintf(name,"AEcastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_AEcastor = new TH1D(name,"A_{Energy}; A_{Energy}; N events",1200,-1.5,1.5);
      m_hVector_AEcastor[j].push_back(histo_AEcastor);

      sprintf(name,"etasignedHF_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_etasignedHF = new TH1D(name,"Signed #eta_{l} Distribution; #tilde{#eta_{l}}; N events",100,-10.,10.);
      m_hVector_etasignedHF[j].push_back(histo_etasignedHF);

      sprintf(name,"etasignedCASTOR_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_etasignedCASTOR = new TH1D(name,"Signed #eta_{l} Distribution; #tilde{#eta_{l}}; N events",100,-10.,10.);
      m_hVector_etasignedCASTOR[j].push_back(histo_etasignedCASTOR);


      // Generator

      sprintf(name,"genProtonMinusXi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genProtonMinusXi = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_genProtonMinusXi[j].push_back(histo_genProtonMinusXi);

      sprintf(name,"genProtonPlusXi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genProtonPlusXi = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_genProtonPlusXi[j].push_back(histo_genProtonPlusXi);

      sprintf(name,"genXiMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiMinus = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_genXiMinus[j].push_back(histo_genXiMinus);

      sprintf(name,"genXiMinusMore_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiMinusMore = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_genXiMinusMore[j].push_back(histo_genXiMinusMore);

      sprintf(name,"genXiPlus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiPlus = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_genXiPlus[j].push_back(histo_genXiPlus);

      sprintf(name,"genXiMinusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiMinusCastor = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_genXiMinusCastor[j].push_back(histo_genXiMinusCastor);

      sprintf(name,"genXiPlusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiPlusCastor = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_genXiPlusCastor[j].push_back(histo_genXiPlusCastor);

      sprintf(name,"genXiMinusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiMinusAll = new TH1D(name,"#xi_{-}; #xi_{-}; N Event",nbinxi,xi_bin);
      m_hVector_genXiMinusAll[j].push_back(histo_genXiMinusAll);

      sprintf(name,"genXiPlusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXiPlusAll = new TH1D(name,"#xi_{+}; #xi_{+}; N Event",nbinxi,xi_bin);
      m_hVector_genXiPlusAll[j].push_back(histo_genXiPlusAll);

      sprintf(name,"genXi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_genXi = new TH1D(name,"#xi; #frac{M_{X}^{2}}{s}; N Event",nbinxi,xi_bin);
      m_hVector_genXi[j].push_back(histo_genXi);

      sprintf(name,"resLeadingLeptonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resLeadingLeptonPt = new TH1D(name,"; #frac{P_{T}^{reco}}{P_{T}^{gen}}; N Event",100,0.6,1.4);
      m_hVector_resLeadingLeptonPt[j].push_back(histo_resLeadingLeptonPt);

      sprintf(name,"resLeadingLeptonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resLeadingLeptonEta = new TH1D(name,"; #frac{#eta^{reco}}{#eta^{gen}}; N Event",100,0.996,1.005);
      m_hVector_resLeadingLeptonEta[j].push_back(histo_resLeadingLeptonEta);

      sprintf(name,"resLeadingLeptonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resLeadingLeptonPhi = new TH1D(name,"; #frac{#phi^{reco}}{#phi^{gen}}; N Event",100,0.996,1.005);
      m_hVector_resLeadingLeptonPhi[j].push_back(histo_resLeadingLeptonPhi);

      sprintf(name,"resMETPt_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resMETPt = new TH1D(name,"; #frac{P_{T}^{reco}}{P_{T}^{gen}}; N Event",200,-2,2);
      m_hVector_resMETPt[j].push_back(histo_resMETPt);

      sprintf(name,"resMETPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resMETPhi = new TH1D(name,"; #frac{#phi^{reco}}{#phi^{gen}}; N Event",200,-2,2);
      m_hVector_resMETPhi[j].push_back(histo_resMETPhi);

      sprintf(name,"resXiPlus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiPlus = new TH1D(name,"; #frac{(#xi^{gen}_{+}-#xi^{reco}_{+})}{#xi^{gen}_{+}}; N Event",300,-6,6);
      m_hVector_resXiPlus[j].push_back(histo_resXiPlus);

      sprintf(name,"resXiMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiMinus = new TH1D(name,"; #frac{(#xi^{gen}_{-}-#xi^{reco}_{-})}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_resXiMinus[j].push_back(histo_resXiMinus);

      sprintf(name,"resXiMinusMore_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiMinusMore = new TH1D(name,"; #frac{(#xi^{gen}_{-}-#xi^{reco}_{-})}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_resXiMinusMore[j].push_back(histo_resXiMinusMore);

      sprintf(name,"resXiPlusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiPlusCastor = new TH1D(name,"; #frac{(#xi^{gen}_{+}-#xi^{reco}_{+})}{#xi^{gen}_{+}}; N Event",300,-6,6);
      m_hVector_resXiPlusCastor[j].push_back(histo_resXiPlusCastor);

      sprintf(name,"resXiMinusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiMinusCastor = new TH1D(name,"; #frac{(#xi^{gen}_{-}-#xi^{reco}_{-})}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_resXiMinusCastor[j].push_back(histo_resXiMinusCastor);

      sprintf(name,"resXiPlusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiPlusAll = new TH1D(name,"; #frac{(#xi^{gen}_{+}-#xi^{reco}_{+})}{#xi^{gen}_{+}}; N Event",300,-6,6);
      m_hVector_resXiPlusAll[j].push_back(histo_resXiPlusAll);

      sprintf(name,"resXiMinusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resXiMinusAll = new TH1D(name,"; #frac{(#xi^{gen}_{-}-#xi^{reco}_{-})}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_resXiMinusAll[j].push_back(histo_resXiMinusAll);

      sprintf(name,"resAbsXiPlus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resAbsXiPlus = new TH1D(name,"; #xi^{gen}_{+}-#xi^{reco}_{+}; N Event",300,-6,6);
      m_hVector_resAbsXiPlus[j].push_back(histo_resAbsXiPlus);

      sprintf(name,"resAbsXiMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resAbsXiMinus = new TH1D(name,"; #xi^{gen}_{-}-#xi^{reco}_{-}; N Event",300,-6,6);
      m_hVector_resAbsXiMinus[j].push_back(histo_resAbsXiMinus);

      sprintf(name,"resAbsXiPlusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resAbsXiPlusCastor = new TH1D(name,"; #xi^{gen}_{+}-#xi^{reco}_{+}; N Event",300,-6,6);
      m_hVector_resAbsXiPlusCastor[j].push_back(histo_resAbsXiPlusCastor);

      sprintf(name,"resAbsXiMinusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resAbsXiMinusCastor = new TH1D(name,"; #xi^{gen}_{-}-#xi^{reco}_{-}; N Event",300,-6,6);
      m_hVector_resAbsXiMinusCastor[j].push_back(histo_resAbsXiMinusCastor);

      sprintf(name,"resAbsXiPlusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resAbsXiPlusAll = new TH1D(name,"; #xi^{gen}_{+}-#xi^{reco}_{+}; N Event",300,-6,6);
      m_hVector_resAbsXiPlusAll[j].push_back(histo_resAbsXiPlusAll);

      sprintf(name,"resAbsXiMinusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resAbsXiMinusAll = new TH1D(name,"; #xi^{gen}_{-}-#xi^{reco}_{-}; N Event",300,-6,6);
      m_hVector_resAbsXiMinusAll[j].push_back(histo_resAbsXiMinusAll);

      sprintf(name,"ratioXiPlus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiPlus = new TH1D(name,"; #frac{#xi^{reco}_{+}}{#xi^{gen}_{+}}; N Event",300,-6,6);
      m_hVector_ratioXiPlus[j].push_back(histo_ratioXiPlus);

      sprintf(name,"ratioXiMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiMinus = new TH1D(name,"; #frac{#xi^{reco}_{-}}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_ratioXiMinus[j].push_back(histo_ratioXiMinus);

      sprintf(name,"ratioXiMinusMore_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiMinusMore = new TH1D(name,"; #frac{#xi^{reco}_{-}}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_ratioXiMinusMore[j].push_back(histo_ratioXiMinusMore);

      sprintf(name,"ratioXiMinusUnique_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiMinusUnique = new TH1D(name,"; #frac{#xi^{reco}_{-}}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_ratioXiMinusUnique[j].push_back(histo_ratioXiMinusUnique);

      sprintf(name,"ratioXiPlusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiPlusCastor = new TH1D(name,"; #frac{#xi^{reco}_{+}}{#xi^{gen}_{+}}; N Event",300,-6,6);
      m_hVector_ratioXiPlusCastor[j].push_back(histo_ratioXiPlusCastor);

      sprintf(name,"ratioXiMinusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiMinusCastor = new TH1D(name,"; #frac{#xi^{reco}_{-}}{#xi^{gen}_{-}}; N Event",1000,-6.,6.);
      m_hVector_ratioXiMinusCastor[j].push_back(histo_ratioXiMinusCastor);

      sprintf(name,"ratioXiPlusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiPlusAll = new TH1D(name,"; #frac{#xi^{reco}_{+}}{#xi^{gen}_{+}}; N Event",300,-6,6);
      m_hVector_ratioXiPlusAll[j].push_back(histo_ratioXiPlusAll);

      sprintf(name,"ratioXiMinusAll_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioXiMinusAll = new TH1D(name,"; #frac{#xi^{reco}_{-}}{#xi^{gen}_{-}}; N Event",300,-6,6);
      m_hVector_ratioXiMinusAll[j].push_back(histo_ratioXiMinusAll);

      sprintf(name,"resHFEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resHFEnergy = new TH1D(name,"; #frac{#sum E_{HF}^{reco}}{#sum E_{HF}^{gen}}; N Event",200,0,20);
      m_hVector_resHFEnergy[j].push_back(histo_resHFEnergy);

      sprintf(name,"resCASTOREnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_resCASTOREnergy = new TH1D(name,"; #frac{#sum E_{CASTOR}^{gen}-#sum E_{CASTOR}^{reco}}{#sum E_{CASTOR}^{gen}}; N Event",200,-10,10);
      m_hVector_resCASTOREnergy[j].push_back(histo_resCASTOREnergy);

      sprintf(name,"ratioCASTOREnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioCASTOREnergy = new TH1D(name,"; #frac{#sum E_{CASTOR}^{reco}}{#sum E_{CASTOR}^{gen}}; N Event",200,-10,10);
      m_hVector_ratioCASTOREnergy[j].push_back(histo_ratioCASTOREnergy);

      sprintf(name,"ratioCASTORPionsEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioCASTORPionsEnergy = new TH1D(name,"; #frac{#sum E_{CASTOR}^{reco}}{#sum E_{CASTOR}^{gen}}; N Event",200,-10,10);
      m_hVector_ratioCASTORPionsEnergy[j].push_back(histo_ratioCASTORPionsEnergy);

      sprintf(name,"ratioCASTOREPhotonEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioCASTOREPhotonEnergy = new TH1D(name,"; #frac{#sum E_{CASTOR}^{reco}}{#sum E_{CASTOR}^{gen}}; N Event",200,-10,10);
      m_hVector_ratioCASTOREPhotonEnergy[j].push_back(histo_ratioCASTOREPhotonEnergy);

      sprintf(name,"ratioCASTOROthersEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ratioCASTOROthersEnergy = new TH1D(name,"; #frac{#sum E_{CASTOR}^{reco}}{#sum E_{CASTOR}^{gen}}; N Event",200,-10,10);
      m_hVector_ratioCASTOROthersEnergy[j].push_back(histo_ratioCASTOROthersEnergy);

      sprintf(name,"correlHFEnergy_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlHFEnergy = new TH2D(name,"; #sum E_{HF, gen} [GeV]; #sum E_{HF} [GeV]",1000,0,1000,1000,0,1000);
      m_hVector_correlHFEnergy[j].push_back(histo_correlHFEnergy);

      sprintf(name,"correlCASTOREnergy_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlCASTOREnergy = new TH2D(name,"; #sum E_{CASTOR, gen} [GeV]; #sum E_{CASTOR} [GeV]",500,0,1500,1000,0,1000);
      m_hVector_correlCASTOREnergy[j].push_back(histo_correlCASTOREnergy);

      sprintf(name,"correlXiPlus_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiPlus = new TH2D(name,"; log_{10}(#xi^{+}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiPlus[j].push_back(histo_correlXiPlus);

      sprintf(name,"correlXiMinus_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinus = new TH2D(name,"; log_{10}(#xi^{-}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinus[j].push_back(histo_correlXiMinus);

      sprintf(name,"correlXiMinusMore_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusMore = new TH2D(name,"; log_{10}(#xi^{-}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusMore[j].push_back(histo_correlXiMinusMore);

      sprintf(name,"correlXiPlusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiPlusCastor = new TH2D(name,"; log_{10}(#xi^{+}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiPlusCastor[j].push_back(histo_correlXiPlusCastor);

      sprintf(name,"correlXiMinusCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusCastor = new TH2D(name,"; log_{10}(#xi^{-}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusCastor[j].push_back(histo_correlXiMinusCastor);

      sprintf(name,"correlXiMinusCastorDef_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusCastorDef = new TH2D(name,"; log_{10}{#sum (E - pz)}/#sqrt{s}, gen at CASTOR; log_{10}{E(1-cos#theta)}/#sqrt{s}, gen at CASTOR",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusCastorDef[j].push_back(histo_correlXiMinusCastorDef);

      sprintf(name,"correlXiPlusAll_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiPlusAll = new TH2D(name,"; log_{10}(#xi^{+}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiPlusAll[j].push_back(histo_correlXiPlusAll);

      sprintf(name,"correlXiMinusAll_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusAll = new TH2D(name,"; log_{10}(#xi^{-}_{gen}); log_{10}(#xi_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusAll[j].push_back(histo_correlXiMinusAll);

      sprintf(name,"correlXiMinusRecoGen_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusRecoGen = new TH2D(name,"; log_{10}(#xi^{-}_{reco}); log_{10}(#xi^{-}_{gen}",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusRecoGen[j].push_back(histo_correlXiMinusRecoGen);

      sprintf(name,"correlXiMinusCastorRecoGen_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusCastorRecoGen = new TH2D(name,"; log_{10}(#xi^{-}_{reco}); log_{10}(#xi^{-}_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusCastorRecoGen[j].push_back(histo_correlXiMinusCastorRecoGen);

      sprintf(name,"correlXiMinusAllRecoGen_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusAllRecoGen = new TH2D(name,"; log_{10}(#xi^{-}_{reco}); log_{10}(#xi^{-}_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusAllRecoGen[j].push_back(histo_correlXiMinusAllRecoGen);

      sprintf(name,"correlXiMinusUniqueRecoGen_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlXiMinusUniqueRecoGen = new TH2D(name,"; log_{10}(#xi^{-}_{reco}); log_{10}(#xi^{-}_{gen})",200,-4,0,200,-4,0);
      m_hVector_correlXiMinusUniqueRecoGen[j].push_back(histo_correlXiMinusUniqueRecoGen);

      sprintf(name,"correlRatioCastor_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlRatioCastor = new TH2D(name,"; #frac{#sum E_{CASTOR}^{gen}}{#sum E_{CASTOR}^{reco}}; #sum E_{CASTOR}^{reco}",500,0,20,500,0,2000);
      m_hVector_correlRatioCastor[j].push_back(histo_correlRatioCastor);

      sprintf(name,"correlRatioCastorPion_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlRatioCastorPion = new TH2D(name,"; #frac{#sum E_{Castor}^{reco}}{#sum E_{Castor}^{gen}}; #sum E_{Castor}^{gen}",500,0,20,500,0,2000);
      m_hVector_correlRatioCastorPion[j].push_back(histo_correlRatioCastorPion);

      sprintf(name,"correlRatioCastorEPhoton_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlRatioCastorEPhoton = new TH2D(name,"; #frac{#sum E_{Castor}^{reco}}{#sum E_{Castor}^{gen}}; #sum E_{Castor}^{gen}",500,0,20,500,0,2000);
      m_hVector_correlRatioCastorEPhoton[j].push_back(histo_correlRatioCastorEPhoton);

      sprintf(name,"correlRatioCastorOthers_%s_%s",tag,Folders.at(j).c_str());
      TH2D *histo_correlRatioCastorOthers = new TH2D(name,"; #frac{#sum E_{Castor}^{reco}}{#sum E_{Castor}^{gen}}; #sum E_{Castor}^{gen}",500,0,20,500,0,2000);
      m_hVector_correlRatioCastorOthers[j].push_back(histo_correlRatioCastorOthers);

      sprintf(name,"gensumEHFplus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_gensumEHFplus = new TH1D(name,"HF^{+} - Sum of Energy; #sum E_{HF^{+}, gen} [GeV]; N events",2000,0,2000);
      m_hVector_gensumEHFplus[j].push_back(histo_gensumEHFplus);

      sprintf(name,"gensumEHFminus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_gensumEHFminus = new TH1D(name,"HF^{-} - Sum of Energy; #sum E_{HF^{-}, gen} [GeV]; N events",2000,0,2000);
      m_hVector_gensumEHFminus[j].push_back(histo_gensumEHFminus);

      sprintf(name,"gensumECastorMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_gensumECastorMinus = new TH1D(name,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500); // 6000,0,3000
      m_hVector_gensumECastorMinus[j].push_back(histo_gensumECastorMinus);

      sprintf(name,"gensumECastorPionsMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_gensumECastorPionsMinus = new TH1D(name,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500);
      m_hVector_gensumECastorPionsMinus[j].push_back(histo_gensumECastorPionsMinus);

      sprintf(name,"gensumECastorEPhotonMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_gensumECastorEPhotonMinus = new TH1D(name,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500);
      m_hVector_gensumECastorEPhotonMinus[j].push_back(histo_gensumECastorEPhotonMinus);

      sprintf(name,"gensumECastorOthersMinus_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_gensumECastorOthersMinus = new TH1D(name,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500);
      m_hVector_gensumECastorOthersMinus[j].push_back(histo_gensumECastorOthersMinus);

      sprintf(name,"genBosonWEta_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_WEtaGen = new TH1D(name, "; #eta_{gen}; N events",50,-5.2,5.2);
      m_hVector_WEtaGen[j].push_back(histo_WEtaGen);

      sprintf(name,"NumberOfProtons_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_NumberOfProtons = new TH1D(name,"; Number of Protons per event; N Event",20,0,20);
      m_hVector_NumberOfProtons[j].push_back(histo_NumberOfProtons);

      sprintf(name,"ProtonsPz_%s_%s",tag,Folders.at(j).c_str());
      TH1D *histo_ProtonPz = new TH1D(name,"; Leading Proton |Pz|; N Event",500,0,7000);
      m_hVector_ProtonPz[j].push_back(histo_ProtonPz);

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
  m_hVector_LeadingLeptonPt[index].at(pileup)->Fill(leptonpt,totalweight);
  m_hVector_LeadingLeptonPhi[index].at(pileup)->Fill(leptonphi,totalweight);
  m_hVector_LeadingLeptonEta[index].at(pileup)->Fill(leptoneta,totalweight);
  m_hVector_DeltaPhi[index].at(pileup)->Fill(fabs(deltaphi),totalweight);
  m_hVector_LeadingLeptonDxy[index].at(pileup)->Fill(fabs(eventdiffW->GetLeadingMuonDxy()),totalweight);
  m_hVector_LeadingLeptonDz[index].at(pileup)->Fill(fabs(eventdiffW->GetLeadingMuonDz()),totalweight);
  m_hVector_SecondLeptonPt[index].at(pileup)->Fill(secondleptonpt,totalweight);
  m_hVector_SecondLeptonPhi[index].at(pileup)->Fill(secondleptonphi,totalweight);
  m_hVector_SecondLeptonEta[index].at(pileup)->Fill(secondleptoneta,totalweight);
  m_hVector_METPt[index].at(pileup)->Fill(eventdiffW->GetMETPt(),totalweight);
  m_hVector_METPhi[index].at(pileup)->Fill(eventdiffW->GetMETPhi(),totalweight);
  m_hVector_METSignificance[index].at(pileup)->Fill(eventdiffW->GetMETsigma(),totalweight);
  m_hVector_DVtxMuon[index].at(pileup)->Fill(eventdiffW->GetDVtxMuon(),totalweight);
  m_hVector_DVtxMuonZ[index].at(pileup)->Fill(eventdiffW->GetDVtxMuonZ(),totalweight);
  m_hVector_DVtxElectron[index].at(pileup)->Fill(eventdiffW->GetDVtxElectron(),totalweight);
  m_hVector_DVtxElectronZ[index].at(pileup)->Fill(eventdiffW->GetDVtxElectronZ(),totalweight);
  m_hVector_DMuonElectron[index].at(pileup)->Fill(eventdiffW->GetDMuonElectron(),totalweight);
  m_hVector_DMuonElectronZ[index].at(pileup)->Fill(eventdiffW->GetDMuonElectronZ(),totalweight);
  m_hVector_DMuons[index].at(pileup)->Fill(eventdiffW->GetDMuons(),totalweight);
  m_hVector_DElectrons[index].at(pileup)->Fill(eventdiffW->GetDElectrons(),totalweight);

  // Detector
  m_hVector_multhf[index].at(pileup)->Fill(eventdiffW->GetMultiplicityHFPlus(),eventdiffW->GetMultiplicityHFMinus(),totalweight);
  m_hVector_minusnhf[index].at(pileup)->Fill(eventdiffW->GetMultiplicityHFMinus(),totalweight);
  m_hVector_plusnhf[index].at(pileup)->Fill(eventdiffW->GetMultiplicityHFPlus(),totalweight);
  m_hVector_sumEHFplus[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus(),totalweight);
  m_hVector_sumEHFminus[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinus(),totalweight);
  m_hVector_sumEHEplus[index].at(pileup)->Fill(eventdiffW->GetSumEHEPlus(),totalweight);
  m_hVector_sumEHEminus[index].at(pileup)->Fill(eventdiffW->GetSumEHEMinus(),totalweight);
  m_hVector_SumEHFplus_S[index].at(pileup)->Fill(eventdiffW->GetSumEHF_SPlus(),totalweight);
  m_hVector_SumEHFminus_S[index].at(pileup)->Fill(eventdiffW->GetSumEHF_SMinus(),totalweight);
  m_hVector_SumEHFplus_L[index].at(pileup)->Fill(eventdiffW->GetSumEHF_LPlus(),totalweight);
  m_hVector_SumEHFminus_L[index].at(pileup)->Fill(eventdiffW->GetSumEHF_LMinus(),totalweight);
  m_hVector_sumEEEminus[index].at(pileup)->Fill(eventdiffW->GetSumEEEMinus(),totalweight);
  m_hVector_sumEEEplus[index].at(pileup)->Fill(eventdiffW->GetSumEEEPlus(),totalweight);
  m_hVector_EnergyHFPlusVsEnergyHFMinus[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus(),eventdiffW->GetSumEHFMinus(),totalweight);
  m_hVector_EnergyEEPlusVsEnergyEEMinus[index].at(pileup)->Fill(eventdiffW->GetSumEEEPlus(),eventdiffW->GetSumEEEMinus(),totalweight);

  for (int k=0; k<eventdiffW->GetEachTowerCounter();k++){
    m_hVector_ECaloVsEta[index].at(pileup)->Fill(eventdiffW->GetEachTowerEta(k),eventdiffW->GetEachTowerEnergy(k),totalweight);
    m_hVector_ECaloVsEtaTProf[index].at(pileup)->Fill(eventdiffW->GetEachTowerEta(k),eventdiffW->GetEachTowerEnergy(k),totalweight);
    m_hVector_EnergyVsEtaBin1D[index].at(pileup)->Fill(eventdiffW->GetEachTowerEta(k),eventdiffW->GetEachTowerEnergy(k)*totalweight);
  }

  if (eventdiffW->GetSumEHFPlus() > eventdiffW->GetSumEHFMinus()){
    m_hVector_SumEHFMax[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus(),totalweight);
    m_hVector_SumEHFMin[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinus(),totalweight);
  }else{
    m_hVector_SumEHFMax[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinus(),totalweight);
    m_hVector_SumEHFMin[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus(),totalweight);
  }

  m_hVector_etcalos_p[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_etcalos_n[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinus(),sumCastorEnergy,totalweight);
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

  m_hVector_EnergyHFMinusVsCastorTProf[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_EnergyHFPlusVsCastorTProf[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_sumECastorAndHFMinus[index].at(pileup)->Fill(sumCastorAndHFMinusEnergy,totalweight);
  m_hVector_CastorMultiplicity[index].at(pileup)->Fill(counterHit,totalweight);
  m_hVector_CastorMultiplicityVsLumi[index].at(pileup)->Fill(eventinfo->GetInstLumiBunch(),counterHit,totalweight);
  m_hVector_SectorVsTotalCastorEnergy[index].at(pileup)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);
  m_hVector_SectorVsTotalCastorEnergyTProf[index].at(pileup)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);
  m_hVector_sumEHFplusPF[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlusPF(),totalweight);
  m_hVector_sumEHFminusPF[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinusPF(),totalweight);
  m_hVector_EPFVsTowerMinus[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinusPF(),eventdiffW->GetSumEHFMinus(),totalweight);
  m_hVector_EPFVsTowerPlus[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlusPF(),eventdiffW->GetSumEHFPlus(),totalweight);

  // Event Information
  m_hVector_lumi[index].at(pileup)->Fill(eventinfo->GetInstLumiBunch(),totalweight);
  m_hVector_tracks[index].at(pileup)->Fill(eventdiff->GetMultiplicityTracks(),totalweight);
  m_hVector_vertex[index].at(pileup)->Fill(eventdiff->GetNVertex(),totalweight);
  m_hVector_pu[index].at(pileup)->Fill(eventinfo->GetNPileUpBx0(),totalweight);

  // Diffractive Variables
  m_hVector_asumE[index].at(pileup)->Fill(aSumE,totalweight);
  if(fabs(etamin)>fabs(etamax)){
    m_hVector_etaminimum[index].at(pileup)->Fill(etamax,totalweight);
  }else{
    m_hVector_etaminimum[index].at(pileup)->Fill(-1*etamin,totalweight);
  }
  m_hVector_etamax[index].at(pileup)->Fill(etamax,totalweight);
  m_hVector_etamin[index].at(pileup)->Fill(etamin,totalweight);
  m_hVector_maxetagap[index].at(pileup)->Fill(fabs(maxLRG),totalweight);
  if(fabs(eventdiffW->GetSumptPFLeft())>fabs(eventdiffW->GetSumptPFRight())){
    m_hVector_SumPTMax[index].at(pileup)->Fill(eventdiffW->GetSumptPFLeft(),totalweight);
    m_hVector_SumPTMin[index].at(pileup)->Fill(eventdiffW->GetSumptPFRight(),totalweight);
  }else{
    m_hVector_SumPTMax[index].at(pileup)->Fill(eventdiffW->GetSumptPFRight(),totalweight);
    m_hVector_SumPTMin[index].at(pileup)->Fill(eventdiffW->GetSumptPFLeft(),totalweight);
  }
  m_hVector_absdeltaEta[index].at(pileup)->Fill(fabs(deltaeta),totalweight);
  m_hVector_deltaEta[index].at(pileup)->Fill(deltaeta,totalweight);
  m_hVector_XiPlus[index].at(pileup)->Fill(xiplus,totalweight);
  m_hVector_XiMinus[index].at(pileup)->Fill(ximinus,totalweight);
  m_hVector_XiPlusCastor[index].at(pileup)->Fill(xipluscastor,totalweight);
  m_hVector_XiMinusCastor[index].at(pileup)->Fill(ximinuscastor,totalweight);
  m_hVector_XiPlusAll[index].at(pileup)->Fill(xiplusall,totalweight);
  m_hVector_XiMinusAll[index].at(pileup)->Fill(ximinusall,totalweight);
  m_hVector_XiMinusUnique[index].at(pileup)->Fill(ximinusunique,totalweight);
  m_hVector_Xi[index].at(pileup)->Fill(xi,totalweight);
  m_hVector_AEcastor[index].at(pileup)->Fill(AEcastor,totalweight);
  m_hVector_etasignedHF[index].at(pileup)->Fill(etasignedHF,totalweight);
  m_hVector_etasignedCASTOR[index].at(pileup)->Fill(etasignedCASTOR,totalweight);

  // Generator
  if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
    m_hVector_genProtonMinusXi[index].at(pileup)->Fill(eventdiff->GetXiGenMinus(),totalweight);
    m_hVector_genProtonPlusXi[index].at(pileup)->Fill(eventdiff->GetXiGenPlus(),totalweight);
    m_hVector_genXiPlus[index].at(pileup)->Fill(xigenplus,totalweight);
    m_hVector_genXiMinus[index].at(pileup)->Fill(xigenminus,totalweight);
    m_hVector_genXiMinusMore[index].at(pileup)->Fill(xigenminusMore,totalweight);
    m_hVector_genXiPlusCastor[index].at(pileup)->Fill(xigenpluscastor,totalweight);
    m_hVector_genXiMinusCastor[index].at(pileup)->Fill(xigenminuscastor,totalweight);
    m_hVector_genXiPlusAll[index].at(pileup)->Fill(xigenplusall,totalweight);
    m_hVector_genXiMinusAll[index].at(pileup)->Fill(xigenminusall,totalweight);
    m_hVector_genXi[index].at(pileup)->Fill(xigen,totalweight);
    m_hVector_resMETPt[index].at(pileup)->Fill(resMETPt);
    m_hVector_resMETPhi[index].at(pileup)->Fill(resMETPhi);
    m_hVector_resLeadingLeptonPt[index].at(pileup)->Fill(resLeadingPt);
    m_hVector_resLeadingLeptonEta[index].at(pileup)->Fill(resLeadingEta);
    m_hVector_resLeadingLeptonPhi[index].at(pileup)->Fill(resLeadingPhi);
    m_hVector_resXiPlus[index].at(pileup)->Fill((xigenplus-xiplus)/xigenplus);
    m_hVector_resXiMinus[index].at(pileup)->Fill((xigenminus-ximinus)/xigenminus);
    m_hVector_resXiMinusMore[index].at(pileup)->Fill((xigenminusMore-ximinusall)/xigenminusMore);
    m_hVector_resXiPlusCastor[index].at(pileup)->Fill((xigenpluscastor-xipluscastor)/xigenpluscastor);
    m_hVector_resXiMinusCastor[index].at(pileup)->Fill((xigenminuscastor-ximinuscastor)/xigenminuscastor);
    m_hVector_resXiPlusAll[index].at(pileup)->Fill((xigenplusall-xiplusall)/xigenplusall);
    m_hVector_resXiMinusAll[index].at(pileup)->Fill((xigenminusall-ximinusall)/xigenminusall);
    m_hVector_resAbsXiPlus[index].at(pileup)->Fill(xigenplus-xiplus);
    m_hVector_resAbsXiMinus[index].at(pileup)->Fill(xigenminus-ximinus);
    m_hVector_resAbsXiPlusCastor[index].at(pileup)->Fill(xigenpluscastor-xipluscastor);
    m_hVector_resAbsXiMinusCastor[index].at(pileup)->Fill(xigenminuscastor-ximinuscastor);
    m_hVector_resAbsXiPlusAll[index].at(pileup)->Fill(xigenplusall-xiplusall);
    m_hVector_resAbsXiMinusAll[index].at(pileup)->Fill(xigenminusall-ximinusall);
    m_hVector_ratioXiPlus[index].at(pileup)->Fill(xiplus/xigenplus);
    m_hVector_ratioXiMinus[index].at(pileup)->Fill(ximinus/xigenminus);
    if(sumCastorEnergy > 0.){
      m_hVector_ratioXiMinusMore[index].at(pileup)->Fill(ximinusall/xigenminusMore);
      m_hVector_ratioXiPlusCastor[index].at(pileup)->Fill(xipluscastor/xigenpluscastor);
      m_hVector_ratioXiMinusCastor[index].at(pileup)->Fill(ximinuscastor/xigenminuscastor);
      m_hVector_ratioXiPlusAll[index].at(pileup)->Fill(xiplusall/xigenplusall);
      m_hVector_ratioXiMinusAll[index].at(pileup)->Fill(ximinusall/xigenminusall);
      m_hVector_ratioXiMinusUnique[index].at(pileup)->Fill(ximinusunique/xigenminusall);
    }
    m_hVector_resHFEnergy[index].at(pileup)->Fill(eventdiffW->GetSumEHFPlus()/eventdiff->GetSumEnergyHFPlusGen());
    m_hVector_resHFEnergy[index].at(pileup)->Fill(eventdiffW->GetSumEHFMinus()/eventdiff->GetSumEnergyHFMinusGen());
    m_hVector_correlHFEnergy[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlusGen(),eventdiffW->GetSumEHFPlus());
    m_hVector_correlHFEnergy[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinusGen(),eventdiffW->GetSumEHFMinus());
    m_hVector_correlXiPlus[index].at(pileup)->Fill(TMath::Log10(xigenplus),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiMinus[index].at(pileup)->Fill(TMath::Log10(xigenminus),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiMinusMore[index].at(pileup)->Fill(TMath::Log10(xigenminusMore),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiPlusCastor[index].at(pileup)->Fill(TMath::Log10(xigenpluscastor),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiMinusCastor[index].at(pileup)->Fill(TMath::Log10(xigenminuscastor),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiPlusAll[index].at(pileup)->Fill(TMath::Log10(xigenplusall),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiMinusAll[index].at(pileup)->Fill(TMath::Log10(xigenminusall),TMath::Log10(xigen),totalweight);
    m_hVector_correlXiMinusCastorDef[index].at(pileup)->Fill(TMath::Log10(xigenminuscastor),TMath::Log10(xigenminuscastorDef),totalweight);
    m_hVector_correlCASTOREnergy[index].at(pileup)->Fill(sumGenCastorEnergy,sumCastorEnergy);
    m_hVector_correlRatioCastor[index].at(pileup)->Fill(sumGenCastorEnergy/sumCastorEnergy,sumCastorEnergy);
    m_hVector_correlRatioCastorPion[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorPions,sumGenCastorPions);
    m_hVector_correlRatioCastorEPhoton[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorEPhoton,sumGenCastorEPhoton);
    m_hVector_correlRatioCastorOthers[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorOthers,sumGenCastorOthers);
    m_hVector_resCASTOREnergy[index].at(pileup)->Fill((sumGenCastorEnergy-sumCastorEnergy)/sumGenCastorEnergy);
    m_hVector_ratioCASTOREnergy[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorEnergy);
    m_hVector_ratioCASTORPionsEnergy[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorPions);
    m_hVector_ratioCASTOREPhotonEnergy[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorEPhoton);
    m_hVector_ratioCASTOROthersEnergy[index].at(pileup)->Fill(sumCastorEnergy/sumGenCastorOthers);
    m_hVector_gensumEHFplus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFPlusGen(),totalweight);
    m_hVector_gensumEHFminus[index].at(pileup)->Fill(eventdiff->GetSumEnergyHFMinusGen(),totalweight);
    m_hVector_gensumECastorMinus[index].at(pileup)->Fill(sumGenCastorEnergy,totalweight);
    m_hVector_gensumECastorPionsMinus[index].at(pileup)->Fill(sumGenCastorPions,totalweight);
    m_hVector_gensumECastorEPhotonMinus[index].at(pileup)->Fill(sumGenCastorEPhoton,totalweight);
    m_hVector_gensumECastorOthersMinus[index].at(pileup)->Fill(sumGenCastorOthers,totalweight);
    m_hVector_NumberOfProtons[index].at(pileup)->Fill(counterproton);
    m_hVector_ProtonPz[index].at(pileup)->Fill(protonpz);
    m_hVector_correlXiMinusRecoGen[index].at(pileup)->Fill(TMath::Log10(ximinus),TMath::Log10(xigenminus),totalweight);
    m_hVector_correlXiMinusCastorRecoGen[index].at(pileup)->Fill(TMath::Log10(ximinuscastor),TMath::Log10(xigenminuscastor),totalweight);
    m_hVector_correlXiMinusAllRecoGen[index].at(pileup)->Fill(TMath::Log10(ximinusall),TMath::Log10(xigenminusall),totalweight);
    m_hVector_correlXiMinusUniqueRecoGen[index].at(pileup)->Fill(TMath::Log10(ximinusunique),TMath::Log10(xigenminusall),totalweight);
    m_hVector_WEtaGen[index].at(pileup)->Fill(bosonWEtaGen,totalweight);
  }

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
      m_hVector_LeadingLeptonPt[j].at(i)->Write();
      m_hVector_LeadingLeptonEta[j].at(i)->Write();
      m_hVector_LeadingLeptonPhi[j].at(i)->Write();
      m_hVector_DeltaPhi[j].at(i)->Write();
      m_hVector_LeadingLeptonDxy[j].at(i)->Write();
      m_hVector_LeadingLeptonDz[j].at(i)->Write();
      m_hVector_SecondLeptonPt[j].at(i)->Write();
      m_hVector_SecondLeptonEta[j].at(i)->Write();
      m_hVector_SecondLeptonPhi[j].at(i)->Write();
      m_hVector_METPt[j].at(i)->Write();
      m_hVector_METSignificance[j].at(i)->Write();
      m_hVector_METPhi[j].at(i)->Write();
      m_hVector_DVtxMuon[j].at(i)->Write();
      m_hVector_DVtxMuonZ[j].at(i)->Write();
      m_hVector_DVtxElectron[j].at(i)->Write();
      m_hVector_DVtxElectronZ[j].at(i)->Write();
      m_hVector_DMuonElectron[j].at(i)->Write();
      m_hVector_DMuonElectronZ[j].at(i)->Write();
      m_hVector_DMuons[j].at(i)->Write();
      m_hVector_DElectrons[j].at(i)->Write();

      // Detector Folder
      foldersFile[1]->cd();
      m_hVector_multhf[j].at(i)->Write();
      m_hVector_minusnhf[j].at(i)->Write();
      m_hVector_plusnhf[j].at(i)->Write();
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
      m_hVector_sumEHFplusPF[j].at(i)->Write();
      m_hVector_sumEHFminusPF[j].at(i)->Write();
      m_hVector_EPFVsTowerMinus[j].at(i)->Write();
      m_hVector_EPFVsTowerPlus[j].at(i)->Write();

      // Event Information
      foldersFile[2]->cd();
      m_hVector_lumi[j].at(i)->Write();
      m_hVector_tracks[j].at(i)->Write();
      m_hVector_vertex[j].at(i)->Write();
      m_hVector_pu[j].at(i)->Write();

      // Diffraction
      foldersFile[3]->cd();
      m_hVector_asumE[j].at(i)->Write();
      m_hVector_etaminimum[j].at(i)->Write();
      m_hVector_etamax[j].at(i)->Write();
      m_hVector_etamin[j].at(i)->Write();
      m_hVector_maxetagap[j].at(i)->Write();
      m_hVector_SumPTMax[j].at(i)->Write();
      m_hVector_SumPTMin[j].at(i)->Write();
      m_hVector_absdeltaEta[j].at(i)->Write();
      m_hVector_deltaEta[j].at(i)->Write();
      m_hVector_XiPlus[j].at(i)->Write();
      m_hVector_XiMinus[j].at(i)->Write();                                                                                                                                
      m_hVector_XiPlusCastor[j].at(i)->Write();
      m_hVector_XiMinusCastor[j].at(i)->Write();
      m_hVector_XiPlusAll[j].at(i)->Write();
      m_hVector_XiMinusAll[j].at(i)->Write();
      m_hVector_XiMinusUnique[j].at(i)->Write();
      m_hVector_Xi[j].at(i)->Write();
      m_hVector_AEcastor[j].at(i)->Write();
      m_hVector_etasignedHF[j].at(i)->Write();
      m_hVector_etasignedCASTOR[j].at(i)->Write();

      // Generator Plots
      if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
	foldersFile[4]->cd();
	m_hVector_genProtonMinusXi[j].at(i)->Write();
	m_hVector_genProtonPlusXi[j].at(i)->Write();
	m_hVector_genXiPlus[j].at(i)->Write();
	m_hVector_genXiMinus[j].at(i)->Write();
	m_hVector_genXiMinusMore[j].at(i)->Write();
	m_hVector_genXiPlusCastor[j].at(i)->Write();
	m_hVector_genXiMinusCastor[j].at(i)->Write();
	m_hVector_genXiPlusAll[j].at(i)->Write();
	m_hVector_genXiMinusAll[j].at(i)->Write();
	m_hVector_genXi[j].at(i)->Write();
	m_hVector_resLeadingLeptonPt[j].at(i)->Write();
	m_hVector_resLeadingLeptonEta[j].at(i)->Write();
	m_hVector_resLeadingLeptonPhi[j].at(i)->Write();
	m_hVector_resMETPt[j].at(i)->Write();
	m_hVector_resMETPhi[j].at(i)->Write();
	m_hVector_resXiPlus[j].at(i)->Write();
	m_hVector_resXiMinus[j].at(i)->Write();
	m_hVector_resXiMinusMore[j].at(i)->Write();
	m_hVector_resXiPlusCastor[j].at(i)->Write();
	m_hVector_resXiMinusCastor[j].at(i)->Write();
	m_hVector_resXiPlusAll[j].at(i)->Write();
	m_hVector_resXiMinusAll[j].at(i)->Write();
	m_hVector_resAbsXiPlus[j].at(i)->Write();
	m_hVector_resAbsXiMinus[j].at(i)->Write();
	m_hVector_resAbsXiPlusCastor[j].at(i)->Write();
	m_hVector_resAbsXiMinusCastor[j].at(i)->Write();
	m_hVector_resAbsXiPlusAll[j].at(i)->Write();
	m_hVector_resAbsXiMinusAll[j].at(i)->Write();
	m_hVector_ratioXiPlus[j].at(i)->Write();
	m_hVector_ratioXiMinus[j].at(i)->Write();
	m_hVector_ratioXiMinusMore[j].at(i)->Write();
	m_hVector_ratioXiMinusUnique[j].at(i)->Write();
	m_hVector_ratioXiPlusCastor[j].at(i)->Write();
	m_hVector_ratioXiMinusCastor[j].at(i)->Write();
	m_hVector_ratioXiPlusAll[j].at(i)->Write();
	m_hVector_ratioXiMinusAll[j].at(i)->Write();
	m_hVector_resHFEnergy[j].at(i)->Write();
	m_hVector_resCASTOREnergy[j].at(i)->Write();
	m_hVector_ratioCASTOREnergy[j].at(i)->Write();
	m_hVector_ratioCASTORPionsEnergy[j].at(i)->Write();
	m_hVector_ratioCASTOREPhotonEnergy[j].at(i)->Write();
	m_hVector_ratioCASTOROthersEnergy[j].at(i)->Write();
	m_hVector_correlHFEnergy[j].at(i)->Write();
	m_hVector_correlCASTOREnergy[j].at(i)->Write();
	m_hVector_correlXiPlus[j].at(i)->Write();
	m_hVector_correlXiMinus[j].at(i)->Write();
	m_hVector_correlXiMinusMore[j].at(i)->Write();
	m_hVector_correlXiPlusCastor[j].at(i)->Write();
	m_hVector_correlXiMinusCastor[j].at(i)->Write();
	m_hVector_correlXiMinusCastorDef[j].at(i)->Write();
	m_hVector_correlXiPlusAll[j].at(i)->Write();
	m_hVector_correlXiMinusAll[j].at(i)->Write();
	m_hVector_correlRatioCastor[j].at(i)->Write();
	m_hVector_correlRatioCastorPion[j].at(i)->Write();
	m_hVector_correlRatioCastorEPhoton[j].at(i)->Write();
	m_hVector_correlRatioCastorOthers[j].at(i)->Write();
	m_hVector_gensumEHFplus[j].at(i)->Write();
	m_hVector_gensumEHFminus[j].at(i)->Write();
	m_hVector_gensumECastorMinus[j].at(i)->Write();
	m_hVector_gensumECastorPionsMinus[j].at(i)->Write();
	m_hVector_gensumECastorEPhotonMinus[j].at(i)->Write();
	m_hVector_gensumECastorOthersMinus[j].at(i)->Write();
	m_hVector_WEtaGen[j].at(i)->Write();
	m_hVector_NumberOfProtons[j].at(i)->Write();
	m_hVector_ProtonPz[j].at(i)->Write();
	m_hVector_correlXiMinusRecoGen[j].at(i)->Write();
	m_hVector_correlXiMinusCastorRecoGen[j].at(i)->Write();
	m_hVector_correlXiMinusAllRecoGen[j].at(i)->Write();
	m_hVector_correlXiMinusUniqueRecoGen[j].at(i)->Write();
      }

    }
  }

}

void DiffractiveW::Run(std::string filein_, std::string processname_, std::string savehistofile_, std::string switchtrigger_, int optTrigger_, double lepton1pt_, double lepton2pt_, int nVertex_, std::string type_, std::string switchlumiweight_, double mcweight_, std::string typesel_, double castorthreshold_, double channelsthreshold_, std::string castorcorrfile_, std::string gapseltype_, std::string pumfile_, std::string pudfile_, std::string corcastor_, std::string corcastorfile_){

  bool debug = false;
  double sqrtS = 7000.;

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
  corcastor = corcastor_;
  corcastorfile = corcastorfile_;

  std::string selStatus;
  std::string TriggerStatus;
  char selCastor[300];

  // Cut Flow
  TH1D *histo_cutflow = new TH1D("CutFlow","Cut Flow; Cuts; Events",30,0,30);

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
  trout->Branch("deltaeta",&bdeltaeta,"bdeltaeta/D");
  trout->Branch("AEcastor",&bAEcastor,"AEcastor/D");
  trout->Branch("etasignedHF",&betasignedHF,"betasignedHF/D");
  trout->Branch("etasignedCASTOR",&betasignedCASTOR,"betasignedCASTOR/D");
  trout->Branch("MaxGap",&bMaxGap,"bMaxGap/D");
  trout->Branch("SumPTMaxLrgPF",&SumPTMaxLrgPF,"SumPTMaxLrgPF/D");
  trout->Branch("SumPTMinLrgPF",&SumPTMinLrgPF,"SumPTMinLrgPF/D");
  trout->Branch("XiPlus",&bXiPlus,"bXiPlus/D");
  trout->Branch("XiMinus",&bXiMinus,"bXiMinus/D");
  trout->Branch("EtaMax",&betamax,"betamax/D");
  trout->Branch("EtaMin",&betamin,"betamin/D");

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
  h_castor_channel = (TH2D*)check2.Get("channelcorrector");

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

  // Fit Castor Correction
  TFile* in_file= new TFile(corcastorfile.c_str());
  TString fname="linearity"; // fit object histogram name

  if (in_file->IsZombie()){
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " There is no the file " << corcastorfile << " or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;
  }                                                                                                                                                                                               

  if (!in_file->FindKey(fname)){
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << " There is no directory/path " << fname << std::endl;
    std::cout << " in the file." << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    return;
  }

  TF1 *CastorReweight = (TF1*)in_file->Get(fname);

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

    double totalASum = eventdiffW->GetSumEHFPlus() + eventdiffW->GetSumEHFMinus();
    if (totalASum > 0.){
      aSumE = (eventdiffW->GetSumEHFPlus() - eventdiffW->GetSumEHFMinus())/(eventdiffW->GetSumEHFPlus() + eventdiffW->GetSumEHFMinus());
    }

    double puweight = -1.;
    double totalcommon = -1.;

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
    bool trigger = false;
    bool vertex = false;
    bool diffseln = false;
    bool diffselp = false;
    bool generator = false;
    bool gencastorlow = false;
    bool gencastorhigh = false;
    bool zeropileup = false;

    bool presel = false;
    bool preselGen = false;
    bool dimass = false;
    bool dimassGen = false;
    bool eleBarrel1 = false;
    bool eleEndCap1 = false;
    bool candSel = false;
    bool isoBarrel1 = false;
    bool isoEndCap1 = false;
    bool isolation = false;
    bool nocosmic = false;
    bool WKinN = false;
    bool WKinP = false;
    bool acceptEvt = false;
    bool xicut = false;

    bool CastorLow = false;

    if (switchtrigger == "trigger_all_electron"){

      if ( (eventdiff->GetRunNumber() >= 146698 && eventdiff->GetRunNumber() <= 148058) && eventdiffW->GetHLTPath(10) > 0) triggerE_a = true;
      if ( (eventdiff->GetRunNumber() >= 148822 && eventdiff->GetRunNumber() <= 149063) && eventdiffW->GetHLTPath(11) > 0) triggerE_b = true;
      if ( (eventdiff->GetRunNumber() >= 149181 && eventdiff->GetRunNumber() <= 149291) && eventdiffW->GetHLTPath(12) > 0) triggerE_c = true;
      if (triggerE_a || triggerE_b || triggerE_c) trigger = true;

      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Electron accepted." << std::endl;
      TriggerStatus = "trigger_all_electron";
    }
    else if (switchtrigger == "trigger_all_muon"){

      if ( (eventdiff->GetRunNumber() >= 146428 && eventdiff->GetRunNumber() <= 147116) && eventdiffW->GetHLTPath(3) > 0) triggerE_a = true;
      if ( (eventdiff->GetRunNumber() >= 147196 && eventdiff->GetRunNumber() <= 149291) && eventdiffW->GetHLTPath(8) > 0) triggerE_b = true;
      if (triggerE_a || triggerE_b) trigger = true;

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
    castorgapGen = false;

    sumGenCastorEnergyExc = 0.;
    sumGenCastorPions = 0.;
    sumGenCastorEPhoton = 0.;
    sumGenCastorOthers = 0.;
    sumGenCastorAll = 0.; 
    sumGenCastorEnergy = 0.;

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

    sumCastorAndHFMinusEnergy = sumCastorEnergy+eventdiffW->GetSumEHFMinus();

    EpluspzGen = 0.;
    EminuspzGen = 0.;
    EminuspzGenMore = 0.;
    EminuspzGenCastor = 0.;

    math::XYZTLorentzVector System(0.,0.,0.,0.);
    math::XYZTLorentzVector SystemX(0.,0.,0.,0.);
    math::XYZTLorentzVector SystemY(0.,0.,0.,0.);
    math::XYZTLorentzVector SystemYMore(0.,0.,0.,0.);
    std::vector< math::XYZTLorentzVector > protonGen;

    if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){

      for (int j=0; j<eventdiffW->GetNParticlesGen();j++){
	if(fabs(eventdiffW->GetParticlesPDGidGen(j))==2212){
	  protonGen.push_back(eventdiffW->GetParticlesP4Gen(j));
	}
      }

      std::sort(protonGen.begin(), protonGen.end(), orderAbsolutPZ());

      if(protonGen.size()>0){
	if(debug){
	  for (unsigned int k=0; k<protonGen.size();k++){
	    std::cout << "Proton: " << protonGen.at(k).pz() <<std::endl;
	  }
	}
	protonpz = fabs(protonGen.at(0).pz());
      }

      for (int k=0; k<eventdiffW->GetNParticlesGen();k++){

	/*
	   if (eventdiffW->GetParticlesEtaGen(k)<-5.2 && eventdiffW->GetParticlesEtaGen(k)>-6.6){
	   sumGenCastorEnergy += eventdiffW->GetParticlesEnergyGen(k);
	   EminuspzGenCastor += eventdiffW->GetParticlesP4Gen(k).energy()-eventdiffW->GetParticlesP4Gen(k).pz(); 
	   }	
	 */

	int part_id = fabs(eventdiffW->GetParticlesPDGidGen(k));

	if (eventdiffW->GetParticlesEtaGen(k)<-5.2 && eventdiffW->GetParticlesEtaGen(k)>-6.6){
	  if(part_id != 12 || part_id != 14 || part_id != 16 || part_id != 13){
	    sumGenCastorEnergyExc += eventdiffW->GetParticlesEnergyGen(k);
	    sumGenCastorEnergy += eventdiffW->GetParticlesEnergyGen(k);
	    EminuspzGenCastor += eventdiffW->GetParticlesP4Gen(k).energy()-eventdiffW->GetParticlesP4Gen(k).pz(); 
	  }
	  if(part_id == 111 || part_id == 211){
	    sumGenCastorPions += eventdiffW->GetParticlesEnergyGen(k);
	  }else if(part_id == 11 || part_id == 22){
	    sumGenCastorEPhoton += eventdiffW->GetParticlesEnergyGen(k);
	  }else{
	    sumGenCastorOthers += eventdiffW->GetParticlesEnergyGen(k);
	  }
	}

	//if( ( fabs(eventdiffW->GetParticlesEtaGen(k)) <= 2.866 && fabs(eventdiffW->GetParticlesEtaGen(k)) > 3.152)){
	System += eventdiffW->GetParticlesP4Gen(k);
	//}
	//if(eventdiffW->GetParticlesEtaGen(k)<4.730 || (eventdiffW->GetParticlesEtaGen(k) <= 2.866 && eventdiffW->GetParticlesEtaGen(k) > 3.152)){
	if(eventdiffW->GetParticlesEtaGen(k)<4.9){ //Proton is already excluded here. Also dissociative system.                                                           
	  EpluspzGen += eventdiffW->GetParticlesP4Gen(k).energy()+eventdiffW->GetParticlesP4Gen(k).pz();
	  SystemX += eventdiffW->GetParticlesP4Gen(k);
	}
	//if(eventdiffW->GetParticlesEtaGen(k)>-6.6 || (eventdiffW->GetParticlesEtaGen(k) >= -2.866 && eventdiffW->GetParticlesEtaGen(k) < -3.152)){
	if(eventdiffW->GetParticlesEtaGen(k)>-4.9){ // Proton is already excluded here. Also dissociative system.
	  EminuspzGen += eventdiffW->GetParticlesP4Gen(k).energy()-eventdiffW->GetParticlesP4Gen(k).pz();
	  SystemY += eventdiffW->GetParticlesP4Gen(k);
	}

	if(eventdiffW->GetParticlesEtaGen(k)>-6.6){ // Proton is already excluded here. Also dissociative system.
	  EminuspzGenMore += eventdiffW->GetParticlesP4Gen(k).energy()-eventdiffW->GetParticlesP4Gen(k).pz();
	  SystemYMore += eventdiffW->GetParticlesP4Gen(k);
	}

      }

      // Excluding leading proton only for diffractive MC, Pomwig.
      size_t pos;
      pos = savehistofile.find("Pomwig");
      if (pos != std::string::npos){
	if(protonGen.size()>0){                                                                                                                    
	  System=System-protonGen.at(0);
	  if(protonGen.at(0).eta() > -6.6 && protonGen.at(0).eta() < 0.){ //proton going to negative side.
	    EminuspzGen = EminuspzGen - (protonGen.at(0).energy()-protonGen.at(0).pz());
	    SystemY = SystemY - protonGen.at(0);
	  }
	}
      }

      // Fill Variables
      counterproton = protonGen.size();
      xigen = System.M2()/(sqrtS*sqrtS);
      //xigenplus = SystemX.M2()/(sqrtS*sqrtS);
      xigenplus = EpluspzGen/sqrtS;
      //xigenminus = SystemY.M2()/(sqrtS*sqrtS);
      xigenminus = EminuspzGen/sqrtS;
      xigenminusMore = EminuspzGenMore/sqrtS;

      if(sumGenCastorEnergy <= 1.) gencastorlow = true;
      if(sumGenCastorEnergy > 1.) gencastorhigh = true;
      if(xigenminus < 0.1) generator = true;

      if(corcastor == "fit_castor"){
	if(debug){
	  std::cout << "Linearity: " << CastorReweight->Eval(sumGenCastorEnergy) << std::endl; //sumGenCastorEnergy 
	  std::cout << "Weight: " << totalcommon << std::endl;
	}
	if(sumGenCastorEnergy < 200.) totalcommon = totalcommon/CastorReweight->Eval(sumGenCastorEnergy); //totalcommon*CastorReweight->Eval(sumGenCastorEnergy);
      }    

      }else{

	xigen = -999.;
	xigenplus = -999.;
	xigenminus = -999.;

	gencastorlow = false;
	gencastorhigh = false;
	generator = false;

      }

      double etacastor = -5.9;
      double theta = 2*std::atan(exp(-etacastor))*180/M_PI;

      ximinuscastor = 2.2768*(sumCastorEnergy*(1-TMath::Cos(M_PI*theta/180.0)))/sqrtS;
      xipluscastor = 2.2768*(sumCastorEnergy*(1+TMath::Cos(M_PI*theta/180.0)))/sqrtS;
      xigenminuscastor = (sumGenCastorEnergy*(1-TMath::Cos(M_PI*theta/180.0)))/sqrtS;
      xigenminuscastorDef = EminuspzGenCastor/sqrtS;
      xigenpluscastor = (sumGenCastorEnergy*(1+TMath::Cos(M_PI*theta/180.0)))/sqrtS;
      xigenminusall = xigenminus + xigenminuscastor;
      xigenplusall = xigenplus + xigenpluscastor;

      if(sumGenCastorEnergy == 0){
	castorgapGen = true;
      }else{
	castorgapGen = false;
      }

      if(debug){
	if(protonGen.size()>0){
	  std::cout << "\nProton: " << protonGen.at(0) << " | eta: " << protonGen.at(0).eta() << std::endl;
	  std::cout << "Theta: " << theta << std::endl;
	  std::cout << "Cos(theta): " << TMath::Cos(M_PI*theta/180.0) << std::endl;
	  std::cout << "1-Cos(theta): " << 1-TMath::Cos(M_PI*theta/180.0) << std::endl;
	  std::cout << "1+Cos(theta): " << 1+TMath::Cos(M_PI*theta/180.0) << std::endl;
	  std::cout << "sqrtS: " << sqrtS << std::endl;
	  std::cout << "<< GEN >>" << std::endl;
	  std::cout << "Gen, xi: " << xigen << std::endl;
	  std::cout << "Gen, xi Castor plus: " << xigenpluscastor << std::endl;
	  std::cout << "Gen, xi Castor minus: " << xigenminuscastor << std::endl;
	  std::cout << "Gen, xi CMS minus: " << xigenminus << std::endl;
	  std::cout << "Gen, xi CMS plus: "  << xigenplus << std::endl;
	  std::cout << "Gen, xi all minus: " << xigenminusall << std::endl;
	  std::cout << "Gen, xi all plus: "  << xigenplusall << std::endl;
	  std::cout << "Gen, Energy Castor: " << sumGenCastorEnergy << std::endl;
	}
      }


      if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
	if(eventinfo->GetNPileUpBx0()<1) zeropileup = true;
      } 

      xi = eventdiffW->GetMx2PF()/TMath::Power(sqrtS,2);

      //if(eventdiff->GetTightNoiseFilter()>0){
      //if (eventdiff->GetNVertex() <= nVertex) vertex = true;
      //}

      if (eventdiff->GetNVertex() <= nVertex && eventdiff->GetMultiplicityTracks()<10.) vertex = true;

      resMETPt = eventdiffW->GetMETPt()/eventdiffW->GetGenNeutrinoPt();
      resMETPhi = eventdiffW->GetMETPhi()/eventdiffW->GetGenNeutrinoPhi();

      if (gapseltype == "PF" || gapseltype == "pf"){
	etamax = eventdiffW->GetEtaMaxPF();
	etamin = eventdiffW->GetEtaMinPF();
	ximinus = 2.1*eventdiffW->GetEminuspzPF()/sqrtS;
	xiplus = 2.1*eventdiffW->GetEpluspzPF()/sqrtS;
	//ximinus = eventdiff->GetXiMinusFromPFCands();
	//xiplus = eventdiff->GetXiPlusFromPFCands();
	maxLRG = eventdiffW->GetLrgPF();
	if (etamax < 3. && etamax > -10.) diffselp = true; // avoid protection -999.
	if (etamin > -3.) diffseln = true; // avoid protection -999.
      }
      else if(gapseltype == "CALO" || gapseltype == "calo"){
	etamax = eventdiffW->GetEtaCaloMax();
	etamin = eventdiffW->GetEtaCaloMin();
	ximinus = 2.1*eventdiffW->GetEminuspzCalo()/sqrtS; 
	xiplus = 2.1*eventdiffW->GetEpluspzCalo()/sqrtS;
	maxLRG = eventdiffW->GetLrgCalo();
	if (eventdiffW->GetSumEHFMinus()==0.) diffseln = true;
	if (eventdiffW->GetSumEHFPlus()==0.) diffselp = true;
      }else{
	std::cout << "" << std::endl;
	std::cout << "Please insert the correct gap type selection: " << std::endl;
	std::cout << ">> HF or hf: for hadron forward calorimeter gap selection." << std::endl;
	std::cout << ">> PF or pf: for etamax, etamin gap selection." << std::endl;
	std::cout << "" << std::endl;
	exit(EXIT_FAILURE);
      }

      xiplusall = xipluscastor + xiplus; 
      ximinusall = ximinuscastor + ximinus; 
      ximinusunique = 2.1593*(ximinuscastor/2.2768 + ximinus/2.1);
      sumGenCastorAll = sumGenCastorPions + sumGenCastorEPhoton + sumGenCastorOthers;

      if(ximinusall < 0.1) xicut = true;

      if(debug){
	//std::cout << "xi: " << xi << std::endl;
	std::cout << "<< RECO >>" << std::endl;
	std::cout << "xi Castor plus: " << xipluscastor << std::endl;
	std::cout << "xi Castor minus: " << ximinuscastor << std::endl;
	std::cout << "xi CMS minus: " << ximinus << std::endl;
	std::cout << "xi CMS plus: "  << xiplus << std::endl;
	std::cout << "xi all minus: " << ximinusall << std::endl;
	std::cout << "xi all plus: "  << xiplusall << std::endl;
	std::cout << "\nEnergy Castor Reco: " << sumCastorEnergy << std::endl;
	std::cout << "Energy Gen Castor Excluded: " << sumGenCastorEnergyExc << std::endl;
	std::cout << "Energy Gen Castor: " << sumGenCastorEnergy << std::endl;
	std::cout << "Energy Gen Castor Pion+Electron+Photon+Others: " << sumGenCastorAll << std::endl; 
	std::cout << "Energy Gen Castor Pion: " << sumGenCastorPions << std::endl;
	std::cout << "Energy Gen Castor Electron+Photon: " << sumGenCastorEPhoton << std::endl;
	std::cout << "Energy Gen Castor Others: " << sumGenCastorOthers << std::endl;
      }

      //
      // Castor
      //
      //if(castoractivity){
      //etamin = -6.;
      //}

      deltaeta = etamax - etamin;

      if(sumCastorEnergy<50.) CastorLow = true;

      if (typesel == "RecoElectron"){
	selStatus = "Reco::Electron";

	//bosonWMass = TMath::Sqrt(2.*eventdiffW->GetLeadingElectronPt()*eventdiffW->GetMETPt()*(1.-TMath::Cos((fabs(eventdiffW->GetLeadingElectronPhi()-eventdiffW->GetMETPhi())))));
	bosonWMass = TransverseMass(eventdiffW->GetLeadingElectronP4(),eventdiffW->GetMETP4()); 
	bosonWEta = eventdiffW->GetLeadingElectronEta();
	bosonWEtaGen = eventdiffW->GetGenLeadingElectronEta();
	bosonWPhi = eventdiffW->GetLeadingElectronPhi();
	bosonWPt = eventdiffW->GetLeadingElectronPt();
	bosonWCharge = eventdiffW->GetLeadingElectronCharge();
	NMuons = eventdiffW->GetMuonsN();
	NElectrons = eventdiffW->GetElectronsN();
	leptoneta = eventdiffW->GetLeadingElectronEta();
	leptonphi = eventdiffW->GetLeadingElectronPhi();
	leptonpt = eventdiffW->GetLeadingElectronPt();
	deltaphi = eventdiffW->GetLeadingElectronPhi()-eventdiffW->GetMETPhi();

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

	resLeadingPt = eventdiffW->GetLeadingElectronPt()/eventdiffW->GetGenLeadingElectronPt();
	resLeadingEta = eventdiffW->GetLeadingElectronEta()/eventdiffW->GetGenLeadingElectronEta();
	resLeadingPhi = eventdiffW->GetLeadingElectronPhi()/eventdiffW->GetGenLeadingElectronPhi();

	nocosmic = true;

	double totalASumCastor = eventdiffW->GetSumEHFMinus() + sumCastorEnergy;

	if(totalASumCastor > 0.){
	  AEcastor = (eventdiffW->GetSumEHFMinus() - sumCastorEnergy)/(eventdiffW->GetSumEHFMinus() + sumCastorEnergy);
	}

	if (eventdiffW->GetLeadingElectronPt() > lepton1pt && eventdiffW->GetMETPt() > lepton2pt) presel = true;
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

	if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
	  math::XYZTLorentzVector DiSystemGen(0.,0.,0.,0.);
	  DiSystemGen += eventdiffW->GetGenLeadingElectronP4();
	  DiSystemGen += eventdiffW->GetGenNeutrinoP4();
	  if (eventdiffW->GetGenLeadingElectronPt() > lepton1pt && eventdiffW->GetGenNeutrinoPt() > lepton2pt) preselGen = true;
	  if (DiSystemGen.M() > 60. && DiSystemGen.M() < 110.) dimassGen = true;
	}else{
	  preselGen = false;
	  dimassGen = false;
	}

      }

      else if (typesel == "RecoMuon"){
	selStatus = "Reco::Muon";

	//bosonWMass = TMath::Sqrt(2.*eventdiffW->GetLeadingMuonPt()*eventdiffW->GetMETPt()*(1.-TMath::Cos((fabs(eventdiffW->GetLeadingMuonPhi()-eventdiffW->GetMETPhi())))));
	bosonWMass = TransverseMass(eventdiffW->GetLeadingMuonP4(),eventdiffW->GetMETP4());
	bosonWEta = eventdiffW->GetLeadingMuonEta();
	bosonWEtaGen = eventdiffW->GetGenLeadingMuonEta();
	bosonWPhi = eventdiffW->GetLeadingMuonPhi();
	bosonWPt = eventdiffW->GetLeadingMuonPt();
	bosonWCharge = eventdiffW->GetLeadingMuonCharge();
	NMuons = eventdiffW->GetMuonsN();
	NElectrons = eventdiffW->GetElectronsN();
	isoRec = eventdiffW->GetLeadingMuonrelIsoDr03();
	leptoneta = eventdiffW->GetLeadingMuonEta();
	leptonphi = eventdiffW->GetLeadingMuonPhi();
	leptonpt = eventdiffW->GetLeadingMuonPt();
	deltaphi = eventdiffW->GetLeadingMuonPhi()-eventdiffW->GetMETPhi();
	secondleptoneta = eventdiffW->GetLeadingElectronEta();
	secondleptonphi = eventdiffW->GetLeadingElectronPhi();
	secondleptonpt = eventdiffW->GetLeadingElectronPt();

	resLeadingPt = eventdiffW->GetLeadingMuonPt()/eventdiffW->GetGenLeadingMuonPt();
	resLeadingEta = eventdiffW->GetLeadingMuonEta()/eventdiffW->GetGenLeadingMuonEta();
	resLeadingPhi = eventdiffW->GetLeadingMuonPhi()/eventdiffW->GetGenLeadingMuonPhi();

	if(NMuons == 1) acceptEvt = true;

	if (eventdiffW->GetLeadingMuonPt() > lepton1pt && eventdiffW->GetMETPt() > lepton2pt && eventdiffW->GetLeadingElectronPt() < 10.) presel = true;
	if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;
	if(isoRec < 0.1 && isoRec > -10.){ 
	  isolation = true;
	}

	bool debugmethod = false;
	if(debugmethod){
	  std::cout << "IsWP80? " << eventdiffW->GetLeadingElectronIsWP80() << std::endl;
	  std::cout << "IsWP95? " << eventdiffW->GetLeadingElectronIsWP95() << std::endl; 
	  std::cout << "PT: " << eventdiffW->GetLeadingMuonPt() << std::endl;
	  std::cout << "Tracker: " << eventdiffW->GetLeadingMuonTrackerHits() << std::endl;
	  std::cout << "Pixel: " << eventdiffW->GetLeadingMuonPixelHits() << std::endl;
	  std::cout << "Chi2: " << eventdiffW->GetLeadingMuonNormalizedChi2() <<std::endl;
	  std::cout << "Matched Stations: " << eventdiffW->GetLeadingMuonMatchedStations() << std::endl;
	  std::cout << "Is Tracker? " << eventdiffW->GetLeadingMuonIsTracker() << std::endl;
	  std::cout << "Is Global? " << eventdiffW->GetLeadingMuonIsGlobal() << std::endl;
	  std::cout << "Dxy: " << eventdiffW->GetLeadingMuonDxy() << std::endl;
	  std::cout << "Dz: " << eventdiffW->GetLeadingMuonDz() << std::endl;
	  std::cout << "Good? " << eventdiffW->GetLeadingMuonIsGood() << std::endl;
	  std::cout << "RelIso: " << eventdiffW->GetLeadingMuonrelIsoDr05() << std::endl;
	} 

	if (fabs(eventdiffW->GetLeadingMuonDxy() < 0.02)) nocosmic = true;
	if (eventdiffW->GetLeadingMuonIsGood()) candSel = true;    

	if (eventdiffW->GetLeadingMuonEta()>0.) WKinP = true;
	if (eventdiffW->GetLeadingMuonEta()<0.) WKinN = true;

	if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetLeadingMuonEta()); }
	else{ etasignedHF = fabs(eventdiffW->GetLeadingMuonEta()); }

	if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetLeadingMuonEta());}
	else {etasignedCASTOR = fabs(eventdiffW->GetLeadingMuonEta());}

	if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
	  math::XYZTLorentzVector DiSystemGen(0.,0.,0.,0.);
	  DiSystemGen += eventdiffW->GetGenLeadingMuonP4();
	  DiSystemGen += eventdiffW->GetGenNeutrinoP4();
	  if (eventdiffW->GetGenLeadingMuonPt() > lepton1pt && eventdiffW->GetGenNeutrinoPt() > lepton2pt) preselGen = true;
	  if (DiSystemGen.M() > 60. && DiSystemGen.M() < 110.) dimassGen = true;
	}else{
	  preselGen = false;
	  dimassGen = false;
	}

      }

      else if (typesel == "PatElectron"){
	selStatus = "Pat::Electron";

	//bosonWMass = TMath::Sqrt(2.*eventdiffW->GetPatElectron1Pt()*eventdiffW->GetMETPt()*(1.-TMath::Cos((fabs(eventdiffW->GetPatElectron1Phi()-eventdiffW->GetMETPhi())))));
	bosonWMass = TransverseMass(eventdiffW->GetPatElectron1P4(),eventdiffW->GetPatMETP4());
	bosonWEta = eventdiffW->GetPatElectron1Eta();
	bosonWEtaGen = eventdiffW->GetGenLeadingElectronEta();
	bosonWPhi = eventdiffW->GetPatElectron1Phi();
	bosonWPt = eventdiffW->GetPatElectron1Pt();
	bosonWCharge = eventdiffW->GetPatMuon1Charge();
	NMuons = eventdiffW->GetPatNMuon();
	NElectrons = eventdiffW->GetPatNElectron();
	leptoneta = eventdiffW->GetPatElectron1Eta();
	leptonphi = eventdiffW->GetPatElectron1Phi();
	leptonpt = eventdiffW->GetPatElectron1Pt();
	deltaphi = eventdiffW->GetPatElectron1Phi()-eventdiffW->GetMETPhi();

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

	resLeadingPt = eventdiffW->GetPatElectron1Pt()/eventdiffW->GetGenLeadingElectronPt();
	resLeadingEta = eventdiffW->GetPatElectron1Eta()/eventdiffW->GetGenLeadingElectronEta();
	resLeadingPhi = eventdiffW->GetPatElectron1Phi()/eventdiffW->GetGenLeadingElectronPhi();

	if (eventdiffW->GetPatElectron1Pt() > lepton1pt && eventdiffW->GetPatMETPt() > lepton2pt) presel = true;
	if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;

	if(NElectrons==1 && NMuons==0) acceptEvt = true;

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

	if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
	  math::XYZTLorentzVector DiSystemGen(0.,0.,0.,0.);
	  DiSystemGen += eventdiffW->GetGenLeadingElectronP4();
	  DiSystemGen += eventdiffW->GetGenNeutrinoP4();
	  if (eventdiffW->GetGenLeadingElectronPt() > lepton1pt && eventdiffW->GetGenNeutrinoPt() > lepton2pt) preselGen = true;
	  if (DiSystemGen.M() > 60. && DiSystemGen.M() < 110.) dimassGen = true;
	}else{
	  preselGen = false;
	  dimassGen = false;
	}


      }

      else if(typesel == "PatMuon"){
	selStatus = "Pat::Muon";

	//bosonWMass = TMath::Sqrt(2.*eventdiffW->GetPatMuon1Pt()*eventdiffW->GetMETPt()*(1.-TMath::Cos((fabs(eventdiffW->GetPatMuon1Phi()-eventdiffW->GetMETPhi())))));
	bosonWMass = TransverseMass(eventdiffW->GetPatMuon1P4(),eventdiffW->GetPatMETP4());
	bosonWEta = eventdiffW->GetPatMuon1Eta();
	bosonWEtaGen = eventdiffW->GetGenLeadingMuonEta();
	bosonWPhi = eventdiffW->GetPatMuon1Phi();
	bosonWPt = eventdiffW->GetPatMuon1Pt();
	bosonWCharge = eventdiffW->GetPatMuon1Charge();
	NMuons = eventdiffW->GetPatNMuon();
	NElectrons = eventdiffW->GetPatNElectron();
	isoRec = eventdiffW->GetLeadingMuonrelIsoDr03();
	leptoneta = eventdiffW->GetPatMuon1Eta();
	leptonphi = eventdiffW->GetPatMuon1Phi();
	leptonpt = eventdiffW->GetPatMuon1Pt();
	deltaphi = eventdiffW->GetPatMuon1Phi()-eventdiffW->GetMETPhi();

	resLeadingPt = eventdiffW->GetPatMuon1Pt()/eventdiffW->GetGenLeadingMuonPt();
	resLeadingEta = eventdiffW->GetPatMuon1Eta()/eventdiffW->GetGenLeadingMuonEta();
	resLeadingPhi = eventdiffW->GetPatMuon1Phi()/eventdiffW->GetGenLeadingMuonPhi();

	if(NMuons ==1 && NElectrons == 0) acceptEvt = true;

	if (eventdiffW->GetPatMuon1Pt() > lepton1pt && eventdiffW->GetPatMETPt() > lepton2pt) presel = true;
	if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;
	if (isoRec < 0.1 && isoRec > -10.) {
	  candSel = true;
	  isolation = true;
	}

	if (eventdiffW->GetPatMuon1Eta()>0.) WKinP = true;
	if (eventdiffW->GetPatMuon1Eta()<0.) WKinN = true;

	if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetPatMuon1Eta()); }
	else{ etasignedHF = fabs(eventdiffW->GetPatMuon1Eta()); }

	if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetPatMuon1Eta());}
	else {etasignedCASTOR = fabs(eventdiffW->GetPatMuon1Eta());}

	if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
	  math::XYZTLorentzVector DiSystemGen(0.,0.,0.,0.);
	  DiSystemGen += eventdiffW->GetGenLeadingMuonP4();
	  DiSystemGen += eventdiffW->GetGenNeutrinoP4();
	  if (eventdiffW->GetGenLeadingMuonPt() > lepton1pt && eventdiffW->GetGenNeutrinoPt() > lepton2pt) preselGen = true;
	  if (DiSystemGen.M() > 60. && DiSystemGen.M() < 110.) dimassGen = true;
	}else{
	  preselGen = false;
	  dimassGen = false;
	}


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
      bSumEnergyHFMinus = eventdiffW->GetSumEHFMinus();
      bSumEnergyHFPlus = eventdiffW->GetSumEHFPlus();
      bsumCastorEnergy = sumCastorEnergy;
      bSectorCastorHit = counterHit;
      bMaxGap = maxLRG;
      if(fabs(eventdiffW->GetSumptPFLeft())>fabs(eventdiffW->GetSumptPFRight())){
	SumPTMaxLrgPF = eventdiffW->GetSumptPFLeft();
	SumPTMinLrgPF = eventdiffW->GetSumptPFRight();
      }else{
	SumPTMaxLrgPF = eventdiffW->GetSumptPFRight();
	SumPTMinLrgPF = eventdiffW->GetSumptPFLeft();
      }
      bXiPlus = xiplus;
      bXiMinus = ximinus;
      betasignedHF = etasignedHF;
      betasignedCASTOR = etasignedCASTOR;
      bAEcastor = AEcastor;
      bdeltaeta = deltaeta;
      betamax = etamax;
      betamin = etamin;

      if(pileup < 21){ // Never comment this line. It is the program defense.

	if(switchtrigger == "trigger" || switchtrigger == "trigger_all_electron" || switchtrigger == "trigger_all_muon"){ 
	  FillHistos(0,pileup,totalcommon); 
	  histo_cutflow->Fill("total",totalcommon); 
	  if(trigger) {
	    ++totalT;
	    FillHistos(1,pileup,totalcommon);
	    histo_cutflow->Fill("trigger",totalcommon);
	  } 
	  if(trigger && vertex) {
	    FillHistos(2,pileup,totalcommon);
	    histo_cutflow->Fill("step2",totalcommon);
	  }
	  if(trigger && vertex && presel){
	    histo_cutflow->Fill("step3",totalcommon);
	    FillHistos(3,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic){
	    histo_cutflow->Fill("step4",totalcommon);
	    FillHistos(4,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel){
	    histo_cutflow->Fill("step5",totalcommon);
	    FillHistos(5,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && acceptEvt){
	    histo_cutflow->Fill("step6",totalcommon);
	    FillHistos(6,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && acceptEvt && dimass) {
	    histo_cutflow->Fill("step7",totalcommon);
	    FillHistos(7,pileup,totalcommon);
	    fOutW->cd();
	    troutW->SetWeight(totalcommon);
	    troutW->Fill();
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln){
	    histo_cutflow->Fill("NGapCMS",totalcommon);
	    FillHistos(8,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp){
	    histo_cutflow->Fill("PGapCMS",totalcommon);
	    FillHistos(9,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln && castorgap){
	    histo_cutflow->Fill("NGapCMSAndCASTOR",totalcommon);
	    FillHistos(10,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity){
	    histo_cutflow->Fill("PGapCMSAndCastorActivity",totalcommon);
	    FillHistos(11,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln && WKinP){
	    histo_cutflow->Fill("NGapCMSAndWKinP",totalcommon);
	    outstring << "HF- Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	    FillHistos(12,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && WKinN){
	    histo_cutflow->Fill("PGapCMSAndWKinN",totalcommon);
	    outstring << "HF+ Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	    FillHistos(13,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln && castorgap && WKinP){
	    histo_cutflow->Fill("NGapCMSAndCASTORAndWKinP",totalcommon);
	    FillHistos(14,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity && WKinN){
	    histo_cutflow->Fill("PGapCMSAndCastorActivityAndWKinN",totalcommon);
	    FillHistos(15,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && castorgap){
	    histo_cutflow->Fill("NGapCASTOR",totalcommon);
	    FillHistos(16,pileup,totalcommon);
	    fOutCASTOR->cd();
	    troutCASTOR->SetWeight(totalcommon);
	    troutCASTOR->Fill();
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && castorgap && WKinP){
	    histo_cutflow->Fill("NGapCASTORAndWKinP",totalcommon);
	    outstring << "CASTOR Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	    FillHistos(17,pileup,totalcommon);
	    fOut->cd();
	    trout->SetWeight(totalcommon);
	    trout->Fill();
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castorgap){
	    histo_cutflow->Fill("PGapCMSAndCASTOR",totalcommon);
	    FillHistos(18,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castorgap && WKinP){
	    histo_cutflow->Fill("PGapCMSAndCASTORAndWKinP",totalcommon);
	    FillHistos(19,pileup,totalcommon);
	  }
	  if(trigger && vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && xicut){
	    histo_cutflow->Fill("#xi",totalcommon);
	    FillHistos(25,pileup,totalcommon);
	  }
	}

	else if (switchtrigger =="no_trigger_nocorrection" || switchtrigger == "no_trigger_correction" ){
	  --totalT;
	  FillHistos(0,pileup,totalcommon);
	  histo_cutflow->Fill("total",totalcommon);
	  if(vertex){
	    histo_cutflow->Fill("step2",totalcommon);
	    FillHistos(2,pileup,totalcommon);
	  }
	  if(vertex && presel){
	    histo_cutflow->Fill("step3",totalcommon);
	    FillHistos(3,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic){
	    histo_cutflow->Fill("step4",totalcommon);
	    FillHistos(4,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel){
	    histo_cutflow->Fill("step5",totalcommon);
	    FillHistos(5,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && acceptEvt){
	    histo_cutflow->Fill("step6",totalcommon);
	    FillHistos(6,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && acceptEvt && dimass) {
	    histo_cutflow->Fill("step7",totalcommon);
	    FillHistos(7,pileup,totalcommon);
	    fOutW->cd();
	    troutW->SetWeight(totalcommon);
	    troutW->Fill();
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln){
	    histo_cutflow->Fill("NGapCMS",totalcommon);
	    FillHistos(8,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp){
	    histo_cutflow->Fill("PGapCMS",totalcommon);
	    FillHistos(9,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln && castorgap){
	    histo_cutflow->Fill("NGapCMSAndCASTOR",totalcommon);
	    FillHistos(10,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity){
	    histo_cutflow->Fill("PGapCMSAndCastorActivity",totalcommon);
	    FillHistos(11,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln && WKinP){
	    histo_cutflow->Fill("NGapCMSAndWKinP",totalcommon);
	    FillHistos(12,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && WKinN){
	    histo_cutflow->Fill("PGapCMSAndWKinN",totalcommon);
	    FillHistos(13,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffseln && castorgap && WKinP){
	    histo_cutflow->Fill("NGapCMSAndCASTORAndWKinP",totalcommon);
	    FillHistos(14,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castoractivity && WKinN){
	    histo_cutflow->Fill("PGapCMSAndCastorActivityAndWKinN",totalcommon);
	    FillHistos(15,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && castorgap){
	    histo_cutflow->Fill("NGapCASTOR",totalcommon);
	    FillHistos(16,pileup,totalcommon);
	    fOutCASTOR->cd();
	    troutCASTOR->SetWeight(totalcommon);
	    troutCASTOR->Fill();
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && castorgap && WKinP){
	    histo_cutflow->Fill("NGapCASTORAndWKinP",totalcommon);
	    FillHistos(17,pileup,totalcommon);
	    fOut->cd();
	    trout->SetWeight(totalcommon);
	    trout->Fill();
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castorgap){
	    histo_cutflow->Fill("PGapCMSAndCASTOR",totalcommon);
	    FillHistos(18,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && diffselp && castorgap && WKinP){
	    histo_cutflow->Fill("PGapCMSAndCASTORAndWKinP",totalcommon);
	    FillHistos(19,pileup,totalcommon);
	  }
	  if(generator){
	    histo_cutflow->Fill("generator",totalcommon);
	    FillHistos(20,pileup,totalcommon);
	  }
	  if(gencastorlow){
	    histo_cutflow->Fill("gencastorlow",totalcommon);
	    FillHistos(21,pileup,totalcommon);
	  }
	  if(gencastorhigh){
	    histo_cutflow->Fill("gencastorhigh",totalcommon);
	    FillHistos(22,pileup,totalcommon);
	  }
	  if(zeropileup && vertex){
	    histo_cutflow->Fill("zeropileup",totalcommon);
	    FillHistos(23,pileup,totalcommon);
	  }
	  if(zeropileup && vertex && castorgap){
	    histo_cutflow->Fill("zeropileup_NGapCASTOR",totalcommon);
	    FillHistos(24,pileup,totalcommon);
	  }
	  if(vertex && presel && nocosmic && isolation && candSel && dimass && acceptEvt && xicut){
	    histo_cutflow->Fill("xi",totalcommon);
	    FillHistos(25,pileup,totalcommon);
	  }
	  if(preselGen && dimassGen){
	    histo_cutflow->Fill("generatorW",totalcommon);
	    FillHistos(26,pileup,totalcommon);
	  }
	  if(preselGen && dimassGen && castorgapGen){
	    histo_cutflow->Fill("generatorWDiff",totalcommon);
	    FillHistos(27,pileup,totalcommon);
	  }

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
      foldersFile[4] = outf->mkdir("Generator");
      histo_cutflow->Write();
      SaveHistos(type,typesel);
      outf->Close();
      std::cout << "\n" << std::endl;

      fOut->cd();
      trout->Write();
      fOut->Close();

      fOutW->cd();
      troutW->Write();
      fOutW->Close();

      fOutCASTOR->cd();
      troutCASTOR->Write();
      fOutCASTOR->Close();

      in_file->cd();
      in_file->Close();

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
      std::string corcastorfile_;
      std::string corcastor_;

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
      if (argc > 19 && strcmp(s1,argv[19]) != 0) corcastor_ = argv[19];
      if (argc > 20 && strcmp(s1,argv[20]) != 0) corcastorfile_ = argv[20];

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
      std::cout << "Castor Correction: " << corcastor_ <<std::endl;
      std::cout << "Castor Correction file: " << corcastorfile_ <<std::endl;

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

	if(corcastor_ == "fit_castor"){ 
	  TFile fitcastorfile(corcastorfile_.c_str());
	  if (fitcastorfile.IsZombie()){
	    std::cout << "-------------------------------------------" << std::endl;
	    std::cout << " There is no the file " << corcastorfile_ << " or the" << std::endl;
	    std::cout << " path is not correct." << std::endl;
	    std::cout << "-------------------------------------------" << std::endl;
	    return 0;
	  }
	}

	DiffractiveW* diffWRun = new DiffractiveW();
	clock_t tStart = clock();
	diffWRun->CreateHistos(type_);
	diffWRun->Run(filein_, processname_, savehistofile_, switchtrigger_, optTrigger_, lepton1pt_, lepton2pt_, nVertex_, type_, switchlumiweight_, mcweight_, typesel_, castorthreshold_, channelsthreshold_, castorcorrfile_, gapseltype_, pumfile_, pudfile_, corcastor_, corcastorfile_);
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
