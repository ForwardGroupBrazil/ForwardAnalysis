//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsCastorAnalysissAnalysis#Macro_Analysis
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
#include "CastorAnalysis.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveZEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace diffractiveAnalysis;
using namespace diffractiveZAnalysis;
using namespace eventInfo;
using namespace reweight;

void CastorAnalysis::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventdiff = new DiffractiveEvent();
  eventCastor = new DiffractiveZEvent();
  eventinfo = new EventInfoEvent();
  diff = tr->GetBranch("DiffractiveAnalysis");
  Castor = tr->GetBranch("DiffractiveZAnalysis");
  info = tr->GetBranch("EventInfo");
  diff->SetAddress(&eventdiff);
  Castor->SetAddress(&eventCastor);
  info->SetAddress(&eventinfo);

}

void CastorAnalysis::CreateHistos(){

  std::string step0 = "without_cuts"; 
  std::string step1 = "with_trigger"; 
  std::string step2 = "step2"; 
  std::string step3 = "zeropileup";

  Folders.push_back(step0);
  Folders.push_back(step1);
  Folders.push_back(step2);
  Folders.push_back(step3);

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name1[300];
    char name2[300];

    sprintf(name1,"CastorCentroid_%s",Folders.at(j).c_str());
    TH2F *histo_castor_centroid = new TH2F(name1,"Castor Centroid Energy; x[cm]; y[cm]",30,-15,15,30,-15,15);
    m_hVector_histo_castor_centroid.push_back(histo_castor_centroid);

    sprintf(name1,"CastorCentroidPhi_%s",Folders.at(j).c_str());
    TH1F *histo_castor_centroid_phi = new TH1F(name1,"Castor Centroid Energy; Sector(#phi,i) [degree]; NEvents",720,0,360);
    m_hVector_histo_castor_centroid_phi.push_back(histo_castor_centroid_phi);

    sprintf(name1,"ECaloVsEta_%s",Folders.at(j).c_str());
    TH2F *histo_ECaloVsEta = new TH2F(name1,"Calorimeter Energy X #eta; #eta; Energy [GeV]", 500, -6, 6, 100, 0., 1000.);
    m_hVector_ECaloVsEta.push_back(histo_ECaloVsEta);

    sprintf(name1,"sumECastorMinus_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorMinus = new TH1F(name1,"Castor Sum of Energy; Energy [GeV]; N events",1500,0,1500); // 6000,0,3000
    m_hVector_sumECastorMinus.push_back(histo_sumECastorMinus);

    sprintf(name1,"gensumECastorMinus_%s",Folders.at(j).c_str());
    TH1D *histo_gensumECastorMinus = new TH1D(name1,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500); // 6000,0,3000
    m_hVector_gensumECastorMinus.push_back(histo_gensumECastorMinus);

    sprintf(name1,"gensumECastorPionsMinus_%s",Folders.at(j).c_str());
    TH1D *histo_gensumECastorPionsMinus = new TH1D(name1,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500);
    m_hVector_gensumECastorPionsMinus.push_back(histo_gensumECastorPionsMinus);

    sprintf(name1,"gensumECastorEPhotonMinuss_%s",Folders.at(j).c_str());
    TH1D *histo_gensumECastorEPhotonMinus = new TH1D(name1,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500);
    m_hVector_gensumECastorEPhotonMinus.push_back(histo_gensumECastorEPhotonMinus);

    sprintf(name1,"gensumECastorOthersMinus_%s",Folders.at(j).c_str());
    TH1D *histo_gensumECastorOthersMinus = new TH1D(name1,"Castor Sum of Energy; Energy, gen [GeV]; N events",1500,0,1500);
    m_hVector_gensumECastorOthersMinus.push_back(histo_gensumECastorOthersMinus);

    sprintf(name1,"ECastorSector_%s",Folders.at(j).c_str());
    TH2F *histo_ECastorSector = new TH2F(name1,"Castor Energy X Sector; Sector; Energy [GeV]", 18, 0, 18, 44, 0., 220.);
    m_hVector_ECastorSector.push_back(histo_ECastorSector);

    sprintf(name1,"ECaloVsEtaTProf_%s",Folders.at(j).c_str());
    TProfile *histo_ECaloVsEtaTProf = new TProfile(name1,"Calorimeter Energy X #eta; #eta; <Energy> [GeV]", 100, -6, 6, 0., 1000.);
    m_hVector_ECaloVsEtaTProf.push_back(histo_ECaloVsEtaTProf);

    sprintf(name1,"ECastorSectorTProf_%s",Folders.at(j).c_str());
    TProfile *histo_ECastorSectorTProf = new TProfile(name1,"Castor Energy X Sector; Sector; <Energy> [GeV]", 100, 0, 18, 0., 7000.);
    m_hVector_ECastorSectorTProf.push_back(histo_ECastorSectorTProf); 

    sprintf(name1,"ECastorSectorBin1D_%s",Folders.at(j).c_str());
    TH1F *histo_ECastorSectorBin1D = new TH1F(name1,"Castor Energy X Sector; Sector; Energy [GeV]", 18, 0, 18);
    m_hVector_ECastorSectorBin1D.push_back(histo_ECastorSectorBin1D);

    sprintf(name1,"EnergyHFMinusVsCastor_%s",Folders.at(j).c_str());
    TH2F *histo_ET_Calos_n = new TH2F(name1,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 6000, 0., 3000. );
    m_hVector_etcalos_n.push_back(histo_ET_Calos_n);

    sprintf(name1,"EnergyHFPlusVsCastorTProf_%s",Folders.at(j).c_str());
    TProfile *histo_EnergyHFPlusVsCastorTProf = new TProfile(name1,"HF^{+} and Castor; #sum Energy HF^{+}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 3000. );
    m_hVector_EnergyHFPlusVsCastorTProf.push_back(histo_EnergyHFPlusVsCastorTProf);

    sprintf(name1,"EnergyHFMinusVsCastorTProf_%s",Folders.at(j).c_str());
    TProfile *histo_EnergyHFMinusVsCastorTProf = new TProfile(name1,"HF^{-} and Castor; #sum Energy HF^{-}; #sum Energy Castor [GeV]; N events", 1000, 0., 1000., 0., 3000. );
    m_hVector_EnergyHFMinusVsCastorTProf.push_back(histo_EnergyHFMinusVsCastorTProf);

    sprintf(name1,"sumECastorAndSumHFMinus_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorAndHFMinus = new TH1F(name1,"HF^{-} and Castor Sum of Energy; Energy [GeV]; N events",6000,0,3000);
    m_hVector_sumECastorAndHFMinus.push_back(histo_sumECastorAndHFMinus);

    sprintf(name1,"CastorMultiplicity_%s",Folders.at(j).c_str());
    TH1F *histo_CastorMultiplicity = new TH1F(name1,"Castor: number of sectors with activity; N Sectors (#phi); N events",18,0,18);
    m_hVector_CastorMultiplicity.push_back(histo_CastorMultiplicity);

    sprintf(name1,"CastorMultiplicityVsLumi_%s",Folders.at(j).c_str());
    TH2F *histo_CastorMultiplicityVsLumi = new TH2F(name1,"CastorMultiplicity Vs Luminosity; Luminosity per Bunch [#mub^{-1}s^{-1}]; Castor Multiplicity",5000,0,2,18,0,18);
    m_hVector_CastorMultiplicityVsLumi.push_back(histo_CastorMultiplicityVsLumi);

    sprintf(name1,"RunNumber_%s",Folders.at(j).c_str());
    TH1F *histo_RunNumber = new TH1F(name1,"Run Number; Run Number; N Event",16000,134000,150000);
    m_hVector_RunNumber.push_back(histo_RunNumber);

    sprintf(name1,"RunNumberZeroCastor_%s",Folders.at(j).c_str());
    TH1F *histo_RunNumberZeroCastor = new TH1F(name1,"Run Number; Run Number; N Event",16000,134000,150000);
    m_hVector_RunNumberZeroCastor.push_back(histo_RunNumberZeroCastor);

    sprintf(name1,"RunNumberHighCastor_%s",Folders.at(j).c_str());
    TH1F *histo_RunNumberHighCastor = new TH1F(name1,"Run Number; Run Number; N Event",16000,134000,150000);
    m_hVector_RunNumberHighCastor.push_back(histo_RunNumberHighCastor);

    sprintf(name1,"SectorVsTotalCastorEnergyTProf_%s",Folders.at(j).c_str());
    TProfile *histo_SectorVsTotalCastorEnergyTProf = new TProfile(name1,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",18,0,18,0,1500);
    m_hVector_SectorVsTotalCastorEnergyTProf.push_back(histo_SectorVsTotalCastorEnergyTProf);

    sprintf(name1,"SectorVsTotalCastorEnergy_%s",Folders.at(j).c_str());
    TH2F *histo_SectorVsTotalCastorEnergy = new TH2F(name1,"Castor Multiplicity Vs CastorEnergy; # Multiplicity; Castor Energy [GeV]",18,0,18,6000,0,3000);
    m_hVector_SectorVsTotalCastorEnergy.push_back(histo_SectorVsTotalCastorEnergy);

    for (int se=1; se<17; se++){
      m_hVector_TotalEnergySectors.push_back( std::vector<TH1F*>() );
      m_hVector_Sector_EnergyVsMultiplicity.push_back( std::vector<TH2F*>() );
      m_hVector_AlongZ_EnergyVsModule.push_back( std::vector<TH2F*>() );
      m_hVector_AlongZ_EnergyVsModuleTProf.push_back( std::vector<TProfile*>() );

      sprintf(name1,"TotalEnergySector%d_%s",se,Folders.at(j).c_str());
      sprintf(name2,"#sum Energy, Castor Sector %d; #sum E_{modules 1,2,3,4,5} [GeV]; N events", se);
      TH1F *histo_TotalEnergySectors = new TH1F(name1,name2,1000,0,500);
      m_hVector_TotalEnergySectors[se-1].push_back(histo_TotalEnergySectors);

      sprintf(name1,"Sector%d_EnergyVsMultiplicity_%s",se,Folders.at(j).c_str());
      sprintf(name2,"Energy for Sector %d Vs Multiplicity; Multiplicity; #sum Energy_{sector} [GeV]", se);
      TH2F *histo_Sector_EnergyVsMultiplicity = new TH2F(name1,name2,18,0,18,1000,0,500);
      m_hVector_Sector_EnergyVsMultiplicity[se-1].push_back(histo_Sector_EnergyVsMultiplicity);

      sprintf(name1,"AlongZ_Sector%d_EnergyVsModule_%s",se,Folders.at(j).c_str());
      sprintf(name2,"Energy Along Z, Sector %d Vs Module; Module Id(Z); Energy [GeV]", se);
      TH2F *histo_AlongZ_EnergyVsModule = new TH2F(name1,name2,7,0,7,4000,0,2000);
      m_hVector_AlongZ_EnergyVsModule[se-1].push_back(histo_AlongZ_EnergyVsModule);

      sprintf(name1,"AlongZ_Sector%d_EnergyVsModuleTProf_%s",se,Folders.at(j).c_str());
      sprintf(name2,"Energy Along Z, Sector %d Vs Module; Module Id (Z); Energy [GeV]", se);
      TProfile *histo_AlongZ_EnergyVsModuleTProf = new TProfile(name1,name2,7,0,7,0,2000);
      m_hVector_AlongZ_EnergyVsModuleTProf[se-1].push_back(histo_AlongZ_EnergyVsModuleTProf);
    }

    for (int bs=1; bs<5; bs++){
      m_hVector_sumEHFplusBinSlice.push_back( std::vector<TH1F*>() );
      m_hVector_sumEHFminusBinSlice.push_back( std::vector<TH1F*>() );
      m_hVector_sumECastorMinusBinSlice.push_back( std::vector<TH1F*>() );

      sprintf(name1,"sumEHFplusBinSlice%d_%s",bs,Folders.at(j).c_str());
      TH1F *histo_sumEHFplusBinSlice = new TH1F(name1,"HF^{+} - Sum of Energy; #sum E_{HF^{+}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHFplusBinSlice[bs-1].push_back(histo_sumEHFplusBinSlice);

      sprintf(name1,"sumEHFminusBinSlice%d_%s",bs,Folders.at(j).c_str());
      TH1F *histo_sumEHFminusBinSlice = new TH1F(name1,"HF^{-} - Sum of Energy; #sum E_{HF^{-}} [GeV]; N events",2000,0,2000);
      m_hVector_sumEHFminusBinSlice[bs-1].push_back(histo_sumEHFminusBinSlice);

      sprintf(name1,"sumECastorMinusBinSlice%d_%s",bs,Folders.at(j).c_str());
      TH1F *histo_sumECastorMinusBinSlice = new TH1F(name1,"Castor Sum of Energy; Energy [GeV]; N events",6000,0,3000);
      m_hVector_sumECastorMinusBinSlice[bs-1].push_back(histo_sumECastorMinusBinSlice);
    }

    sprintf(name1,"CastorBadChannelVsRun_%s",Folders.at(j).c_str());
    TH2F *histo_CastorBadChannelVsRun = new TH2F(name1,"Castor Bad Channel Vs Run; # Castor Bad Channel(id); Run",226,0,226,16000,134000,150000);
    m_hVector_CastorBadChannelVsRun.push_back(histo_CastorBadChannelVsRun);

    sprintf(name1,"CastorBadRuns_%s",Folders.at(j).c_str());
    TH1F *histo_CastorBadRuns = new TH1F(name1,"Castor Bad Run Number; Bad Run Number; N Event",16000,134000,150000);
    m_hVector_CastorBadRuns.push_back(histo_CastorBadRuns);

    sprintf(name1,"CastorBadChannels_%s",Folders.at(j).c_str());
    TH1F *histo_CastorBadChannels = new TH1F(name1,"Castor Bad Channels Number; Channel(Id); N Event",225,0,225);
    m_hVector_CastorBadChannels.push_back(histo_CastorBadChannels);

    sprintf(name1,"MultiplicityPerModule_%s",Folders.at(j).c_str());
    TH2F *histo_CastorMultiplicityPerModule = new TH2F(name1,"Multiplicity per Module; Module Id (Z); Multiplicity",7,0,7,17,0,17);
    m_hVector_CastorMultiplicityPerModule.push_back(histo_CastorMultiplicityPerModule);

    sprintf(name1,"MultiplicityPerModuleTProf_%s",Folders.at(j).c_str());
    TProfile *histo_CastorMultiplicityPerModuleTProf = new TProfile(name1,"Multiplicity per Module; Module Id (Z); Multiplicity",7,0,7,0,17);
    m_hVector_CastorMultiplicityPerModuleTProf.push_back(histo_CastorMultiplicityPerModuleTProf);

    for (int i=1;i<6;i++){
      m_hVector_CastorMultiplicityModule.push_back( std::vector<TH1F*>() );
      m_hVector_CastorMultiplicityModuleAll.push_back( std::vector<TH2F*>() );

      sprintf(name1,"CastorMultiplicityModule%d_%s",i,Folders.at(j).c_str());
      sprintf(name2,"Castor Module %d: number of channels with activity; N Sectors (#phi); N events", i);
      TH1F *histo_CastorMultiplicityModule = new TH1F(name1,name2,17,0,17);
      m_hVector_CastorMultiplicityModule[i-1].push_back(histo_CastorMultiplicityModule);

      sprintf(name1,"CastorMultiplicityModule%dVsTotal_%s",i,Folders.at(j).c_str());
      sprintf(name2,"Castor Module %d Vs All Modules Multiplicities; N Castor Module %d; N Castor Module 1,2,3,4 and 5", i, i);
      TH2F *histo_CastorMultiplicityModuleAll = new TH2F(name1,name2,17,0,17,17,0,17);
      m_hVector_CastorMultiplicityModuleAll[i-1].push_back(histo_CastorMultiplicityModuleAll);

    }

    sprintf(name1,"EnergyPerModule_%s",Folders.at(j).c_str());
    TH2F *histo_CastorEnergyPerModule = new TH2F(name1,"Energy per Module; Module Id (Z); Energy [GeV]",6,0,6,4000,0,2000);
    m_hVector_CastorEnergyPerModule.push_back(histo_CastorEnergyPerModule);

    sprintf(name1,"EnergyPerModuleTProf_%s",Folders.at(j).c_str());
    TProfile *histo_CastorEnergyPerModuleTProf = new TProfile(name1,"Energy per Module; Module Id (Z); Energy [GeV]",6,0,6,0,2000);
    m_hVector_CastorEnergyPerModuleTProf.push_back(histo_CastorEnergyPerModuleTProf);

    for (int i=1;i<6;i++){
      m_hVector_CastorModuleFraction.push_back( std::vector<TH1F*>() );
      sprintf(name1,"CastorFraction_Module_%d_%s",i,Folders.at(j).c_str());
      sprintf(name2,"Castor Module %d; #frac{N_module %d}{N_module(1,2,3,4,5)}; N events", i,i);
      TH1F *histo_CastorModuleFraction = new TH1F(name1,name2,100,0,1);
      m_hVector_CastorModuleFraction[i-1].push_back(histo_CastorModuleFraction);
    }

    sprintf(name1,"CastorMappingMultiplicities3D_%s",Folders.at(j).c_str());
    TH2F *histo_CastorMappingMultiplicity3D = new TH2F(name1,"Mapping Multiplicities; Sector (#phi); Module Id (Z)",16,0,16,5,1,6);
    m_hVector_CastorMappingMultiplicity3D.push_back(histo_CastorMappingMultiplicity3D);

    sprintf(name1,"CastorMappingEnergy3D_%s",Folders.at(j).c_str());
    TH2F *histo_CastorMappingEnergy3D = new TH2F(name1,"Mapping Energy [GeV]; Sector (#phi); Module Id (Z)",16,0,16,5,1,6);
    m_hVector_CastorMappingEnergy3D.push_back(histo_CastorMappingEnergy3D);

    for (int i=1;i<5;i++){
      m_hVector_CastorMappingMultiplicity.push_back( std::vector<TH2F*>() );
      m_hVector_CastorMappingEnergy.push_back( std::vector<TH2F*>() );

      sprintf(name1,"CastorMappingMultiplicitiesSelection%d_%s",i-1,Folders.at(j).c_str());
      TH2F *histo_CastorMappingMultiplicity = new TH2F(name1,"Mapping Multiplicities; Module Id (Z); Sector (#phi)",5,1,6,16,1,17);
      m_hVector_CastorMappingMultiplicity[i-1].push_back(histo_CastorMappingMultiplicity);

      sprintf(name1,"CastorMappingEnergySelection%d_%s",i-1,Folders.at(j).c_str());
      TH2F *histo_CastorMappingEnergy = new TH2F(name1,"Mapping Energy [GeV]; Module Id (Z); Sector (#phi)",5,1,6,16,1,17);
      m_hVector_CastorMappingEnergy[i-1].push_back(histo_CastorMappingEnergy);
    }

    sprintf(name1,"FirstModuleHitCastor_%s",Folders.at(j).c_str());
    TH1F *histo_FirstModuleHitCastor = new TH1F(name1,"First Module Hit in CASTOR; Module Id (Z); N Events",7,0,7);
    m_hVector_FirstModuleHitCastor.push_back(histo_FirstModuleHitCastor);

    sprintf(name1,"CastorMultiplicityChannels_%s",Folders.at(j).c_str());
    TH1F *histo_CastorMultiplicityChannels = new TH1F(name1,"Castor: number of channels with activity; N channels; N events",80,0,80);
    m_hVector_CastorMultiplicityChannels.push_back(histo_CastorMultiplicityChannels);

  }

  for (int i=1;i<6;i++){
    for (int j=1;j<6;j++){
      m_hVector_CastorMappingEnergySnapshot.push_back( std::vector<TH2F*>() );
      m_hVector_CastorMappingMultiplicitySnapshot.push_back( std::vector<TH2F*>() );

      char namea[300];
      char nameb[300];

      sprintf(namea,"CastorMappingMultiplicities_module_%d_snapshot_%d",i,j);
      sprintf(nameb,"Castor Mapping Module %d, Snapshot %d. Multiplicities; Sector (#phi); Module Id (Z)", i,j);
      TH2F *histo_CastorMappingMultiplicitySnapshot = new TH2F(namea,nameb,16,1,17,5,1,6);
      m_hVector_CastorMappingMultiplicitySnapshot[i-1].push_back(histo_CastorMappingMultiplicitySnapshot);

      sprintf(namea,"CastorMappingEnergy_module_%d_snapshot_%d",i,j);
      sprintf(nameb,"Castor Mapping Module %d, Snapshot %d. Energy [GeV]; Sector (#phi); Module Id (Z)", i, j);
      TH2F *histo_CastorMappingEnergySnapshot = new TH2F(namea,nameb,16,1,17,5,1,6);
      m_hVector_CastorMappingEnergySnapshot[i-1].push_back(histo_CastorMappingEnergySnapshot);
    }
  }

}

void CastorAnalysis::FillHistos(int index, double totalweight){

  bool debug = false;
  bool debugmult = false;
  bool debug_correction = false;
  bool debug_factor = false;
  bool debub_error = false;

  sumCastorEnergy = 0.; 
  sumCastorAndHFMinusEnergy = 0.;
  SectorZeroCastorCounter = 0.;
  num_phi = 0.;
  num_x_centroid = 0;
  num_y_centroid = 0.;
  x_temp = 0.;
  y_temp = 0.;
  x_centroid = 0.;
  y_centroid = 0.;
  phi_average = 0.;
  SectorCastorHit = 0;

  int multicounter[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double energymodule[5]={0.,0.,0.,0.,0.};
  double castorId[16] = {11.25,33.75,56.25,78.75,101.25,123.75,146.25,168.75,191.25,213.75,236.25,258.75,281.25,303.75,326.75,348.75};
  int mapping[5][16];
  double energymapping[5][16];
  double energycorr[5][16];
  double energymappingunc[5][16];

  for (int i=0; i<5; i++){
    for (int j=0; j<16; j++){
      mapping[i][j]=0;
      energymapping[i][j]=0;
      energycorr[i][j]=0;
      energymappingunc[i][j]=0;
    }
  }

  for(int i=1; i<6;i++){
    for(int j=1; j<17; j++){
      if (switchtrigger == "no_trigger_correction"){
	energycorr[i-1][j-1]=h_castor_channel->GetBinContent(i,j);
	if (debug_correction) std::cout << "GetBinContent(" << i << "," << j << "): " << h_castor_channel->GetBinContent(i,j) << std::endl;
	if (debug_factor) std::cout << "Factor: "<< energycorr[i-1][j-1] << std::endl;
      }else{
	energycorr[i-1][j-1]=1.;
	if (debug_factor) std::cout << "Factor: "<< energycorr[i-1][j-1] << std::endl;
      }
    }
  }

  if (ChannelThreshold > 0. && SectorThreshold < 0.){
    sprintf(selCastor,">> Castor Channel Threshold: %0.2f GeV",ChannelThreshold);
    for (int i=0; i < 16; i++){
      CastorEnergySector[i]=0.;
      if (i==4 || i==5){ // removing channel 5 and 6, module1
	if (eventCastor->GetCastorModule2Energy(i)*energycorr[1][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule2Energy(i)*energycorr[1][i];
	if (eventCastor->GetCastorModule3Energy(i)*energycorr[2][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule3Energy(i)*energycorr[2][i];
	if (eventCastor->GetCastorModule4Energy(i)*energycorr[3][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule4Energy(i)*energycorr[3][i];
	if (eventCastor->GetCastorModule5Energy(i)*energycorr[4][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule5Energy(i)*energycorr[4][i];
	if (debug) std::cout << "Removing Channel " << i << "."<< std::endl;
      }else{
	if (eventCastor->GetCastorModule1Energy(i)*energycorr[0][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule1Energy(i)*energycorr[0][i];
	if (eventCastor->GetCastorModule2Energy(i)*energycorr[1][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule2Energy(i)*energycorr[1][i];
	if (eventCastor->GetCastorModule3Energy(i)*energycorr[2][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule3Energy(i)*energycorr[2][i];
	if (eventCastor->GetCastorModule4Energy(i)*energycorr[3][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule4Energy(i)*energycorr[3][i];
	if (eventCastor->GetCastorModule5Energy(i)*energycorr[4][i] > ChannelThreshold) CastorEnergySector[i] += eventCastor->GetCastorModule5Energy(i)*energycorr[4][i];
      }

      if (eventCastor->GetCastorModule1Energy(i)*energycorr[0][i] > ChannelThreshold){
	if (i==0 || i==1 || i==2 || i==3 || i==6 || i==7 || i==8 || i==9 || i==10 || i==11 || i==12 || i==13 || i==14 || i==15){ //remove channel 5 and 6, module1
	  multicounter[0]++;
	  energymodule[0]+=eventCastor->GetCastorModule1Energy(i)*energycorr[0][i];
	  m_hVector_AlongZ_EnergyVsModule[i].at(index)->Fill(1,eventCastor->GetCastorModule1Energy(i)*energycorr[0][i],totalweight);
	  m_hVector_AlongZ_EnergyVsModuleTProf[i].at(index)->Fill(1,eventCastor->GetCastorModule1Energy(i)*energycorr[0][i],totalweight);
	  m_hVector_CastorMappingEnergy3D.at(index)->Fill(i,1,eventCastor->GetCastorModule1Energy(i)*energycorr[0][i]);
	  m_hVector_CastorMappingMultiplicity3D.at(index)->Fill(i,1,totalweight);
	  mapping[0][i]=1;
	  energymapping[0][i]=eventCastor->GetCastorModule1Energy(i)*energycorr[0][i];
	  // Parameters from http://indico.cern.ch/getFile.py/access?contribId=73&sessionId=14&resId=0&materialId=slides&confId=96930
	  energymappingunc[0][i]=energymapping[0][i]*sqrt((0.18/energymapping[0][i])+0.002);
	  energymappinguncsqr[0][0][i]+=pow(energymappingunc[0][i],2);
	  if (debub_error) std::cout << "Mapping(1,"<<i<<"), Energy: " << energymapping[0][i] << ", error: " << energymappingunc[0][i] << std::endl;
	  if (debug) std::cout << "Using channel " << i << ", module 1." << std::endl;
	}
      }
      if (eventCastor->GetCastorModule2Energy(i)*energycorr[1][i] > ChannelThreshold){
	multicounter[1]++;
	energymodule[1]+=eventCastor->GetCastorModule2Energy(i)*energycorr[1][i];
	m_hVector_AlongZ_EnergyVsModule[i].at(index)->Fill(2,eventCastor->GetCastorModule2Energy(i)*energycorr[1][i],totalweight);
	m_hVector_AlongZ_EnergyVsModuleTProf[i].at(index)->Fill(2,eventCastor->GetCastorModule2Energy(i)*energycorr[1][i],totalweight);
	m_hVector_CastorMappingEnergy3D.at(index)->Fill(i,2,eventCastor->GetCastorModule2Energy(i)*energycorr[1][i]);
	m_hVector_CastorMappingMultiplicity3D.at(index)->Fill(i,2,totalweight);
	mapping[1][i]=1;
	energymapping[1][i]=eventCastor->GetCastorModule2Energy(i)*energycorr[1][i];
	// Parameters from http://indico.cern.ch/getFile.py/access?contribId=73&sessionId=14&resId=0&materialId=slides&confId=96930
	energymappingunc[1][i]=energymapping[1][i]*sqrt((0.18/energymapping[1][i])+0.002);
	energymappinguncsqr[0][1][i]+=pow(energymappingunc[1][i],2);
	if (debub_error) std::cout << "Mapping(2,"<<i<<"), Energy: " << energymapping[1][i] << ", error: " << energymappingunc[1][i] << std::endl;
      }
      if (eventCastor->GetCastorModule3Energy(i)*energycorr[2][i] > ChannelThreshold){
	multicounter[2]++;
	energymodule[2]+=eventCastor->GetCastorModule3Energy(i)*energycorr[2][i];
	m_hVector_AlongZ_EnergyVsModule[i].at(index)->Fill(3,eventCastor->GetCastorModule3Energy(i)*energycorr[2][i],totalweight);
	m_hVector_AlongZ_EnergyVsModuleTProf[i].at(index)->Fill(3,eventCastor->GetCastorModule3Energy(i)*energycorr[2][i],totalweight);
	m_hVector_CastorMappingEnergy3D.at(index)->Fill(i,3,eventCastor->GetCastorModule3Energy(3)*energycorr[2][i]);
	m_hVector_CastorMappingMultiplicity3D.at(index)->Fill(i,3,totalweight);
	mapping[2][i]=1;
	energymapping[2][i]=eventCastor->GetCastorModule3Energy(i)*energycorr[2][i];
	// Parameters from http://indico.cern.ch/getFile.py/access?contribId=73&sessionId=14&resId=0&materialId=slides&confId=96930
	energymappingunc[2][i]=energymapping[2][i]*sqrt((3.49/energymapping[2][i])+0.033);
	energymappinguncsqr[0][2][i]+=pow(energymappingunc[2][i],2);
	if (debub_error) std::cout << "Mapping(3,"<<i<<"), Energy: " << energymapping[2][i] << ", error: " << energymappingunc[2][i] << std::endl;
      }
      if (eventCastor->GetCastorModule4Energy(i)*energycorr[3][i] > ChannelThreshold){
	multicounter[3]++;
	energymodule[3]+=eventCastor->GetCastorModule4Energy(i)*energycorr[3][i];
	m_hVector_AlongZ_EnergyVsModule[i].at(index)->Fill(4,eventCastor->GetCastorModule4Energy(i)*energycorr[3][i],totalweight);
	m_hVector_AlongZ_EnergyVsModuleTProf[i].at(index)->Fill(4,eventCastor->GetCastorModule4Energy(i)*energycorr[3][i],totalweight);
	m_hVector_CastorMappingEnergy3D.at(index)->Fill(i,4,eventCastor->GetCastorModule4Energy(i)*energycorr[3][i]);
	m_hVector_CastorMappingMultiplicity3D.at(index)->Fill(i,4,totalweight);
	mapping[3][i]=1;
	energymapping[3][i]=eventCastor->GetCastorModule4Energy(i)*energycorr[3][i];
	// Parameters from http://indico.cern.ch/getFile.py/access?contribId=73&sessionId=14&resId=0&materialId=slides&confId=96930
	energymappingunc[3][i]=energymapping[3][i]*sqrt((3.49/energymapping[3][i])+0.033);
	energymappinguncsqr[0][3][i]+=pow(energymappingunc[3][i],2);
	if (debub_error) std::cout << "Mapping(4,"<<i<<"), Energy: " << energymapping[3][i] << ", error: " << energymappingunc[3][i] << std::endl;
      }
      if (eventCastor->GetCastorModule5Energy(i)*energycorr[4][i] > ChannelThreshold){
	multicounter[4]++;
	energymodule[4]+=eventCastor->GetCastorModule5Energy(i)*energycorr[4][i];
	m_hVector_AlongZ_EnergyVsModule[i].at(index)->Fill(5,eventCastor->GetCastorModule5Energy(i)*energycorr[4][i],totalweight);
	m_hVector_AlongZ_EnergyVsModuleTProf[i].at(index)->Fill(5,eventCastor->GetCastorModule5Energy(i)*energycorr[4][i],totalweight);
	m_hVector_CastorMappingEnergy3D.at(index)->Fill(i,5,eventCastor->GetCastorModule5Energy(i)*energycorr[4][i]);
	m_hVector_CastorMappingMultiplicity3D.at(index)->Fill(i,5,totalweight);
	mapping[4][i]=1;
	energymapping[4][i]=eventCastor->GetCastorModule5Energy(i)*energycorr[4][i];
	// Parameters from http://indico.cern.ch/getFile.py/access?contribId=73&sessionId=14&resId=0&materialId=slides&confId=96930
	energymappingunc[4][i]=energymapping[4][i]*sqrt((3.49/energymapping[4][i])+0.033);
	energymappinguncsqr[0][4][i]+=pow(energymappingunc[4][i],2);
	if (debub_error) std::cout << "Mapping(5,"<<i<<"), Energy: " << energymapping[4][i] << ", error: " << energymappingunc[4][i] << std::endl;
      }
      if (CastorEnergySector[i] > ChannelThreshold){
	SectorCastorHit++;
	sumCastorEnergy+=CastorEnergySector[i];
	x_temp = 15*TMath::Cos(castorId[i]*M_PI/180.0);
	y_temp = 15*TMath::Sin(castorId[i]*M_PI/180.0);
	num_phi += castorId[i]*CastorEnergySector[i];
	num_x_centroid += x_temp*CastorEnergySector[i];
	num_y_centroid += y_temp*CastorEnergySector[i];
	m_hVector_ECastorSector.at(index)->Fill(i+1,CastorEnergySector[i],totalweight);
	m_hVector_ECastorSectorTProf.at(index)->Fill(i+1,CastorEnergySector[i],totalweight);
	m_hVector_ECastorSectorBin1D.at(index)->Fill(i+1,CastorEnergySector[i]);
      }
      else{
	sumCastorEnergy+=0;
	SectorZeroCastorCounter++;
	num_x_centroid += 0;
	num_y_centroid += 0;
	num_phi += 0;
	m_hVector_ECastorSector.at(index)->Fill(i+1,0);
      }
    }

    for (int j=0; j<5; j++){
      if (debug){ std::cout << "\nMultiplicity Module " << j+1 << ":" << multicounter[j] << std::endl;
	std::cout << "Castor Multiplicity: " << SectorCastorHit << std::endl;
      }
      m_hVector_CastorMultiplicityModule[j].at(index)->Fill(multicounter[j],totalweight);
      m_hVector_CastorMultiplicityPerModule.at(index)->Fill(j+1,multicounter[j],totalweight);
      m_hVector_CastorMultiplicityPerModuleTProf.at(index)->Fill(j+1,multicounter[j],totalweight);
      m_hVector_CastorEnergyPerModule.at(index)->Fill(j+1,energymodule[j],totalweight);
      m_hVector_CastorEnergyPerModuleTProf.at(index)->Fill(j+1,energymodule[j],totalweight);
      m_hVector_CastorMultiplicityModuleAll[j].at(index)->Fill(multicounter[j],SectorCastorHit,totalweight);
    }

    if (sumCastorEnergy > 0.){
      x_centroid = num_x_centroid/sumCastorEnergy;
      y_centroid = num_y_centroid/sumCastorEnergy;
      phi_average = num_phi/sumCastorEnergy;
      m_hVector_histo_castor_centroid.at(index)->Fill(x_centroid,y_centroid);
      m_hVector_histo_castor_centroid_phi.at(index)->Fill(phi_average);
    }

    if (debug){
      std::cout << "SectorCastorHit: " << SectorCastorHit << std::endl;
      std::cout << "Module 1: " << multicounter[0] <<  std::endl;
      std::cout << "Module 2: " << multicounter[1] <<  std::endl;
      std::cout << "Module 3: " << multicounter[2] <<  std::endl;
      std::cout << "Module 4: " << multicounter[3] <<  std::endl;
      std::cout << "Module 5: " << multicounter[4] <<  std::endl;
    }

    int sumMulti = multicounter[0]+multicounter[1]+multicounter[2]+multicounter[3]+multicounter[4];
    double percentage[5];

    for (int l=0;l<5;l++){
      percentage[l] = 0.;
      if (sumMulti > 0) {percentage[l] = multicounter[l]/(double)sumMulti;
	m_hVector_CastorModuleFraction[l].at(index)->Fill(percentage[l]);
      }if (debug){
	std::cout << "Percentage: " << percentage[l] << std::endl;
	std::cout << "Module["<<l+1<<"]: " << multicounter[l] << std::endl;
	std::cout << "Castor total multiplicity: " << sumMulti << std::endl;
      }
    }

    if (index==7){
      for (int l=0; l<5; l++){
	if (percentage[l] >= 0.4 && chk[l]<5) {
	  for (int j=0; j<5; j++){
	    for (int k=0;k<16;k++){
	      if (mapping[j][k]==1){
		m_hVector_CastorMappingEnergySnapshot[l].at(chk[l])->Fill(k+1,j+1,energymapping[j][k]);
		m_hVector_CastorMappingMultiplicitySnapshot[l].at(chk[l])->Fill(k+1,j+1);
		if (debug) std::cout << "Mapping[" << j << "," << k << "]" << " | chk[l]: "<< chk[l] << " | Module: " << l <<  std::endl;
	      }
	    }
	  }
	  chk[l]++;
	  if (debug) std::cout << "Counter["<< l << "]: "<< chk[l] << std::endl;
	}
      }
    }

    for(int id=0; id < 16; id++){
      if (CastorEnergySector[id] > ChannelThreshold){
	m_hVector_TotalEnergySectors[id].at(index)->Fill(CastorEnergySector[id],totalweight);
	m_hVector_Sector_EnergyVsMultiplicity[id].at(index)->Fill(SectorCastorHit,CastorEnergySector[id],totalweight);
      }else{
	m_hVector_TotalEnergySectors[id].at(index)->Fill(0.);
      }
    }

    for(int h=0;h<5;h++){ // Loop Modules
      for (int i=1;i<17;i++){ //loop sectors
	m_hVector_CastorMappingEnergy[0].at(index)->Fill(h+1,i,energymapping[h][i-1]);
	m_hVector_CastorMappingMultiplicity[0].at(index)->Fill(h+1,i,mapping[h][i-1]);

	/* In case, multiple cuts.
	//for(int l=0;l<6;l++){ // Loop Cuts
	if(energymodule[x] > l*50 && energymodule[x] < l*50+50){
	m_hVector_CastorMappingEnergy[l+1].at(index)->Fill(x+1,i,energymapping[x][i-1]);
	m_hVector_CastorMappingMultiplicity[l+1].at(index)->Fill(x+1,i,mapping[x][i-1]);
	energymappinguncsqr[l+1][x][i-1]+= pow(energymappingunc[x][i-1],2); 
	}

	}*/

	if(energymodule[h] > 0 && energymodule[h] < 50){
	  m_hVector_CastorMappingEnergy[1].at(index)->Fill(h+1,i,energymapping[h][i-1]);
	  m_hVector_CastorMappingMultiplicity[1].at(index)->Fill(h+1,i,mapping[h][i-1]);
	  energymappinguncsqr[1][h][i-1]+= pow(energymappingunc[h][i-1],2);
	}

	if(energymodule[h] > 50 && energymodule[h] < 100){
	  m_hVector_CastorMappingEnergy[2].at(index)->Fill(h+1,i,energymapping[h][i-1]);
	  m_hVector_CastorMappingMultiplicity[2].at(index)->Fill(h+1,i,mapping[h][i-1]);
	  energymappinguncsqr[2][h][i-1]+= pow(energymappingunc[h][i-1],2);
	}

	if(energymodule[h] > 100 && energymodule[h] < 300){
	  m_hVector_CastorMappingEnergy[3].at(index)->Fill(h+1,i,energymapping[h][i-1]);
	  m_hVector_CastorMappingMultiplicity[3].at(index)->Fill(h+1,i,mapping[h][i-1]);
	  energymappinguncsqr[3][h][i-1]+= pow(energymappingunc[h][i-1],2);
	}

      }
    }

    int checkcounter[5]={0};
    for(int i=0;i<16;i++){
      checkcounter[0] += mapping[0][i];
      checkcounter[1] += mapping[1][i];
      checkcounter[2] += mapping[2][i];
      checkcounter[3] += mapping[3][i];
      checkcounter[4] += mapping[4][i];
    }

    if(checkcounter[0]>0){
      m_hVector_FirstModuleHitCastor.at(index)->Fill(1,totalweight);
    }
    if(checkcounter[0]==0 && checkcounter[1]>0){
      m_hVector_FirstModuleHitCastor.at(index)->Fill(2,totalweight);
    }
    if(checkcounter[0]==0 && checkcounter[1]==0 && checkcounter[2]>0){
      m_hVector_FirstModuleHitCastor.at(index)->Fill(3,totalweight);
    }
    if(checkcounter[0]==0 && checkcounter[1]==0 && checkcounter[2]==0 && checkcounter[3]>0){
      m_hVector_FirstModuleHitCastor.at(index)->Fill(4,totalweight);
    }
    if(checkcounter[0]==0 && checkcounter[1]==0 && checkcounter[2]==0 && checkcounter[3]==0 && checkcounter[4]>0){
      m_hVector_FirstModuleHitCastor.at(index)->Fill(5,totalweight);
    }

    int totalchannels = checkcounter[0] + checkcounter[1] + checkcounter[2] + checkcounter[3] + checkcounter[4];
    m_hVector_CastorMultiplicityChannels.at(index)->Fill(totalchannels,totalweight);

  }

  if (ChannelThreshold < 0. && SectorThreshold > 0.){
    sprintf(selCastor,">> Castor Sector Threshold: %0.2f GeV",SectorThreshold);
    for (int i=0; i < 16; i++){
      CastorEnergySector[i]=0.;
      if (i==4 || i==5){ // removing channel 5 and 6, module1
	CastorEnergySector[i] += eventCastor->GetCastorModule2Energy(i)*energycorr[1][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule3Energy(i)*energycorr[2][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule4Energy(i)*energycorr[3][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule5Energy(i)*energycorr[4][i];
      }else{
	CastorEnergySector[i] += eventCastor->GetCastorModule1Energy(i)*energycorr[0][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule2Energy(i)*energycorr[1][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule3Energy(i)*energycorr[2][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule4Energy(i)*energycorr[3][i];
	CastorEnergySector[i] += eventCastor->GetCastorModule5Energy(i)*energycorr[4][i];
      }
      if (CastorEnergySector[i] > SectorThreshold){
	SectorCastorHit++;
	sumCastorEnergy+=CastorEnergySector[i];
	x_temp = 15*TMath::Cos(castorId[i]*M_PI/180.0);
	y_temp = 15*TMath::Sin(castorId[i]*M_PI/180.0);
	num_phi += castorId[i]*CastorEnergySector[i];
	num_x_centroid += x_temp*CastorEnergySector[i];
	num_y_centroid += y_temp*CastorEnergySector[i];
	m_hVector_ECastorSector.at(index)->Fill(i+1,CastorEnergySector[i],totalweight);
	m_hVector_ECastorSectorTProf.at(index)->Fill(i+1,CastorEnergySector[i],totalweight);
	m_hVector_ECastorSectorBin1D.at(index)->Fill(i+1,CastorEnergySector[i]);
      }
      else{
	sumCastorEnergy+=0;
	SectorZeroCastorCounter++;
	num_x_centroid += 0;
	num_y_centroid += 0;
	num_phi += 0;
	m_hVector_ECastorSector.at(index)->Fill(i+1,0);
      }
    }

    if (sumCastorEnergy > 0.){
      x_centroid = num_x_centroid/sumCastorEnergy;
      y_centroid = num_y_centroid/sumCastorEnergy;
      phi_average = num_phi/sumCastorEnergy;
      m_hVector_histo_castor_centroid.at(index)->Fill(x_centroid,y_centroid);
      m_hVector_histo_castor_centroid_phi.at(index)->Fill(phi_average);
    }

    for(int id=0; id < 16; id++){
      if (CastorEnergySector[id] > SectorThreshold){
	m_hVector_TotalEnergySectors[id].at(index)->Fill(CastorEnergySector[id],totalweight);
	m_hVector_Sector_EnergyVsMultiplicity[id].at(index)->Fill(SectorCastorHit,CastorEnergySector[id],totalweight);
      }else{
	m_hVector_TotalEnergySectors[id].at(index)->Fill(0.);
      }
    }

  }

  if (debugmult && (multicounter[0] > SectorCastorHit)) {
    std::cout << selCastor << std::endl;
    std::cout << "\nBug: multicounter[" << 0 << "]:" << multicounter[0] << " | SectorCounter: " << SectorCastorHit << std::endl;
    std::cout << "Sum Castor Energy: " << sumCastorEnergy << std::endl;
    std::cout << "Sector Threshold: " << SectorThreshold << " | Channel Threshold: " << ChannelThreshold << std::endl;
    for (int i=0; i < 16; i++){
      std::cout << "Energy per Sector " << i+1 << " (m1 + m2 + m3 + m4 + m5): "<< CastorEnergySector[i] << std::endl;
      std::cout << "Energy per Module 1: " << eventCastor->GetCastorModule1Energy(i)*energycorr[0][i] << std::endl;
      std::cout << "Energy per Module 2: " << eventCastor->GetCastorModule2Energy(i)*energycorr[1][i] << std::endl;
      std::cout << "Energy per Module 3: " << eventCastor->GetCastorModule3Energy(i)*energycorr[2][i] << std::endl;
      std::cout << "Energy per Module 4: " << eventCastor->GetCastorModule4Energy(i)*energycorr[3][i] << std::endl;
      std::cout << "Energy per Module 5: " << eventCastor->GetCastorModule5Energy(i)*energycorr[4][i] << std::endl;
    }
  }

  for (int k=0; k<eventCastor->GetEachTowerCounter();k++){
    m_hVector_ECaloVsEta.at(index)->Fill(eventCastor->GetEachTowerEta(k),eventCastor->GetEachTowerEnergy(k),totalweight);
    m_hVector_ECaloVsEtaTProf.at(index)->Fill(eventCastor->GetEachTowerEta(k),eventCastor->GetEachTowerEnergy(k),totalweight);
  }

  sumCastorAndHFMinusEnergy = sumCastorEnergy+eventdiff->GetSumEnergyHFMinus();

  m_hVector_sumECastorMinus.at(index)->Fill(sumCastorEnergy);
  m_hVector_gensumECastorMinus.at(index)->Fill(sumGenCastorEnergy);
  m_hVector_gensumECastorPionsMinus.at(index)->Fill(sumGenCastorPions);
  m_hVector_gensumECastorEPhotonMinus.at(index)->Fill(sumGenCastorEPhoton);
  m_hVector_gensumECastorOthersMinus.at(index)->Fill(sumGenCastorOthers);

  m_hVector_etcalos_n.at(index)->Fill(eventCastor->GetSumEHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_EnergyHFPlusVsCastorTProf.at(index)->Fill(eventCastor->GetSumEHFPlus(),sumCastorEnergy,totalweight);
  m_hVector_EnergyHFMinusVsCastorTProf.at(index)->Fill(eventCastor->GetSumEHFMinus(),sumCastorEnergy,totalweight);
  m_hVector_sumECastorAndHFMinus.at(index)->Fill(sumCastorAndHFMinusEnergy,totalweight);
  m_hVector_CastorMultiplicity.at(index)->Fill(SectorCastorHit,totalweight);
  m_hVector_CastorMultiplicityVsLumi.at(index)->Fill(eventinfo->GetInstLumiBunch(),SectorCastorHit,totalweight);
  if (SectorCastorHit < 1) m_hVector_RunNumberZeroCastor.at(index)->Fill(eventdiff->GetRunNumber());
  if (SectorCastorHit > 15) m_hVector_RunNumberHighCastor.at(index)->Fill(eventdiff->GetRunNumber());
  m_hVector_RunNumber.at(index)->Fill(eventdiff->GetRunNumber());
  m_hVector_SectorVsTotalCastorEnergy.at(index)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);
  m_hVector_SectorVsTotalCastorEnergyTProf.at(index)->Fill(SectorCastorHit,sumCastorEnergy,totalweight);

  for (int l=0; l<eventCastor->GetCastorNumberBadChannels(); l++){
    m_hVector_CastorBadChannelVsRun.at(index)->Fill(eventCastor->GetCastorBadChannels(l),eventdiff->GetRunNumber());
    m_hVector_CastorBadChannels.at(index)->Fill(eventCastor->GetCastorBadChannels(l));
    if (eventCastor->GetCastorBadChannels(l) >=1 && eventCastor->GetCastorBadChannels(l) <= 80) m_hVector_CastorBadRuns.at(index)->Fill(eventdiff->GetRunNumber());
  }

  if (SectorCastorHit >= 0 && SectorCastorHit <= 4){
    m_hVector_sumECastorMinusBinSlice[0].at(index)->Fill(sumCastorEnergy,totalweight);
    m_hVector_sumEHFplusBinSlice[0].at(index)->Fill(eventCastor->GetSumEHFPlus(),totalweight);
    m_hVector_sumEHFminusBinSlice[0].at(index)->Fill(eventCastor->GetSumEHFMinus(),totalweight);
  }

  if (SectorCastorHit >= 5 && SectorCastorHit <= 8){
    m_hVector_sumECastorMinusBinSlice[1].at(index)->Fill(sumCastorEnergy,totalweight);
    m_hVector_sumEHFplusBinSlice[1].at(index)->Fill(eventCastor->GetSumEHFPlus(),totalweight);
    m_hVector_sumEHFminusBinSlice[1].at(index)->Fill(eventCastor->GetSumEHFMinus(),totalweight);
  }

  if (SectorCastorHit >= 9 && SectorCastorHit <= 12){
    m_hVector_sumECastorMinusBinSlice[2].at(index)->Fill(sumCastorEnergy,totalweight);
    m_hVector_sumEHFplusBinSlice[2].at(index)->Fill(eventCastor->GetSumEHFPlus(),totalweight);
    m_hVector_sumEHFminusBinSlice[2].at(index)->Fill(eventCastor->GetSumEHFMinus(),totalweight);
  }

  if (SectorCastorHit >= 13 && SectorCastorHit <= 16){
    m_hVector_sumECastorMinusBinSlice[3].at(index)->Fill(sumCastorEnergy,totalweight);
    m_hVector_sumEHFplusBinSlice[3].at(index)->Fill(eventCastor->GetSumEHFPlus(),totalweight);
    m_hVector_sumEHFminusBinSlice[3].at(index)->Fill(eventCastor->GetSumEHFMinus(),totalweight);
  }

}

void CastorAnalysis::SaveHistos(){


  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    m_hVector_histo_castor_centroid[j]->Write();
    m_hVector_histo_castor_centroid_phi[j]->Write();
    m_hVector_ECastorSector[j]->Write();
    m_hVector_ECastorSectorTProf[j]->Write();
    m_hVector_ECastorSectorBin1D[j]->Write();
    m_hVector_ECaloVsEta[j]->Write();
    m_hVector_ECaloVsEtaTProf[j]->Write();
    m_hVector_sumECastorMinus[j]->Write();
    m_hVector_gensumECastorMinus[j]->Write();
    m_hVector_gensumECastorPionsMinus[j]->Write();
    m_hVector_gensumECastorEPhotonMinus[j]->Write();
    m_hVector_gensumECastorOthersMinus[j]->Write();
    m_hVector_etcalos_n[j]->Write();
    m_hVector_EnergyHFPlusVsCastorTProf[j]->Write();
    m_hVector_EnergyHFMinusVsCastorTProf[j]->Write();
    m_hVector_sumECastorAndHFMinus[j]->Write();
    m_hVector_FirstModuleHitCastor[j]->Write();
    m_hVector_CastorMultiplicity[j]->Write();
    m_hVector_CastorMultiplicityChannels[j]->Write();
    m_hVector_CastorMultiplicityVsLumi[j]->Write();
    m_hVector_RunNumberZeroCastor[j]->Write();
    m_hVector_RunNumberHighCastor[j]->Write();
    m_hVector_RunNumber[j]->Write();
    m_hVector_SectorVsTotalCastorEnergy[j]->Write();
    m_hVector_SectorVsTotalCastorEnergyTProf[j]->Write();
    m_hVector_CastorBadChannelVsRun[j]->Write();
    m_hVector_CastorBadRuns[j]->Write();
    m_hVector_CastorBadChannels[j]->Write();  
    m_hVector_CastorMultiplicityPerModule[j]->Write();
    m_hVector_CastorMultiplicityPerModuleTProf[j]->Write();
    m_hVector_CastorEnergyPerModule[j]->Write();
    m_hVector_CastorEnergyPerModuleTProf[j]->Write();
    m_hVector_CastorMappingMultiplicity3D[j]->Write();
    m_hVector_CastorMappingEnergy3D[j]->Write();

    for (int id=0; id<16; id++){
      m_hVector_TotalEnergySectors[id].at(j)->Write();
      m_hVector_Sector_EnergyVsMultiplicity[id].at(j)->Write();
      m_hVector_AlongZ_EnergyVsModule[id].at(j)->Write();
      m_hVector_AlongZ_EnergyVsModuleTProf[id].at(j)->Write();
    }

    for (int id=0; id<4; id++){
      m_hVector_sumECastorMinusBinSlice[id].at(j)->Write();
      m_hVector_sumEHFplusBinSlice[id].at(j)->Write();
      m_hVector_sumEHFminusBinSlice[id].at(j)->Write();
    }

    for (int id=0;id<5;id++){
      m_hVector_CastorMultiplicityModule[id].at(j)->Write();
      m_hVector_CastorMultiplicityModuleAll[id].at(j)->Write();
      m_hVector_CastorModuleFraction[id].at(j)->Write();
    }

    for (int id=0;id<4;id++){
      for (int i=0;i<16;i++){
	m_hVector_CastorMappingEnergy[id].at(j)->SetBinError(1,i+1,sqrt(energymappinguncsqr[id][0][i]));
	m_hVector_CastorMappingEnergy[id].at(j)->SetBinError(2,i+1,sqrt(energymappinguncsqr[id][1][i]));
	m_hVector_CastorMappingEnergy[id].at(j)->SetBinError(3,i+1,sqrt(energymappinguncsqr[id][2][i]));
	m_hVector_CastorMappingEnergy[id].at(j)->SetBinError(4,i+1,sqrt(energymappinguncsqr[id][3][i]));
	m_hVector_CastorMappingEnergy[id].at(j)->SetBinError(5,i+1,sqrt(energymappinguncsqr[id][4][i]));
      }
    }

    for (int id=0;id<4;id++){
      m_hVector_CastorMappingMultiplicity[id].at(j)->Write();
      m_hVector_CastorMappingEnergy[id].at(j)->Write();
    }

  }

  for (int i=0;i<5;i++){
    for (int j=0;j<5;j++){
      m_hVector_CastorMappingMultiplicitySnapshot[i].at(j)->Write();
      m_hVector_CastorMappingEnergySnapshot[i].at(j)->Write();
    }
  }

}

void CastorAnalysis::Run(std::string filein_, std::string processname_, std::string savehistofile_, std::string switchtrigger_, int optTrigger_, double lepton1pt_, double lepton2pt_, int nVertex_, std::string typesel_, std::string switchlumiweight_, double mcweight_, double SectorThreshold_, double ChannelThreshold_, std::string channelcorrfile_, std::string pumfile_, std::string pudfile_){

  bool debug = false;

  for (int i=0;i<5;i++){
    chk[i] = 0;
  }

  for (int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      for(int k=0;k<16;k++){
	energymappinguncsqr[i][j][k]=0;
      }
    }
  }

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  switchtrigger = switchtrigger_;
  nVertex = nVertex_;
  optTrigger = optTrigger_;
  lepton1pt = lepton1pt_;
  lepton2pt = lepton2pt_;
  typesel = typesel_;
  mcweight = mcweight_;
  switchlumiweight = switchlumiweight_;
  SectorThreshold = SectorThreshold_;
  ChannelThreshold = ChannelThreshold_;
  channelcorrfile = channelcorrfile_;
  pumfile = pumfile_;
  pudfile = pudfile_;

  std::string selStatus;
  std::string TriggerStatus;

  TFile check1(filein.c_str());
  TFile check2(channelcorrfile.c_str());

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

  if (check2.IsZombie()){
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " There is no the file " << channelcorrfile << " or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;
  }
  else if (!check2.GetListOfKeys()->Contains("channelcorrector")){
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " There is no CASTOR channel corrector histogram" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;
  }

  h_castor_channel = (TH2F*)check2.Get("channelcorrector");

  TFile *outf = new TFile(savehistofile.c_str(),"RECREATE");
  TString outtxt = savehistofile;
  outtxt.ReplaceAll("root","txt");
  std::ofstream outstring(outtxt);

  outstring << "" << std::endl;
  outstring << "<< Gold Events >>" << std::endl;
  outstring << "" << std::endl;
  outstring << "Please, insert this events in another text file to be used by PickEvent Tool. " << std::endl;
  outstring << "Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents " << std::endl;
  outstring << ">>---------------------------------------------------------------------- " << std::endl;

  int triggercounter[20]={0};
  int totalT=0;

  isoTk1 = -999.;
  isoTk2 = -999.;
  isoEcal1 = -999.;
  isoEcal2 = -999.;
  isoHcal1 = -999.;
  isoHcal2 = -999.;  
  innerHits1 = -999;
  Dcot1 = -999.;
  Dist1 = -999.;
  DeltaEtaTkClu1 = -999.;
  DeltaPhiTkClu1 = -999.;
  sigmaIeIe1 = -999.;
  HE1 = -999.;
  innerHits2 = -999;
  Dcot2 = -999.;
  Dist2 = -999.;
  DeltaEtaTkClu2 = -999.;
  DeltaPhiTkClu2 = -999.;
  sigmaIeIe2 = -999.;
  HE2 = -999.;

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  NEntries = tr->GetEntries();
  std::cout << "" << std::endl;
  std::cout<< "Reading Tree: "<< NEntries << " events"<<std::endl;
  std::cout << "" << std::endl;

  std::string status;  

  // Lumiweight PU
  edm::LumiReWeighting LumiWeights_(pumfile.c_str(),pudfile.c_str(),"pileUpBx0_complete_without_cuts","pileup");

  for(int i=0;i<NEntries;i++){

    tr->GetEntry(i);

    for (int nt=0;nt<20;nt++){
      if(eventCastor->GetHLTPath(nt)>0){
	triggercounter[nt]++;
      }
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

    bool zeropileup = false;

    if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
      if(eventinfo->GetNPileUpBx0()<1) zeropileup = true;


      sumGenCastorPions = 0.;
      sumGenCastorEPhoton = 0.;
      sumGenCastorOthers = 0.;
      sumGenCastorEnergy = 0.;

      for (int k=0; k<eventCastor->GetNParticlesGen();k++){

	int part_id = fabs(eventCastor->GetParticlesPDGidGen(k));

	if (eventCastor->GetParticlesEtaGen(k)<-5.2 && eventCastor->GetParticlesEtaGen(k)>-6.6){
	  if(part_id != 12 || part_id != 14 || part_id != 16 || part_id != 13){
	    sumGenCastorEnergy += eventCastor->GetParticlesEnergyGen(k);
	  }
	  if(part_id == 111 || part_id == 211){
	    sumGenCastorPions += eventCastor->GetParticlesEnergyGen(k);
	  }else if(part_id == 11 || part_id == 22){
	    sumGenCastorEPhoton += eventCastor->GetParticlesEnergyGen(k);
	  }else{
	    sumGenCastorOthers += eventCastor->GetParticlesEnergyGen(k);
	  }
	}
      }

    }

    if (!debug){
      if (i==0) {
	std::cout << "" << std::endl;
	std::cout<< "Status Bar" << std::endl;
	std::cout << "" << std::endl;
      }
      loadBar(i,NEntries);
    }

    if (debug){
      std::cout << " " << std::endl;
      std::cout << "# of active towers: " << eventCastor->GetEachTowerCounter() << std::endl;
      std::cout << " " <<std::endl;
      for (int ic=0; ic<eventCastor->GetEachTowerCounter(); ic++){
	std::cout << "Each Tower Energy: " << eventCastor->GetEachTowerEnergy(ic) << std::endl;
      }
      std::cout << " " <<std::endl;
    }

    if (debug){
      if( i % 1000 == 0 ){
	std::cout << "\nEvent " << i << std::endl;
	std::cout << "Nr. events Bx 0: " << eventinfo->GetNPileUpBx0() << std::endl;
	std::cout << "Lumi per Bunch: " << eventinfo->GetInstLumiBunch() << std::endl;
      }
    }

    bool triggerE_a = false;
    bool triggerE_b = false;
    bool triggerE_c = false;
    bool trigger = false;
    bool vertex = false;

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

    if (switchtrigger == "trigger_all_electron"){

      if ( (eventdiff->GetRunNumber() >= 146698 && eventdiff->GetRunNumber() <= 148058) && eventCastor->GetHLTPath(10) > 0) triggerE_a = true;
      if ( (eventdiff->GetRunNumber() >= 148822 && eventdiff->GetRunNumber() <= 149063) && eventCastor->GetHLTPath(11) > 0) triggerE_b = true;
      if ( (eventdiff->GetRunNumber() >= 149181 && eventdiff->GetRunNumber() <= 149291) && eventCastor->GetHLTPath(12) > 0) triggerE_c = true;
      if (triggerE_a || triggerE_b || triggerE_c) trigger = true;

      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Electron accepted." << std::endl;
      TriggerStatus = "trigger_all_electron";
    }
    else if (switchtrigger == "trigger_all_muon"){

      if ( (eventdiff->GetRunNumber() >= 146428 && eventdiff->GetRunNumber() <= 147116) && eventCastor->GetHLTPath(3) > 0) triggerE_a = true;
      if ( (eventdiff->GetRunNumber() >= 147196 && eventdiff->GetRunNumber() <= 149291) && eventCastor->GetHLTPath(8) > 0) triggerE_b = true;
      if (triggerE_a || triggerE_b) trigger = true;

      if (debug) std::cout << "\nTrigger Status: " << trigger << ", Trigger All Muon accepted." << std::endl;
      TriggerStatus = "trigger_all_muon";
    }
    else if (switchtrigger == "trigger"){
      if (eventCastor->GetHLTPath(optTrigger) > 0) trigger = true;
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

    //if(eventdiff->GetTightNoiseFilter()>0){
    //if (eventdiff->GetNVertex() <= nVertex) vertex = true;
    //}                                                                                                                                                                             
    if (eventdiff->GetNVertex() <= nVertex) vertex = true;

    if (switchlumiweight =="mc_lumi_weight" || switchlumiweight == "mc_lumi_pu_weight"){
      if(eventinfo->GetNPileUpBx0()<1) zeropileup = true;
    }

    if (typesel == "RecoElectron"){
      selStatus = "Reco::Electron";
      isoTk1 = eventCastor->GetLeadingElectronTkDr03()/eventCastor->GetLeadingElectronPt();
      isoEcal1 = eventCastor->GetLeadingElectronEcalDr03()/eventCastor->GetLeadingElectronPt();
      isoHcal1 = eventCastor->GetLeadingElectronHcalDr03()/eventCastor->GetLeadingElectronPt();
      isoTk2 = eventCastor->GetSecondElectronTkDr03()/eventCastor->GetSecondElectronPt();
      isoEcal2 = eventCastor->GetSecondElectronEcalDr03()/eventCastor->GetSecondElectronPt();
      isoHcal2 = eventCastor->GetSecondElectronHcalDr03()/eventCastor->GetSecondElectronPt();
      innerHits1 = eventCastor->GetLeadingElectronInnerHits();
      Dcot1 = eventCastor->GetLeadingElectronDCot();
      Dist1 = eventCastor->GetLeadingElectronDist();
      DeltaEtaTkClu1 = eventCastor->GetLeadingElectronDeltaEtaTkClu();
      DeltaPhiTkClu1 = eventCastor->GetLeadingElectronDeltaPhiTkClu();
      sigmaIeIe1 = eventCastor->GetLeadingElectronSigmaIeIe();
      HE1 = eventCastor->GetLeadingElectronHE();
      innerHits2 = eventCastor->GetSecondElectronInnerHits();
      Dcot2 = eventCastor->GetSecondElectronDCot();
      Dist2 = eventCastor->GetSecondElectronDist();
      DeltaEtaTkClu2 = eventCastor->GetSecondElectronDeltaEtaTkClu();
      DeltaPhiTkClu2 = eventCastor->GetSecondElectronDeltaPhiTkClu();
      sigmaIeIe2 = eventCastor->GetSecondElectronSigmaIeIe();
      HE2 = eventCastor->GetSecondElectronHE();

      if (eventCastor->GetLeadingElectronPt() > lepton1pt && eventCastor->GetSecondElectronPt() > lepton2pt) presel = true;
      if (eventCastor->GetLeadingElectronCharge()*eventCastor->GetSecondElectronCharge()==-1) charge = true;
      if (eventCastor->GetDiElectronMass() > 60. && eventCastor->GetDiElectronMass() < 110.) dimass = true;
      if (eventCastor->GetElectronsN() > 1) nSel = true;

      //Isolation Electron
      if ((fabs (eventCastor->GetLeadingElectronEta()) <= 1.4442) ){
	if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1 = true;
      }

      if ((fabs (eventCastor->GetLeadingElectronEta()) >= 1.5660) && (fabs (eventCastor->GetLeadingElectronEta()) <= 2.5)){
	if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1 = true;
      }

      if ((fabs (eventCastor->GetSecondElectronEta()) <= 1.4442) ){
	if (isoTk2<0.09 && isoEcal2<0.07 && isoHcal2<0.10) isoBarrel2 = true;
      }

      if ((fabs (eventCastor->GetSecondElectronEta()) >= 1.5660) && (fabs (eventCastor->GetSecondElectronEta()) <= 2.5)){
	if (isoTk2<0.04 && isoEcal2<0.05 && isoHcal2<0.025) isoEndCap2 = true;
      }

      if ((isoEndCap1 || isoBarrel1) && (isoEndCap2 || isoBarrel2)) isolation = true;

      //Electron Candidate Selection
      if ((fabs (eventCastor->GetLeadingElectronEta()) <= 1.4442) ){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.004 && fabs (DeltaPhiTkClu1) < 0.06 && sigmaIeIe1 < 0.01 && HE1 < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (eventCastor->GetLeadingElectronEta()) >= 1.5660) && (fabs (eventCastor->GetLeadingElectronEta()) <= 2.5)){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.007 && fabs (DeltaPhiTkClu1) < 0.03 && sigmaIeIe1 < 0.03 && HE1 < 0.025) eleEndCap1 = true;
      }

      if ((fabs (eventCastor->GetSecondElectronEta()) <= 1.4442) ){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.004 && fabs (DeltaPhiTkClu2) < 0.06 && sigmaIeIe2 < 0.01 && HE2 < 0.04 ) eleBarrel2 = true;
      }

      if ((fabs (eventCastor->GetSecondElectronEta()) >= 1.5660) && (fabs (eventCastor->GetSecondElectronEta()) <= 2.5)){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.007 && fabs (DeltaPhiTkClu2) < 0.03 && sigmaIeIe2 < 0.03 && HE2 < 0.025) eleEndCap2 = true;
      }

      if ((eleEndCap1 || eleBarrel1) && (eleEndCap2 || eleBarrel2)) candSel = true;
    }

    else if (typesel == "RecoMuon"){
      selStatus = "Reco::Muon";
      if (eventCastor->GetLeadingMuonPt() > lepton1pt && eventCastor->GetSecondMuonPt() > lepton2pt) presel = true;
      if (eventCastor->GetLeadingMuonCharge()*eventCastor->GetSecondMuonCharge()==-1) charge = true;
      if (eventCastor->GetDiMuonMass() > 60. && eventCastor->GetDiMuonMass() < 110.) dimass = true;
      if (eventCastor->GetMuonsN() > 1) nSel = true;
      if (eventCastor->GetLeadingMuonSumPtR03() < 3 && eventCastor->GetSecondMuonSumPtR03() < 3 ) { 
	isolation = true;
	candSel = true;
      }

    }

    else if (typesel == "PatElectron"){
      selStatus = "Pat::Electron";
      isoTk1 = eventCastor->GetPatElectron1TkDr03()/eventCastor->GetPatElectron1Pt();
      isoEcal1 = eventCastor->GetPatElectron1EcalDr03()/eventCastor->GetPatElectron1Pt();
      isoHcal1 = eventCastor->GetPatElectron1HcalDr03()/eventCastor->GetPatElectron1Pt();
      isoTk2 = eventCastor->GetPatElectron2TkDr03()/eventCastor->GetPatElectron2Pt();
      isoEcal2 = eventCastor->GetPatElectron2EcalDr03()/eventCastor->GetPatElectron2Pt();
      isoHcal2 = eventCastor->GetPatElectron2HcalDr03()/eventCastor->GetPatElectron2Pt();
      innerHits1 = eventCastor->GetPatElectron1InnerHits();
      Dcot1 = eventCastor->GetPatElectron1DCot();
      Dist1 = eventCastor->GetPatElectron1Dist();
      DeltaEtaTkClu1 = eventCastor->GetPatElectron1DeltaEtaTkClu();
      DeltaPhiTkClu1 = eventCastor->GetPatElectron1DeltaPhiTkClu();
      sigmaIeIe1 = eventCastor->GetPatElectron1SigmaIeIe();
      HE1 = eventCastor->GetPatElectron1HE();
      innerHits2 = eventCastor->GetPatElectron2InnerHits();
      Dcot2 = eventCastor->GetPatElectron2DCot();
      Dist2 = eventCastor->GetPatElectron2Dist();
      DeltaEtaTkClu2 = eventCastor->GetPatElectron2DeltaEtaTkClu();
      DeltaPhiTkClu2 = eventCastor->GetPatElectron2DeltaPhiTkClu();
      sigmaIeIe2 = eventCastor->GetPatElectron2SigmaIeIe();
      HE2 = eventCastor->GetPatElectron2HE();

      if (eventCastor->GetPatElectron1Pt() > lepton1pt && eventCastor->GetPatElectron2Pt() > lepton2pt) presel = true;
      if (eventCastor->GetPatElectron1Charge()*eventCastor->GetPatElectron2Charge()==-1) charge = true;
      if (eventCastor->GetPatDiElectronMass() > 60. && eventCastor->GetPatDiElectronMass() < 110.) dimass = true;
      if (eventCastor->GetPatNElectron() > 1) nSel = true;

      //Isolation Electron
      if ((fabs (eventCastor->GetPatElectron1Eta()) <= 1.4442) ){
	if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1 = true;
      }

      if ((fabs (eventCastor->GetPatElectron1Eta()) >= 1.5660) && (fabs (eventCastor->GetPatElectron1Eta()) <= 2.5)){
	if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1 = true;
      }

      if ((fabs (eventCastor->GetPatElectron2Eta()) <= 1.4442) ){
	if (isoTk2<0.09 && isoEcal2<0.07 && isoHcal2<0.10) isoBarrel2 = true;
      }

      if ((fabs (eventCastor->GetPatElectron2Eta()) >= 1.5660) && (fabs (eventCastor->GetPatElectron2Eta()) <= 2.5)){
	if (isoTk2<0.04 && isoEcal2<0.05 && isoHcal2<0.025) isoEndCap2 = true;
      }

      if ((isoEndCap1 || isoBarrel1) && (isoEndCap2 || isoBarrel2)) isolation = true;

      //Electron Candidate Selection
      if ((fabs (eventCastor->GetPatElectron1Eta()) <= 1.4442) ){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.004 && fabs (DeltaPhiTkClu1) < 0.06 && sigmaIeIe1 < 0.01 && HE1 < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (eventCastor->GetPatElectron1Eta()) >= 1.5660) && (fabs (eventCastor->GetPatElectron1Eta()) <= 2.5)){
	if (innerHits1 == 0 && (fabs (Dcot1) >= 0.02 || fabs (Dist1) >= 0.02 ) && fabs (DeltaEtaTkClu1) < 0.007 && fabs (DeltaPhiTkClu1) < 0.03 && sigmaIeIe1 < 0.03 && HE1 < 0.025) eleEndCap1 = true;
      }

      if ((fabs (eventCastor->GetPatElectron2Eta()) <= 1.4442) ){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.004 && fabs (DeltaPhiTkClu2) < 0.06 && sigmaIeIe2 < 0.01 && HE2 < 0.04 ) eleBarrel2 = true;
      }

      if ((fabs (eventCastor->GetPatElectron2Eta()) >= 1.5660) && (fabs (eventCastor->GetPatElectron2Eta()) <= 2.5)){
	if (innerHits2 == 0 && (fabs (Dcot2) >= 0.02 || fabs (Dist2) >= 0.02 ) && fabs (DeltaEtaTkClu2) < 0.007 && fabs (DeltaPhiTkClu2) < 0.03 && sigmaIeIe2 < 0.03 && HE2 < 0.025) eleEndCap2 = true;
      }

      if ((eleEndCap1 || eleBarrel1) && (eleEndCap2 || eleBarrel2)) candSel = true;
    }

    else if(typesel == "PatMuon"){
      selStatus = "Pat::Muon";
      if (eventCastor->GetPatMuon1Pt() > lepton1pt && eventCastor->GetPatMuon2Pt() > lepton2pt) presel = true;
      if (eventCastor->GetPatMuon1Charge()*eventCastor->GetPatMuon2Charge()==-1) charge = true;
      if (eventCastor->GetPatDiMuonMass() > 60. && eventCastor->GetPatDiMuonMass() < 110.) dimass = true;
      if (eventCastor->GetPatNMuon() > 1) nSel = true; 
      if (eventCastor->GetPatMuon1SumPtR03() < 3 && eventCastor->GetPatMuon2SumPtR03() < 3 ) {
	candSel = true;
	isolation = true;
      }
    }
    else{
      exit(EXIT_FAILURE);
    }

    if(switchtrigger == "trigger" || switchtrigger == "trigger_all_electron" || switchtrigger == "trigger_all_muon"){ 
      FillHistos(0,totalcommon); 
      if(trigger) {
	++totalT;
	FillHistos(1,totalcommon);
      } 
      if(trigger && vertex && presel && nSel && charge && dimass && isolation && candSel) FillHistos(2,totalcommon);
    }

    else if (switchtrigger=="no_trigger_nocorrection" || switchtrigger =="no_trigger_correction"){
      --totalT;
      FillHistos(0,totalcommon);
      if(vertex && presel && nSel && charge && dimass && isolation && candSel) FillHistos(2,totalcommon);
      if(zeropileup && vertex) FillHistos(3,totalcommon);
    }
    else{
      exit(EXIT_FAILURE);
    }
  }

  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << ">> Input file: " << filein << std::endl;
  outstring << ">> Output file: " << savehistofile << std::endl;
  outstring << ">> Channel Correction file: " << channelcorrfile << std::endl;
  outstring << ">> TTree Name: " << processname << std::endl;
  outstring << " " << std::endl;
  outstring << "<< OPTIONS >>" << std::endl; 
  outstring << " " << std::endl;
  outstring << ">> Trigger: " << TriggerStatus << std::endl;
  outstring << ">> # Vertex: " << nVertex << std::endl;
  outstring << ">> Lepton1(pT) > " << lepton1pt <<std::endl;
  outstring << ">> Lepton2(pT) > " << lepton2pt <<std::endl;
  outstring << ">> Type of Selection: " << selStatus << std::endl;
  outstring << ">> Event Weight: " << switchlumiweight << " | MC Weight:" << mcweight << std::endl;
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
  SaveHistos();
  outf->Close();
  std::cout << "\n" << std::endl;
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
  std::string typesel_;
  std::string pumfile_;
  std::string pudfile_;
  double mcweight_;
  std::string switchlumiweight_;
  double SectorThreshold_;
  double ChannelThreshold_;
  std::string channelcorrfile_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0) filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0) processname_ = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0) savehistofile_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0) switchtrigger_ = argv[4];
  if (argc > 5 && strcmp(s1,argv[5]) != 0) optTrigger_   = atoi(argv[5]);
  if (argc > 6 && strcmp(s1,argv[6]) != 0) lepton1pt_ = atof(argv[6]);
  if (argc > 7 && strcmp(s1,argv[7]) != 0) lepton2pt_ = atof(argv[7]);
  if (argc > 8 && strcmp(s1,argv[8]) != 0) nVertex_ = atoi(argv[8]);
  if (argc > 9 && strcmp(s1,argv[9]) != 0) switchlumiweight_ = argv[9];
  if (argc > 10 && strcmp(s1,argv[10]) != 0) mcweight_ = atof(argv[10]);
  if (argc > 11 && strcmp(s1,argv[11]) != 0) typesel_ = argv[11];
  if (argc > 12 && strcmp(s1,argv[12]) != 0) SectorThreshold_ = atof(argv[12]);
  if (argc > 13 && strcmp(s1,argv[13]) != 0) ChannelThreshold_ = atof(argv[13]);
  if (argc > 14 && strcmp(s1,argv[14]) != 0) channelcorrfile_ = argv[14];
  if (argc > 15 && strcmp(s1,argv[15]) != 0) pumfile_ = argv[15];
  if (argc > 16 && strcmp(s1,argv[16]) != 0) pudfile_ = argv[16];

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
  std::cout << "Type of Selection: " << typesel_ << std::endl;
  std::cout << "Luminosity Weight Option: " << switchlumiweight_ << std::endl;
  std::cout << "MC Weight: " << mcweight_ << std::endl;
  std::cout << "Sector Castor Threshold [GeV]: " << SectorThreshold_ <<  std::endl;
  std::cout << "Channel Castor Threshold [GeV]: " << ChannelThreshold_ <<  std::endl;
  std::cout << "Channel Correction File: " << channelcorrfile_ <<  std::endl;
  std::cout << "" << std::endl;

  if (switchtrigger_=="trigger" || switchtrigger_=="no_trigger_nocorrection" || switchtrigger_=="no_trigger_correction" || switchtrigger_=="trigger_all_electron" || switchtrigger_ == "trigger_all_muon") {}
  else{
    std::cout << "Please Insert type of Swithtrigger: " << std::endl;
    std::cout << "1) trigger: run with trigger. Need optTrigger >=0;" << std::endl;
    std::cout << "2) trigger_all_electron: all trigger electron path will be accepted. Do not require optTrigger." << std::endl;
    std::cout << "3) trigger_all_muon: all trigger muon path will be accepted. Do not require optTrigger." << std::endl;
    std::cout << "4) no_trigger: run without trigger." << std::endl;
    return 0;
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

  if (nVertex_ <= 0 || optTrigger_ < 0 || lepton1pt_ < 0 || lepton2pt_ < 0 || mcweight_ <= 0 ){
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << " Pay attention on the input numbers parameters" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << ">> Requirements:                             " << std::endl;
    std::cout << "I)   nVertex_ > 0 " << std::endl;
    std::cout << "II)  optTrigger >= 0" << std::endl;
    std::cout << "III) SectorThreshold_ > 0" << std::endl;
    std::cout << "IV)  ChannelThreshold_ >= 0" << std::endl;  
    std::cout << "V )  Lepton1pt_ and Lepton2pt_ >= 0" << std::endl;
    std::cout << "VI)  mcweight_ > 0" << std::endl;
    std::cout << "----------------------------------------------" << std::endl; 
    return 0;
  }

  if ( (SectorThreshold_ < 0 && ChannelThreshold_ < 0) || (SectorThreshold_ > 0 && ChannelThreshold_ > 0) ){
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "   Pay attention on the input numbers parameters" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << ">> Requirements:                             " << std::endl;
    std::cout << "If castorthreshold_ and channelsthreshold are both" << std::endl;
    std::cout << "positive or negative, you can not run program. "           << std::endl;
    std::cout << "" << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;
    std::cout << ">> Instructions: " << std::endl;
    std::cout << "I) SectorThreshold_ > 0 and ChannelThreshold_ < 0, allow Sector threshold." << std::endl;
    std::cout << "I) SectorThreshold_ < 0 and ChannelThreshold_ > 0, allow Channel threshold." << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;
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

  clock_t tStart = clock();


  CastorAnalysis* CastorRun = new CastorAnalysis();
  CastorRun->CreateHistos();
  CastorRun->Run(filein_, processname_, savehistofile_, switchtrigger_, optTrigger_, lepton1pt_, lepton2pt_, nVertex_, typesel_, switchlumiweight_, mcweight_, SectorThreshold_, ChannelThreshold_, channelcorrfile_, pumfile_, pudfile_);
  std::cout<< "Time taken: " << (double)(clock() - tStart)/(60*CLOCKS_PER_SEC) << " min" << std::endl; 
  return 0;
}

#endif
