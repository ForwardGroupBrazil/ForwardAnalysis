//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsDiffractiveWsAnalysis#Macro_Analysis
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

void DiffractiveW::CreateHistos(std::string type){

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
  std::string step10 = "NGapCMSAndZKinP";
  std::string step11 = "PGapCMSAndZKinN";
  std::string step12 = "NGapCMSAndCASTORAndZKinP";
  std::string step13 = "PGapCMSAndCastorActivityAndZKinN";
  std::string step14 = "NGapCASTOR";
  std::string step15 = "NGapCASTORAndZKinP";
  std::string step16 = "PGapCMSAndCASTOR";
  std::string step17 = "PGapCMSAndCASTORAndZKinP";

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

    m_hVector_WMuonMass.push_back( std::vector<TH1F*>() );
    m_hVector_WMuonEta.push_back( std::vector<TH1F*>() );
    m_hVector_WMuonPt.push_back( std::vector<TH1F*>() );
    m_hVector_WMuonPhi.push_back( std::vector<TH1F*>() );
    m_hVector_WElectronMass.push_back( std::vector<TH1F*>() );
    m_hVector_WElectronEta.push_back( std::vector<TH1F*>() );
    m_hVector_WElectronPt.push_back( std::vector<TH1F*>() );
    m_hVector_WElectronPhi.push_back( std::vector<TH1F*>() );

    for (int k=0;k<nloop;k++){

      if (type=="multiple_pileup"){
	sprintf(tag,"multiple_pileup_%i",k);
      }
      else{
	sprintf(tag,"single");
      }


      char name1[300];
      sprintf(name1,"WMuonMass_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WMuonMass = new TH1F(name1,"Boson W Invariant Mass Distribution; M_{#mu#nu} [GeV]; N events",500,0,500);
      m_hVector_WMuonMass[j].push_back(histo_WMuonMass);

      char name2[300];
      sprintf(name2,"WMuonPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WMuonPt = new TH1F(name2,"Boson W Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_WMuonPt[j].push_back(histo_WMuonPt);

      char name3[300];
      sprintf(name3,"WMuonEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WMuonEta = new TH1F(name3,"Boson W #eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_WMuonEta[j].push_back(histo_WMuonEta);

      char name4[300];
      sprintf(name4,"WMuonPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WMuonPhi = new TH1F(name4,"Boson W #phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_WMuonPhi[j].push_back(histo_WMuonPhi);

      char name5[300];
      sprintf(name5,"WElectronMass_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WElectronMass = new TH1F(name5,"Boson W Invariant Mass Distribution; M_{#mu#nu} [GeV]; N events",500,0,500);
      m_hVector_WElectronMass[j].push_back(histo_WElectronMass);

      char name6[300];
      sprintf(name6,"WElectronPt_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WElectronPt = new TH1F(name6,"Boson W Pt Distribution; P_{T} [GeV.c^{-1}]; N events",200,0,1000);
      m_hVector_WElectronPt[j].push_back(histo_WElectronPt);

      char name7[300];
      sprintf(name7,"WElectronEta_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WElectronEta = new TH1F(name7,"Boson W #eta Distribution; #eta; N events",50,-5.2,5.2);
      m_hVector_WElectronEta[j].push_back(histo_WElectronEta);

      char name8[300];
      sprintf(name8,"WElectronPhi_%s_%s",tag,Folders.at(j).c_str());
      TH1F *histo_WElectronPhi = new TH1F(name8,"Boson W #phi Distribution; #phi [rad]; N events",60,-3.3,3.3);
      m_hVector_WElectronPhi[j].push_back(histo_WElectronPhi);

    }
  }
}

void DiffractiveW::FillHistos(int index, int pileup, double totalweight){

  m_hVector_WMuonMass[index].at(pileup)->Fill(bosonWMass,totalweight);
  m_hVector_WMuonEta[index].at(pileup)->Fill(eventdiffW->GetLeadingMuonEta(),totalweight);
  m_hVector_WMuonPt[index].at(pileup)->Fill(eventdiffW->GetLeadingMuonPt(),totalweight);
  m_hVector_WMuonPhi[index].at(pileup)->Fill(eventdiffW->GetLeadingMuonPhi(),totalweight);
  m_hVector_WElectronMass[index].at(pileup)->Fill(bosonWMass,totalweight);
  m_hVector_WElectronEta[index].at(pileup)->Fill(eventdiffW->GetLeadingElectronEta(),totalweight);
  m_hVector_WElectronPt[index].at(pileup)->Fill(eventdiffW->GetLeadingElectronPt(),totalweight);
  m_hVector_WElectronPhi[index].at(pileup)->Fill(eventdiffW->GetLeadingElectronPhi(),totalweight);

}

void DiffractiveW::SaveHistos(std::string type,std::string typesel){

  // Creating Correlation Histograms

  int ipileup;

  if (type=="multiple_pileup") ipileup=21;
  else ipileup=1;

  for (int i = 0; i < ipileup; i++){
    for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

      if(typesel=="RecoMuon" || typesel=="RecoElectron"){
	// Lepton Kinematics Folder
	foldersFile[0]->cd();
	m_hVector_WMuonMass[j].at(i)->Write();
	m_hVector_WMuonEta[j].at(i)->Write();
	m_hVector_WMuonPt[j].at(i)->Write();
	m_hVector_WMuonPhi[j].at(i)->Write();
	m_hVector_WElectronMass[j].at(i)->Write();
	m_hVector_WElectronEta[j].at(i)->Write();
	m_hVector_WElectronPt[j].at(i)->Write();
	m_hVector_WElectronPhi[j].at(i)->Write();
      }

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
      if(eventdiffW->GetHLTPath(nt) > 0){
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
    bool ZKinN = false;
    bool ZKinP = false;

    double bosonWMass = -999.;
    if (typesel == "RecoMuon"){
      bosonWMass = TMath::Sqrt(2.*(eventdiffW->GetLeadingMuonPt()*eventdiffW->GetMETPt() - (eventdiffW->GetLeadingMuonPt()*TMath::Cos(eventdiffW->GetLeadingMuonPhi())*eventdiffW->GetMETPt()*TMath::Cos(eventdiffW->GetMETPhi()) + eventdiffW->GetLeadingMuonPt()*TMath::Sin(eventdiffW->GetLeadingMuonPhi())*eventdiffW->GetMETPt()*TMath::Sin(eventdiffW->GetMETPhi()) ) ) );
    }
    if (typesel == "PatMuon"){
      bosonWMass = TMath::Sqrt(2.*(eventdiffW->GetPatMuon1Pt()*eventdiffW->GetPatMETPt() - (eventdiffW->GetPatMuon1Pt()*TMath::Cos(eventdiffW->GetPatMuon1Phi())*eventdiffW->GetPatMETPt()*TMath::Cos(eventdiffW->GetPatMETPhi()) + eventdiffW->GetPatMuon1Pt()*TMath::Sin(eventdiffW->GetPatMuon1Phi())*eventdiffW->GetPatMETPt()*TMath::Sin(eventdiffW->GetPatMETPhi()) ) ) );
    }
    if (typesel == "RecoElectron"){
      bosonWMass = TMath::Sqrt(2.*(eventdiffW->GetLeadingElectronPt()*eventdiffW->GetMETPt() - (eventdiffW->GetLeadingElectronPt()*TMath::Cos(eventdiffW->GetLeadingElectronPhi())*eventdiffW->GetMETPt()*TMath::Cos(eventdiffW->GetMETPhi()) + eventdiffW->GetLeadingElectronPt()*TMath::Sin(eventdiffW->GetLeadingElectronPhi())*eventdiffW->GetMETPt()*TMath::Sin(eventdiffW->GetMETPhi()) ) ) );
    }
    if (typesel == "PatElectron"){
      bosonWMass = TMath::Sqrt(2.*(eventdiffW->GetPatElectron1Pt()*eventdiffW->GetPatMETPt() - (eventdiffW->GetPatElectron1Pt()*TMath::Cos(eventdiffW->GetPatElectron1Phi())*eventdiffW->GetPatMETPt()*TMath::Cos(eventdiffW->GetPatMETPhi()) + eventdiffW->GetPatElectron1Pt()*TMath::Sin(eventdiffW->GetPatElectron1Phi())*eventdiffW->GetPatMETPt()*TMath::Sin(eventdiffW->GetPatMETPhi()) ) ) );
    }

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

    if (typesel == "RecoElectron"){
      selStatus = "Reco::Electron";
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
      else{
	AEcastor = 999.;
      }

      if (eventdiffW->GetLeadingElectronPt() > lepton1pt && eventdiffW->GetMETPt() > lepton2pt) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;

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

      if (eventdiffW->GetLeadingElectronEta()>0.) ZKinP = true;
      if (eventdiffW->GetLeadingElectronEta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetLeadingElectronEta()); }
      else{ etasignedHF = fabs(eventdiffW->GetLeadingElectronEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetLeadingElectronEta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetLeadingElectronEta());}

    }

    else if (typesel == "RecoMuon"){
      selStatus = "Reco::Muon";
      if (eventdiffW->GetLeadingMuonPt() > lepton1pt && eventdiffW->GetMETPt() > lepton2pt) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;
      if (eventdiffW->GetLeadingMuonSumPtR03() < 3) { 
	isolation = true;
	candSel = true;
      }

      if (eventdiffW->GetLeadingMuonEta()>0.) ZKinP = true;
      if (eventdiffW->GetLeadingMuonEta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetLeadingMuonEta()); }
      else{ etasignedHF = fabs(eventdiffW->GetLeadingMuonEta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetLeadingMuonEta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetLeadingMuonEta());}

    }

    else if (typesel == "PatElectron"){
      selStatus = "Pat::Electron";
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

      if (eventdiffW->GetPatElectron1Pt() > lepton1pt && eventdiffW->GetPatMETPt() > lepton2pt) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;

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

      if (eventdiffW->GetPatElectron1Eta()>0.) ZKinP = true;
      if (eventdiffW->GetPatElectron1Eta()<0.) ZKinN = true;

      if (diffselp || diffseln){ etasignedHF = -1*fabs(eventdiffW->GetPatElectron1Eta()); }
      else{ etasignedHF = fabs(eventdiffW->GetPatElectron1Eta()); }

      if (castorgap){ etasignedCASTOR = -1*fabs(eventdiffW->GetPatElectron1Eta());}
      else {etasignedCASTOR = fabs(eventdiffW->GetPatElectron1Eta());}

    }

    else if(typesel == "PatMuon"){
      selStatus = "Pat::Muon";
      if (eventdiffW->GetPatMuon1Pt() > lepton1pt && eventdiffW->GetPatMETPt() > lepton2pt) presel = true;
      if (bosonWMass > 60. && bosonWMass < 110.) dimass = true;
      if (eventdiffW->GetPatMuon1SumPtR03() < 3) {
	candSel = true;
	isolation = true;
      }

      if (eventdiffW->GetPatMuon1Eta()>0.) ZKinP = true;
      if (eventdiffW->GetPatMuon1Eta()<0.) ZKinN = true;

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
    bMaxGapPF = eventdiffW->GetMaxGapPF();
    bPTMaxGapMaxPF = eventdiffW->GetPTMaxGapMaxPF();
    bPTMinGapMaxPF = eventdiffW->GetPTMinGapMaxPF();
    bXiPlusFromPFCands = eventdiff->GetXiPlusFromPFCands();
    bXiMinusFromPFCands = eventdiff->GetXiMinusFromPFCands();
    betasignedHF = etasignedHF;
    betasignedCASTOR = etasignedCASTOR;
    bAEcastor = AEcastor;
    betamax = eventdiff->GetEtaMaxFromPFCands();
    betamin = etamin_;
    betalimmin = eventdiffW->GetLimMinusGapPF();
    betalimmax = eventdiffW->GetLimPlusGapPF();

    if(pileup < 21){ // Never comment this line. It is the program defense.

      if(switchtrigger == "trigger" || switchtrigger == "trigger_all_electron" || switchtrigger == "trigger_all_muon"){ 
	FillHistos(0,pileup,totalcommon); 
	if(trigger) {
	  ++totalT;
	  FillHistos(1,pileup,totalcommon);
	} 
	if(trigger && vertex && presel) FillHistos(2,pileup,totalcommon);
	if(trigger && vertex && presel && dimass) FillHistos(3,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation) FillHistos(4,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel) {
	  FillHistos(5,pileup,totalcommon);
	  fOutZ->cd();
	  troutZ->Fill();
	}
	if(trigger && vertex && presel && dimass && isolation && candSel && diffseln) FillHistos(6,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && diffselp) FillHistos(7,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && diffseln && castorgap) FillHistos(8,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && diffselp && castoractivity) FillHistos(9,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && diffseln && ZKinP){
	  outstring << "HF- Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(10,pileup,totalcommon);
	}
	if(trigger && vertex && presel && dimass && isolation && candSel && diffselp && ZKinN){
	  outstring << "HF+ Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(11,pileup,totalcommon);
	}
	if(trigger && vertex && presel && dimass && isolation && candSel && diffseln && castorgap && ZKinP) FillHistos(12,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && diffselp && castoractivity && ZKinN) FillHistos(13,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && castorgap){
	  FillHistos(14,pileup,totalcommon);
	  fOutCASTOR->cd();
	  troutCASTOR->Fill();
	}
	if(trigger && vertex && presel && dimass && isolation && candSel && castorgap && ZKinP){
	  outstring << "CASTOR Gap, W Candidate: " << eventdiff->GetRunNumber() << ":" << eventdiff->GetLumiSection() << ":" << eventdiff->GetEventNumber() << std::endl;
	  FillHistos(15,pileup,totalcommon);
	  fOut->cd();
	  trout->Fill();
	}

	if(trigger && vertex && presel && dimass && isolation && candSel && diffselp && castorgap) FillHistos(16,pileup,totalcommon);
	if(trigger && vertex && presel && dimass && isolation && candSel && diffselp && castorgap && ZKinP) FillHistos(17,pileup,totalcommon);

      }

      else if (switchtrigger =="no_trigger_nocorrection" || switchtrigger == "no_trigger_correction" ){
	--totalT;
	FillHistos(0,pileup,totalcommon);
	if(vertex && presel) FillHistos(2,pileup,totalcommon);
	if(vertex && presel && dimass) FillHistos(3,pileup,totalcommon);
	if(vertex && presel && dimass && isolation) FillHistos(4,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel) {
	  FillHistos(5,pileup,totalcommon);
	  fOutZ->cd();
	  troutZ->Fill();
	}
	if(vertex && presel && dimass && isolation && candSel && diffseln) FillHistos(6,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffselp) FillHistos(7,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffseln && castorgap) FillHistos(8,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffselp && castoractivity) FillHistos(9,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffseln && ZKinP) FillHistos(10,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffselp && ZKinN) FillHistos(11,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffseln && castorgap && ZKinP) FillHistos(12,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffselp && castoractivity && ZKinN) FillHistos(13,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && castorgap){ 
	  FillHistos(14,pileup,totalcommon);
	  fOutCASTOR->cd();
	  troutCASTOR->Fill();
	}
	if(vertex && presel && dimass && isolation && candSel && castorgap && ZKinP){
	  FillHistos(15,pileup,totalcommon);
	  fOut->cd();
	  trout->Fill();
	}
	if(vertex && presel && dimass && isolation && candSel && diffselp && castorgap) FillHistos(16,pileup,totalcommon);
	if(vertex && presel && dimass && isolation && candSel && diffselp && castorgap && ZKinP) FillHistos(17,pileup,totalcommon);
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

    DiffractiveW* diffWRun = new DiffractiveW();
    diffWRun->CreateHistos(type_);
    diffWRun->Run(filein_, processname_, savehistofile_, switchtrigger_, optTrigger_, lepton1pt_, lepton2pt_, nVertex_, type_, switchlumiweight_, mcweight_, typesel_, castorthreshold_, channelsthreshold_, castorcorrfile_, gapseltype_, pumfile_, pudfile_);
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
