//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
// Project: Exclusive Dijets Analysis
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsExclusiveDijetsAnalysis
//

//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TMath.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "EffMacro.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ExclusiveDijetsEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace diffractiveAnalysis;
using namespace exclusiveDijetsAnalysis;
using namespace eventInfo;
using namespace reweight;

void EffMacro::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  tr = (TTree*)inf->Get(processinput.c_str());
  eventdiff = new DiffractiveEvent();
  eventexcl = new ExclusiveDijetsEvent();
  eventinfo = new EventInfoEvent();
  diff = tr->GetBranch("DiffractiveAnalysis");
  excl = tr->GetBranch("ExclusiveDijetsAnalysis");
  info = tr->GetBranch("EventInfo");
  diff->SetAddress(&eventdiff);
  excl->SetAddress(&eventexcl);
  info->SetAddress(&eventinfo);

}

void EffMacro::Run(std::string filein_, std::string savehistofile_, std::string processname_, int optnVertex_, int optTrigger_, bool switchPreSel_, bool switchVertex_, bool switchTrigger_, double channelsthreshold_){

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  filein = filein_;
  optnVertex = optnVertex_;
  optTrigger = optTrigger_;
  switchPreSel = switchPreSel_;
  switchVertex = switchVertex_;
  switchTrigger = switchTrigger_;
  channelsthreshold = channelsthreshold_;

  std::cout << "" << std::endl;
  std::cout << "Running..." << std::endl;
  std::cout << "" << std::endl;
  std::cout << "<< INPUTS >>" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Input file: " << filein << std::endl;
  std::cout << " " << std::cout;
  std::cout << "Output file: " << savehistofile << std::endl;
  std::cout << " " << std::cout; 
  std::cout << "# Vertex: " << optnVertex << std::endl;
  std::cout << "Trigger Option: " << optTrigger << std::endl;
  std::cout << " " << std::endl;
  std::cout << "--> TRUE = 1 FALSE = 0" << std::endl;
  std::cout << "Vertex Switch: " << switchVertex << std::endl;
  std::cout << "Trigger Switch: " << switchTrigger << std::endl;
  std::cout << "Pre-Selection Switch: " << switchPreSel << std::endl;
  std::cout << "CASTOR Threshold: " << channelsthreshold << std::endl;
  std::cout << " " << std::endl;


  // Code Protection
  if (optnVertex == 0){

    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "Please, restart your setup. Respect the Condition # Vertex > 0)" << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
    return;

  }

  TFile check1(filein.c_str());

  if (check1.IsZombie()){

    std::cout << "----------------------------------------------" << std::endl;
    std::cout << " There is no PatTuple file or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << " Edit the source and recompile." << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    return;

  }
  //--------------------------------------------------------------------------------------------------------------------------


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

  // Root file with histograms
  TFile *outf = new TFile(savehistofile.c_str(),"RECREATE");

  // File with Number of Events
  TString outtxt = savehistofile;
  outtxt.ReplaceAll("root","txt");  
  std::ofstream outstring(outtxt); 

  int NEVENTS = tr->GetEntries();

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  unsigned NEntries = tr->GetEntries();
  std::cout << "" << std::endl;
  std::cout<< "Reading Tree: "<< NEntries << " events"<<std::endl;
  std::cout << "" << std::endl;

  int decade = 0;

  int TotalE = 0;
  int counterTrigger = 0;
  int counterPreSel=0;
  int counterVertex = 0;
  int counterVertex_castorgap = 0;
  int counterAllstep4_4 = 0;
  int counterAllstep4_3 = 0;
  int counterAllstep4_2 = 0;
  int counterAllstep4_1 = 0;
  int counterAllstep4_4_castorgap = 0;
  int counterAllstep4_3_castorgap = 0;
  int counterAllstep4_2_castorgap = 0;
  int counterAllstep4_1_castorgap = 0;

  std::vector <std::string> Folders;
  Folders.push_back("without_cuts");
  Folders.push_back("with_trigger");
  Folders.push_back("with_trigger_presel");
  Folders.push_back("with_trigger_presel_vertex");
  Folders.push_back("All_step4_4");
  Folders.push_back("All_step4_3");
  Folders.push_back("All_step4_2");
  Folders.push_back("All_step4_1");
  Folders.push_back("with_trigger_presel_vertex_castorgap");
  Folders.push_back("All_step4_4_castorgap");
  Folders.push_back("All_step4_3_castorgap");
  Folders.push_back("All_step4_2_castorgap");
  Folders.push_back("All_step4_1_castorgap");

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++)
  {

    char name1[300];
    sprintf(name1,"Events_%s",Folders.at(j).c_str());
    TH1D *histo_m_Evt_lumis = new TH1D(name1,"; Lumis; N events",100,0,2.0);
    m_hVector_Evt_lumis.push_back(histo_m_Evt_lumis);

    char name2[300];
    sprintf(name2,"Eff_%s",Folders.at(j).c_str());
    TH1D *histo_m_Eff_lumis = new TH1D(name2,"; Lumis; Efficiency",100,0,2.0);
    m_hVector_Eff_lumis.push_back(histo_m_Eff_lumis);

  }

  for(int i=0;i<NEVENTS;i++) {

    ++TotalE;

    double progress = 10.0*i/(1.0*NEVENTS);
    int l = TMath::FloorNint(progress); 
    if (l > decade){
      std::cout <<"\n<<<<<< STATUS >>>>>>" << std::endl; 
      std::cout<<10*l<<" % completed." << std::endl;
      std::cout <<"<<<<<<<<<<>>>>>>>>>>\n" << std::endl;
    }
    decade = l;          

    tr->GetEntry(i);

    if( i % 1000 == 0 ){
      std::cout << "\nEvent " << i << std::endl;
    }

    bool trigger = false;
    bool presel = false;
    bool castorgap = false;
    bool castoractivity = false;
    bool vertex = false;
    bool eta4 = false;
    bool eta3 = false;
    bool eta2 = false;
    bool eta1 = false;
    int SectorCastorHit = 0;

    for (int i=0; i < 16; i++){
      CastorEnergySector[i]=0.;
      if (i==4 || i==5){
	if (eventdiff->GetCastorModule2Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule2Energy(i);
	if (eventdiff->GetCastorModule3Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule3Energy(i);
	if (eventdiff->GetCastorModule4Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule4Energy(i);
	if (eventdiff->GetCastorModule5Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule5Energy(i);
      }else{
	if (eventdiff->GetCastorModule1Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule1Energy(i);
	if (eventdiff->GetCastorModule2Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule2Energy(i);
	if (eventdiff->GetCastorModule3Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule3Energy(i);
	if (eventdiff->GetCastorModule4Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule4Energy(i);
	if (eventdiff->GetCastorModule5Energy(i) > channelsthreshold) CastorEnergySector[i]+=eventdiff->GetCastorModule5Energy(i);
      }
    }

    for (l=0; l<16;l++){
      if (CastorEnergySector[l] >= channelsthreshold){
	++SectorCastorHit;
      }
    }

    if (SectorCastorHit >= 1){
      castoractivity = true;
    }else{
      castorgap = true;
    }

    if (!switchTrigger || (switchTrigger && eventexcl->GetHLTPath(optTrigger))) trigger = true;
    if (!switchPreSel || (switchPreSel && ( (eventdiff->GetSumEnergyHFMinus() < 30 && eventdiff->GetSumEnergyHFPlus() < 30) || (eventdiff->GetEtaMinFromPFCands() < -990 && eventdiff->GetEtaMaxFromPFCands() < -990 && castorgap) ))) presel = true;
    if (!switchVertex || (switchVertex && eventexcl->GetNVertex()<= optnVertex)) vertex = true;
    if ((eventdiff->GetEtaMinFromPFCands() > -4. && eventdiff->GetEtaMaxFromPFCands() < 4.) || (eventdiff->GetEtaMinFromPFCands() < -990 && eventdiff->GetEtaMaxFromPFCands() < -990 && castorgap)) eta4 = true;
    if ((eventdiff->GetEtaMinFromPFCands() > -3. && eventdiff->GetEtaMaxFromPFCands() < 3.) || (eventdiff->GetEtaMinFromPFCands() < -990 && eventdiff->GetEtaMaxFromPFCands() < -990 && castorgap)) eta3 = true;
    if ((eventdiff->GetEtaMinFromPFCands() > -2. && eventdiff->GetEtaMaxFromPFCands() < 2.) || (eventdiff->GetEtaMinFromPFCands() < -990 && eventdiff->GetEtaMaxFromPFCands() < -990 && castorgap)) eta2 = true;
    if ((eventdiff->GetEtaMinFromPFCands() > -1. && eventdiff->GetEtaMaxFromPFCands() < 1.) || (eventdiff->GetEtaMinFromPFCands() < -990 && eventdiff->GetEtaMaxFromPFCands() < -990 && castorgap)) eta1 = true;

    m_hVector_Evt_lumis.at(0)->Fill(eventinfo->GetInstLumiBunch());
    m_hVector_Eff_lumis.at(0)->Fill(eventinfo->GetInstLumiBunch());

    if(trigger){
      ++counterTrigger;     
      m_hVector_Evt_lumis.at(1)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(1)->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel){
      ++counterPreSel;
      m_hVector_Evt_lumis.at(2)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(2)->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel && vertex){
      ++counterVertex;
      m_hVector_Evt_lumis.at(3)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(3)->Fill(eventinfo->GetInstLumiBunch());
    }
    
    if(trigger && presel && vertex && eta4){
      ++counterAllstep4_4;
      m_hVector_Evt_lumis.at(4)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(4)->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel && vertex && eta3){
      ++counterAllstep4_3;
      m_hVector_Evt_lumis[5]->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis[5]->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel && vertex && eta2){
      ++counterAllstep4_2;
      m_hVector_Evt_lumis.at(6)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(6)->Fill(eventinfo->GetInstLumiBunch());
    }
    
    if(trigger && presel && vertex && eta1){
      ++counterAllstep4_1;
      m_hVector_Evt_lumis.at(7)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(7)->Fill(eventinfo->GetInstLumiBunch());
    }
    
    if(trigger && presel && vertex && castorgap){
      ++counterVertex_castorgap;
      m_hVector_Evt_lumis.at(8)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(8)->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel && vertex && eta4 && castorgap){
      ++counterAllstep4_4_castorgap;
      m_hVector_Evt_lumis.at(9)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(9)->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel && vertex && eta3 && castorgap){
      ++counterAllstep4_3_castorgap;
      m_hVector_Evt_lumis[10]->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis[10]->Fill(eventinfo->GetInstLumiBunch());
    }

    if(trigger && presel && vertex && eta2 && castorgap){
      ++counterAllstep4_2_castorgap;
      m_hVector_Evt_lumis.at(11)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(11)->Fill(eventinfo->GetInstLumiBunch());
    }
    
    if(trigger && presel && vertex && eta1 && castorgap){
      ++counterAllstep4_1_castorgap;
      m_hVector_Evt_lumis.at(12)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(12)->Fill(eventinfo->GetInstLumiBunch());
    }
    
  }

  //Scalling Plots
  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){
    m_hVector_Eff_lumis.at(j)->Scale(1./TotalE);
  }

  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << "Input file: " << filein << std::endl;
  outstring << "Output file: " << savehistofile << std::endl;
  outstring << " " << std::endl;
  outstring << "# Vertex: " << optnVertex << std::endl;
  outstring << "Trigger Option: " << optTrigger << std::endl;
  outstring << " " << std::endl;
  outstring << "--> TRUE = 1 FALSE = 0" << std::endl;
  outstring << "Trigger Switch: " << switchTrigger << std::endl;
  outstring << "Vertex  Switch: " << switchVertex << std::endl;
  outstring << "Pre-Selection Switch: " << switchPreSel << std::endl;
  outstring << "" << std::endl;
  outstring << "<< EVENT INFO >> " << std::endl;
  outstring << " " << std::endl;
  outstring << "Total # of Events without Weight: " << TotalE << std::endl;
  outstring << " " << std::endl;
  outstring << "Trigger: " << counterTrigger << std::endl;
  outstring << "Trigger + Pre Sel.: " << counterPreSel << std::endl;
  outstring << "Trigger + Pre Sel. + Vertex: " << counterVertex << std::endl;
  outstring << "STEP4_4 (CMS Acceptance): " << counterAllstep4_4 << std::endl;
  outstring << "STEP4_3 (CMS Acceptance): " << counterAllstep4_3 << std::endl;
  outstring << "STEP4_2 (CMS Acceptance): " << counterAllstep4_2 << std::endl;
  outstring << "STEP4_1 (CMS Acceptance): " << counterAllstep4_1 << std::endl;
  outstring << "STEP4_4_Castor (CMS Acceptance): " << counterAllstep4_4_castorgap << std::endl;
  outstring << "STEP4_3_Castor (CMS Acceptance): " << counterAllstep4_3_castorgap << std::endl;
  outstring << "STEP4_2_Castor (CMS Acceptance): " << counterAllstep4_2_castorgap << std::endl;
  outstring << "STEP4_1_Castor (CMS Acceptance): " << counterAllstep4_1_castorgap << std::endl;
  outstring << " " << std::endl;
  outstring << "<< LEGEND >> " << std::endl;
  outstring << "STEP4_X: Trigger + Pre Sel. + Vertex + Eta_max < X and Eta_min > X." << std::endl;
  outstring << " " << std::endl;

  outf->Write();
  outf->Close();
  outstring.close();

}

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
int main(int argc, char **argv)
{
  AutoLibraryLoader::enable();

  const char *s1="*";
  std::string filein_;
  std::string savehistofile_;
  std::string processname_;
  int optnVertex_;
  int optTrigger_;
  bool switchPreSel_;
  bool switchVertex_;
  bool switchTrigger_;
  double channelsthreshold_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0)  filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0)  savehistofile_  = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0)  processname_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0)  optnVertex_ = atoi(argv[4]);
  if (argc > 5 && strcmp(s1,argv[5]) != 0)  optTrigger_   = atoi(argv[5]);
  if (argc > 6 && strcmp(s1,argv[6]) != 0)  switchPreSel_   = atoi(argv[6]);
  if (argc > 7 && strcmp(s1,argv[7]) != 0)  switchVertex_   = atoi(argv[7]);
  if (argc > 8 && strcmp(s1,argv[8]) != 0)  switchTrigger_   = atoi(argv[8]);
  if (argc > 9 && strcmp(s1,argv[9]) != 0)  channelsthreshold_   = atof(argv[9]);

 if (channelsthreshold_ < 0){
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << " Pay attention on the input numbers parameters" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << ">> Requirements: " << std::endl;
      std::cout << ">> channelthreshold must be >= 0 GeV" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      return 0;
    }

  EffMacro* exclDijets = new EffMacro();   
  exclDijets->Run(filein_, savehistofile_, processname_, optnVertex_, optTrigger_, switchPreSel_, switchVertex_, switchTrigger_, channelsthreshold_);



  return 0;
}
#endif
