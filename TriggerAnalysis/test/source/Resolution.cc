//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
// Project: Exclusive Dijets Analysis
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsDiffractiveWAnalysis
//
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

#include "statusbar.h"
#include "Resolution.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveWEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace diffractiveAnalysis;
using namespace diffractiveWAnalysis;
using namespace eventInfo;
using namespace reweight;

void Resolution::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventdiff = new DiffractiveEvent();
  eventBoson = new DiffractiveWEvent();
  eventinfo = new EventInfoEvent();
  diff = tr->GetBranch("DiffractiveAnalysis");
  boson = tr->GetBranch("DiffractiveWAnalysis");
  info = tr->GetBranch("EventInfo");
  diff->SetAddress(&eventdiff);
  boson->SetAddress(&eventBoson);
  info->SetAddress(&eventinfo);

}

void Resolution::Run(std::string filein_, std::string savehistofile_, std::string processname_, int optTriggerRef_, int optTriggerRefOR_, int optTrigger_, int optTriggerOR_, int bin_, double channelsthreshold_){

  bool debug = false;

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  optTriggerRef = optTriggerRef_;
  optTriggerRefOR = optTriggerRefOR_;
  optTrigger = optTrigger_;
  optTriggerOR = optTriggerOR_;
  bin = bin_;
  channelsthreshold = channelsthreshold_;

  TFile check1(filein.c_str());

  if (check1.IsZombie()){

    std::cout << "----------------------------------------------" << std::endl;
    std::cout << " There is no PatTuple file or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << " Edit the source and recompile." << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    return;

  }

  if (check1.GetDirectory(processname.c_str())){
    LoadFile(filein,processname);
  }else{
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << " There is no directory/path " << processname << std::endl;
    std::cout << " in the file." << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    return;
  }

  TFile *outf = new TFile(savehistofile.c_str(),"RECREATE");

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

  int triggercounter[20]={0};
  double TotalE = 0.;
  double counter[18] = {0};
  double deltaphi_ = 0.;

  TH1D *histo_checktrigger = new TH1D("checktrigger","Check Trigger; Trigger; N events",21,0,21);
  TH1D *histo_checkcuts = new TH1D("checkcuts","Check Cuts; Cut Order; N events",19,0,19);

  std::vector <std::string> Folders;
  Folders.push_back("without_cuts");
  Folders.push_back("with_RefTrigger");
  Folders.push_back("W_selection");

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name1[300];
    sprintf(name1,"Events_%s",Folders.at(j).c_str());
    TH1D *histo_m_Evt_lumis = new TH1D(name1,"; Lumis; N events",bin,0,2.0);
    m_hVector_Evt_lumis.push_back(histo_m_Evt_lumis);

    char name2[300];
    sprintf(name2,"Eff_%s",Folders.at(j).c_str());
    TH1D *histo_m_Eff_lumis = new TH1D(name2,"; Lumis; Efficiency",bin,0,2.0);
    m_hVector_Eff_lumis.push_back(histo_m_Eff_lumis);

    char name3[300];
    sprintf(name3,"Events_PFEtamax_%s",Folders.at(j).c_str());
    TH1D *histo_m_Evt_PFEtamax = new TH1D(name3,";Particle Flow #eta_{max}; N events",bin,0,5.5);
    m_hVector_Evt_pfetamax.push_back(histo_m_Evt_PFEtamax);

    char name4[300];
    sprintf(name4,"Events_PFEtamin_%s",Folders.at(j).c_str());
    TH1D *histo_m_Evt_PFEtamin = new TH1D(name4,";Particle Flow #eta_{min}; N events",bin,-5.5,0);
    m_hVector_Evt_pfetamin.push_back(histo_m_Evt_PFEtamin);

  }

  for(int i=0;i<NEVENTS;i++) {

    tr->GetEntry(i);

    ++TotalE;

    if (!debug){
      if (i==0) {
	std::cout << "" << std::endl;
	std::cout<< "Status Bar" << std::endl;
	std::cout << "" << std::endl;
      }
      loadBar(i,NEVENTS,100,100);
    }

    for (int nt=0;nt<20;nt++){
      if(eventBoson->GetHLTPath(nt)>0){
	if(debug) std::cout << "Trigger: " << eventBoson->GetHLTPath(nt) << std::endl;
	histo_checktrigger->Fill(nt,eventBoson->GetHLTPath(nt));
	triggercounter[nt]++;
      }
    }

    counter[0]++;

    m_hVector_Evt_lumis.at(0)->Fill(eventinfo->GetInstLumiBunch());
    m_hVector_Eff_lumis.at(0)->Fill(eventinfo->GetInstLumiBunch());
    m_hVector_Evt_pfetamax.at(0)->Fill(eventdiff->GetEtaMaxFromPFCands());
    m_hVector_Evt_pfetamin.at(0)->Fill(eventdiff->GetEtaMinFromPFCands());

    if(eventBoson->GetHLTPath(optTriggerRef)>0 || eventBoson->GetHLTPath(optTriggerRefOR)>0 ){

      counter[1]++;
      m_hVector_Evt_lumis.at(1)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Eff_lumis.at(1)->Fill(eventinfo->GetInstLumiBunch());
      m_hVector_Evt_pfetamax.at(1)->Fill(eventdiff->GetEtaMaxFromPFCands());
      m_hVector_Evt_pfetamin.at(1)->Fill(eventdiff->GetEtaMinFromPFCands());

    }
  }

  //Scalling Plots
  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){
    m_hVector_Eff_lumis.at(j)->Scale(1./counter[j]);
  }

  histo_checkcuts->Fill("No Cuts",1);
  histo_checkcuts->Fill("Ref. Trigger",1);
  histo_checkcuts->Fill("W Boson Selection",1);

  for(int i=0; i<18; i++){
    histo_checkcuts->SetBinContent(i+1,counter[i]);
  }

  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << "Input file: " << filein << std::endl;
  outstring << "Output file: " << savehistofile << std::endl;
  outstring << " " << std::endl;
  outstring << "Trigger Ref Option: " << optTriggerRef << std::endl;
  outstring << "Trigger Ref Option OR: " << optTriggerRefOR << std::endl;
  outstring << "Trigger Option: " << optTrigger << std::endl;
  outstring << "Trigger Option OR: " << optTriggerOR << std::endl;
  outstring << "Bin: " << bin << std::endl;
  outstring << "CASTOR Threshold: " << channelsthreshold << " GeV" << std::endl; 
  outstring << " " << std::endl;
  outstring << "Number of Events: " << TotalE << std::endl;
  outstring << "Number of Events Without Cuts: " << counter[0] << std::endl;
  outstring << "Number of Events With Ref. Trigger: " << counter[1] << std::endl;
  outstring << "Number of W Boson Events: " << counter[2] << std::endl;
  outstring << "" << std::endl;
  outstring << "Total Trigger Fired: " <<  std::endl;
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
  int optTrigger_;
  int optTriggerOR_;
  int optTriggerRef_;
  int optTriggerRefOR_;
  int bin_;
  double channelsthreshold_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0)  filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0)  savehistofile_  = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0)  processname_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0)  optTriggerRef_   = atoi(argv[4]);
  if (argc > 5 && strcmp(s1,argv[5]) != 0)  optTriggerRefOR_   = atoi(argv[5]);
  if (argc > 6 && strcmp(s1,argv[6]) != 0)  optTrigger_   = atoi(argv[6]);
  if (argc > 7 && strcmp(s1,argv[7]) != 0)  optTriggerOR_   = atoi(argv[7]);
  if (argc > 8 && strcmp(s1,argv[8]) != 0)  bin_   = atoi(argv[8]);
  if (argc > 9 && strcmp(s1,argv[9]) != 0)  channelsthreshold_ = atof(argv[9]);

  std::cout << "" << std::endl;
  std::cout << "Running..." << std::endl;
  std::cout << "" << std::endl;
  std::cout << "<< INPUTS >>" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Input file: " << filein_ << std::endl;
  std::cout << "Output file: " << savehistofile_ << std::endl;
  std::cout << "Reference Trigger Option: " << optTriggerRef_ << std::endl;
  std::cout << "Reference Trigger Option OR: " << optTriggerRefOR_ << std::endl;
  std::cout << "Trigger Option: " << optTrigger_ << std::endl;
  std::cout << "Trigger Option OR: " << optTriggerOR_ << std::endl;
  std::cout << "Bin: " << bin_ << std::endl;
  std::cout << " " << std::endl;

  if (channelsthreshold_ < 0){
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << " Pay attention on the input numbers parameters" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << ">> Requirements: " << std::endl;
    std::cout << ">> channelthreshold must be >= 0 GeV" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    return 0;
  }

  Resolution* triggereff = new Resolution();   
  triggereff->Run(filein_, savehistofile_, processname_, optTriggerRef_, optTriggerRefOR_, optTrigger_, optTriggerOR_, bin_, channelsthreshold_);

  return 0;

}
#endif

