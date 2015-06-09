//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Autor: Diego Figueiredo 
// dmf@cern.ch

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
#include "TriggerEfficiency.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/TriggerEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace triggerAnalysis;
using namespace eventInfo;
using namespace reweight;

void TriggerEfficiency::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventtrigger = new TriggerEvent();
  eventinfo = new EventInfoEvent();
  trigger = tr->GetBranch("TriggerAnalysis");
  info = tr->GetBranch("EventInfo");
  trigger->SetAddress(&eventtrigger);
  info->SetAddress(&eventinfo);

}

void TriggerEfficiency::Run(std::string filein_, std::string savehistofile_, std::string processname_, std::string typesel_, int optTriggerRef_, int optTrigger_){

  bool debug = false;

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  typesel = typesel_;
  optTriggerRef = optTriggerRef_;
  optTrigger = optTrigger_;

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

  std::cout << "" << std::endl;
  std::cout<< "Reading Tree: "<< NEVENTS << " events"<<std::endl;
  std::cout << "" << std::endl;

  TH1D *histo_checktrigger = new TH1D("checktrigger","Check Trigger; Trigger; N events",21,0,21);
  TH1D *histo_checkcuts = new TH1D("checkcuts","Check Cuts; Cut Order; N events",19,0,19);

  std::vector <std::string> Folders;
  Folders.push_back("without_cuts");
  Folders.push_back("with_RefTrigger");
  Folders.push_back("with_Trigger");
  Folders.push_back("step1");
  Folders.push_back("step2");
  Folders.push_back("step3");

  int triggercounter[20]={0};
  double counter[6] = {0};

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name[300];

    sprintf(name,"leading_lepton_pt_%s",Folders.at(j).c_str());
    TH1D *histo_leading_pt = new TH1D(name,"; P_{T} [GeV.c^{-1}]; N events",500,0,1000);
    m_hVector_leading_pt.push_back(histo_leading_pt);

    sprintf(name,"second_lepton_pt_%s",Folders.at(j).c_str());
    TH1D *histo_second_pt = new TH1D(name,"; P_{T} [GeV.c^{-1}]; N events",500,0,1000);
    m_hVector_second_pt.push_back(histo_second_pt);

    sprintf(name,"dilepton_mass_%s",Folders.at(j).c_str());
    TH1D *histo_dilepton_mass = new TH1D(name,"; M_{ll} [GeV]; N events",500,0,500);
    m_hVector_dilepton_mass.push_back(histo_dilepton_mass);

  }

  for(int i=0;i<NEVENTS;i++) {

    tr->GetEntry(i);

    if (!debug){
      if (i==0) {
	std::cout << "" << std::endl;
	std::cout<< "Status Bar" << std::endl;
	std::cout << "" << std::endl;
      }
      loadBar(i,NEVENTS,100,100);
    }

    for (int i=0;i<20;i++){
      if(eventtrigger->GetHLTPath(i)>0){
	if(debug) std::cout << "Trigger: " << eventtrigger->GetHLTPath(i) << std::endl;
	histo_checktrigger->Fill(i,eventtrigger->GetHLTPath(i));
	triggercounter[i]++;
      }
    }

    math::XYZTLorentzVector DiSystem(0.,0.,0.,0.);
    double leadingleptonpt = -999.;
    double secondleptonpt = -999.;
    double ptcut = 999.;
    bool tagcheck = false;

    if(typesel=="Electron" || typesel=="ELECTRON"){
      leadingleptonpt = eventtrigger->GetLeadingElectronPt();
      secondleptonpt = eventtrigger->GetSecondElectronPt();
      DiSystem += eventtrigger->GetLeadingElectronP4();
      DiSystem += eventtrigger->GetSecondElectronP4();
      ptcut = 25.;
      tagcheck = eventtrigger->GetLeadingElectronIsWP80();
    }
    else if(typesel=="Muon" || typesel=="MUON"){
      leadingleptonpt = eventtrigger->GetLeadingMuonPt();
      secondleptonpt = eventtrigger->GetSecondMuonPt();
      DiSystem += eventtrigger->GetLeadingMuonP4();
      DiSystem += eventtrigger->GetSecondMuonP4();
      ptcut = 20.;
      tagcheck = eventtrigger->GetLeadingMuonIsGood(); 
    }else{
      std::cout << "Please, put type of particle!" << std::endl;
      exit(EXIT_FAILURE);
    }

    bool refTrigger = false;
    bool Trigger = false;
    bool presel = false;
    bool dimass = false;

    if(eventtrigger->GetHLTPath(optTriggerRef)>0) refTrigger = true;
    if(eventtrigger->GetHLTPath(optTrigger)>0) Trigger = true;
    if(leadingleptonpt > ptcut) presel = true;
    if(secondleptonpt > ptcut && DiSystem.M() > 60. && DiSystem.M() < 110.) dimass = true;

    counter[0]++;
    m_hVector_leading_pt.at(0)->Fill(leadingleptonpt);
    m_hVector_second_pt.at(0)->Fill(secondleptonpt);
    m_hVector_dilepton_mass.at(0)->Fill(DiSystem.M());

    if(refTrigger){
      counter[1]++;
      m_hVector_leading_pt.at(1)->Fill(leadingleptonpt);
      m_hVector_second_pt.at(1)->Fill(secondleptonpt);
      m_hVector_dilepton_mass.at(1)->Fill(DiSystem.M());
    }

    if(refTrigger && Trigger){
      counter[2]++;
      m_hVector_leading_pt.at(2)->Fill(leadingleptonpt);
      m_hVector_second_pt.at(2)->Fill(secondleptonpt);
      m_hVector_dilepton_mass.at(2)->Fill(DiSystem.M());
    }

    if(refTrigger && Trigger && presel){
      counter[3]++;
      m_hVector_leading_pt.at(3)->Fill(leadingleptonpt);
      m_hVector_second_pt.at(3)->Fill(secondleptonpt);
      m_hVector_dilepton_mass.at(3)->Fill(DiSystem.M());
    }

    if(refTrigger && Trigger && presel && tagcheck){
      counter[4]++;
      m_hVector_leading_pt.at(4)->Fill(leadingleptonpt);
      m_hVector_second_pt.at(4)->Fill(secondleptonpt);
      m_hVector_dilepton_mass.at(4)->Fill(DiSystem.M());
    }

    if(refTrigger && Trigger && presel && tagcheck && dimass){
      counter[5]++;
      m_hVector_leading_pt.at(5)->Fill(leadingleptonpt);
      m_hVector_second_pt.at(5)->Fill(secondleptonpt);
      m_hVector_dilepton_mass.at(5)->Fill(DiSystem.M());
    }

  }

  //Scalling Plots
  histo_checkcuts->Fill("No Cuts",1);
  histo_checkcuts->Fill("Ref. Trigger",1);
  histo_checkcuts->Fill("Trigger",1);
  histo_checkcuts->Fill("P_{T} cut",1);
  histo_checkcuts->Fill("Quality",1);
  histo_checkcuts->Fill("Dimass",1);

  for(std::vector<std::string>::size_type j=0; j<Folders.size(); j++){
    histo_checkcuts->SetBinContent(j+1,counter[j]);
  }

  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << "Input file: " << filein << std::endl;
  outstring << "Output file: " << savehistofile << std::endl;
  outstring << " " << std::endl;
  outstring << "Trigger Ref Option: " << optTriggerRef << std::endl;
  outstring << "Trigger Option: " << optTrigger << std::endl;
  outstring << " " << std::endl;
  outstring << "Number of Events Without Cuts: " << counter[0] << std::endl;
  outstring << "Number of Events With Ref. Trigger: " << counter[1] << std::endl;
  outstring << "" << std::endl;
  outstring << "Total Trigger Fired: " <<  std::endl;
  for(int j=0; j<20; j++){
    outstring << "Trigger " << j << ": " << triggercounter[j] << std::endl;
  }
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
  std::string typesel_;
  int optTrigger_;
  int optTriggerRef_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0)  filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0)  savehistofile_  = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0)  processname_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0)  typesel_ = argv[4];
  if (argc > 5 && strcmp(s1,argv[5]) != 0)  optTriggerRef_   = atoi(argv[5]);
  if (argc > 6 && strcmp(s1,argv[6]) != 0)  optTrigger_   = atoi(argv[6]);

  std::cout << "" << std::endl;
  std::cout << "Running..." << std::endl;
  std::cout << "" << std::endl;
  std::cout << "<< INPUTS >>" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Input file: " << filein_ << std::endl;
  std::cout << "Output file: " << savehistofile_ << std::endl;
  std::cout << "Type Selection: " << typesel_ << std::endl;
  std::cout << "Reference Trigger Option: " << optTriggerRef_ << std::endl;
  std::cout << "Trigger Option: " << optTrigger_ << std::endl;
  std::cout << " " << std::endl;

  if (typesel_=="Muon" || typesel_=="Electron" || typesel_=="ELECTRON" || typesel_=="MUON") {}
  else{
    std::cout << "Please Insert type of Selections: " << std::endl;
    std::cout << "1) Muon/MUON: selections with Reco::Muon." << std::endl;
    std::cout << "2) Electron/ELECTRON: selections with Reco::Electron." << std::endl;
    return 0;
  }

  TriggerEfficiency* triggereff = new TriggerEfficiency();   
  triggereff->Run(filein_, savehistofile_, processname_, typesel_, optTriggerRef_, optTrigger_);

  return 0;

}
#endif

