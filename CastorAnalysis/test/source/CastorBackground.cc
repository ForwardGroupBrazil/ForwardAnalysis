//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/FwdPhysicsExclusiveDijetsAnalysis#Example_Analysis_Macro
//

//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom2.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "statusbar.h"
#include "CastorBackground.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ZeroBiasEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveWEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace zerobiasAnalysis;
using namespace diffractiveWAnalysis;
using namespace diffractiveAnalysis;
using namespace eventInfo;
using namespace reweight;

void CastorBackground::LoadFile(std::string fileinput, std::string processinput, std::string fileinputrandom, std::string processinputrandom){

  inf = NULL;
  tr  = NULL;

  infrandom = NULL;
  trrandom = NULL;

  std::string fdirectory = processinput + "/ProcessedTree";
  std::string fdirectoryrandom = processinputrandom + "/ProcessedTree";

  inf = TFile::Open(fileinput.c_str(),"read");
  tr = (TTree*)inf->Get(fdirectory.c_str());

  infrandom = TFile::Open(fileinputrandom.c_str(),"read");
  trrandom = (TTree*)infrandom->Get(fdirectoryrandom.c_str());

  eventCastor = new DiffractiveWEvent();
  eventInfo = new EventInfoEvent();
  eventRandom = new ZeroBiasEvent();
  eventDiff = new DiffractiveEvent();

  Castor = tr->GetBranch("DiffractiveWAnalysis");
  info = tr->GetBranch("EventInfo");
  Diff = tr->GetBranch("DiffractiveAnalysis");
  random = trrandom->GetBranch("ZeroBiasAnalysis");

  Castor->SetAddress(&eventCastor);
  Diff->SetAddress(&eventDiff);
  info->SetAddress(&eventInfo);
  random->SetAddress(&eventRandom);

}

void CastorBackground::CreateHistos(std::string type){

  std::string step0 = "without_cuts_" + type;
  std::string step1 = "step1_" + type;
  std::string step2 = "step2_" + type;
  std::string step3 = "step3_" + type;

  if(type == "mixture"){
    Folders.push_back(step0);
    Folders.push_back(step1);
    Folders.push_back(step2);
    Folders.push_back(step3);
  }
  if(type=="mc"){
    Folders.push_back(step0);
    Folders.push_back(step1);
    Folders.push_back(step3);
  }
  if(type=="data"){
    Folders.push_back(step0);
    Folders.push_back(step3);
  }

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name[300];
    sprintf(name,"Tracks_%s",Folders.at(j).c_str());
    TH1F *histo_Tracks = new TH1F(name,"Tracks Multiplicity; n Tracks; N events",50,0,150);
    m_hVector_tracks.push_back(histo_Tracks);

    sprintf(name,"vertex_%s",Folders.at(j).c_str());
    TH1F *histo_vertex = new TH1F(name,"Number of Vertex; # Vertex; N events",25,0,25);
    m_hVector_vertex.push_back(histo_vertex);

    sprintf(name,"sumECastorMinusBySector_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorMinusBySector = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_TotalCastorEnergyBySector.push_back(histo_sumECastorMinusBySector);

    sprintf(name,"sumECastorMinusBySector20PercUp_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorMinusBySector20PercUp = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_TotalCastorEnergyBySector20PercUp.push_back(histo_sumECastorMinusBySector20PercUp);

    sprintf(name,"sumECastorMinusBySector20PercDow_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorMinusBySector20PercDow = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_TotalCastorEnergyBySector20PercDow.push_back(histo_sumECastorMinusBySector20PercDow);

    sprintf(name,"sumECastorMinusBySector10PercUp_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorMinusBySector10PercUp = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_TotalCastorEnergyBySector10PercUp.push_back(histo_sumECastorMinusBySector10PercUp);

    sprintf(name,"sumECastorMinusBySector10PercDow_%s",Folders.at(j).c_str());
    TH1F *histo_sumECastorMinusBySector10PercDow = new TH1F(name,"Castor Sum of Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_TotalCastorEnergyBySector10PercDow.push_back(histo_sumECastorMinusBySector10PercDow);

    for (int bs=0; bs<16; bs++){
      m_hVector_SectorCastorEnergy.push_back( std::vector<TH1F*>() );

      char namesector[300];
      char namesectort[300];
      sprintf(namesector,"Sector%d_CastorSumEnergy_%s",bs+1,Folders.at(j).c_str());
      sprintf(namesectort,"#sum Energy, Castor Sector %d; #sum E_{modules 1,2,3,4,5} [GeV]; N events", bs+1);
      TH1F *histo_SectorCastorEnergy = new TH1F(namesector,namesectort,2000,0,1000);
      m_hVector_SectorCastorEnergy[bs].push_back(histo_SectorCastorEnergy);

    }
  }

}

void CastorBackground::FillHistos(int index){

  m_hVector_tracks[index]->Fill(eventDiff->GetMultiplicityTracks());
  m_hVector_vertex[index]->Fill(eventDiff->GetNVertex());
  m_hVector_TotalCastorEnergyBySector[index]->Fill(sumCastorEnergy);
  m_hVector_TotalCastorEnergyBySector20PercUp[index]->Fill(sumCastorEnergy20Up);
  m_hVector_TotalCastorEnergyBySector20PercDow[index]->Fill(sumCastorEnergy20Dw);
  m_hVector_TotalCastorEnergyBySector10PercUp[index]->Fill(sumCastorEnergy10Up);
  m_hVector_TotalCastorEnergyBySector10PercDow[index]->Fill(sumCastorEnergy10Dw);

  for (int i=0; i < 16; i++){
    m_hVector_SectorCastorEnergy[i].at(index)->Fill(CastorEnergySector[i]);
  }                                                                                                                                                        

}

void CastorBackground::SaveHistos(){

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){
    m_hVector_tracks[j]->Write();
    m_hVector_vertex[j]->Write();
    m_hVector_TotalCastorEnergyBySector[j]->Write();
    m_hVector_TotalCastorEnergyBySector20PercUp[j]->Write();
    m_hVector_TotalCastorEnergyBySector20PercDow[j]->Write();
    m_hVector_TotalCastorEnergyBySector10PercUp[j]->Write();
    m_hVector_TotalCastorEnergyBySector10PercDow[j]->Write();
    for (int id=0; id<16; id++){
      m_hVector_SectorCastorEnergy[id].at(j)->Write();
    }
  }

}

void CastorBackground::Run(std::string filein_, std::string processname_, std::string fileinrandom_, std::string processnamerandom_, std::string savehistofile_, std::string type_, double castorthreshold_){

  bool debug = false;
  sumCastorEnergy = 0.;                                                                                                                                      
  sumCastorEnergy20Up = 0.;
  sumCastorEnergy20Dw = 0.;
  sumCastorEnergy10Up = 0.;
  sumCastorEnergy10Dw = 0.;

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  filein = filein_;
  fileinrandom = fileinrandom_;
  savehistofile = savehistofile_;
  processname = processname_;
  processnamerandom = processnamerandom_;
  type = type_;
  castorthreshold = castorthreshold_;

  TFile check1(filein.c_str());
  TFile check2(fileinrandom.c_str());
  if (check1.IsZombie()){
    std::cout << "\n----------------------------------------------" << std::endl;
    std::cout << " There is no the file " << filein << " or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;
  }
  if (check2.IsZombie()){
    std::cout << "\n----------------------------------------------" << std::endl;
    std::cout << " There is no the file " << fileinrandom << " or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;
  }

  if (check1.GetDirectory(processname.c_str()) && check2.GetDirectory(processnamerandom.c_str())){
    LoadFile(filein, processname, fileinrandom, processnamerandom);
  }
  else {
    std::cout << "\n-------------------------------------------------" << std::endl;
    std::cout << " There is no directory/path " << processname << std::endl;
    std::cout << " or " << processnamerandom << " in the file." << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    return;
  }

  TFile *outf = new TFile(savehistofile.c_str(),"RECREATE");
  TString outtxt = savehistofile;
  outtxt.ReplaceAll("root","txt");
  std::ofstream outstring(outtxt);

  int NEVENTS = tr->GetEntries();
  int NEVENTSRnd = trrandom->GetEntries();

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
  std::string status;

  TRandom r;
  int e1 = r.Uniform(0,NEVENTS);
  //int sec_count = 0;

  for(int i=0;i<NEVENTS;i++){

    tr->GetEntry(i);
    //++sec_count;

    //if(sec_count > NEVENTSRnd) break;
    trrandom->GetEntry(e1);

    if(debug){
      std::cout << "\n Random Castor Tower Energy: " << eventRandom->GetCastorTowerEnergy(0) << std::endl; 
      std::cout << "\n Castor Tower Energy: " << eventCastor->GetCastorTowerEnergy(0) << " | Event: " << i << std::endl; 
    }

    if (!debug){
      if (i==0) {
	std::cout << "" << std::endl;
	std::cout<< "Status Bar" << std::endl;
	std::cout << "" << std::endl;
      }
      loadBar(i,NEVENTS);
    }

    bool zeropileup = false;
    bool randomNoVtx = false;
    bool vtx = false;

    if(type=="mc" || type == "mixture"){
      if(eventInfo->GetNPileUpBx0()<1) zeropileup = true;
    }
    if(eventRandom->GetNVertex()<1) randomNoVtx = true;
    if(eventDiff->GetNVertex()==1) vtx = true;

    if (type == "mixture"){

      status = "mixture";

      for (int i=0; i < 16; i++){

	CastorEnergySector[i]=0;
	CastorEnergySector20Up[i]=0;
	CastorEnergySector20Dw[i]=0;
	CastorEnergySector10Up[i]=0;
	CastorEnergySector10Dw[i]=0;

	if (i==4 || i==5){
	  CastorEnergySector[i]=eventCastor->GetCastorModule2Energy(i) + eventRandom->GetCastorModule2Energy(i) + eventCastor->GetCastorModule3Energy(i) + eventRandom->GetCastorModule3Energy(i) + eventCastor->GetCastorModule4Energy(i) + eventRandom->GetCastorModule4Energy(i) + eventCastor->GetCastorModule5Energy(i) + eventRandom->GetCastorModule5Energy(i);
	  CastorEnergySector20Up[i]=CastorEnergySector[i]+0.2*CastorEnergySector[i];
	  CastorEnergySector20Dw[i]=CastorEnergySector[i]-0.2*CastorEnergySector[i];
	  CastorEnergySector10Up[i]=CastorEnergySector[i]+0.1*CastorEnergySector[i];
	  CastorEnergySector10Dw[i]=CastorEnergySector[i]-0.1*CastorEnergySector[i];
	}else{
	  CastorEnergySector[i]=eventCastor->GetCastorModule1Energy(i) + eventRandom->GetCastorModule1Energy(i) + eventCastor->GetCastorModule2Energy(i) + eventRandom->GetCastorModule2Energy(i) + eventCastor->GetCastorModule3Energy(i) + eventRandom->GetCastorModule3Energy(i) + eventCastor->GetCastorModule4Energy(i) + eventRandom->GetCastorModule4Energy(i) + eventCastor->GetCastorModule5Energy(i) + eventRandom->GetCastorModule5Energy(i);
	  CastorEnergySector20Up[i]=CastorEnergySector[i]+0.2*CastorEnergySector[i];
	  CastorEnergySector20Dw[i]=CastorEnergySector[i]-0.2*CastorEnergySector[i];
	  CastorEnergySector10Up[i]=CastorEnergySector[i]+0.1*CastorEnergySector[i];
	  CastorEnergySector10Dw[i]=CastorEnergySector[i]-0.1*CastorEnergySector[i]; 
	}

	if (CastorEnergySector[i] >= castorthreshold){
	  sumCastorEnergy+=CastorEnergySector[i];
	}

	if (CastorEnergySector20Up[i] >= castorthreshold){
	  sumCastorEnergy20Up+=CastorEnergySector20Up[i];
	}

	if (CastorEnergySector20Dw[i] >= castorthreshold){
	  sumCastorEnergy20Dw+=CastorEnergySector20Dw[i];
	}

	if (CastorEnergySector10Up[i] >= castorthreshold){
	  sumCastorEnergy10Up+=CastorEnergySector10Up[i];
	}

	if (CastorEnergySector10Dw[i] >= castorthreshold){
	  sumCastorEnergy10Dw+=CastorEnergySector10Dw[i];
	}

      }

      if (debug){
	for (int i=0; i<16; i++){
	  std::cout << "\nType: " << type << std::endl;  
	  std::cout << "Castor Sector(" << i+1 << "): " << CastorEnergySector[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), more 20%: " << CastorEnergySector20Up[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), less 20%: " << CastorEnergySector20Dw[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), more 10%: " << CastorEnergySector10Up[i] << std::endl;                                               
	  std::cout << "Castor Sector(" << i+1 << "), less 10%: " << CastorEnergySector10Dw[i] << std::endl;
	}
      }

      FillHistos(0);
      if (zeropileup) FillHistos(1);
      if (zeropileup && randomNoVtx) FillHistos(2); 
      if (zeropileup && randomNoVtx && vtx) FillHistos(3);

    }else{

      for (int i=0; i < 16; i++){

	CastorEnergySector[i]=0;
	CastorEnergySector20Up[i]=0;
	CastorEnergySector20Dw[i]=0;
	CastorEnergySector10Up[i]=0;
	CastorEnergySector10Dw[i]=0;

	if (i==4 || i==5){
	  CastorEnergySector[i]=eventCastor->GetCastorModule2Energy(i) + eventCastor->GetCastorModule3Energy(i) + eventCastor->GetCastorModule4Energy(i) + eventCastor->GetCastorModule5Energy(i);
	  CastorEnergySector20Up[i]=CastorEnergySector[i]+0.2*CastorEnergySector[i];
	  CastorEnergySector20Dw[i]=CastorEnergySector[i]-0.2*CastorEnergySector[i];
	  CastorEnergySector10Up[i]=CastorEnergySector[i]+0.1*CastorEnergySector[i];
	  CastorEnergySector10Dw[i]=CastorEnergySector[i]-0.1*CastorEnergySector[i];
	}else{
	  CastorEnergySector[i]=eventCastor->GetCastorModule1Energy(i) + eventCastor->GetCastorModule2Energy(i) + eventCastor->GetCastorModule3Energy(i) + eventCastor->GetCastorModule4Energy(i) + eventCastor->GetCastorModule5Energy(i);
	  CastorEnergySector20Up[i]=CastorEnergySector[i]+0.2*CastorEnergySector[i];
	  CastorEnergySector20Dw[i]=CastorEnergySector[i]-0.2*CastorEnergySector[i];
	  CastorEnergySector10Up[i]=CastorEnergySector[i]+0.1*CastorEnergySector[i];
	  CastorEnergySector10Dw[i]=CastorEnergySector[i]-0.1*CastorEnergySector[i]; 
	}

	if (CastorEnergySector[i] >= castorthreshold){
	  sumCastorEnergy+=CastorEnergySector[i];
	}

	if (CastorEnergySector20Up[i] >= castorthreshold){
	  sumCastorEnergy20Up+=CastorEnergySector20Up[i];
	}

	if (CastorEnergySector20Dw[i] >= castorthreshold){
	  sumCastorEnergy20Dw+=CastorEnergySector20Dw[i];
	}

	if (CastorEnergySector10Up[i] >= castorthreshold){
	  sumCastorEnergy10Up+=CastorEnergySector10Up[i];
	}

	if (CastorEnergySector10Dw[i] >= castorthreshold){
	  sumCastorEnergy10Dw+=CastorEnergySector10Dw[i];
	}

      }

      if (debug){
	for (int i=0; i<16; i++){
	  std::cout << "\nType: " << type << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "): " << CastorEnergySector[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), more 20%: " << CastorEnergySector20Up[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), less 20%: " << CastorEnergySector20Dw[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), more 10%: " << CastorEnergySector10Up[i] << std::endl;
	  std::cout << "Castor Sector(" << i+1 << "), less 10%: " << CastorEnergySector10Dw[i] << std::endl;
	}
      }

      if (type == "mc"){
	status = "mc";
	FillHistos(0);
	if (zeropileup) FillHistos(1);
	if (vtx) FillHistos(2);
      }
      else if (type == "data"){
	status = "data";
	FillHistos(0);
	if (vtx) FillHistos(1);
      }
      else {
	std::cout << "\n Unrecognized Type of Selection." << std::endl;
	return;
      }

    }

  }   

  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << "Input file: " << filein << std::endl;
  outstring << "Output file: " << savehistofile << std::endl;
  outstring << "Sector Threshold: " << castorthreshold << std::endl;
  outstring << " " << std::endl;
  outstring << "Type: " << status << std::endl;
  outstring << "" << std::endl;

  SaveHistos();
  outf->Close();

}

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
int main(int argc, char **argv)
{
  AutoLibraryLoader::enable();

  const char *s1="*";
  std::string filein_;
  std::string fileinrandom_;
  std::string savehistofile_;
  std::string processname_;
  std::string processnamerandom_;
  std::string type_;
  double castorthreshold_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0) filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0) processname_  = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0) fileinrandom_ = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0) processnamerandom_  = argv[4];
  if (argc > 5 && strcmp(s1,argv[5]) != 0) savehistofile_  = argv[5];
  if (argc > 6 && strcmp(s1,argv[6]) != 0) type_  = argv[6];
  if (argc > 7 && strcmp(s1,argv[7]) != 0) castorthreshold_ = atof(argv[7]);

  if (type_=="mc" || type_=="mixture" || type_=="data") {}
  else{
    std::cout << "Please Insert type of selection: " << std::endl;
    std::cout << "1) mc: only Monte Carlo." << std::endl;
    std::cout << "2) mixture: random ZeroBias event + MC(zeropileup)." << std::endl;
    std::cout << "3) data: only data." << std::endl;
    return 0;
  }

  if ( castorthreshold_ < 0){
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "   Pay attention on the input numbers parameters" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    return 0;
  }

  CastorBackground* castor = new CastorBackground();   
  castor->CreateHistos(type_);
  castor->Run(filein_, processname_, fileinrandom_, processnamerandom_, savehistofile_, type_, castorthreshold_);

  std::cout << "\n" << std::endl;

  return 0;
}
#endif
