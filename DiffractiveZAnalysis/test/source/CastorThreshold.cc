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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "statusbar.h"
#include "CastorThreshold.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ExclusiveDijetsEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveZEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace diffractiveAnalysis;
using namespace diffractiveZAnalysis;
using namespace exclusiveDijetsAnalysis;
using namespace eventInfo;
using namespace reweight;

#define PI 3.14159265

void CastorThreshold::LoadFile(std::string fileinput, std::string processinput){

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


void CastorThreshold::CreateHistos(std::string type){

  std::string step0 = "without_cuts_" + type;
  std::string step1 = "with_type_" + type;  

  Folders.push_back(step0);
  Folders.push_back(step1);

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name1[300];
    sprintf(name1,"Tracks_%s",Folders.at(j).c_str());
    TH1F *histo_Tracks = new TH1F(name1,"Tracks Multiplicity; n Tracks; N events",50,0,150);
    m_hVector_tracks.push_back(histo_Tracks);

    char name2[300];
    sprintf(name2,"vertex_%s",Folders.at(j).c_str());
    TH1F *histo_vertex = new TH1F(name2,"Number of Vertex; # Vertex; N events",25,0,25);
    m_hVector_vertex.push_back(histo_vertex);

    char name3[300];
    sprintf(name3,"RunNumber_%s",Folders.at(j).c_str());
    TH1F *histo_RunNumber = new TH1F(name3,"Run Number; Run Number; N Event",16000,134000,150000);
    m_hVector_RunNumber.push_back(histo_RunNumber);

    for (int bs=0; bs<16; bs++){
      m_hVector_SectorCastorEnergy.push_back( std::vector<TH1F*>() );

      char name4[300];
      char name4t[300];
      sprintf(name4,"Sector%d_CastorSumEnergy_%s",bs+1,Folders.at(j).c_str());
      sprintf(name4t,"#sum Energy, Castor Sector %d; #sum E_{modules 1,2,3,4,5} [GeV]; N events", bs+1);
      TH1F *histo_SectorCastorEnergy = new TH1F(name4,name4t,200,-10,10);
      m_hVector_SectorCastorEnergy[bs].push_back(histo_SectorCastorEnergy);

    }

    char name5[300];
    sprintf(name5,"CastorSumOfEnergyForAllSectors_%s",Folders.at(j).c_str());
    TH1F *histo_AllSectorsCastorEnergy = new TH1F(name5,"Castor Sector: #sum Energy; #sum E_{modules 1,2,3,4,5} per sector [GeV]; N events",200,-10,10);
    m_hVector_AllSectorsCastorEnergy.push_back(histo_AllSectorsCastorEnergy);

    for (int bs=0; bs<80; bs++){
      m_hVector_ChannelCastorEnergy.push_back( std::vector<TH1F*>() );

      char name6[300];
      char name6t[300];
      sprintf(name6,"Channel%d_Energy_%s",bs+1,Folders.at(j).c_str());
      sprintf(name6t,"#Castor Channel %d; #E_{Channel %d} [GeV]; N events", bs+1, bs+1);
      TH1F *histo_ChannelCastorEnergy = new TH1F(name6,name6t,200,-1,1);
      m_hVector_ChannelCastorEnergy[bs].push_back(histo_ChannelCastorEnergy);
    }

  }

}

void CastorThreshold::FillHistos(int index){

  bool debug = false;

  m_hVector_tracks[index]->Fill(eventdiff->GetMultiplicityTracks());
  m_hVector_vertex[index]->Fill(eventdiff->GetNVertex());
  m_hVector_RunNumber[index]->Fill(eventdiff->GetRunNumber());

  for (int i=0; i < 16; i++){
    CastorEnergySector[i]=0;
    if (i==4 || i==5){
      CastorEnergySector[i]=eventCastor->GetCastorModule2Energy(i)+eventCastor->GetCastorModule3Energy(i)+eventCastor->GetCastorModule4Energy(i)+eventCastor->GetCastorModule5Energy(i);
    }else{
      CastorEnergySector[i]=eventCastor->GetCastorModule1Energy(i)+eventCastor->GetCastorModule2Energy(i)+eventCastor->GetCastorModule3Energy(i)+eventCastor->GetCastorModule4Energy(i)+eventCastor->GetCastorModule5Energy(i);
    }
    m_hVector_SectorCastorEnergy[i].at(index)->Fill(CastorEnergySector[i]);
    m_hVector_AllSectorsCastorEnergy[index]->Fill(CastorEnergySector[i]);
  }

  for (int i=0; i < 16; i++){
    m_hVector_ChannelCastorEnergy[i].at(index)->Fill(eventCastor->GetCastorModule1Energy(i));
    m_hVector_ChannelCastorEnergy[16+i].at(index)->Fill(eventCastor->GetCastorModule2Energy(i));
    m_hVector_ChannelCastorEnergy[32+i].at(index)->Fill(eventCastor->GetCastorModule3Energy(i));
    m_hVector_ChannelCastorEnergy[48+i].at(index)->Fill(eventCastor->GetCastorModule4Energy(i));
    m_hVector_ChannelCastorEnergy[64+i].at(index)->Fill(eventCastor->GetCastorModule5Energy(i));   
  }


  if (debug){
    for (int i=0; i<16; i++){
      std::cout << "\nCastor Sector(" << i+1 << "): " << CastorEnergySector[i] << std::endl;
    }
  } 

}

void CastorThreshold::SaveHistos(){

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){
    m_hVector_tracks[j]->Write();
    m_hVector_vertex[j]->Write();
    m_hVector_RunNumber[j]->Write();
    m_hVector_AllSectorsCastorEnergy[j]->Write();
    for (int id=0; id<16; id++){
      m_hVector_SectorCastorEnergy[id].at(j)->Write();
    }

    for (int id=0; id<80; id++){
      m_hVector_ChannelCastorEnergy[id].at(j)->Write();
    }
  }

}

void CastorThreshold::Run(std::string filein_, std::string savehistofile_, std::string processname_, std::string type_, int runmin_, int runmax_){

  bool debug = false;

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  type = type_;
  runmin = runmin_;
  runmax = runmax_;

  TFile check1(filein.c_str());

  if (check1.IsZombie()){

    std::cout << "\n----------------------------------------------" << std::endl;
    std::cout << " There is no the file " << filein << " or the"   << std::endl;
    std::cout << " path is not correct." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    return;

  }

  if (check1.GetDirectory(processname.c_str())){
    LoadFile(filein,processname);
  }

  else {
    std::cout << "\n-------------------------------------------------" << std::endl;
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
  int triggercounter[20]={0};

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  std::string status;  

  for(int i=0;i<NEVENTS;i++){

    tr->GetEntry(i);

    if (!debug){
      if (i==0) {
	std::cout << "" << std::endl;
	std::cout<< "Status Bar" << std::endl;
	std::cout << "" << std::endl;
      }
      loadBar(i,NEVENTS,100,100);
    }

    for (int nt=0;nt<20;nt++){
      if(eventCastor->GetHLTPath(nt)>0){
	triggercounter[nt]++;
      }
    }

    bool triggerZeroBias = false;
    bool triggerHLTPlus = false;
    bool triggerHLTMinus = false;
    bool vertex = false;
    bool tracks = false;
    bool collisions = false;
    bool nocollisions = false;
    bool unpaired = false;
    bool runselection = false;

    if (eventdiff->GetRunNumber() >= runmin && eventdiff->GetRunNumber() <= runmax) runselection = true;
    if (eventCastor->GetHLTPath(0)>0) triggerZeroBias = true;
    if (eventCastor->GetHLTPath(1)>0) triggerHLTPlus = true;
    if (eventCastor->GetHLTPath(2)>0) triggerHLTMinus = true;
    if (eventdiff->GetMultiplicityTracks() > 0) tracks = true;  
    if (eventdiff->GetNVertex() > 0) vertex = true;

    if (type == "collisions"){
      if (triggerZeroBias && vertex && tracks) collisions = true;
      status = "collisions";
      if (runselection) FillHistos(0); 
      if (runselection && collisions) FillHistos(1);
    }

    else if (type == "no_collisions"){
      if(triggerZeroBias && !vertex && !tracks) nocollisions = true; 
      status = "no collisions";
      if (runselection) FillHistos(0);
      if (runselection && nocollisions) FillHistos(1);
    }

    else if (type == "unpaired"){
      if((triggerHLTPlus || triggerHLTMinus) && !vertex && !tracks) unpaired = true;
      status = "unpaired";
      if (runselection) FillHistos(0);
      if (runselection && unpaired) FillHistos(1);
    }

    else if (type == "collisionsmc"){
      if (vertex && tracks) collisions = true;
      status = "collisionsmc";
      FillHistos(0);
      if (collisions) FillHistos(1);
    }

    else if (type == "unpairedmc"){
      if(!vertex && !tracks) unpaired = true;
      status = "unpairedmc";
      FillHistos(0);
      if (unpaired) FillHistos(1);
    }

    else {
      std::cout << "\n Unrecognized Type of Selection." << std::endl;
      return;
    }

  }   


  outstring << "" << std::endl;
  outstring << "<< INPUTS >>" << std::endl;
  outstring << " " << std::endl;
  outstring << "Input file: " << filein << std::endl;
  outstring << "Output file: " << savehistofile << std::endl;
  outstring << " " << std::endl;
  outstring << "Type: " << status << std::endl;
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
  std::string savehistofile_;
  std::string processname_;
  std::string type_;
  int runmin_;
  int runmax_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0)  filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0)  savehistofile_  = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0)  processname_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0)  type_  = argv[4];
  if (argc > 5 && strcmp(s1,argv[5]) != 0)  runmin_  = atoi(argv[5]);
  if (argc > 6 && strcmp(s1,argv[6]) != 0)  runmax_  = atoi(argv[6]);

  if (type_=="collisions" || type_=="no_collisions" || type_=="unpaired" || type_=="collisionsmc" || type_=="unpairedmc") {}
  else{
    std::cout << "Please Insert type of selection: " << std::endl;
    std::cout << "1) collisions: ZeroBias trigger. Tracks and Vertex > 0." << std::endl;
    std::cout << "2) no_collisions: ZeroBias trigger. No tracks and no vertex." << std::endl;
    std::cout << "3) unpaired: HLT_BPTX Minus or Plus. No tracks and no vertex." << std::endl;
    std::cout << "4) collisionsmc: Tracks and Vertex > 0." << std::endl;
    std::cout << "5) unpairedmc: No tracks and no vertex." << std::endl;
    return 0;
  }

  if (runmin_ < 0 || runmax_ < 0) {
    std::cout << "Please Insert Run Min. or Run Max. > 0." << std::endl;
    return 0;
  }

  CastorThreshold* castor = new CastorThreshold();   
  castor->CreateHistos(type_);
  castor->Run(filein_, savehistofile_, processname_, type_, runmin_, runmax_);

  std::cout << "\n" << std::endl;

  return 0;
}
#endif
