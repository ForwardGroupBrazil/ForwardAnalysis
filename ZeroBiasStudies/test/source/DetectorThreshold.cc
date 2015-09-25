//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
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
#include "DetectorThreshold.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ZeroBiasEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace zerobiasAnalysis;
using namespace eventInfo;
using namespace reweight;

void DetectorThreshold::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventZeroBias = new ZeroBiasEvent();
  eventinfo = new EventInfoEvent();
  zeroBias = tr->GetBranch("ZeroBiasAnalysis");
  info = tr->GetBranch("EventInfo");
  zeroBias->SetAddress(&eventZeroBias);
  info->SetAddress(&eventinfo);

}


void DetectorThreshold::CreateHistos(std::string type){

  std::string step0 = "without_cuts_" + type;
  std::string step1 = "with_type_" + type;  

  Folders.push_back(step0);
  Folders.push_back(step1);

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name1[600];
    sprintf(name1,"Tracks_%s",Folders.at(j).c_str());
    TH1D *histo_Tracks = new TH1D(name1,"Tracks Multiplicity; n Tracks; N events",50,0,150);
    m_hVector_tracks.push_back(histo_Tracks);

    char name2[600];
    sprintf(name2,"vertex_%s",Folders.at(j).c_str());
    TH1D *histo_vertex = new TH1D(name2,"Number of Vertex; # Vertex; N events",25,0,25);
    m_hVector_vertex.push_back(histo_vertex);

    char name3[600];
    sprintf(name3,"ETowerHFPlus_%s",Folders.at(j).c_str());
    TH1D *histo_HFCaloPlus = new TH1D(name3,"HF^{+} Tower Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_HFCaloPlus.push_back(histo_HFCaloPlus);

    char name4[600];
    sprintf(name4,"ETowerHFMinus_%s",Folders.at(j).c_str());
    TH1D *histo_HFCaloMinus = new TH1D(name4,"HF^{-} Tower Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_HFCaloMinus.push_back(histo_HFCaloMinus);

    char name5[600];
    sprintf(name5,"EpfHFPlus_%s",Folders.at(j).c_str());
    TH1D *histo_HFpfPlus = new TH1D(name5,"HF^{+} PF Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_HFpfPlus.push_back(histo_HFpfPlus);

    char name6[600];
    sprintf(name6,"EpfHFMinus_%s",Folders.at(j).c_str());
    TH1D *histo_HFpfMinus = new TH1D(name6,"HF^{-} PF Energy; Energy [GeV]; N events",2000,0,1000);
    m_hVector_HFpfMinus.push_back(histo_HFpfMinus);

    for (int bs=0; bs<16; bs++){
      m_hVector_SectorCastorEnergy.push_back( std::vector<TH1D*>() );
      char castorname[600];
      char subtitle[600];
      sprintf(castorname,"Sector%d_CastorSumEnergy_%s",bs+1,Folders.at(j).c_str());
      sprintf(subtitle,"#sum Energy, Castor Sector %d; #sum E_{modules 1,2,3,4,5} [GeV]; N events", bs+1);
      TH1D *histo_SectorCastorEnergy = new TH1D(castorname,subtitle,2000,-100,100);
      m_hVector_SectorCastorEnergy[bs].push_back(histo_SectorCastorEnergy);

    }

    char name7[600];
    sprintf(name7,"CastorSumOfEnergyForAllSectors_%s",Folders.at(j).c_str());
    TH1D *histo_AllSectorsCastorEnergy = new TH1D(name7,"Castor Sector: #sum Energy; #sum E_{modules 1,2,3,4,5} all sectors [GeV]; N events",2000,-100,100);
    m_hVector_AllSectorsCastorEnergy.push_back(histo_AllSectorsCastorEnergy);

    for (int bs=0; bs<80; bs++){
      m_hVector_ChannelCastorEnergy.push_back( std::vector<TH1D*>() );
      char castornamecha[600];
      char subtitlecha[600];
      sprintf(castornamecha,"Channel%d_Energy_%s",bs+1,Folders.at(j).c_str());
      sprintf(subtitlecha,"#Castor Channel %d; #E_{Channel %d} [GeV]; N events", bs+1, bs+1);
      TH1D *histo_ChannelCastorEnergy = new TH1D(castornamecha,subtitlecha,200,-1,1);
      m_hVector_ChannelCastorEnergy[bs].push_back(histo_ChannelCastorEnergy);
    }

  }

}

void DetectorThreshold::FillHistos(int index){

  bool debug = false;

  m_hVector_tracks[index]->Fill(eventZeroBias->GetMultiplicityTracks());
  m_hVector_vertex[index]->Fill(eventZeroBias->GetNVertex());

  for (int i=0; i < eventZeroBias->GetNHFPlus(); i++){
    m_hVector_HFCaloPlus[index]->Fill(eventZeroBias->GetEHFPlus(i));
  }

  for (int i=0; i < eventZeroBias->GetNHFMinus(); i++){
    m_hVector_HFCaloMinus[index]->Fill(eventZeroBias->GetEHFMinus(i));
  }

  for (int i=0; i < eventZeroBias->GetPFNHFPlus(); i++){
    m_hVector_HFpfPlus[index]->Fill(eventZeroBias->GetPFHFPlus(i));
  }

  for (int i=0; i < eventZeroBias->GetPFNHFMinus(); i++){
    m_hVector_HFpfMinus[index]->Fill(eventZeroBias->GetPFHFMinus(i));
  }

  for (int i=0; i < 16; i++){
    CastorEnergySector[i]=0;
    if (i==4 || i==5){
      CastorEnergySector[i]=eventZeroBias->GetCastorModule2Energy(i)+eventZeroBias->GetCastorModule3Energy(i)+eventZeroBias->GetCastorModule4Energy(i)+eventZeroBias->GetCastorModule5Energy(i);
    }else{
      CastorEnergySector[i]=eventZeroBias->GetCastorModule1Energy(i)+eventZeroBias->GetCastorModule2Energy(i)+eventZeroBias->GetCastorModule3Energy(i)+eventZeroBias->GetCastorModule4Energy(i)+eventZeroBias->GetCastorModule5Energy(i);
    }
    m_hVector_SectorCastorEnergy[i].at(index)->Fill(CastorEnergySector[i]);
    EnergyAllCastor+=CastorEnergySector[i];
  }

  for (int i=0; i < 16; i++){
    m_hVector_ChannelCastorEnergy[i].at(index)->Fill(eventZeroBias->GetCastorModule1Energy(i));
    m_hVector_ChannelCastorEnergy[16+i].at(index)->Fill(eventZeroBias->GetCastorModule2Energy(i));
    m_hVector_ChannelCastorEnergy[32+i].at(index)->Fill(eventZeroBias->GetCastorModule3Energy(i));
    m_hVector_ChannelCastorEnergy[48+i].at(index)->Fill(eventZeroBias->GetCastorModule4Energy(i));
    m_hVector_ChannelCastorEnergy[64+i].at(index)->Fill(eventZeroBias->GetCastorModule5Energy(i));   
  }

  if (debug){
    for (int i=0; i<16; i++){
      std::cout << "\nCastor Sector(" << i+1 << "): " << CastorEnergySector[i] << std::endl;
    }
  } 

  m_hVector_AllSectorsCastorEnergy[index]->Fill(EnergyAllCastor);

}

void DetectorThreshold::SaveHistos(){

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    m_hVector_tracks[j]->Write();
    m_hVector_vertex[j]->Write();
    m_hVector_HFCaloPlus[j]->Write();
    m_hVector_HFCaloMinus[j]->Write();
    m_hVector_HFpfPlus[j]->Write();
    m_hVector_HFpfMinus[j]->Write();
    m_hVector_AllSectorsCastorEnergy[j]->Write();

    for (int id=0; id<16; id++){
      m_hVector_SectorCastorEnergy[id].at(j)->Write();
    }

    for (int id=0; id<80; id++){
      m_hVector_ChannelCastorEnergy[id].at(j)->Write();
    }
  }

}

void DetectorThreshold::Run(std::string filein_, std::string savehistofile_, std::string processname_, std::string type_, int runmin_, int runmax_){

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
      loadBar(i,NEVENTS);
    }

    for (int nt=0;nt<20;nt++){
      if(eventZeroBias->GetHLTPath(nt)>0){
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

    if (eventinfo->GetRunNumber() >= runmin && eventinfo->GetRunNumber() <= runmax) runselection = true;
    if (eventZeroBias->GetHLTPath(0)>0) triggerZeroBias = true;
    if (eventZeroBias->GetHLTPath(1)>0) triggerHLTPlus = true;
    if (eventZeroBias->GetHLTPath(2)>0) triggerHLTMinus = true;
    if (eventZeroBias->GetMultiplicityTracks() > 0) tracks = true;  
    if (eventZeroBias->GetNVertex() > 0) vertex = true;

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
  for (int i = 0; i<20; i++){
    outstring << "Trigger " << i << ": " << triggercounter[i] << std::endl;
  }
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

  DetectorThreshold* castor = new DetectorThreshold();   
  castor->CreateHistos(type_);
  castor->Run(filein_, savehistofile_, processname_, type_, runmin_, runmax_);

  std::cout << "\n" << std::endl;

  return 0;
}
#endif
