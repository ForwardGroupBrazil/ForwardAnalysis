//-------------------------------------------------------------------------------------------------------->>
// UNIVERSIDADE DO ESTADO DO RIO DE JANEIRO - CMS/Brasil
//-------------------------------------------------------------------------------------------------------->>
//
// Instructions to run the code: (a) using script ro run multiple parameters or (b) bash command line. 
//
// First Step: compile the code.
//             ?> scram b -r
//
// Second Step: copy from $CMSSW_BASE/test/slc.../Binary file to your own area with pileup root files.
//
// (A) SCRIPT TO RUN MULTIPLE PARAMETERS
// -------------------------------------
// $> python CreateMCPileUp.py
//
//
// (B) COMMAND LINE
// ----------------
//
// $> ./mcPU_Distributions "Inputfile.root" "outputfile.root" "CMSSW Process_Name" "TTree_name" <Number of Bins>
//
// Example: 
//
// $> ./mcPU_Distributions "In.root" "Out.root" "RECO" "ttree" 3
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

#include "mcPU_Distributions.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/EventInfoEvent.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace eventInfo;

void mcPU_Distributions::LoadFile(std::string fileinput, std::string processinput){

  inf = NULL;
  tr  = NULL;
  inf = TFile::Open(fileinput.c_str(),"read");
  std::string fdirectory = processinput + "/ProcessedTree";
  tr = (TTree*)inf->Get(fdirectory.c_str());
  eventinfo = new EventInfoEvent();
  info = tr->GetBranch("EventInfo");
  info->SetAddress(&eventinfo);

}

void mcPU_Distributions::Run(std::string filein_, std::string savehistofile_, std::string processname_, int nbins_){

  filein = filein_;
  savehistofile = savehistofile_;
  processname = processname_;
  nbins = nbins_;

  std::cout << "Running..." << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Input file: " << filein << std::endl;
  std::cout << "Output file: " << savehistofile << std::endl;
  std::cout << "Process Name: " << processname_ << std::endl;
  std::cout << "N Bins PU Distributions: " << nbins << std::endl;
  std::cout << "" << std::endl;

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
  }
  else{
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << " There is no directory/path " << processname << std::endl;
    std::cout << " in the file." << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    return;
  }

  TFile *outf = new TFile(savehistofile.c_str(),"RECREATE");

  int NEVENTS = tr->GetEntries();

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  for(int m=0;m<2;m++) {
    tr->GetEntry(m);
    if (eventinfo->GetNPileUpBx0()==-1 && eventinfo->GetNPileUpBxm1()==-1 && eventinfo->GetNPileUpBxp1()==-1 ){
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << " There is no Pile Up TTree information in your PATTuplefile."   << std::endl;
      std::cout << " Please, use another PATTuple with PU information to run mul- " << std::endl;
      std::cout << " tiple PU option." << std::endl;
      std::cout << "--------------------------------------------------------------" << std::endl;
      return;
    }
  }

  // Event by Event Analysis
  //////////////////////////

  std::vector <std::string> Folders;
  Folders.push_back("without_cuts");

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char name1[300];
    sprintf(name1,"pileUpBx0_complete_%s",Folders.at(j).c_str());
    TH1D *histo_pumcbx0 = new TH1D(name1,"PileUp Monte Carlo; # Pile Up; N events",nbins,0,nbins);
    m_hVector_pumcbx0.push_back(histo_pumcbx0);

    char name2[300];
    sprintf(name2,"pileUpBxm1_complete_%s",Folders.at(j).c_str());
    TH1D *histo_pumcbxm1 = new TH1D(name2,"PileUp Monte Carlo; # Pile Up; N events",nbins,0,nbins);
    m_hVector_pumcbxm1.push_back(histo_pumcbxm1);

    char name3[300];
    sprintf(name3,"pileUpBxp1_complete_%s",Folders.at(j).c_str());
    TH1D *histo_pumcbxp1 = new TH1D(name3,"PileUp Monte Carlo; # Pile Up; N events",nbins,0,nbins);
    m_hVector_pumcbxp1.push_back(histo_pumcbxp1);

  }

  for(int i=0;i<NEVENTS;i++) {

    tr->GetEntry(i);

    // Without Cuts          
    ////////////////////////////////////////////////
    m_hVector_pumcbx0.at(0)->Fill(eventinfo->GetNPileUpBx0());
    m_hVector_pumcbxm1.at(0)->Fill(eventinfo->GetNPileUpBxm1());
    m_hVector_pumcbxp1.at(0)->Fill(eventinfo->GetNPileUpBxp1());
    //////////////////////////////////////////////////

  }// Run All Events

  outf->Write();
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
  int nbins_;

  if (argc > 1 && strcmp(s1,argv[1]) != 0)  filein_ = argv[1];
  if (argc > 2 && strcmp(s1,argv[2]) != 0)  savehistofile_  = argv[2];
  if (argc > 3 && strcmp(s1,argv[3]) != 0)  processname_  = argv[3];
  if (argc > 4 && strcmp(s1,argv[4]) != 0)  nbins_ = atoi(argv[4]);

  mcPU_Distributions* puObj = new mcPU_Distributions();   
  puObj->Run(filein_, savehistofile_, processname_, nbins_);
  return 0;
}
#endif
