#ifndef mcPU_Distributions_h
#define mcPU_Distributions_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class EventInfoEvent;

class mcPU_Distributions {

   TFile* inf;
   TTree* tr;
   TBranch *info;
   EventInfoEvent *eventinfo;

   std::string filein;
   std::string fileinput;
   std::string processname;
   std::string savehistofile;
   int nbins;
   std::vector<TH1D*> m_hVector_pumcbx0;
   std::vector<TH1D*> m_hVector_pumcbxm1;
   std::vector<TH1D*> m_hVector_pumcbxp1;

   public :
   mcPU_Distributions() {}
   ~mcPU_Distributions() { inf->Close(); }
   
   void Run(std::string, std::string, std::string, int);
   void LoadFile(std::string,std::string);
   void FillHistograms();

};
#endif
