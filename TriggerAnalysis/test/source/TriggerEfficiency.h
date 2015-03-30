#ifndef TriggerEfficiency_h
#define TriggerEfficiency_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class TriggerEvent;
class EventInfoEvent;

class TriggerEfficiency {

  TFile* inf;
  TTree* tr;
  TBranch *trigger;
  TBranch *info;
  TriggerEvent *eventtrigger;
  EventInfoEvent *eventinfo;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  std::string typesel;
  int optTrigger;
  int optTriggerRef;
  std::vector<TH1D*> m_hVector_leading_pt;
  std::vector<TH1D*> m_hVector_second_pt;
  std::vector<TH1D*> m_hVector_dilepton_mass;

  public :
  TriggerEfficiency() {}
  ~TriggerEfficiency() { inf->Close(); }

  void Run(std::string, std::string, std::string, std::string, int, int);
  void LoadFile(std::string,std::string);
  void FillHistograms();

};
#endif

