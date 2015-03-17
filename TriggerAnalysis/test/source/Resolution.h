#ifndef Resolution_h
#define Resolution_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class DiffractiveEvent;
class DiffractiveWEvent;
class EventInfoEvent;

class Resolution {

  TFile* inf;
  TTree* tr;
  TBranch *diff;
  TBranch *excl;
  TBranch *info;
  DiffractiveEvent *eventdiff;
  DiffractiveWEvent *eventexcl;
  EventInfoEvent *eventinfo;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  int optTrigger;
  int optTriggerOR;
  int optTriggerRef;
  int optTriggerRefOR;
  int bin;
  double channelsthreshold;
  double CastorEnergySector[16];
  std::vector<TH1D*> m_hVector_Evt_lumis;
  std::vector<TH1D*> m_hVector_Eff_lumis;
  std::vector<TH1D*> m_hVector_Evt_pfetamax;
  std::vector<TH1D*> m_hVector_Evt_pfetamin;
  public :
  Resolution() {}
  ~Resolution() { inf->Close(); }

  void Run(std::string, std::string, std::string, int, int, int, int, int, double);
  void LoadFile(std::string,std::string);
  void FillHistograms();

};
#endif

