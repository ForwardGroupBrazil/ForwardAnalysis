#ifndef CastorThreshold_h
#define CastorThreshold_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class DiffractiveEvent;
class ExclusiveDijetsEvent;
class EventInfoEvent;
class DiffractiveZEvent;

class CastorThreshold {

  TFile* inf;
  TTree* tr;
  TBranch *diff;
  TBranch *Castor;
  TBranch *info;
  DiffractiveEvent *eventdiff;
  ExclusiveDijetsEvent *eventexcl;
  EventInfoEvent *eventinfo;
  DiffractiveZEvent *eventCastor;

  std::string fileinput;
  std::string processinput;
  int index;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  std::string type;
  int runmin;
  int runmax;

  double CastorEnergySector[16];

  std::vector<TH1F*> m_hVector_tracks;
  std::vector<TH1F*> m_hVector_vertex;
  std::vector<TH1F*> m_hVector_RunNumber;
  std::vector<TH1F*> m_hVector_AllSectorsCastorEnergy;

  std::vector<std::vector<TH1F*> > m_hVector_SectorCastorEnergy;
  std::vector<std::vector<TH1F*> > m_hVector_ChannelCastorEnergy;

  std::vector <std::string> Folders;

  //
  public :
  CastorThreshold() {}
  ~CastorThreshold() { inf->Close(); }

  void Run(std::string, std::string, std::string, std::string, int, int);
  void LoadFile(std::string,std::string);
  void CreateHistos(std::string);
  void FillHistos(int);
  void SaveHistos();

};
#endif
