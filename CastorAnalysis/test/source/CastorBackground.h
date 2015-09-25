#ifndef CastorBackground_h
#define CastorBackground_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class EventInfoEvent;
class ZeroBiasEvent;
class DiffractiveWEvent;
class DiffractiveEvent;

class CastorBackground {

  TFile* inf;
  TTree* tr;

  TFile* infrandom;
  TTree* trrandom;

  TBranch *Castor;
  TBranch *info;
  TBranch *random;
  TBranch *Diff;

  EventInfoEvent *eventInfo;
  DiffractiveWEvent *eventCastor;
  ZeroBiasEvent *eventRandom;
  DiffractiveEvent *eventDiff;

  std::string fileinput;
  std::string fileinputrandom;

  std::string processinput;
  std::string processinputrandom;
  int index;

  std::string filein;
  std::string fileinrandom;
  std::string processname;
  std::string processnamerandom;
  std::string savehistofile;
  std::string type;

  double castorthreshold;
  double sumCastorEnergy;
  double sumCastorEnergy20Up;
  double sumCastorEnergy20Dw;
  double sumCastorEnergy10Up;
  double sumCastorEnergy10Dw;

  double CastorEnergySector[16];
  double CastorEnergySector20Up[16];
  double CastorEnergySector20Dw[16];
  double CastorEnergySector10Up[16];
  double CastorEnergySector10Dw[16];

  std::vector<TH1F*> m_hVector_tracks;
  std::vector<TH1F*> m_hVector_vertex;
  std::vector<TH1F*> m_hVector_RunNumber;
  std::vector<std::vector<TH1F*> > m_hVector_SectorCastorEnergy;
  std::vector<TH1F*> m_hVector_TotalCastorEnergyBySector;
  std::vector<TH1F*> m_hVector_TotalCastorEnergyBySector20PercUp;
  std::vector<TH1F*> m_hVector_TotalCastorEnergyBySector20PercDow;
  std::vector<TH1F*> m_hVector_TotalCastorEnergyBySector10PercUp;
  std::vector<TH1F*> m_hVector_TotalCastorEnergyBySector10PercDow;

  std::vector <std::string> Folders;

  //
  public :
  CastorBackground() {}
  ~CastorBackground() { inf->Close(); }

  void Run(std::string, std::string, std::string, std::string, std::string, std::string, double);
  void LoadFile(std::string,std::string,std::string,std::string);
  void CreateHistos(std::string);
  void FillHistos(int);
  void SaveHistos();

};
#endif
