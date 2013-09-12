#ifndef CastorAnalysis_h
#define CastorAnalysis_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class DiffractiveEvent;
class DiffractiveZEvent;
class EventInfoEvent;

class CastorAnalysis {

  TFile* effcut;
  TFile* efftrigger;
  TFile* inf;
  TFile* pudata;
  TFile* pumc;
  TTree* tr;
  TBranch *diff;
  TBranch *Castor;
  TBranch *info;
  DiffractiveEvent *eventdiff;
  DiffractiveZEvent *eventCastor;
  EventInfoEvent *eventinfo;

  std::string fileinput;
  std::string processinput;
  int index;
  int l, k;

  double isoTk1;
  double isoTk2;
  double isoEcal1;
  double isoEcal2;
  double isoHcal1;
  double isoHcal2;
  int innerHits1;
  double Dcot1;
  double Dist1;
  double DeltaEtaTkClu1;
  double DeltaPhiTkClu1;
  double sigmaIeIe1;
  double HE1;
  int innerHits2;
  double Dcot2;
  double Dist2;
  double DeltaEtaTkClu2;
  double DeltaPhiTkClu2;
  double sigmaIeIe2;
  double HE2;
  double sumCastorAndHFMinusEnergy;
  int SectorCastorHit;
  double sumCastorEnergy;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  std::string switchtrigger;
  std::string typesel;
  int nVertex;
  int optTrigger;
  double lepton1pt;
  double lepton2pt;
  int SectorZeroCastorCounter;

  double x_centroid;
  double y_centroid;
  double x_temp;
  double y_temp;
  double num_x_centroid;
  double num_y_centroid;
  double num_phi;
  double phi_average;

  std::vector<TH2F*> m_hVector_histo_castor_centroid;
  std::vector<TH1F*> m_hVector_histo_castor_centroid_phi;

  std::vector <std::string> Folders;

  public :
  CastorAnalysis() {}
  ~CastorAnalysis() {
    inf->Close();
  }

  void Run(std::string, std::string, std::string, std::string, int, double, double, int, std::string);
  void LoadFile(std::string,std::string);
  void CreateHistos();
  void FillHistos(int);
  void SaveHistos();

};
#endif
