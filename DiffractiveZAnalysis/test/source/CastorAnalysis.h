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
  int SectorCastorHitTh;
  double sumCastorEnergyTh;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  std::string switchtrigger;
  std::string typesel;
  int nVertex;
  int optTrigger;
  double lepton1pt;
  double lepton2pt;
  double SectorThreshold;
  double ChannelThreshold;
  int SectorZeroCastorCounter;
  int NEntries;

  int chk[5];
  double x_centroid;
  double y_centroid;
  double x_temp;
  double y_temp;
  double num_x_centroid;
  double num_y_centroid;
  double num_phi;
  double phi_average;
  double CastorEnergySector[16];
  double CastorEnergySectorTh[16];
  char selCastor[300];

  std::vector<TH2F*> m_hVector_histo_castor_centroid;
  std::vector<TH1F*> m_hVector_histo_castor_centroid_phi;
  std::vector<TProfile*> m_hVector_SectorVsTotalCastorEnergyTProf;
  std::vector<TH2F*> m_hVector_SectorVsTotalCastorEnergy;
  std::vector<TH1F*> m_hVector_RunNumber;
  std::vector<TH1F*> m_hVector_RunNumberZeroCastor;
  std::vector<TH1F*> m_hVector_RunNumberHighCastor;
  std::vector<TH2F*> m_hVector_CastorMultiplicityVsLumi;
  std::vector<TH1F*> m_hVector_CastorMultiplicity;
  std::vector<TH1F*> m_hVector_sumECastorAndHFMinus;
  std::vector<TProfile*> m_hVector_EnergyHFMinusVsCastorTProf;
  std::vector<TProfile*> m_hVector_EnergyHFPlusVsCastorTProf;
  std::vector<TH2F*> m_hVector_etcalos_n;
  std::vector<TH2F*> m_hVector_etcalos_p;
  std::vector<TH1F*> m_hVector_ECastorSectorBin1D;
  std::vector<TProfile*> m_hVector_ECastorSectorTProf;
  std::vector<TProfile*> m_hVector_ECaloVsEtaTProf;
  std::vector<TH2F*> m_hVector_ECastorSector;
  std::vector<TH1F*> m_hVector_sumECastorMinus;
  std::vector<TH2F*> m_hVector_ECaloVsEta;
  std::vector<TH2F*> m_hVector_CastorBadChannelVsRun;
  std::vector<TH1F*> m_hVector_CastorBadRuns;
  std::vector<TH1F*> m_hVector_CastorBadChannels;
  std::vector<TProfile*> m_hVector_CastorMultiplicityPerModuleTProf;
  std::vector<TH2F*> m_hVector_CastorMultiplicityPerModule;
  std::vector<TProfile*> m_hVector_CastorEnergyPerModuleTProf;
  std::vector<TH2F*> m_hVector_CastorEnergyPerModule;

  std::vector<std::vector<TH1F*> > m_hVector_TotalEnergySectors;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFplusBinSlice;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFminusBinSlice; 
  std::vector<std::vector<TH2F*> > m_hVector_Sector_EnergyVsMultiplicity; 
  std::vector<std::vector<TH1F*> > m_hVector_sumECastorMinusBinSlice; 
  std::vector<std::vector<TH1F*> > m_hVector_CastorMultiplicityModule;
  std::vector<std::vector<TH2F*> > m_hVector_AlongZ_EnergyVsModule;
  std::vector<std::vector<TProfile*> > m_hVector_AlongZ_EnergyVsModuleTProf;
  std::vector<std::vector<TH2F*> > m_hVector_CastorMultiplicityModuleAll;

  std::vector<TH2F*> m_hVector_CastorMappingMultiplicity;
  std::vector<TH2F*> m_hVector_CastorMappingEnergy;

  std::vector<std::vector<TH2F*> > m_hVector_CastorMappingMultiplicitySnapshot;
  std::vector<std::vector<TH2F*> > m_hVector_CastorMappingEnergySnapshot;

  std::vector <std::string> Folders;

  public :
  CastorAnalysis() {}
  ~CastorAnalysis() {
    inf->Close();
  }

  void Run(std::string, std::string, std::string, std::string, int, double, double, int, std::string, double, double);
  void LoadFile(std::string,std::string);
  void CreateHistos();
  void FillHistos(int);
  void SaveHistos();

};
#endif
