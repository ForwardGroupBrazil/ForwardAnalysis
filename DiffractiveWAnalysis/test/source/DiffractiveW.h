#ifndef DiffractiveW_h
#define DiffractiveW_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

class DiffractiveEvent;
class DiffractiveWEvent;
class EventInfoEvent;

class DiffractiveW {

  TFile* effcut;
  TFile* efftrigger;
  TFile* inf;
  TFile* pudata;
  TFile* pumc;
  TFile* fOut;
  TFile* fOutZ;
  TFile* fOutCASTOR;
  TTree* tr;
  TBranch *diff;
  TBranch *diffW;
  TBranch *info;
  TTree* trout;
  TTree* troutZ;
  TTree* troutCASTOR;

  TH1* h_castor_channel;
  DiffractiveEvent *eventdiff;
  DiffractiveWEvent *eventdiffW;
  EventInfoEvent *eventinfo;

  std::string fileinput;
  std::string processinput;
  std::string gapseltype;
  int index;
  int pileup;
  int totalweight;
  double aSumE;
  int l, k;

  double bosonWMass;
  double sumCastorEnergy;
  double isoTk1;
  double isoTk2;
  double isoEcal1;
  double isoHcal1;
  int innerHits1;
  double Dcot1;
  double Dist1;
  double DeltaEtaTkClu1;
  double DeltaPhiTkClu1;
  double sigmaIeIe1;
  double HE1;
  double sumCastorAndHFMinusEnergy;
  int SectorCastorHit;
  double castorthreshold;
  double channelsthreshold;
  double etamin_;
  bool castoractivity;
  bool castorgap;
  double AEcastor;
  double etasignedHF;
  double etasignedCASTOR;
  int counterHit;

  int bRunNumber;
  int bLumiSection;
  int bEventNumber;
  double bDiBosonPt;
  double bDiBosonEta;
  double bDiBosonPhi;
  double bDiBosonMass;
  int bMultiplicityTracks;
  double bSumEEEMinus;
  double bSumEEEPlus;
  double bSumEnergyHFMinus;
  double bSumEnergyHFPlus;
  double bsumCastorEnergy;
  double bSectorCastorHit;
  double bAEcastor;
  double betasignedHF;
  double betasignedCASTOR;
  double bMaxGapPF;
  double bPTMinGapMaxPF;
  double bPTMaxGapMaxPF;
  double bXiPlusFromPFCands;
  double bXiMinusFromPFCands;
  double betamax;
  double betamin;
  double betalimmin;
  double betalimmax;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  std::string switchtrigger;
  std::string type;
  std::string typesel;
  std::string switchlumiweight;
  std::string castorcorrfile;
  std::string pumfile;
  std::string pudfile;

  double mcweight;
  int nVertex;
  int optTrigger;
  double lepton1pt;
  double lepton2pt;
  int SectorZeroCastorCounter;
  double CastorEnergySector[16];

  std::vector<std::vector<TH1F*> >    m_hVector_WMuonMass;
  std::vector<std::vector<TH1F*> >    m_hVector_WMuonEta;
  std::vector<std::vector<TH1F*> >    m_hVector_WMuonPt;
  std::vector<std::vector<TH1F*> >    m_hVector_WMuonPhi;
  std::vector<std::vector<TH1F*> >    m_hVector_WElectronMass;
  std::vector<std::vector<TH1F*> >    m_hVector_WElectronEta;
  std::vector<std::vector<TH1F*> >    m_hVector_WElectronPt;
  std::vector<std::vector<TH1F*> >    m_hVector_WElectronPhi;

  std::vector <std::string> Folders;

  TDirectory *foldersFile[1];

  public :
  DiffractiveW() {}
  ~DiffractiveW() {
    inf->Close();
  }

  void Run(std::string, std::string, std::string, std::string, int, double, double, int, std::string, std::string, double, std::string, double, double, std::string, std::string, std::string, std::string);
  void LoadFile(std::string,std::string);
  void CreateHistos(std::string);
  void FillHistos(int, int, double);
  void SaveHistos(std::string, std::string);

};
#endif
