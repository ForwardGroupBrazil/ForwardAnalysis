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
  TFile* fOutW;
  TFile* fOutCASTOR;
  TTree* tr;
  TBranch *diff;
  TBranch *diffW;
  TBranch *info;
  TTree* trout;
  TTree* troutW;
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

  double bosonWMass;
  double bosonWEta;
  double bosonWPhi;
  double bosonWPt;
  int bosonWCharge;
  int NElectrons;
  int NMuons;
  double isoTk1;
  double isoEcal1;
  double isoHcal1;
  double isoRec;
  int innerHits1;
  double Dcot1;
  double Dist1;
  double DeltaEtaTkClu1;
  double DeltaPhiTkClu1;
  double sigmaIeIe1;
  double HE1;
  double aSumE;
  double deltaetapf;
  double deltaetapfcastor;

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
  double sumCastorEnergy;

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
  double SumPTMinLrgPF;
  double SumPTMaxLrgPF;
  double bXiPlusFromPFCands;
  double bXiMinusFromPFCands;
  double betamax;
  double betamin;
  double betamincastor;
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

  // Kinematics
  std::vector<std::vector<TH1F*> > m_hVector_WMass;
  std::vector<std::vector<TH1F*> > m_hVector_WEta;
  std::vector<std::vector<TH1F*> > m_hVector_WPt;
  std::vector<std::vector<TH1F*> > m_hVector_WPhi;
  std::vector<std::vector<TH1F*> > m_hVector_WCharge;
  std::vector<std::vector<TH1F*> > m_hVector_NElectrons;
  std::vector<std::vector<TH1F*> > m_hVector_NMuons;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonTkDr03;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonEcalDr03;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonHcalDr03;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonIsolation;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonInnerHits;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonDCot;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonDist;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonDeltaEtaTkClu;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonDeltaPhiTkClu;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonSigmaIeIe;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonHE;

  // Detector
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFplus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFminus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHEplus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHEminus;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFplus_S;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFminus_S;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFplus_L;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFminus_L;
  std::vector<std::vector<TH1F*> > m_hVector_sumEEEplus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEEEminus;
  std::vector<std::vector<TH2F*> > m_hVector_EnergyHFPlusVsEnergyHFMinus;
  std::vector<std::vector<TH2F*> > m_hVector_EnergyEEPlusVsEnergyEEMinus;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFMax;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFMin;
  std::vector<std::vector<TH2F*> > m_hVector_etcalos_p;
  std::vector<std::vector<TH2F*> > m_hVector_etcalos_n;
  std::vector<std::vector<TH2F*> > m_hVector_ECaloVsEta;
  std::vector<std::vector<TProfile*> > m_hVector_ECaloVsEtaTProf;
  std::vector<std::vector<TH1F*> > m_hVector_EnergyVsEtaBin1D;
  std::vector<std::vector<TH1F*> > m_hVector_sumECastorMinus;
  std::vector<std::vector<TH2F*> > m_hVector_ECastorSector;
  std::vector<std::vector<TProfile*> > m_hVector_ECastorSectorTProf;
  std::vector<std::vector<TH1F*> > m_hVector_ECastorSectorBin1D;
  std::vector<std::vector<TProfile*> > m_hVector_EnergyHFMinusVsCastorTProf;
  std::vector<std::vector<TProfile*> > m_hVector_EnergyHFPlusVsCastorTProf;
  std::vector<std::vector<TH1F*> > m_hVector_sumECastorAndHFMinus;
  std::vector<std::vector<TH1F*> > m_hVector_CastorMultiplicity;
  std::vector<std::vector<TH2F*> > m_hVector_CastorMultiplicityVsLumi;
  std::vector<std::vector<TH2F*> > m_hVector_SectorVsTotalCastorEnergy;
  std::vector<std::vector<TProfile*> > m_hVector_SectorVsTotalCastorEnergyTProf;

  //Event Info
  std::vector<std::vector<TH1F*> > m_hVector_lumi;
  std::vector<std::vector<TH1F*> > m_hVector_tracks;
  std::vector<std::vector<TH1F*> > m_hVector_vertex;

  // Diffraction
  std::vector<std::vector<TH1F*> > m_hVector_asumE;
  std::vector<std::vector<TH2F*> > m_hVector_multhf;
  std::vector<std::vector<TH1F*> > m_hVector_pfetamax;
  std::vector<std::vector<TH1F*> > m_hVector_pfetamin;
  std::vector<std::vector<TH1F*> > m_hVector_pfetamincastor;
  std::vector<std::vector<TH1F*> > m_hVector_maxetagap;
  std::vector<std::vector<TH1F*> > m_hVector_SumPTMax;
  std::vector<std::vector<TH1F*> > m_hVector_SumPTMin;
  std::vector<std::vector<TH1F*> > m_hVector_absdeltaEtaPF;
  std::vector<std::vector<TH1F*> > m_hVector_deltaEtaPF;
  std::vector<std::vector<TH1F*> > m_hVector_absdeltaEtaPFCastor;
  std::vector<std::vector<TH1F*> > m_hVector_deltaEtaPFCastor;
  std::vector<std::vector<TH1F*> > m_hVector_XiPlusPF;
  std::vector<std::vector<TH1F*> > m_hVector_XiMinusPF;
  std::vector<std::vector<TH1F*> > m_hVector_XiPF;
  std::vector<std::vector<TH1F*> > m_hVector_AEcastor;
  std::vector<std::vector<TH1F*> > m_hVector_etasignedHF;
  std::vector<std::vector<TH1F*> > m_hVector_etasignedCASTOR;
 
  std::vector <std::string> Folders;
  TDirectory *foldersFile[4];

  public :
  DiffractiveW() {}
  ~DiffractiveW() {
    inf->Close();
  }

  void Run(std::string, std::string, std::string, std::string, int, double, double, int, std::string, std::string, double, std::string, double, double, std::string, std::string, std::string, std::string);
  void LoadFile(std::string,std::string);
  void CreateHistos(std::string);
  void CleanVariables();
  void FillHistos(int, int, double);
  void SaveHistos(std::string, std::string);

};
#endif
