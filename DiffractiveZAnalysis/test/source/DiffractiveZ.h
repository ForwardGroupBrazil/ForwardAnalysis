#ifndef DiffractiveZ_h
#define DiffractiveZ_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

class DiffractiveEvent;
class DiffractiveZEvent;
class EventInfoEvent;

class DiffractiveZ {

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
  TBranch *diffZ;
  TBranch *info;
  TTree* trout;
  TTree* troutZ;
  TTree* troutCASTOR;

  TH1* h_castor_channel;
  DiffractiveEvent *eventdiff;
  DiffractiveZEvent *eventdiffZ;
  EventInfoEvent *eventinfo;

  std::string fileinput;
  std::string processinput;
  std::string gapseltype;
  int index;
  int pileup;
  int totalweight;

  double nMuons;
  double nElectrons;
  double dileptonMass;
  double dileptonEta;
  double dileptonPhi;
  double dileptonPt;
  double lepton1Pt;
  double lepton1Eta;
  double lepton1Phi;
  double lepton2Pt;
  double lepton2Eta;
  double lepton2Phi;
  double isoTk1;
  double isoEcal1;
  double isoHcal1;
  double innerHits1;
  double Dcot1;
  double Dist1;
  double DeltaEtaTkClu1;
  double DeltaPhiTkClu1;
  double sigmaIeIe1;
  double HE1;
  double isoTk2;
  double isoEcal2;
  double isoHcal2;
  double innerHits2;
  double Dcot2;
  double Dist2;
  double DeltaEtaTkClu2;
  double DeltaPhiTkClu2;
  double sigmaIeIe2;
  double HE2;
  double deltaetaleptons;
  double deltaphileptons;
  double deltaptleptons;
  double cone03tracks;
  double cone04tracks;
  double cone05tracks;
  double isoRec1;
  double isoRec2;
  double deltaeta;
  double deltaetacastor;
  double etamin_;
  double etamin;
  double etamax;
  double xiplus;
  double ximinus;
  double maxLRG;

  double aSumE;
  double AEcastor;
  double etasignedHF;
  double etasignedCASTOR;

  double sumCastorEnergy;
  double sumCastorAndHFMinusEnergy;
  int SectorCastorHit;
  double castorthreshold;
  double channelsthreshold;
  bool castoractivity;
  bool castorgap;
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
  double bdeltaeta;
  double bAEcastor;
  double betasignedHF;
  double betasignedCASTOR;
  double bMaxGap;
  double SumPTMinLrgPF;
  double SumPTMaxLrgPF;
  double bXiPlus;
  double bXiMinus;
  double betamax;
  double betamin;
  double betamincastor;

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
  std::vector<std::vector<TH1F*> > m_hVector_ElectronsN;
  std::vector<std::vector<TH1F*> > m_hVector_MuonsN;
  std::vector<std::vector<TH1F*> > m_hVector_BosonZPt;
  std::vector<std::vector<TH1F*> > m_hVector_BosonZEta;
  std::vector<std::vector<TH1F*> > m_hVector_BosonZPhi;
  std::vector<std::vector<TH1F*> > m_hVector_BosonZMass;
  std::vector<std::vector<TH1F*> > m_hVector_LeptonsPt;
  std::vector<std::vector<TH1F*> > m_hVector_LeptonsEta;
  std::vector<std::vector<TH1F*> > m_hVector_LeptonsPhi;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonPt;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonEta;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonPhi;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonCharge;
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
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonPt;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonEta;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonPhi;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonCharge;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonTkDr03;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonEcalDr03;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonHcalDr03;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonIsolation;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonInnerHits;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonDCot;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonDist;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonDeltaEtaTkClu;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonDeltaPhiTkClu;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonSigmaIeIe;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonHE;
  std::vector<std::vector<TH1F*> > m_hVector_deltaphiLeptons;
  std::vector<std::vector<TH1F*> > m_hVector_deltapTLeptons;
  std::vector<std::vector<TH1F*> > m_hVector_deltaetaLeptons;
  std::vector<std::vector<TH1F*> > m_hVector_tracksOutLeptonsCone03;
  std::vector<std::vector<TH1F*> > m_hVector_tracksOutLeptonsCone04;
  std::vector<std::vector<TH1F*> > m_hVector_tracksOutLeptonsCone05;

  // Event Info
  std::vector<std::vector<TH1F*> > m_hVector_RunNumber;
  std::vector<std::vector<TH1F*> > m_hVector_RunNumberZeroCastor;
  std::vector<std::vector<TH1F*> > m_hVector_RunNumberHighCastor;
  std::vector<std::vector<TH1F*> > m_hVector_vertex;
  std::vector<std::vector<TH1F*> > m_hVector_lumi;
  std::vector<std::vector<TH2F*> > m_hVector_vertexvslumi;
  std::vector<std::vector<TH1F*> > m_hVector_tracks;
  std::vector<std::vector<TH1F*> > m_hVector_tracksLow;

  // Detector
  std::vector<std::vector<TH2F*> > m_hVector_ECaloVsEta;
  std::vector<std::vector<TProfile*> > m_hVector_ECaloVsEtaTProf;
  std::vector<std::vector<TH1F*> > m_hVector_EnergyVsEtaBin1D;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFplus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFminus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHEplus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHEminus;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFplus_S;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFminus_S;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFplus_L;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFminus_L;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFMax;
  std::vector<std::vector<TH1F*> > m_hVector_SumEHFMin;
  std::vector<std::vector<TH2F*> > m_hVector_EnergyHFPlusVsEnergyHFMinus;
  std::vector<std::vector<TH2F*> > m_hVector_EnergyEEPlusVsEnergyEEMinus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEEEminus;
  std::vector<std::vector<TH1F*> > m_hVector_sumEEEplus;
  std::vector<std::vector<TH2F*> > m_hVector_multhf;
  std::vector<std::vector<TH2F*> > m_hVector_etcalos_p;
  std::vector<std::vector<TH2F*> > m_hVector_etcalos_n;
  std::vector<std::vector<TH2F*> > m_hVector_ECastorSector;
  std::vector<std::vector<TProfile*> > m_hVector_ECastorSectorTProf;
  std::vector<std::vector<TH1F*> > m_hVector_CastorMultiplicity;
  std::vector<std::vector<TH2F*> > m_hVector_CastorMultiplicityVsLumi;
  std::vector<std::vector<TH2F*> > m_hVector_SectorVsTotalCastorEnergy;
  std::vector<std::vector<TProfile*> > m_hVector_SectorVsTotalCastorEnergyTProf;
  std::vector<std::vector<TH1F*> > m_hVector_ECastorSectorBin1D;
  std::vector<std::vector<TH1F*> > m_hVector_sumECastorMinus;
  std::vector<std::vector<TH1F*> > m_hVector_sumECastorMinusLow;
  std::vector<std::vector<TH1F*> > m_hVector_sumECastorAndHFMinus;
  std::vector<std::vector<TProfile*> > m_hVector_EnergyHFPlusVsCastorTProf;
  std::vector<std::vector<TProfile*> > m_hVector_EnergyHFMinusVsCastorTProf;

  // Diffraction
  std::vector<std::vector<TH1F*> > m_hVector_asumE;
  std::vector<std::vector<TH1F*> > m_hVector_AEcastor;
  std::vector<std::vector<TH1F*> > m_hVector_etasignedHF;
  std::vector<std::vector<TH1F*> > m_hVector_etasignedCASTOR;
  std::vector<std::vector<TH1F*> > m_hVector_XiPlus;
  std::vector<std::vector<TH1F*> > m_hVector_XiMinus;
  std::vector<std::vector<TH1F*> > m_hVector_Xi;
  std::vector<std::vector<TH1F*> > m_hVector_etamincastor;
  std::vector<std::vector<TH1F*> > m_hVector_absdeltaEta;
  std::vector<std::vector<TH1F*> > m_hVector_deltaEta;
  std::vector<std::vector<TH1F*> > m_hVector_absdeltaEtaCastor;
  std::vector<std::vector<TH1F*> > m_hVector_deltaEtaCastor;
  std::vector<std::vector<TH1F*> > m_hVector_maxetagap;
  std::vector<std::vector<TH1F*> > m_hVector_LimPlusgap;
  std::vector<std::vector<TH1F*> > m_hVector_LimMinusgap;
  std::vector<std::vector<TH1F*> > m_hVector_SumPTMax;
  std::vector<std::vector<TH1F*> > m_hVector_SumPTMin;
  std::vector<std::vector<TH1F*> > m_hVector_etamax;
  std::vector<std::vector<TH1F*> > m_hVector_etamin;

  std::vector <std::string> Folders;
  TDirectory *foldersFile[4];

  public :
  DiffractiveZ() {}
  ~DiffractiveZ() {
    inf->Close();
  }

  void Run(std::string, std::string, std::string, std::string, int, double, double, int, std::string, std::string, double, std::string, double, double, std::string, std::string, std::string, std::string);
  void LoadFile(std::string,std::string);
  void CreateHistos(std::string);
  void FillHistos(int, int, double);
  void SaveHistos(std::string, std::string);
  void CleanVariables();

};
#endif
