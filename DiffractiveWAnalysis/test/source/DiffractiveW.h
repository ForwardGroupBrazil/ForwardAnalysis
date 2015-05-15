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
  double totalweight;

  double leptonpt;
  double leptoneta;
  double leptonphi;
  double deltaphi;
  double secondleptonpt;
  double secondleptonphi;
  double secondleptoneta;

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
  double deltaeta;
  double etamin;
  double etamax;
  double xi;
  double xiplus;
  double ximinus;
  double xigen;
  double xigenplus;
  double xigenminus;
  double maxLRG;
  double resLeadingPt;
  double resLeadingEta;
  double resLeadingPhi;
  double resMETPt;
  double resMETPhi;

  double sumCastorAndHFMinusEnergy;
  int SectorCastorHit;
  double castorthreshold;
  double channelsthreshold;
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
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonPt;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonEta;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonPhi;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonPt;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonEta;
  std::vector<std::vector<TH1F*> > m_hVector_SecondLeptonPhi;
  std::vector<std::vector<TH1F*> > m_hVector_DeltaPhi;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonDxy;
  std::vector<std::vector<TH1F*> > m_hVector_LeadingLeptonDz;
  std::vector<std::vector<TH1F*> > m_hVector_METPt;
  std::vector<std::vector<TH1F*> > m_hVector_METSignificance;
  std::vector<std::vector<TH1F*> > m_hVector_METPhi;
  std::vector<std::vector<TH1F*> > m_hVector_DVtxMuon;
  std::vector<std::vector<TH1F*> > m_hVector_DVtxMuonZ;
  std::vector<std::vector<TH1F*> > m_hVector_DVtxElectron;
  std::vector<std::vector<TH1F*> > m_hVector_DVtxElectronZ;
  std::vector<std::vector<TH1F*> > m_hVector_DMuonElectron;
  std::vector<std::vector<TH1F*> > m_hVector_DMuonElectronZ;
  std::vector<std::vector<TH1F*> > m_hVector_DMuons;
  std::vector<std::vector<TH1F*> > m_hVector_DElectrons;

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
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFplusPF;
  std::vector<std::vector<TH1F*> > m_hVector_sumEHFminusPF;
  std::vector<std::vector<TH2F*> > m_hVector_EPFVsTowerMinus;
  std::vector<std::vector<TH2F*> > m_hVector_EPFVsTowerPlus;

  //Event Info
  std::vector<std::vector<TH1F*> > m_hVector_lumi;
  std::vector<std::vector<TH1F*> > m_hVector_tracks;
  std::vector<std::vector<TH1F*> > m_hVector_vertex;
  std::vector<std::vector<TH1F*> > m_hVector_pu;


  // Diffraction
  std::vector<std::vector<TH1F*> > m_hVector_asumE;
  std::vector<std::vector<TH2F*> > m_hVector_multhf;
  std::vector<std::vector<TH1F*> > m_hVector_minusnhf;
  std::vector<std::vector<TH1F*> > m_hVector_plusnhf;
  std::vector<std::vector<TH1F*> > m_hVector_etamax;
  std::vector<std::vector<TH1F*> > m_hVector_etamin;
  std::vector<std::vector<TH1F*> > m_hVector_maxetagap;
  std::vector<std::vector<TH1F*> > m_hVector_SumPTMax;
  std::vector<std::vector<TH1F*> > m_hVector_SumPTMin;
  std::vector<std::vector<TH1F*> > m_hVector_absdeltaEta;
  std::vector<std::vector<TH1F*> > m_hVector_deltaEta;
  std::vector<std::vector<TH1F*> > m_hVector_XiPlus;
  std::vector<std::vector<TH1F*> > m_hVector_XiMinus;
  std::vector<std::vector<TH1F*> > m_hVector_Xi;
  std::vector<std::vector<TH1F*> > m_hVector_AEcastor;
  std::vector<std::vector<TH1F*> > m_hVector_etasignedHF;
  std::vector<std::vector<TH1F*> > m_hVector_etasignedCASTOR;

  // Generator
  std::vector<std::vector<TH1F*> > m_hVector_genProtonMinusXi;
  std::vector<std::vector<TH1F*> > m_hVector_genProtonPlusXi;
  std::vector<std::vector<TH1F*> > m_hVector_genXiPlus;
  std::vector<std::vector<TH1F*> > m_hVector_genXiMinus;
  std::vector<std::vector<TH1F*> > m_hVector_genXi;
  std::vector<std::vector<TH1F*> > m_hVector_resLeadingLeptonPt;
  std::vector<std::vector<TH1F*> > m_hVector_resLeadingLeptonEta;
  std::vector<std::vector<TH1F*> > m_hVector_resLeadingLeptonPhi;
  std::vector<std::vector<TH1F*> > m_hVector_resMETPt;
  std::vector<std::vector<TH1F*> > m_hVector_resMETPhi;
  std::vector<std::vector<TH1F*> > m_hVector_resXiPlus;
  std::vector<std::vector<TH1F*> > m_hVector_resXiMinus;
  std::vector<std::vector<TH1F*> > m_hVector_resHFEnergy;
  std::vector<std::vector<TH1F*> > m_hVector_resCASTOREnergy;
  std::vector<std::vector<TH2F*> > m_hVector_correlHFEnergy;
  std::vector<std::vector<TH2F*> > m_hVector_correlCASTOREnergy;
  std::vector<std::vector<TH2F*> > m_hVector_correlXiPlus;
  std::vector<std::vector<TH2F*> > m_hVector_correlXiMinus;

  std::vector <std::string> Folders;
  TDirectory *foldersFile[5];

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
