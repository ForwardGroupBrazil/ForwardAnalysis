#ifndef ExclusiveDijet_h
#define ExclusiveDijet_h

#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

class DiffractiveEvent;
class ExclusiveDijetsEvent;
class EventInfoEvent;

class ExclusiveDijet {

  TFile* effcut;
  TFile* efftrigger;
  TFile* inf;
  TFile* pudata;
  TFile* pumc;
  TFile* fOut;
  TFile* fOutAll;
  TTree* tr;
  TTree* trout;
  TTree* troutAll;
  TH1* h_castor_channel;
  TBranch *diff;
  TBranch *excl;
  TBranch *info;
  DiffractiveEvent *eventdiff;
  ExclusiveDijetsEvent *eventexcl;
  EventInfoEvent *eventinfo;

  std::string fileinput;
  std::string processinput;
  int index;
  int pileup;
  int totalweight;
  double deltaphi, aSumE, absdeltaetapf, deltaetapf, ptJet1, ptJet2;

  int counterinfcut;
  int counterinftrigger;

  std::string filein;
  std::string processname;
  std::string savehistofile;
  std::string switchtrigger;
  std::string type;
  std::string jetunc;
  std::string switchpucorr;
  std::string pudatafile;
  std::string pumcfile;
  std::string cutcorrfile;
  std::string triggercorrfile;
  std::string switchcutcorr;
  std::string switchtriggercorr;
  std::string switchlumiweight;
  std::string switchcastor;
  std::string castorcorrfile;
  std::string switchpresel;
  double lumiweight;
  std::string switchmceventweight;
  int optnVertex;
  int optTrigger;
  double jet1pT;
  double jet2pT;
  double castorthreshold;
  double sumCastorEnergy;
  double sumCastorAndHFMinusEnergy;
  int SectorCastorHit;
  int SectorZeroCastorCounter;
  double etamin_;
  double energycorr[5][16];

  int bRunNumber;
  int bLumiSection;
  int bEventNumber;
  double bptJet1;
  double bptJet2;
  double bUnc1;
  double bUnc2;
  double bLeadingJetEta;
  double bSecondJetEta;
  double bLeadingJetPhi;
  double bSecondJetPhi;
  double bJetsDeltaEta;
  double bJetsDeltaphi;
  double bJetsDeltaPt;
  double bMassDijets;
  double bSumEnergyHFPlus;
  double bSumEnergyHFMinus;
  double bSumEnergyHEPlus;
  double bSumEnergyHEMinus;
  double bSumEHFPFlowPlus;
  double bSumEHFPFlowMinus;
  double bsumCastorEnergy;
  double bsumCastorAndHFMinusEnergy;
  double bMultiplicityHFPlus;
  double bMultiplicityHFMinus;
  double bMultiplicityTracks;
  double bTracksNonCone;
  double bTracksTransverse;
  double bTracksOutsideJets;
  double bRjjFromJets;
  double bEtaMaxFromPFCands;
  double betamin;
  double baSumE;
  double bdeltaetapf;
  double bXiPlusFromPFCands;
  double bXiMinusFromPFCands;

  std::vector<std::vector<TH1D*> > m_hVector_rjj;
  std::vector<std::vector<TH1D*> > m_hVector_detagen;
  std::vector<std::vector<TH1D*> > m_hVector_mxGen;
  std::vector<std::vector<TH2F*> > m_hVector_multhf;
  std::vector<std::vector<TH2F*> > m_hVector_etcalos;
  std::vector<std::vector<TH1D*> > m_hVector_tracks;
  std::vector<std::vector<TH1D*> > m_hVector_pfetamax;
  std::vector<std::vector<TH1D*> > m_hVector_pfetamin;
  std::vector<std::vector<TH1D*> > m_hVector_asumE;
  std::vector<std::vector<TH1D*> > m_hVector_deltaetajets;
  std::vector<std::vector<TH1D*> > m_hVector_deltaphijets;
  std::vector<std::vector<TH1D*> > m_hVector_deltaptjets;
  std::vector<std::vector<TH1D*> > m_hVector_dijetmass;
  std::vector<std::vector<TH1D*> > m_hVector_ptjet1;
  std::vector<std::vector<TH1D*> > m_hVector_ptjet2;
  std::vector<std::vector<TH1D*> > m_hVector_etajet1;
  std::vector<std::vector<TH1D*> > m_hVector_etajet2;
  std::vector<std::vector<TH1D*> > m_hVector_phijet1;
  std::vector<std::vector<TH1D*> > m_hVector_phijet2;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHFplus;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHFminus;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHEplus;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHEminus;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHFpfplus;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHFpfminus;
  std::vector<std::vector<TH1D*> > m_hVector_sumECastor;
  std::vector<std::vector<TH1D*> > m_hVector_deltaEtaPF;
  std::vector<std::vector<TH1D*> > m_hVector_absdeltaEtaPF;
  std::vector<std::vector<TH1D*> > m_hVector_vertex;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHFplusVsiEta;
  std::vector<std::vector<TH1D*> > m_hVector_sumEHFminusVsiEta;
  std::vector<std::vector<TH1D*> > m_hVector_lumi;
  std::vector<std::vector<TH2D*> > m_hVector_sumEHFplusVsetaMax;
  std::vector<std::vector<TH2D*> > m_hVector_sumEHFminusVsetaMin;
  std::vector<std::vector<TH2D*> > m_hVector_sumEHFpfplusVsetaMax;
  std::vector<std::vector<TH2D*> > m_hVector_sumEHFpfminusVsetaMin;
  std::vector<std::vector<TH1D*> > m_hVector_uncJet1;
  std::vector<std::vector<TH1D*> > m_hVector_uncJet2;
  std::vector<std::vector<TH1D*> > m_hVector_sumECastorHFMinus;
  std::vector<std::vector<TH1D*> > m_hVector_TracksNonCone;
  std::vector<std::vector<TH1D*> > m_hVector_TracksTransverse;
  std::vector<std::vector<TH1D*> > m_hVector_TracksOutsideJets;
  std::vector<std::vector<TH1D*> > m_hVector_XiPlusPF;
  std::vector<std::vector<TH1D*> > m_hVector_XiMinusPF;
  std::vector<std::vector<TH1D*> > m_hVector_XiPF;
  std::vector<std::vector<TH1D*> > m_hVector_ptjets;
  std::vector<std::vector<TH1D*> > m_hVector_etajets;
  std::vector<std::vector<TH1D*> > m_hVector_phijets;

  std::vector <std::string> Folders;

  public :
  ExclusiveDijet() {}
  ~ExclusiveDijet() {
    inf->Close();
  }

  void Run(std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, std::string, int, int, double, double, std::string, std::string, double, std::string);
  void LoadFile(std::string,std::string);
  void CreateHistos(std::string, std::string);
  void FillHistos(int, int, double);
  void SaveHistos(std::string);
  double* cutCorrection();
  double* triggerCorrection();

};
#endif
