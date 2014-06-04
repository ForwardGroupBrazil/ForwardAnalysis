#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveWEvent.h"
#include <cstdio>

const char* DiffractiveWEvent::name = "DiffractiveWEvent";

DiffractiveWEvent::DiffractiveWEvent() {}

DiffractiveWEvent::~DiffractiveWEvent() {}

void DiffractiveWEvent::reset(){

  size_t len_hltTrigResults = sizeof(hltTrigResults_)/sizeof(int);
  for(size_t k = 0; k < len_hltTrigResults; ++k) hltTrigResults_[k] = 0;

  LeadingElectronPt_=-999.;
  LeadingElectronEta_=-999.;
  LeadingElectronPhi_=-999.;
  LeadingElectronCharge_=-999.;
  LeadingElectronEt_=-999.;
  ElectronsN_= -1;

  LeadingMuonPt_=-999.;
  LeadingMuonEta_=-999.;
  LeadingMuonPhi_=-999.;
  LeadingMuonCharge_=-999.;
  LeadingMuonEt_=-999.;
  MuonsN_= -1;

  VertexMultiplicity_.clear();
  VertexChiNorm_.clear();
  VertexNDOF_.clear();
  Vz_.clear();
  Vx_.clear();
  Vy_.clear();
  TracksPt_.clear();
  EachTowerEnergy_.clear();
  EachTowerEta_.clear();
  CastorTowerEnergy_.clear();
  CastorModule1Energy_.clear();
  CastorModule2Energy_.clear();
  CastorModule3Energy_.clear();
  CastorModule4Energy_.clear();
  CastorModule5Energy_.clear();
  CastorBadChannels_.clear();
  CastorNumberBadChannels_ = -999;
  EachTowerCounter_ = -1;

  PrimaryGapMaxGen_=-999.;
  SecondGapMaxGen_=-999.;
  TracksPtGen_.clear();
  EtaOfTracksPtGen_.clear();
  NTracksGen_=-999;
  Mx2PlusGen_=-999.;
  Mx2MinusGen_=-999.;
  Mx2Gen_=-999.;
  Mx2ZGen_=-999.;
  NMx2PlusGen_=-999;
  NMx2MinusGen_=-999;
  EtaGaplimPlusGen_=-999.;
  EtaGaplimMinusGen_=-999.;
  NParticlesGen_=-999;
  sumECastorMinusGen_=-999.;
  sumECastorPlusGen_=-999.;
  sumEZDCMinusGen_=-999.;
  sumEZDCPlusGen_=-999.;
  EtaOutcomingProtonGen_=-999.;
  xLGen_=-999.;
  xLMostEnergeticGen_=-999.;
  xiZMinusGen_=-999.;
  xiZPlusGen_=-999.;
  EtaZGen_=-999.;
  EnergyZGen_=-999.;
  pDissMassGen_=-999.;
  xLpDissMass_=-999.;

  SumEHFPlus_ = -999.; 
  SumEHF_SPlus_ = -999.; 
  SumEHF_LPlus_ = -999.; 
  SumEtHFPlus_ = -999.; 
  SumEHFMinus_ = -999.; 
  SumEHF_SMinus_ = -999.; 
  SumEHF_LMinus_ = -999.; 
  SumEtHFMinus_ = -999.; 
  SumEHEPlus_ = -999.; 
  SumEtHEPlus_ = -999.; 
  SumEHEMinus_ = -999.; 
  SumEtHEMinus_ = -999.; 
  SumEHBPlus_ = -999.; 
  SumEtHBPlus_ = -999.;   
  SumEHBMinus_ = -999.; 
  SumEtHBMinus_ = -999.; 
  SumEEEPlus_ = -999.;    
  SumEtEEPlus_ = -999.; 
  SumEEEMinus_ = -999.; 
  SumEtEEMinus_ = -999.; 
  SumEEBPlus_ = -999.; 
  SumEtEBPlus_ = -999.; 
  SumEEBMinus_ = -999.; 
  SumEtEBMinus_ = -999.; 
  EPZCaloPlus_ = -999.; 
  EPZCaloMinus_ = -999.; 
  XiCaloPlus_ = -999.; 
  XiCaloMinus_ = -999.; 
  EtaCaloMax_ = -999.; 
  EtaCaloMin_ = -999.; 

  Vertex_ = -999;
  Xi_PF_minus_ = -999.;
  Xi_PF_plus_ = -999.;
  Xi_mass_ = -999.;
  Xi_massnow_ = -999.;

  Epz_PF_minus_ = -999.;
  Epz_PF_plus_ = -999.;
  MultiplicityPF_ = -999;
  SumEtaTimesEnergyPF_ = -999.;
  SumpxModulePF_ = -999.;
  SumpyModulePF_ = -999.;
  SumpzModulePF_ = -999.;
  SumpxPF_ = -999.;
  SumpyPF_ = -999.;
  SumpzPF_ = -999.;
  SumEnergyPF_ = -999.;
  MuEnergyPF_ = -999.;
  ElectronEnergyPF_ = -999.;
  MaxGapPF_ = -999.;
  SecondMaxGapPF_ = -999.;
  LimPlusGapPF_ = -999.;
  LimMinusGapPF_ = -999.;
  PTMaxGapMaxPF_ = -999.;
  PTMinGapMaxPF_ = -999.;
  MultiplicityGapPlusPF_ = -999;
  MultiplicityGapMinusPF_ = -999;

  patNMuon_ = -1;

  patMuon1Pt_=-999.;
  patMuon1Charge_=-999;
  patMuon1Phi_=-999.;
  patMuon1Eta_=-999.;
  patMuon1Et_=-999.;

  patMuon1SumPtR03_=-999.;
  patMuon1EmEtR03_=-999.;
  patMuon1HadEtR03_=-999.;   
  patMuon1SumPtR05_=-999.;
  patMuon1EmEtR05_=-999.;
  patMuon1HadEtR05_=-999.;   

  patMuon1relIsoDr03_=-999.;
  patMuon1relIsoDr05_=-999.;

  patMuon1relIso_=-999.;

  patNElectron_ = -1;

  patElectron1Pt_=-999.;
  patElectron1Charge_=-999;
  patElectron1Phi_=-999.;
  patElectron1Eta_=-999.;
  patElectron1Et_=-999.;

  patElectron1TkDr03_=-999.;    
  patElectron1EcalDr03_=-999.;
  patElectron1HcalDr03_=-999.;

  patElectron1TkDr04_=-999.;
  patElectron1EcalDr04_=-999.;
  patElectron1HcalDr04_=-999.;

  patElectron1relIsoDr03_=-999.;
  patElectron1relIsoDr03_=-999.;

  LeadingElectronTkDr03_=-999.;
  LeadingElectronEcalDr03_=-999.;
  LeadingElectronHcalDr03_=-999.;

  LeadingElectronTkDr04_=-999.;
  LeadingElectronEcalDr04_=-999.;
  LeadingElectronHcalDr04_=-999.;

  LeadingElectronrelIsoDr03_=-999.;
  LeadingElectronrelIsoDr03_=-999.;

  LeadingMuonSumPtR03_=-999.;
  LeadingMuonEmEtR03_=-999.;
  LeadingMuonHadEtR03_=-999.;
  LeadingMuonSumPtR05_=-999.;
  LeadingMuonEmEtR05_=-999.;
  LeadingMuonHadEtR05_=-999.;

  LeadingMuonrelIsoDr03_=-999.;
  LeadingMuonrelIsoDr05_=-999.;

  LeadingElectronDeltaPhiTkClu_ = -999.;
  LeadingElectronDeltaEtaTkClu_ = -999.;
  LeadingElectronSigmaIeIe_ = -999.;
  LeadingElectronDCot_ = -999.;
  LeadingElectronDist_ = -999.;
  LeadingElectronInnerHits_ = -999.;
  LeadingElectronHE_ = -999.;
  patElectron1DeltaPhiTkClu_ = -999.;
  patElectron1DeltaEtaTkClu_ = -999.;
  patElectron1SigmaIeIe_ = -999.;
  patElectron1DCot_ = -999.;
  patElectron1Dist_ = -999.;
  patElectron1InnerHits_ = -999.;
  patElectron1HE_ = -999.;

  fmetPt_ = -999.;
  fmetPhi_ = -999.;
  fmetEt_ = -999.;
  fmetSumEt_ = -999.;
  fmetpx_ = -999.;
  fmetpy_ = -999.;

  fpatmetPt_ = -999.;
  fpatmetPhi_ = -999.;
  fpatmetEt_ = -999.;
  fpatmetSumEt_ = -999.;
  fpatmetpx_ = -999.;
  fpatmetpy_ = -999.;

  Mwenu_ = -999.;
  Mwmunu_ = -999.;
  Mpatwenu_ = -999.;
  Mpatmunu_ = -999.;

  //ZDCdigifC_.clear();

}
