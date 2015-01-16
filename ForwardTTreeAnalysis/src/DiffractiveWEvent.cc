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
  BosonElectronMass_=-999.;
  BosonElectronPt_ = -999.;
  BosonElectronEta_ = -999.;
  BosonElectronPhi_ = -999.;
  ElectronsN_= -1;

  LeadingMuonPt_=-999.;
  LeadingMuonEta_=-999.;
  LeadingMuonPhi_=-999.;
  LeadingMuonCharge_=-999.;
  BosonMuonMass_=-999.;
  BosonMuonPt_ = -999.;
  BosonMuonEta_ = -999.;
  BosonMuonPhi_ = -999.; 
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

  xigenminus_=-999.;
  xigenplus_=-999.;

  mxgenminus_=-999.;
  mxgenplus_=-999.;
  mx2genminus_=-999.;
  mx2genplus_=-999.;
  mxgenleft_=-999.;
  mx2genright_=-999.;
  mx2genleft_=-999.;
  mxgenright_=-999.;
  lrggen_=-999.;
  etamaxgen_=-999.;
  etamingen_=-999.;
  epluspzgen_=-999.;
  eminuspzgen_=-999.;
  etexpoplusgen_=-999.;
  etexpominusgen_=-999.;
  sumECastorMinusGen_=-999.;
  sumptgenleft_=-999.;
  sumptgenright_=-999.;

  mxgenminusCMS_=-999.;
  mxgenplusCMS_=-999.;
  mx2genminusCMS_=-999.;
  mx2genplusCMS_=-999.;
  mxgenleftCMS_=-999.;
  mxgenrightCMS_=-999.;
  mx2genleftCMS_=-999.;
  mx2genrightCMS_=-999.;
  lrggenCMS_=-999.;
  etamaxgenCMS_=-999.;
  etamingenCMS_=-999.;
  epluspzgenCMS_=-999.;
  eminuspzgenCMS_=-999.;
  etexpoplusgenCMS_=-999.;
  etexpominusgenCMS_=-999.;
  sumECastorMinusGenCMS_=-999.;
  sumptgenleftCMS_=-999.;
  sumptgenrightCMS_=-999.;

  BosonElectronMassPF_=-999.;
  BosonMuonMassPF_=-999.;

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
  EtEtaCaloPlus_ = -999.; 
  EtEtaCaloMinus_ = -999.; 
  EtaCaloMax_ = -999.; 
  EtaCaloMin_ = -999.; 
  lrgCalo_ = -999.;

  Vertex_ = -999;
  tracketamax_=-999.;
  tracketamin_=-999.;
  mxpfminus_ =-999.;
  mxpfplus_ = -999.;
  mx2pfminus_ = -999.;
  mx2pfplus_ = -999.;
  etamaxpf_ = -999.;
  etaminpf_ = -999.;
  epluspzpf_= -999.;
  eminuspzpf_=-999.;
  etexpopluspf_=-999.;
  etexpominuspf_=-999.;
  lrgPF_=-999.;
  mxpfleft_ =-999.;
  mxpfright_=-999.;
  mx2pfleft_=-999.;
  mx2pfright_=-999.;
  sumptpfleft_=-999.;
  sumptpfright_=-999.;

  mxpfnowminus_ =-999.;
  mxpfnowplus_ = -999.;
  mx2pfnowminus_ = -999.;
  mx2pfnowplus_ = -999.;
  etamaxpfnow_ = -999.;
  etaminpfnow_ = -999.;
  epluspzpfnow_= -999.;
  eminuspzpfnow_=-999.;
  etexpopluspfnow_=-999.;
  etexpominuspfnow_=-999.;
  lrgPFnow_=-999.;
  mxpfnowleft_ =-999.;
  mxpfnowright_=-999.;
  mx2pfnowleft_=-999.;
  mx2pfnowright_=-999.;
  sumptpfnowleft_=-999.;
  sumptpfnowright_=-999.;

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

  patMuon1TrackerHits_ = -999.;
  patMuon1PixelHits_ = -999.;
  patMuon1NormalizedChi2_ = -999.;
  patMuon1MatchedStations_ = -999.;
  patMuon1Dxy_ = -999.;
  patMuon1IsGlobal_ = false;
  patMuon1IsTracker_ = false;

  patBosonMuonMass_ = -999.;
  patBosonMuonPt_ = -999.;
  patBosonMuonEta_ = -999.;
  patBosonMuonPhi_ = -999.;

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

  patBosonElectronMass_ = -999.;
  patBosonElectronPt_ = -999.;
  patBosonElectronPhi_ = -999.;
  patBosonElectronEta_ = -999.;

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
  LeadingMuonTrackerHits_ = -999.;
  LeadingMuonPixelHits_ = -999.;
  LeadingMuonNormalizedChi2_ = -999.;
  LeadingMuonMatchedStations_ = -999.;
  LeadingMuonDxy_ = -999.;
  LeadingMuonIsGlobal_ = false;
  LeadingMuonIsTracker_ = false;

  TracksNonConeMuon03_ = -1;
  TracksNonConeElectron03_ = -1;
  TracksNonConepatMuon03_ = -1;
  TracksNonConepatElectron03_ = -1;

  TracksNonConeMuon04_ = -1;
  TracksNonConeElectron04_ = -1;
  TracksNonConepatMuon04_ = -1;
  TracksNonConepatElectron04_ = -1;

  TracksNonConeMuon05_ = -1;
  TracksNonConeElectron05_ = -1;
  TracksNonConepatMuon05_ = -1;
  TracksNonConepatElectron05_ = -1;

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

  //ZDCdigifC_.clear();

}
