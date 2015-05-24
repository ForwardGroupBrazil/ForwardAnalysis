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

  SecondElectronPt_ = -999.;
  SecondElectronEta_ = -999.;
  SecondElectronPhi_ = -999.;
  SecondElectronCharge_ = -999;


  SecondElectronTkDr03_ = -999.;
  SecondElectronEcalDr03_ = -999.;
  SecondElectronHcalDr03_= -999.;

  SecondElectronTkDr04_= -999.;
  SecondElectronEcalDr04_= -999.;
  SecondElectronHcalDr04_= -999.;

  SecondElectronrelIsoDr03_= -999.;
  SecondElectronrelIsoDr04_= -999.;

  SecondElectronDeltaPhiTkClu_= -999.;
  SecondElectronDeltaEtaTkClu_= -999.;
  SecondElectronSigmaIeIe_= -999.;
  SecondElectronDCot_= -999.;
  SecondElectronDist_= -999.;
  SecondElectronInnerHits_= -999.;
  SecondElectronHE_= -999.;
  SecondElectronIsWP95_= false;
  SecondElectronIsWP80_= false;

  SecondMuonPt_= -999.;
  SecondMuonEta_= -999.;
  SecondMuonPhi_= -999.;
  SecondMuonCharge_= -999;

  SecondMuonSumPtR03_= -999.;
  SecondMuonEmEtR03_= -999.;
  SecondMuonHadEtR03_= -999.;
  SecondMuonSumPtR05_= -999.;
  SecondMuonEmEtR05_= -999.;
  SecondMuonHadEtR05_= -999.;

  SecondMuonrelIsoDr03_= -999.;
  SecondMuonrelIsoDr05_= -999.;
  SecondMuonTrackerHits_= -999.;
  SecondMuonPixelHits_= -999.;
  SecondMuonNormalizedChi2_= -999.;
  SecondMuonMatchedStations_= -999.;
  SecondMuonDxy_= -999.;
  SecondMuonDz_= -999.;
  SecondMuonIsGlobal_= false;
  SecondMuonIsTracker_ = false;
  SecondMuonIsGood_ = false;

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
  castorparticlespdgidgen_.clear();
  castorparticlesenergygen_.clear();
  castorparticlesetagen_.clear();
  castornparticlesgen_ = -1;
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

  mxGen_ = -999.;
  mx2Gen_ = -999.;
  epluspzGenLim_ = -999.;
  eminuspzGenLim_ = -999.;
  etexpoplusGenLim_ = -999.;
  etexpominusGenLim_ = -999.;

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
  patMuon1Dz_= -999.;
  patMuon1IsGlobal_ = false;
  patMuon1IsTracker_ = false;
  patMuon1IsGood_ = false;

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

  patMuon2Pt_= -999.;
  patMuon2Charge_= -999;
  patMuon2Phi_= -999.;
  patMuon2Eta_= -999.;
  patMuon2Et_= -999.;

  patMuon2SumPtR03_= -999.;
  patMuon2EmEtR03_= -999.;
  patMuon2HadEtR03_= -999.;
  patMuon2SumPtR05_= -999.;
  patMuon2EmEtR05_= -999.;
  patMuon2HadEtR05_= -999.;

  patMuon2relIsoDr03_= -999.;
  patMuon2relIsoDr05_= -999.;
  patMuon2relIso_= -999.;

  patMuon2TrackerHits_= -999.;
  patMuon2PixelHits_= -999.;
  patMuon2NormalizedChi2_= -999.;
  patMuon2MatchedStations_= -999.;
  patMuon2Dxy_ = -999.;
  patMuon2Dz_= -999.;
  patMuon2IsGlobal_ = false;
  patMuon2IsTracker_ = false;
  patMuon2IsGood_ = false;

  patElectron2Pt_= -999.;
  patElectron2Charge_= -999;
  patElectron2Phi_= -999.;
  patElectron2Eta_= -999.;
  patElectron2Et_= -999.;

  patElectron2TkDr03_= -999.;
  patElectron2EcalDr03_= -999.;
  patElectron2HcalDr03_= -999.;

  patElectron2TkDr04_= -999.;
  patElectron2EcalDr04_= -999.;
  patElectron2HcalDr04_= -999.;

  patElectron2relIsoDr03_= -999.;
  patElectron2relIsoDr04_= -999.;

  patElectron2DeltaPhiTkClu_= -999.;
  patElectron2DeltaEtaTkClu_= -999.;
  patElectron2SigmaIeIe_= -999.;
  patElectron2DCot_= -999.;
  patElectron2Dist_= -999.;
  patElectron2InnerHits_= -999.;
  patElectron2HE_= -999.;
  patElectron2IsWP95_= false;
  patElectron2IsWP80_ = false;

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
  LeadingMuonDz_ = -999.;
  LeadingMuonIsGlobal_ = false;
  LeadingMuonIsTracker_ = false;
  LeadingMuonIsGood_ = false;

  TracksNonConeSecondMuonR03_=-1;
  TracksNonConeSecondElectronR03_=-1;
  TracksNonConePatMuon2R03_=-1;
  TracksNonConePatElectron2R03_=-1;

  TracksNonConeSecondMuonR04_=-1;
  TracksNonConeSecondElectronR04_=-1;
  TracksNonConePatMuon2R04_=-1;
  TracksNonConePatElectron2R04_=-1;

  TracksNonConeSecondMuonR05_=-1;
  TracksNonConeSecondElectronR05_=-1;
  TracksNonConePatMuon2R05_=-1;
  TracksNonConePatElectron2R05_=-1;

  TracksNonConeSecondMuonR03_=-1;
  TracksNonConeSecondElectronR03_=-1;
  TracksNonConePatMuon2R03_=-1;
  TracksNonConePatElectron2R03_=-1;

  TracksNonConeSecondMuonR04_=-1;
  TracksNonConeSecondElectronR04_=-1;
  TracksNonConePatMuon2R04_=-1;
  TracksNonConePatElectron2R04_=-1;

  TracksNonConeSecondMuonR05_=-1;
  TracksNonConeSecondElectronR05_=-1;
  TracksNonConePatMuon2R05_=-1;
  TracksNonConePatElectron2R05_=-1;

  LeadingElectronDeltaPhiTkClu_ = -999.;
  LeadingElectronDeltaEtaTkClu_ = -999.;
  LeadingElectronSigmaIeIe_ = -999.;
  LeadingElectronDCot_ = -999.;
  LeadingElectronDist_ = -999.;
  LeadingElectronInnerHits_ = -999.;
  LeadingElectronHE_ = -999.;
  LeadingElectronIsWP95_ = false;
  LeadingElectronIsWP80_ = false;

  patElectron1DeltaPhiTkClu_ = -999.;
  patElectron1DeltaEtaTkClu_ = -999.;
  patElectron1SigmaIeIe_ = -999.;
  patElectron1DCot_ = -999.;
  patElectron1Dist_ = -999.;
  patElectron1InnerHits_ = -999.;
  patElectron1HE_ = -999.;
  patElectron1IsWP95_ = false;
  patElectron1IsWP80_ = false;

  metPt_ = -999.;
  metPhi_ = -999.;
  metEt_ = -999.;
  metSumEt_ = -999.;
  metpx_ = -999.;
  metpy_ = -999.;
  metsigma_ = -999.;

  patmetPt_ = -999.;
  patmetPhi_ = -999.;
  patmetEt_ = -999.;
  patmetSumEt_ = -999.;
  patmetpx_ = -999.;
  patmetpy_ = -999.;
  patmetsigma_ = -999.;

  dvtxmuon_ = -999.;
  dvtxmuonZ_ = -999.;
  dvtxelectron_ = -999.;
  dvtxelectronZ_ = -999.;
  dmuonelectron_ = -999.;
  dmuonelectronZ_ = -999.;
  dmuons_ = -999.;
  dmuonsZ_ = -999.;
  delectrons_ = -999.;
  delectronsZ_ = -999.;

  GenLeadingElectronPt_ = -999.;
  GenLeadingElectronEta_ = -999.;
  GenLeadingElectronPhi_ = -999.;

  GenLeadingMuonPt_ = -999.;
  GenLeadingMuonEta_ = -999.;
  GenLeadingMuonPhi_ = -999.;

  GenNeutrinoPt_ = -999.;
  GenNeutrinoPhi_ = -999.;
  GenNeutrinoPx_ = -999.;
  GenNeutrinoPy_ = -999.;

  SumEHFPluspf_ = -999.;
  SumEHFMinuspf_ = -999.;
  mxpf_ = -999.;
  mx2pf_ = -999.;
  mxpfnow_ = -999.;
  mx2pfnow_ = -999.;

  //ZDCdigifC_.clear();

}
