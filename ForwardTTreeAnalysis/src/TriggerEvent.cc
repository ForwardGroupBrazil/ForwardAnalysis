#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/TriggerEvent.h"
#include <cstdio>

const char* TriggerEvent::name = "TriggerEvent";

TriggerEvent::TriggerEvent() {}

TriggerEvent::~TriggerEvent() {}

void TriggerEvent::reset(){

  size_t len_hltTrigResults = sizeof(hltTrigResults_)/sizeof(int);
  for(size_t k = 0; k < len_hltTrigResults; ++k) hltTrigResults_[k] = 0;

  vtx_ = -999;
  ntrk_ = -999;

  LeadingElectronPt_=-999.;
  LeadingElectronEta_=-999.;
  LeadingElectronPhi_=-999.;
  LeadingElectronCharge_=-999.;
  ElectronsN_= -1;

  LeadingMuonPt_=-999.;
  LeadingMuonEta_=-999.;
  LeadingMuonPhi_=-999.;
  LeadingMuonCharge_=-999.;
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

  TracksNonConeSecondMuonR04_=-1;
  TracksNonConeSecondElectronR04_=-1;

  TracksNonConeSecondMuonR05_=-1;
  TracksNonConeSecondElectronR05_=-1;

  TracksNonConeSecondMuonR03_=-1;
  TracksNonConeSecondElectronR03_=-1;

  TracksNonConeSecondMuonR04_=-1;
  TracksNonConeSecondElectronR04_=-1;

  TracksNonConeSecondMuonR05_=-1;
  TracksNonConeSecondElectronR05_=-1;

  LeadingElectronDeltaPhiTkClu_ = -999.;
  LeadingElectronDeltaEtaTkClu_ = -999.;
  LeadingElectronSigmaIeIe_ = -999.;
  LeadingElectronDCot_ = -999.;
  LeadingElectronDist_ = -999.;
  LeadingElectronInnerHits_ = -999.;
  LeadingElectronHE_ = -999.;
  LeadingElectronIsWP95_ = false;
  LeadingElectronIsWP80_ = false;

  metPt_ = -999.;
  metPhi_ = -999.;
  metEt_ = -999.;
  metSumEt_ = -999.;
  metpx_ = -999.;
  metpy_ = -999.;
  metsigma_ = -999.;

}
