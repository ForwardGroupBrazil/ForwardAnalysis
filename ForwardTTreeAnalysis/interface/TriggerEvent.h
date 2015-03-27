#ifndef TriggerEvent_h
#define TriggerEvent_h

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/DetSetVector.h"

namespace triggerAnalysis {
  class TriggerAnalysis;
}

class TriggerEvent {
  public:
    typedef triggerAnalysis::TriggerAnalysis analysis_type;
    static const char* name;

    typedef reco::Particle::LorentzVector LorentzVector;

    TriggerEvent();
    ~TriggerEvent();

    void SetNVertex(int fvtx)    {vtx_ = fvtx;}
    void SetMultiplicityTracks(int fntrk) {ntrk_ = fntrk;}
    void SetHLTPath(int idx, int fHLTBit)         { hltTrigResults_[idx] = fHLTBit;}

    void SetMETPt(double fmetPt)    {metPt_ = fmetPt;}
    void SetMETPhi(double fmetPhi)    {metPhi_ = fmetPhi;}
    void SetMETEt(double fmetEt)    {metEt_ = fmetEt;}
    void SetMETSumEt(double fmetSumEt)    {metSumEt_ = fmetSumEt;}
    void SetMETpx(double fmetpx)    {metpx_ = fmetpx;}
    void SetMETpy(double fmetpy)    {metpy_ = fmetpy;}
    void SetMETP4(LorentzVector fmetp4)    {metp4_     = fmetp4;}
    void SetMETsigma(double fmetsigma)    {metsigma_     = fmetsigma;}

    void SetLeadingElectronPt(double fLeadingElectronPt)    { LeadingElectronPt_     = fLeadingElectronPt;}
    void SetLeadingElectronEta(double fLeadingElectronEta)  { LeadingElectronEta_     = fLeadingElectronEta;}
    void SetLeadingElectronPhi(double fLeadingElectronPhi)  { LeadingElectronPhi_    = fLeadingElectronPhi;}
    void SetLeadingElectronP4(LorentzVector fLeadingElectronP4)    { LeadingElectronP4_     = fLeadingElectronP4;}
    void SetLeadingElectronCharge(int fLeadingElectronCharge)  { LeadingElectronCharge_     = fLeadingElectronCharge;}
    void SetElectronsN(int fElectronsN)  { ElectronsN_    = fElectronsN;}

    void SetLeadingElectronTkDr03(double fLeadingElectronTkDr03)    {LeadingElectronTkDr03_ = fLeadingElectronTkDr03;}
    void SetLeadingElectronEcalDr03(double fLeadingElectronEcalDr03)    {LeadingElectronEcalDr03_ = fLeadingElectronEcalDr03;}
    void SetLeadingElectronHcalDr03(double fLeadingElectronHcalDr03)    {LeadingElectronHcalDr03_ = fLeadingElectronHcalDr03;}

    void SetLeadingElectronTkDr04(double fLeadingElectronTkDr04)    {LeadingElectronTkDr04_ = fLeadingElectronTkDr04;}
    void SetLeadingElectronEcalDr04(double fLeadingElectronEcalDr04)    {LeadingElectronEcalDr04_ = fLeadingElectronEcalDr04;}
    void SetLeadingElectronHcalDr04(double fLeadingElectronHcalDr04)    {LeadingElectronHcalDr04_ = fLeadingElectronHcalDr04;}

    void SetLeadingElectronrelIsoDr03(double fLeadingElectronrelIsoDr03)    {LeadingElectronrelIsoDr03_ = fLeadingElectronrelIsoDr03;}
    void SetLeadingElectronrelIsoDr04(double fLeadingElectronrelIsoDr04)    {LeadingElectronrelIsoDr04_ = fLeadingElectronrelIsoDr04;}

    void SetLeadingMuonPt(double fLeadingMuonPt)    { LeadingMuonPt_     = fLeadingMuonPt;}
    void SetLeadingMuonEta(double fLeadingMuonEta)  { LeadingMuonEta_     = fLeadingMuonEta;}
    void SetLeadingMuonPhi(double fLeadingMuonPhi)  { LeadingMuonPhi_    = fLeadingMuonPhi;}
    void SetLeadingMuonP4(LorentzVector fLeadingMuonP4)    { LeadingMuonP4_     = fLeadingMuonP4;}
    void SetLeadingMuonCharge(int fLeadingMuonCharge)  { LeadingMuonCharge_     = fLeadingMuonCharge;}
    void SetMuonsN(int fMuonsN)  { MuonsN_    = fMuonsN;}

    void SetLeadingMuonSumPtR03(double fLeadingMuonSumPtR03)    {LeadingMuonSumPtR03_ = fLeadingMuonSumPtR03;}
    void SetLeadingMuonEmEtR03(double fLeadingMuonEmEtR03)    {LeadingMuonEmEtR03_ = fLeadingMuonEmEtR03;}
    void SetLeadingMuonHadEtR03(double fLeadingMuonHadEtR03)    {LeadingMuonHadEtR03_ = fLeadingMuonHadEtR03;}
    void SetLeadingMuonSumPtR05(double fLeadingMuonSumPtR05)    {LeadingMuonSumPtR05_ = fLeadingMuonSumPtR05;}
    void SetLeadingMuonEmEtR05(double fLeadingMuonEmEtR05)    {LeadingMuonEmEtR05_ = fLeadingMuonEmEtR05;}
    void SetLeadingMuonHadEtR05(double fLeadingMuonHadEtR05)    {LeadingMuonHadEtR05_ = fLeadingMuonHadEtR05;}

    void SetLeadingMuonrelIsoDr03(double fLeadingMuonrelIsoDr03)    {LeadingMuonrelIsoDr03_ = fLeadingMuonrelIsoDr03;}
    void SetLeadingMuonrelIsoDr05(double fLeadingMuonrelIsoDr05)    {LeadingMuonrelIsoDr05_ = fLeadingMuonrelIsoDr05;}
    void SetLeadingMuonTrackerHits(double fLeadingMuonTrackerHits)    {LeadingMuonTrackerHits_ = fLeadingMuonTrackerHits;}
    void SetLeadingMuonPixelHits(double fLeadingMuonPixelHits)    {LeadingMuonPixelHits_ = fLeadingMuonPixelHits;}
    void SetLeadingMuonNormalizedChi2(double fLeadingMuonNormalizedChi2)    {LeadingMuonNormalizedChi2_ = fLeadingMuonNormalizedChi2;}
    void SetLeadingMuonMatchedStations(double fLeadingMuonMatchedStations)    {LeadingMuonMatchedStations_ = fLeadingMuonMatchedStations;}
    void SetLeadingMuonDxy(double fLeadingMuonDxy)    {LeadingMuonDxy_ = fLeadingMuonDxy;}
    void SetLeadingMuonDz(double fLeadingMuonDz)    {LeadingMuonDz_ = fLeadingMuonDz;}
    void SetLeadingMuonIsGlobal(bool fLeadingMuonIsGlobal)    {LeadingMuonIsGlobal_ = fLeadingMuonIsGlobal;}
    void SetLeadingMuonIsTracker(bool fLeadingMuonIsTracker)    {LeadingMuonIsTracker_ = fLeadingMuonIsTracker;}
    void SetLeadingMuonIsGood(bool fLeadingMuonIsGood)    {LeadingMuonIsGood_ = fLeadingMuonIsGood;}

    void SetSecondElectronPt(double fSecondElectronPt)    { SecondElectronPt_     = fSecondElectronPt;}
    void SetSecondElectronEta(double fSecondElectronEta)  { SecondElectronEta_     = fSecondElectronEta;}
    void SetSecondElectronPhi(double fSecondElectronPhi)  { SecondElectronPhi_    = fSecondElectronPhi;}
    void SetSecondElectronP4(LorentzVector fSecondElectronP4)    { SecondElectronP4_     = fSecondElectronP4;}
    void SetSecondElectronCharge(int fSecondElectronCharge)  { SecondElectronCharge_     = fSecondElectronCharge;}

    void SetSecondElectronTkDr03(double fSecondElectronTkDr03)    {SecondElectronTkDr03_ = fSecondElectronTkDr03;}
    void SetSecondElectronEcalDr03(double fSecondElectronEcalDr03)    {SecondElectronEcalDr03_ = fSecondElectronEcalDr03;}
    void SetSecondElectronHcalDr03(double fSecondElectronHcalDr03)    {SecondElectronHcalDr03_ = fSecondElectronHcalDr03;}

    void SetSecondElectronTkDr04(double fSecondElectronTkDr04)    {SecondElectronTkDr04_ = fSecondElectronTkDr04;}
    void SetSecondElectronEcalDr04(double fSecondElectronEcalDr04)    {SecondElectronEcalDr04_ = fSecondElectronEcalDr04;}
    void SetSecondElectronHcalDr04(double fSecondElectronHcalDr04)    {SecondElectronHcalDr04_ = fSecondElectronHcalDr04;}

    void SetSecondElectronrelIsoDr03(double fSecondElectronrelIsoDr03)    {SecondElectronrelIsoDr03_ = fSecondElectronrelIsoDr03;}
    void SetSecondElectronrelIsoDr04(double fSecondElectronrelIsoDr04)    {SecondElectronrelIsoDr04_ = fSecondElectronrelIsoDr04;}

    void SetSecondElectronDeltaPhiTkClu(double fSecondElectronDeltaPhiTkClu)    {SecondElectronDeltaPhiTkClu_ = fSecondElectronDeltaPhiTkClu;}
    void SetSecondElectronDeltaEtaTkClu(double fSecondElectronDeltaEtaTkClu)    {SecondElectronDeltaEtaTkClu_ = fSecondElectronDeltaEtaTkClu;}
    void SetSecondElectronSigmaIeIe(double fSecondElectronSigmaIeIe)    {SecondElectronSigmaIeIe_ = fSecondElectronSigmaIeIe;}
    void SetSecondElectronDCot(double fSecondElectronDCot)    {SecondElectronDCot_ = fSecondElectronDCot;}
    void SetSecondElectronDist(double fSecondElectronDist)    {SecondElectronDist_ = fSecondElectronDist;}
    void SetSecondElectronInnerHits(double fSecondElectronInnerHits)    {SecondElectronInnerHits_ = fSecondElectronInnerHits;}
    void SetSecondElectronHE(double fSecondElectronHE)    {SecondElectronHE_ = fSecondElectronHE;}
    void SetSecondElectronIsWP95(bool fSecondElectronIsWP95)    {SecondElectronIsWP95_ = fSecondElectronIsWP95;}
    void SetSecondElectronIsWP80(bool fSecondElectronIsWP80)    {SecondElectronIsWP80_ = fSecondElectronIsWP80;}

    void SetSecondMuonPt(double fSecondMuonPt)    { SecondMuonPt_     = fSecondMuonPt;}
    void SetSecondMuonEta(double fSecondMuonEta)  { SecondMuonEta_     = fSecondMuonEta;}
    void SetSecondMuonPhi(double fSecondMuonPhi)  { SecondMuonPhi_    = fSecondMuonPhi;}
    void SetSecondMuonP4(LorentzVector fSecondMuonP4)    { SecondMuonP4_     = fSecondMuonP4;}
    void SetSecondMuonCharge(int fSecondMuonCharge)  { SecondMuonCharge_     = fSecondMuonCharge;}

    void SetSecondMuonSumPtR03(double fSecondMuonSumPtR03)    {SecondMuonSumPtR03_ = fSecondMuonSumPtR03;}
    void SetSecondMuonEmEtR03(double fSecondMuonEmEtR03)    {SecondMuonEmEtR03_ = fSecondMuonEmEtR03;}
    void SetSecondMuonHadEtR03(double fSecondMuonHadEtR03)    {SecondMuonHadEtR03_ = fSecondMuonHadEtR03;}
    void SetSecondMuonSumPtR05(double fSecondMuonSumPtR05)    {SecondMuonSumPtR05_ = fSecondMuonSumPtR05;}
    void SetSecondMuonEmEtR05(double fSecondMuonEmEtR05)    {SecondMuonEmEtR05_ = fSecondMuonEmEtR05;}
    void SetSecondMuonHadEtR05(double fSecondMuonHadEtR05)    {SecondMuonHadEtR05_ = fSecondMuonHadEtR05;}

    void SetSecondMuonrelIsoDr03(double fSecondMuonrelIsoDr03)    {SecondMuonrelIsoDr03_ = fSecondMuonrelIsoDr03;}
    void SetSecondMuonrelIsoDr05(double fSecondMuonrelIsoDr05)    {SecondMuonrelIsoDr05_ = fSecondMuonrelIsoDr05;}
    void SetSecondMuonTrackerHits(double fSecondMuonTrackerHits)    {SecondMuonTrackerHits_ = fSecondMuonTrackerHits;}
    void SetSecondMuonPixelHits(double fSecondMuonPixelHits)    {SecondMuonPixelHits_ = fSecondMuonPixelHits;}
    void SetSecondMuonNormalizedChi2(double fSecondMuonNormalizedChi2)    {SecondMuonNormalizedChi2_ = fSecondMuonNormalizedChi2;}
    void SetSecondMuonMatchedStations(double fSecondMuonMatchedStations)    {SecondMuonMatchedStations_ = fSecondMuonMatchedStations;}
    void SetSecondMuonDxy(double fSecondMuonDxy)    {SecondMuonDxy_ = fSecondMuonDxy;}
    void SetSecondMuonDz(double fSecondMuonDz)    {SecondMuonDz_ = fSecondMuonDz;}
    void SetSecondMuonIsGlobal(bool fSecondMuonIsGlobal)    {SecondMuonIsGlobal_ = fSecondMuonIsGlobal;}
    void SetSecondMuonIsTracker(bool fSecondMuonIsTracker)    {SecondMuonIsTracker_ = fSecondMuonIsTracker;}
    void SetSecondMuonIsGood(bool fSecondMuonIsGood)    {SecondMuonIsGood_ = fSecondMuonIsGood;}

    void SetTracksNonConeLeadingMuonR03(int fTracksNonConeLeadingMuonR03)    {TracksNonConeLeadingMuonR03_ = fTracksNonConeLeadingMuonR03;}
    void SetTracksNonConeLeadingElectronR03(int fTracksNonConeLeadingElectronR03)    {TracksNonConeLeadingElectronR03_ = fTracksNonConeLeadingElectronR03;}

    void SetTracksNonConeLeadingMuonR04(int fTracksNonConeLeadingMuonR04)    {TracksNonConeLeadingMuonR04_ = fTracksNonConeLeadingMuonR04;}
    void SetTracksNonConeLeadingElectronR04(int fTracksNonConeLeadingElectronR04)    {TracksNonConeLeadingElectronR04_ = fTracksNonConeLeadingElectronR04;}

    void SetTracksNonConeLeadingMuonR05(int fTracksNonConeLeadingMuonR05)    {TracksNonConeLeadingMuonR05_ = fTracksNonConeLeadingMuonR05;}
    void SetTracksNonConeLeadingElectronR05(int fTracksNonConeLeadingElectronR05)    {TracksNonConeLeadingElectronR05_ = fTracksNonConeLeadingElectronR05;}

    void SetTracksNonConeSecondMuonR03(int fTracksNonConeSecondMuonR03)    {TracksNonConeSecondMuonR03_ = fTracksNonConeSecondMuonR03;}
    void SetTracksNonConeSecondElectronR03(int fTracksNonConeSecondElectronR03)    {TracksNonConeSecondElectronR03_ = fTracksNonConeSecondElectronR03;}

    void SetTracksNonConeSecondMuonR04(int fTracksNonConeSecondMuonR04)    {TracksNonConeSecondMuonR04_ = fTracksNonConeSecondMuonR04;}
    void SetTracksNonConeSecondElectronR04(int fTracksNonConeSecondElectronR04)    {TracksNonConeSecondElectronR04_ = fTracksNonConeSecondElectronR04;}

    void SetTracksNonConeSecondMuonR05(int fTracksNonConeSecondMuonR05)    {TracksNonConeSecondMuonR05_ = fTracksNonConeSecondMuonR05;}
    void SetTracksNonConeSecondElectronR05(int fTracksNonConeSecondElectronR05)    {TracksNonConeSecondElectronR05_ = fTracksNonConeSecondElectronR05;}

    void SetLeadingElectronDeltaPhiTkClu(double fLeadingElectronDeltaPhiTkClu)    {LeadingElectronDeltaPhiTkClu_ = fLeadingElectronDeltaPhiTkClu;}
    void SetLeadingElectronDeltaEtaTkClu(double fLeadingElectronDeltaEtaTkClu)    {LeadingElectronDeltaEtaTkClu_ = fLeadingElectronDeltaEtaTkClu;}
    void SetLeadingElectronSigmaIeIe(double fLeadingElectronSigmaIeIe)    {LeadingElectronSigmaIeIe_ = fLeadingElectronSigmaIeIe;}
    void SetLeadingElectronDCot(double fLeadingElectronDCot)    {LeadingElectronDCot_ = fLeadingElectronDCot;}
    void SetLeadingElectronDist(double fLeadingElectronDist)    {LeadingElectronDist_ = fLeadingElectronDist;}
    void SetLeadingElectronInnerHits(double fLeadingElectronInnerHits)    {LeadingElectronInnerHits_ = fLeadingElectronInnerHits;}
    void SetLeadingElectronHE(double fLeadingElectronHE)    {LeadingElectronHE_ = fLeadingElectronHE;}
    void SetLeadingElectronIsWP95(bool fLeadingElectronIsWP95)    {LeadingElectronIsWP95_ = fLeadingElectronIsWP95;}
    void SetLeadingElectronIsWP80(bool fLeadingElectronIsWP80)    {LeadingElectronIsWP80_ = fLeadingElectronIsWP80;}

    int GetNVertex() const {return vtx_;}
    int GetMultiplicityTracks() const {return ntrk_;}
    int GetHLTPath(int idx) const { return hltTrigResults_[idx]; }

    double GetLeadingElectronPt() const {return LeadingElectronPt_;}
    double GetLeadingElectronEta() const {return LeadingElectronEta_;}
    double GetLeadingElectronPhi() const {return LeadingElectronPhi_;}
    const LorentzVector& GetLeadingElectronP4() const {return LeadingElectronP4_;}
    int GetLeadingElectronCharge() const {return LeadingElectronCharge_;}
    int GetElectronsN() const {return ElectronsN_;}
    double GetLeadingElectronTkDr03() const {return LeadingElectronTkDr03_;}
    double GetLeadingElectronEcalDr03() const {return LeadingElectronEcalDr03_;}
    double GetLeadingElectronHcalDr03() const {return LeadingElectronHcalDr03_;}

    double GetLeadingElectronTkDr04() const {return LeadingElectronTkDr04_;}
    double GetLeadingElectronEcalDr04() const {return LeadingElectronEcalDr04_;}
    double GetLeadingElectronHcalDr04() const {return LeadingElectronHcalDr04_;}

    double GetLeadingElectronrelIsoDr03() const {return LeadingElectronrelIsoDr03_;}
    double GetLeadingElectronrelIsoDr04() const {return LeadingElectronrelIsoDr04_;}

    bool GetLeadingElectronIsWP95() const {return LeadingElectronIsWP95_;}
    bool GetLeadingElectronIsWP80() const {return LeadingElectronIsWP80_;}

    double GetLeadingMuonPt() const {return LeadingMuonPt_;}
    double GetLeadingMuonEta() const {return LeadingMuonEta_;}
    double GetLeadingMuonPhi() const {return LeadingMuonPhi_;}
    const LorentzVector& GetLeadingMuonP4() const {return LeadingMuonP4_;}
    int GetLeadingMuonCharge() const {return LeadingMuonCharge_;}

    double GetLeadingMuonSumPtR03() const {return LeadingMuonSumPtR03_;}
    double GetLeadingMuonEmEtR03() const {return LeadingMuonEmEtR03_;}
    double GetLeadingMuonHadEtR03() const {return LeadingMuonHadEtR03_;}
    double GetLeadingMuonSumPtR05() const {return LeadingMuonSumPtR05_;}
    double GetLeadingMuonEmEtR05() const {return LeadingMuonEmEtR05_;}
    double GetLeadingMuonHadEtR05() const {return LeadingMuonHadEtR05_;}

    double GetLeadingMuonrelIsoDr03() const {return LeadingMuonrelIsoDr03_;}
    double GetLeadingMuonrelIsoDr05() const {return LeadingMuonrelIsoDr05_;}

    int GetMuonsN() const {return MuonsN_;}
    double GetLeadingMuonTrackerHits() const {return LeadingMuonTrackerHits_;}
    double GetLeadingMuonPixelHits() const {return LeadingMuonPixelHits_;}
    double GetLeadingMuonNormalizedChi2() const {return LeadingMuonNormalizedChi2_;}
    double GetLeadingMuonMatchedStations() const {return LeadingMuonMatchedStations_;}
    double GetLeadingMuonDxy() const {return LeadingMuonDxy_;}
    double GetLeadingMuonDz() const {return LeadingMuonDz_;}
    bool GetLeadingMuonIsGlobal() const {return LeadingMuonIsGlobal_;}
    bool GetLeadingMuonIsTracker() const {return LeadingMuonIsTracker_;}
    bool GetLeadingMuonIsGood() const {return LeadingMuonIsGood_;}

    double GetSecondElectronPt() const {return SecondElectronPt_;}
    double GetSecondElectronEta() const {return SecondElectronEta_;}
    double GetSecondElectronPhi() const {return SecondElectronPhi_;}
    const LorentzVector& GetSecondElectronP4() const {return SecondElectronP4_;}
    int GetSecondElectronCharge() const {return SecondElectronCharge_;}
    double GetSecondElectronTkDr03() const {return SecondElectronTkDr03_;}
    double GetSecondElectronEcalDr03() const {return SecondElectronEcalDr03_;}
    double GetSecondElectronHcalDr03() const {return SecondElectronHcalDr03_;}

    double GetSecondElectronTkDr04() const {return SecondElectronTkDr04_;}
    double GetSecondElectronEcalDr04() const {return SecondElectronEcalDr04_;}
    double GetSecondElectronHcalDr04() const {return SecondElectronHcalDr04_;}

    double GetSecondElectronrelIsoDr03() const {return SecondElectronrelIsoDr03_;}
    double GetSecondElectronrelIsoDr04() const {return SecondElectronrelIsoDr04_;}

    bool GetSecondElectronIsWP95() const {return SecondElectronIsWP95_;}
    bool GetSecondElectronIsWP80() const {return SecondElectronIsWP80_;}

    double GetSecondElectronDeltaPhiTkClu() const {return SecondElectronDeltaPhiTkClu_;}
    double GetSecondElectronDeltaEtaTkClu() const {return SecondElectronDeltaEtaTkClu_;}
    double GetSecondElectronSigmaIeIe() const {return SecondElectronSigmaIeIe_;}
    double GetSecondElectronDCot() const {return SecondElectronDCot_;}
    double GetSecondElectronDist() const {return SecondElectronDist_;}
    double GetSecondElectronInnerHits() const {return SecondElectronInnerHits_;}
    double GetSecondElectronHE() const {return SecondElectronHE_;}

    double GetSecondMuonPt() const {return SecondMuonPt_;}
    double GetSecondMuonEta() const {return SecondMuonEta_;}
    double GetSecondMuonPhi() const {return SecondMuonPhi_;}
    const LorentzVector& GetSecondMuonP4() const {return SecondMuonP4_;}
    int GetSecondMuonCharge() const {return SecondMuonCharge_;}

    double GetSecondMuonSumPtR03() const {return SecondMuonSumPtR03_;}
    double GetSecondMuonEmEtR03() const {return SecondMuonEmEtR03_;}
    double GetSecondMuonHadEtR03() const {return SecondMuonHadEtR03_;}
    double GetSecondMuonSumPtR05() const {return SecondMuonSumPtR05_;}
    double GetSecondMuonEmEtR05() const {return SecondMuonEmEtR05_;}
    double GetSecondMuonHadEtR05() const {return SecondMuonHadEtR05_;}

    double GetSecondMuonrelIsoDr03() const {return SecondMuonrelIsoDr03_;}
    double GetSecondMuonrelIsoDr05() const {return SecondMuonrelIsoDr05_;}

    double GetSecondMuonTrackerHits() const {return SecondMuonTrackerHits_;}
    double GetSecondMuonPixelHits() const {return SecondMuonPixelHits_;}
    double GetSecondMuonNormalizedChi2() const {return SecondMuonNormalizedChi2_;}
    double GetSecondMuonMatchedStations() const {return SecondMuonMatchedStations_;}
    double GetSecondMuonDxy() const {return SecondMuonDxy_;}
    double GetSecondMuonDz() const {return SecondMuonDz_;}
    bool GetSecondMuonIsGlobal() const {return SecondMuonIsGlobal_;}
    bool GetSecondMuonIsTracker() const {return SecondMuonIsTracker_;}
    bool GetSecondMuonIsGood() const {return SecondMuonIsGood_;}

    int GetTracksNonConeLeadingMuonR03() const {return TracksNonConeLeadingMuonR03_;}
    int GetTracksNonConeLeadingElectronR03() const {return TracksNonConeLeadingElectronR03_;}

    int GetTracksNonConeLeadingMuonR04() const {return TracksNonConeLeadingMuonR04_;}
    int GetTracksNonConeLeadingElectronR04() const {return TracksNonConeLeadingElectronR04_;}

    int GetTracksNonConeLeadingMuonR05() const {return TracksNonConeLeadingMuonR05_;}
    int GetTracksNonConeLeadingElectronR05() const {return TracksNonConeLeadingElectronR05_;}

    int GetTracksNonConeSecondMuonR03() const {return TracksNonConeSecondMuonR03_;}
    int GetTracksNonConeSecondElectronR03() const {return TracksNonConeSecondElectronR03_;}

    int GetTracksNonConeSecondMuonR04() const {return TracksNonConeSecondMuonR04_;}
    int GetTracksNonConeSecondElectronR04() const {return TracksNonConeSecondElectronR04_;}

    int GetTracksNonConeSecondMuonR05() const {return TracksNonConeSecondMuonR05_;}
    int GetTracksNonConeSecondElectronR05() const {return TracksNonConeSecondElectronR05_;}

    double GetLeadingElectronDeltaPhiTkClu() const {return LeadingElectronDeltaPhiTkClu_;}
    double GetLeadingElectronDeltaEtaTkClu() const {return LeadingElectronDeltaEtaTkClu_;}
    double GetLeadingElectronSigmaIeIe() const {return LeadingElectronSigmaIeIe_;}
    double GetLeadingElectronDCot() const {return LeadingElectronDCot_;}
    double GetLeadingElectronDist() const {return LeadingElectronDist_;}
    double GetLeadingElectronInnerHits() const {return LeadingElectronInnerHits_;}
    double GetLeadingElectronHE() const {return LeadingElectronHE_;}

    double GetMETPt() const {return metPt_;}
    double GetMETPhi() const {return metPhi_;}
    double GetMETEt() const {return metEt_;}
    double GetMETSumEt() const {return metSumEt_;}
    double GetMETpx() const {return metpx_;}
    double GetMETpy() const {return metpy_;}
    const LorentzVector& GetMETP4() const {return metp4_;}
    double GetMETsigma() const {return metsigma_;}

  private:
    friend class triggerAnalysis::TriggerAnalysis;

    void reset();

    int vtx_;
    int ntrk_;
    int hltTrigResults_[20];

    double LeadingElectronPt_;
    double LeadingElectronEta_;
    double LeadingElectronPhi_;
    LorentzVector LeadingElectronP4_;
    int LeadingElectronCharge_;
    int ElectronsN_;

    double LeadingElectronTkDr03_;
    double LeadingElectronEcalDr03_;
    double LeadingElectronHcalDr03_;

    double LeadingElectronTkDr04_;
    double LeadingElectronEcalDr04_;
    double LeadingElectronHcalDr04_;

    double LeadingElectronrelIsoDr03_;
    double LeadingElectronrelIsoDr04_;

    double LeadingMuonPt_;
    double LeadingMuonEta_;
    double LeadingMuonPhi_;
    LorentzVector LeadingMuonP4_;
    int LeadingMuonCharge_;
    int MuonsN_;

    double LeadingMuonSumPtR03_;
    double LeadingMuonEmEtR03_;
    double LeadingMuonHadEtR03_;
    double LeadingMuonSumPtR05_;
    double LeadingMuonEmEtR05_;
    double LeadingMuonHadEtR05_;

    double LeadingMuonrelIsoDr03_;
    double LeadingMuonrelIsoDr05_;
    double LeadingMuonTrackerHits_;
    double LeadingMuonPixelHits_;
    double LeadingMuonNormalizedChi2_;
    double LeadingMuonMatchedStations_;
    double LeadingMuonDxy_;
    double LeadingMuonDz_;
    bool LeadingMuonIsGlobal_;
    bool LeadingMuonIsTracker_;
    bool LeadingMuonIsGood_;

    double SecondElectronPt_;
    double SecondElectronEta_;
    double SecondElectronPhi_;
    LorentzVector SecondElectronP4_;
    int SecondElectronCharge_;

    double SecondElectronTkDr03_;
    double SecondElectronEcalDr03_;
    double SecondElectronHcalDr03_;

    double SecondElectronTkDr04_;
    double SecondElectronEcalDr04_;
    double SecondElectronHcalDr04_;

    double SecondElectronrelIsoDr03_;
    double SecondElectronrelIsoDr04_;

    double SecondElectronDeltaPhiTkClu_;
    double SecondElectronDeltaEtaTkClu_;
    double SecondElectronSigmaIeIe_;
    double SecondElectronDCot_;
    double SecondElectronDist_;
    double SecondElectronInnerHits_;
    double SecondElectronHE_;
    bool SecondElectronIsWP95_;
    bool SecondElectronIsWP80_;

    double SecondMuonPt_;
    double SecondMuonEta_;
    double SecondMuonPhi_;
    LorentzVector SecondMuonP4_;
    int SecondMuonCharge_;

    double SecondMuonSumPtR03_;
    double SecondMuonEmEtR03_;
    double SecondMuonHadEtR03_;
    double SecondMuonSumPtR05_;
    double SecondMuonEmEtR05_;
    double SecondMuonHadEtR05_;

    double SecondMuonrelIsoDr03_;
    double SecondMuonrelIsoDr05_;
    double SecondMuonTrackerHits_;
    double SecondMuonPixelHits_;
    double SecondMuonNormalizedChi2_;
    double SecondMuonMatchedStations_;
    double SecondMuonDxy_;
    double SecondMuonDz_;
    bool SecondMuonIsGlobal_;
    bool SecondMuonIsTracker_;
    bool SecondMuonIsGood_;

    int TracksNonConeLeadingMuonR03_;
    int TracksNonConeLeadingElectronR03_;

    int TracksNonConeLeadingMuonR04_;
    int TracksNonConeLeadingElectronR04_;

    int TracksNonConeLeadingMuonR05_;
    int TracksNonConeLeadingElectronR05_;

    int TracksNonConeSecondMuonR03_;
    int TracksNonConeSecondElectronR03_;

    int TracksNonConeSecondMuonR04_;
    int TracksNonConeSecondElectronR04_;

    int TracksNonConeSecondMuonR05_;
    int TracksNonConeSecondElectronR05_;

    double LeadingElectronDeltaPhiTkClu_;
    double LeadingElectronDeltaEtaTkClu_;
    double LeadingElectronSigmaIeIe_;
    double LeadingElectronDCot_;
    double LeadingElectronDist_;
    double LeadingElectronInnerHits_;
    double LeadingElectronHE_;
    bool LeadingElectronIsWP95_;
    bool LeadingElectronIsWP80_;

    double metPt_;
    double metPhi_;
    double metEt_;
    double metSumEt_;
    double metpx_;
    double metpy_;
    LorentzVector metp4_;
    double metsigma_;

};

#endif    
