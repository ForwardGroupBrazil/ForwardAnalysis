#ifndef DiffractiveWEvent_h
#define DiffractiveWEvent_h

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/DetSetVector.h"

namespace diffractiveWAnalysis {
  class DiffractiveWAnalysis;
}

class DiffractiveWEvent {
  public:
    typedef diffractiveWAnalysis::DiffractiveWAnalysis analysis_type;
    static const char* name;

    typedef reco::Particle::LorentzVector LorentzVector;

    DiffractiveWEvent();
    ~DiffractiveWEvent();

    void SetHLTPath(int idx, int fHLTBit)         { hltTrigResults_[idx] = fHLTBit;}
    void SetBosonElectronMass(double fBosonElectronMass) { BosonElectronMass_ = fBosonElectronMass;}
    void SetBosonElectronPt(double fBosonElectronPt) { BosonElectronPt_ = fBosonElectronPt;}
    void SetBosonElectronEta(double fBosonElectronEta) { BosonElectronEta_ = fBosonElectronEta;}
    void SetBosonElectronPhi(double fBosonElectronPhi) { BosonElectronPhi_ = fBosonElectronPhi;}
    void SetBosonMuonMass(double fBosonMuonMass) { BosonMuonMass_ = fBosonMuonMass;}
    void SetBosonMuonPt(double fBosonMuonPt) { BosonMuonPt_ = fBosonMuonPt;}
    void SetBosonMuonEta(double fBosonMuonEta) { BosonMuonEta_ = fBosonMuonEta;}
    void SetBosonMuonPhi(double fBosonMuonPhi) { BosonMuonPhi_ = fBosonMuonPhi;}

    void SetMETPt(double fmetPt)    {fmetPt_ = fmetPt;}
    void SetMETPhi(double fmetPhi)    {fmetPhi_ = fmetPhi;}
    void SetMETEt(double fmetEt)    {fmetEt_ = fmetEt;}
    void SetMETSumEt(double fmetSumEt)    {fmetSumEt_ = fmetSumEt;}
    void SetMETpx(double fmetpx)    {fmetpx_ = fmetpx;}
    void SetMETpy(double fmetpy)    {fmetpy_ = fmetpy;}
    void SetMETP4(LorentzVector fmetp4)    {fmetp4_     = fmetp4;}

    void SetPatMETPt(double fpatmetPt)    {fpatmetPt_ = fpatmetPt;}
    void SetPatMETPhi(double fpatmetPhi)    {fpatmetPhi_ = fpatmetPhi;}
    void SetPatMETEt(double fpatmetEt)    {fpatmetEt_ = fpatmetEt;}
    void SetPatMETSumEt(double fpatmetSumEt)    {fpatmetSumEt_ = fpatmetSumEt;}
    void SetPatMETpx(double fpatmetpx)    {fpatmetpx_ = fpatmetpx;}
    void SetPatMETpy(double fpatmetpy)    {fpatmetpy_ = fpatmetpy;}
    void SetPatMETP4(LorentzVector fpatmetp4)    {fpatmetp4_     = fpatmetp4;}

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
    void SetLeadingMuonIsGlobal(bool fLeadingMuonIsGlobal)    {LeadingMuonIsGlobal_ = fLeadingMuonIsGlobal;}
    void SetLeadingMuonIsTracker(bool fLeadingMuonIsTracker)    {LeadingMuonIsTracker_ = fLeadingMuonIsTracker;}

    void SetVertexMultiplicity(const std::vector<double>& fVertexMultiplicity) { VertexMultiplicity_ = fVertexMultiplicity; }
    void SetVertexChiNorm(const std::vector<double>& fVertexChiNorm) { VertexChiNorm_ = fVertexChiNorm; }
    void SetVertexNDOF(const std::vector<double>& fVertexNDOF) { VertexNDOF_ = fVertexNDOF; }
    void SetVz(const std::vector<double>& fVz) { Vz_ = fVz; }
    void SetVx(const std::vector<double>& fVx) { Vx_ = fVx; }
    void SetVy(const std::vector<double>& fVy) { Vy_ = fVy; }
    void SetTracksPt(const std::vector<std::vector<double> >& fTracksPt) { TracksPt_ = fTracksPt; }
    void SetTrackEtaMin(double ftracketamin)    { tracketamin_     = ftracketamin;}
    void SetTrackEtaMax(double ftracketamax)    { tracketamax_     = ftracketamax;}
    void SetZDCdigifC(const std::vector<std::vector<double> >& fZDCdigifC) { ZDCdigifC_ = fZDCdigifC; }

    void SetXiGenMinus(double fxigenminus)    { xigenminus_     = fxigenminus;}
    void SetXiGenPlus(double fxigenplus)    { xigenplus_     = fxigenplus;}

    void SetMxGenMinus(double fmxgenminus)    { mxgenminus_     = fmxgenminus;}
    void SetMxGenPlus(double fmxgenplus)    { mxgenplus_     = fmxgenplus;}
    void SetMx2GenMinus(double fmx2genminus)    { mx2genminus_     = fmx2genminus;}
    void SetMx2GenPlus(double fmx2genplus)    { mx2genplus_     = fmx2genplus;}
    void SetMxGenLeft(double fmxgenleft)    { mxgenleft_     = fmxgenleft;}
    void SetMxGenRight(double fmxgenright)    { mxgenright_     = fmxgenright;}
    void SetMx2GenLeft(double fmx2genleft)    { mx2genleft_     = fmx2genleft;}
    void SetMx2GenRight(double fmx2genright)    { mx2genright_     = fmx2genright;}
    void SetLrgGen(double flrggen)    { lrggen_     = flrggen;}
    void SetEtaMaxGen(double fetamaxgen)    { etamaxgen_     = fetamaxgen;}
    void SetEtaMinGen(double fetamingen)    { etamingen_     = fetamingen;}
    void SetEpluspzGen(double fepluspzgen)    { epluspzgen_     = fepluspzgen;}
    void SetEminuspzGen(double feminuspzgen)    { eminuspzgen_     = feminuspzgen;}
    void SetEtExpoPlusGen(double fetexpoplusgen)    { etexpoplusgen_     = fetexpoplusgen;}
    void SetEtExpoMinusGen(double fetexpominusgen)    { etexpominusgen_     = fetexpominusgen;}
    void SetsumECastorMinusGen(double fsumECastorMinusGen)    { sumECastorMinusGen_     = fsumECastorMinusGen;}
    void SetSumptGenLeft(double fsumptgenleft)    { sumptgenleft_     = fsumptgenleft;}
    void SetSumptGenRight(double fsumptgenright)    { sumptgenright_     = fsumptgenright;}

    void SetMxGenMinusCMS(double fmxgenminusCMS)    { mxgenminusCMS_     = fmxgenminusCMS;}
    void SetMxGenPlusCMS(double fmxgenplusCMS)    { mxgenplusCMS_     = fmxgenplusCMS;}
    void SetMx2GenMinusCMS(double fmx2genminusCMS)    { mx2genminusCMS_     = fmx2genminusCMS;}
    void SetMx2GenPlusCMS(double fmx2genplusCMS)    { mx2genplusCMS_     = fmx2genplusCMS;}
    void SetMxGenLeftCMS(double fmxgenleftCMS)    { mxgenleftCMS_     = fmxgenleftCMS;}
    void SetMxGenRightCMS(double fmxgenrightCMS)    { mxgenrightCMS_     = fmxgenrightCMS;}
    void SetMx2GenLeftCMS(double fmx2genleftCMS)    { mx2genleftCMS_     = fmx2genleftCMS;}
    void SetMx2GenRightCMS(double fmx2genrightCMS)    { mx2genrightCMS_     = fmx2genrightCMS;}
    void SetLrgGenCMS(double flrggenCMS)    { lrggenCMS_     = flrggenCMS;}
    void SetEtaMaxGenCMS(double fetamaxgenCMS)    { etamaxgenCMS_     = fetamaxgenCMS;}
    void SetEtaMinGenCMS(double fetamingenCMS)    { etamingenCMS_     = fetamingenCMS;}
    void SetEpluspzGenCMS(double fepluspzgenCMS)    { epluspzgenCMS_     = fepluspzgenCMS;}
    void SetEminuspzGenCMS(double feminuspzgenCMS)    { eminuspzgenCMS_     = feminuspzgenCMS;}
    void SetEtExpoPlusGenCMS(double fetexpoplusgenCMS)    { etexpoplusgenCMS_     = fetexpoplusgenCMS;}
    void SetEtExpoMinusGenCMS(double fetexpominusgenCMS)    { etexpominusgenCMS_     = fetexpominusgenCMS;}
    void SetsumECastorMinusGenCMS(double fsumECastorMinusGenCMS)    { sumECastorMinusGenCMS_     = fsumECastorMinusGenCMS;}
    void SetSumptGenLeftCMS(double fsumptgenleftCMS)    { sumptgenleftCMS_     = fsumptgenleftCMS;}
    void SetSumptGenRightCMS(double fsumptgenrightCMS)    { sumptgenrightCMS_     = fsumptgenrightCMS;}

    void SetBosonElectronMassPF(double fBosonElectronMassPF) { BosonElectronMassPF_ = fBosonElectronMassPF;}
    void SetBosonMuonMassPF(double fBosonMuonMassPF) { BosonMuonMassPF_ = fBosonMuonMassPF;}

    void SetSumEHFPlus(double fSumEHFPlus)    { SumEHFPlus_ = fSumEHFPlus;}
    void SetSumEHF_SPlus(double fSumEHF_SPlus)    { SumEHF_SPlus_ = fSumEHF_SPlus;}
    void SetSumEHF_LPlus(double fSumEHF_LPlus)    { SumEHF_LPlus_ = fSumEHF_LPlus;}
    void SetSumEtHFPlus(double fSumEtHFPlus)    { SumEtHFPlus_ = fSumEtHFPlus;}

    void SetSumEHFMinus(double fSumEHFMinus)    { SumEHFMinus_ = fSumEHFMinus;}
    void SetSumEHF_SMinus(double fSumEHF_SMinus)    { SumEHF_SMinus_ = fSumEHF_SMinus;}
    void SetSumEHF_LMinus(double fSumEHF_LMinus)    { SumEHF_LMinus_ = fSumEHF_LMinus;}
    void SetSumEtHFMinus(double fSumEtHFMinus)    { SumEtHFMinus_ = fSumEtHFMinus;}

    void SetSumEHEPlus(double fSumEHEPlus)    { SumEHEPlus_ = fSumEHEPlus;}
    void SetSumEtHEPlus(double fSumEtHEPlus)    { SumEtHEPlus_ = fSumEtHEPlus;}
    void SetSumEHEMinus(double fSumEHEMinus)    { SumEHEMinus_ = fSumEHEMinus;}
    void SetSumEtHEMinus(double fSumEtHEMinus)    { SumEtHEMinus_ = fSumEtHEMinus;}

    void SetSumEHBPlus(double fSumEHBPlus)    { SumEHBPlus_ = fSumEHBPlus;}
    void SetSumEtHBPlus(double fSumEtHBPlus)    { SumEtHBPlus_ = fSumEtHBPlus;}
    void SetSumEHBMinus(double fSumEHBMinus)    { SumEHBMinus_ = fSumEHBMinus;}
    void SetSumEtHBMinus(double fSumEtHBMinus)    { SumEtHBMinus_ = fSumEtHBMinus;}

    void SetSumEEEPlus(double fSumEEEPlus)    { SumEEEPlus_ = fSumEEEPlus;}
    void SetSumEtEEPlus(double fSumEtEEPlus)    { SumEtEEPlus_ = fSumEtEEPlus;}
    void SetSumEEEMinus(double fSumEEEMinus)    { SumEEEMinus_ = fSumEEEMinus;}
    void SetSumEtEEMinus(double fSumEtEEMinus)    { SumEtEEMinus_ = fSumEtEEMinus;}

    void SetSumEEBPlus(double fSumEEBPlus)    { SumEEBPlus_ = fSumEEBPlus;}
    void SetSumEtEBPlus(double fSumEtEBPlus)    { SumEtEBPlus_ = fSumEtEBPlus;}
    void SetSumEEBMinus(double fSumEEBMinus)    { SumEEBMinus_ = fSumEEBMinus;}
    void SetSumEtEBMinus(double fSumEtEBMinus)    { SumEtEBMinus_ = fSumEtEBMinus;}

    void SetEpluspzCalo(double fEPZCaloPlus)    { EPZCaloPlus_ = fEPZCaloPlus;}
    void SetEminuspzCalo(double fEPZCaloMinus)    { EPZCaloMinus_ = fEPZCaloMinus;}
    void SetEtExpoPlusCalo(double fEtEtaCaloPlus)    { EtEtaCaloPlus_ = fEtEtaCaloPlus;}
    void SetEtExpoMinusCalo(double fEtEtaCaloMinus)    { EtEtaCaloMinus_ = fEtEtaCaloMinus;}

    void SetEtaCaloMax(double fEtaCaloMax)    { EtaCaloMax_ = fEtaCaloMax;}
    void SetEtaCaloMin(double fEtaCaloMin)    { EtaCaloMin_ = fEtaCaloMin;}
    void SetLrgCalo(double flrgCalo)    { lrgCalo_     = flrgCalo;}

    void SetMultiplicityHFPlus(int fMultiplicityHFPlus)    { MultiplicityHFPlus_ = fMultiplicityHFPlus;}
    void SetMultiplicityHEPlus(int fMultiplicityHEPlus)    { MultiplicityHEPlus_ = fMultiplicityHEPlus;}
    void SetMultiplicityEEPlus(int fMultiplicityEEPlus)    { MultiplicityEEPlus_ = fMultiplicityEEPlus;}
    void SetMultiplicityHFMinus(int fMultiplicityHFMinus)    { MultiplicityHFMinus_ = fMultiplicityHFMinus;}
    void SetMultiplicityHEMinus(int fMultiplicityHEMinus)    { MultiplicityHEMinus_ = fMultiplicityHEMinus;}
    void SetMultiplicityEEMinus(int fMultiplicityEEMinus)    { MultiplicityEEMinus_ = fMultiplicityEEMinus;}

    void SetVertex(int fVertex)    { Vertex_ = fVertex;}
    void SetMxPFMinus(double fmxpfminus)    { mxpfminus_     = fmxpfminus;}
    void SetMxPFPlus(double fmxpfplus)    { mxpfplus_     = fmxpfplus;}  
    void SetMx2PFMinus(double fmx2pfminus)    { mx2pfminus_     = fmx2pfminus;}
    void SetMx2PFPlus(double fmx2pfplus)    { mx2pfplus_     = fmx2pfplus;}
    void SetEtaMaxPF(double fetamaxpf)    { etamaxpf_     = fetamaxpf;}
    void SetEtaMinPF(double fetaminpf)    { etaminpf_     = fetaminpf;}
    void SetEpluspzPF(double fepluspzpf)    { epluspzpf_     = fepluspzpf;}
    void SetEminuspzPF(double feminuspzpf)    { eminuspzpf_     = feminuspzpf;}
    void SetEtExpoPlusPF(double fetexpopluspf)    { etexpopluspf_     = fetexpopluspf;}
    void SetEtExpoMinusPF(double fetexpominuspf)    { etexpominuspf_     = fetexpominuspf;}
    void SetLrgPF(double flrgPF)    { lrgPF_     = flrgPF;}
    void SetMxPFLeft(double fmxpfleft)    { mxpfleft_     = fmxpfleft;}
    void SetMxPFRight(double fmxpfright)    { mxpfright_     = fmxpfright;}
    void SetMx2PFLeft(double fmx2pfleft)    { mx2pfleft_     = fmx2pfleft;}
    void SetMx2PFRight(double fmx2pfright)    { mx2pfright_     = fmx2pfright;}
    void SetSumptPFLeft(double fsumptpfleft)    { sumptpfleft_     = fsumptpfleft;}
    void SetSumptPFRight(double fsumptpfright)    { sumptpfright_     = fsumptpfright;}

    void SetMxPFNoWMinus(double fmxpfnowminus)    { mxpfnowminus_     = fmxpfnowminus;} 
    void SetMxPFNoWPlus(double fmxpfnowplus)    { mxpfnowplus_     = fmxpfnowplus;}  
    void SetMx2PFNoWMinus(double fmx2pfnowminus)    { mx2pfnowminus_     = fmx2pfnowminus;}
    void SetMx2PFNoWPlus(double fmx2pfnowplus)    { mx2pfnowplus_     = fmx2pfnowplus;}
    void SetEtaMaxPFNoW(double fetamaxpfnow)    { etamaxpfnow_     = fetamaxpfnow;}
    void SetEtaMinPFNoW(double fetaminpfnow)    { etaminpfnow_     = fetaminpfnow;}
    void SetEpluspzPFNoW(double fepluspzpfnow)    { epluspzpfnow_     = fepluspzpfnow;}
    void SetEminuspzPFNoW(double feminuspzpfnow)    { eminuspzpfnow_     = feminuspzpfnow;}
    void SetEtExpoPlusPFNoW(double fetexpopluspfnow)    { etexpopluspfnow_     = fetexpopluspfnow;}
    void SetEtExpoMinusPFNoW(double fetexpominuspfnow)    { etexpominuspfnow_     = fetexpominuspfnow;}
    void SetLrgPFNoW(double flrgPFnow)    { lrgPFnow_     = flrgPFnow;}
    void SetMxPFNoWLeft(double fmxpfnowleft)    { mxpfnowleft_     = fmxpfnowleft;}
    void SetMxPFNoWRight(double fmxpfnowright)    { mxpfnowright_     = fmxpfnowright;}
    void SetMx2PFNoWLeft(double fmx2pfnowleft)    { mx2pfnowleft_     = fmx2pfnowleft;}
    void SetMx2PFNoWRight(double fmx2pfnowright)    { mx2pfnowright_     = fmx2pfnowright;}
    void SetSumptPFNoWLeft(double fsumptpfnowleft)    { sumptpfnowleft_     = fsumptpfnowleft;}
    void SetSumptPFNoWRight(double fsumptpfnowright)    { sumptpfnowright_     = fsumptpfnowright;}

    void SetPatNMuon(int fpatNMuon)    {patNMuon_ = fpatNMuon;}
    void SetPatMuon1Pt(double fpatMuon1Pt)    {patMuon1Pt_ = fpatMuon1Pt;}
    void SetPatMuon1Charge(int fpatMuon1Charge)    {patMuon1Charge_ = fpatMuon1Charge;}
    void SetPatMuon1Phi(double fpatMuon1Phi)    {patMuon1Phi_ = fpatMuon1Phi;}
    void SetPatMuon1Eta(double fpatMuon1Eta)    {patMuon1Eta_ = fpatMuon1Eta;}
    void SetPatMuon1Et(double fpatMuon1Et)    {patMuon1Et_ = fpatMuon1Et;}
    void SetPatMuon1P4(LorentzVector fpatMuon1P4)    {patMuon1P4_ = fpatMuon1P4;}

    void SetPatMuon1SumPtR03(double fpatMuon1SumPtR03)    {patMuon1SumPtR03_ = fpatMuon1SumPtR03;}
    void SetPatMuon1EmEtR03(double fpatMuon1EmEtR03)    {patMuon1EmEtR03_ = fpatMuon1EmEtR03;}
    void SetPatMuon1HadEtR03(double fpatMuon1HadEtR03)    {patMuon1HadEtR03_ = fpatMuon1HadEtR03;}    
    void SetPatMuon1SumPtR05(double fpatMuon1SumPtR05)    {patMuon1SumPtR05_ = fpatMuon1SumPtR05;}
    void SetPatMuon1EmEtR05(double fpatMuon1EmEtR05)    {patMuon1EmEtR05_ = fpatMuon1EmEtR05;}
    void SetPatMuon1HadEtR05(double fpatMuon1HadEtR05)    {patMuon1HadEtR05_ = fpatMuon1HadEtR05;}    

    void SetPatMuon1relIsoDr03(double fpatMuon1relIsoDr03)    {patMuon1relIsoDr03_ = fpatMuon1relIsoDr03;}
    void SetPatMuon1relIsoDr05(double fpatMuon1relIsoDr05)    {patMuon1relIsoDr05_ = fpatMuon1relIsoDr05;}
    void SetPatMuon1relIso(double fpatMuon1relIso)    {patMuon1relIso_ = fpatMuon1relIso;}

    void SetPatMuon1TrackerHits(double fpatMuon1TrackerHits)    {patMuon1TrackerHits_ = fpatMuon1TrackerHits;}
    void SetPatMuon1PixelHits(double fpatMuon1PixelHits)    {patMuon1PixelHits_ = fpatMuon1PixelHits;}
    void SetPatMuon1NormalizedChi2(double fpatMuon1NormalizedChi2)    {patMuon1NormalizedChi2_ = fpatMuon1NormalizedChi2;}
    void SetPatMuon1MatchedStations(double fpatMuon1MatchedStations)    {patMuon1MatchedStations_ = fpatMuon1MatchedStations;}
    void SetPatMuon1Dxy(double fpatMuon1Dxy)    {patMuon1Dxy_ = fpatMuon1Dxy;}
    void SetPatMuon1IsGlobal(bool fpatMuon1IsGlobal)    {patMuon1IsGlobal_ = fpatMuon1IsGlobal;}
    void SetPatMuon1IsTracker(bool fpatMuon1IsTracker)    {patMuon1IsTracker_ = fpatMuon1IsTracker;}

    void SetPatBosonMuonMass(double fpatBosonMuonMass) { patBosonMuonMass_ = fpatBosonMuonMass;}
    void SetPatBosonMuonPt(double fpatBosonMuonPt) { patBosonMuonPt_ = fpatBosonMuonPt;}
    void SetPatBosonMuonEta(double fpatBosonMuonEta) { patBosonMuonEta_ = fpatBosonMuonEta;}
    void SetPatBosonMuonPhi(double fpatBosonMuonPhi) { patBosonMuonPhi_ = fpatBosonMuonPhi;}

    void SetPatNElectron(int fpatNElectron)    {patNElectron_ = fpatNElectron;}
    void SetPatElectron1Pt(double fpatElectron1Pt)    {patElectron1Pt_ = fpatElectron1Pt;}
    void SetPatElectron1Charge(int fpatElectron1Charge)    {patElectron1Charge_ = fpatElectron1Charge;}
    void SetPatElectron1Phi(double fpatElectron1Phi)    {patElectron1Phi_ = fpatElectron1Phi;}
    void SetPatElectron1Eta(double fpatElectron1Eta)    {patElectron1Eta_ = fpatElectron1Eta;}
    void SetPatElectron1Et(double fpatElectron1Et)    {patElectron1Et_ = fpatElectron1Et;}
    void SetPatElectron1P4(LorentzVector fpatElectron1P4)    {patElectron1P4_ = fpatElectron1P4;}

    void SetPatElectron1TkDr03(double fpatElectron1TkDr03)    {patElectron1TkDr03_ = fpatElectron1TkDr03;}
    void SetPatElectron1EcalDr03(double fpatElectron1EcalDr03)    {patElectron1EcalDr03_ = fpatElectron1EcalDr03;}
    void SetPatElectron1HcalDr03(double fpatElectron1HcalDr03)    {patElectron1HcalDr03_ = fpatElectron1HcalDr03;}

    void SetPatElectron1TkDr04(double fpatElectron1TkDr04)    {patElectron1TkDr04_ = fpatElectron1TkDr04;}
    void SetPatElectron1EcalDr04(double fpatElectron1EcalDr04)    {patElectron1EcalDr04_ = fpatElectron1EcalDr04;}
    void SetPatElectron1HcalDr04(double fpatElectron1HcalDr04)    {patElectron1HcalDr04_ = fpatElectron1HcalDr04;}

    void SetPatElectron1relIsoDr03(double fpatElectron1relIsoDr03)    {patElectron1relIsoDr03_ = fpatElectron1relIsoDr03;}
    void SetPatElectron1relIsoDr04(double fpatElectron1relIsoDr04)    {patElectron1relIsoDr04_ = fpatElectron1relIsoDr04;}

    void SetPatBosonElectronMass(double fpatBosonElectronMass) { patBosonElectronMass_ = fpatBosonElectronMass;}
    void SetPatBosonElectronPt(double fpatBosonElectronPt) { patBosonElectronPt_ = fpatBosonElectronPt;}
    void SetPatBosonElectronEta(double fpatBosonElectronEta) { patBosonElectronEta_ = fpatBosonElectronEta;}
    void SetPatBosonElectronPhi(double fpatBosonElectronPhi) { patBosonElectronPhi_ = fpatBosonElectronPhi;}

    void SetCastorTowerEnergy(const std::vector<double>& fCastorTowerEnergy) { CastorTowerEnergy_ = fCastorTowerEnergy; }
    void SetCastorModule1Energy(const std::vector<double>& fCastorModule1Energy) { CastorModule1Energy_ = fCastorModule1Energy; }
    void SetCastorModule2Energy(const std::vector<double>& fCastorModule2Energy) { CastorModule2Energy_ = fCastorModule2Energy; }
    void SetCastorModule3Energy(const std::vector<double>& fCastorModule3Energy) { CastorModule3Energy_ = fCastorModule3Energy; }
    void SetCastorModule4Energy(const std::vector<double>& fCastorModule4Energy) { CastorModule4Energy_ = fCastorModule4Energy; }
    void SetCastorModule5Energy(const std::vector<double>& fCastorModule5Energy) { CastorModule5Energy_ = fCastorModule5Energy; }
    void SetCastorBadChannels(const std::vector<int>& fCastorBadChannels) { CastorBadChannels_ = fCastorBadChannels; }
    void SetCastorNumberBadChannels(int fCastorNumberBadChannels) { CastorNumberBadChannels_ = fCastorNumberBadChannels;}

    void SetEachTowerEta(const std::vector<double>& fEachTowerEta) { EachTowerEta_ = fEachTowerEta; }
    void SetEachTowerEnergy(const std::vector<double>& fEachTowerEnergy) { EachTowerEnergy_ = fEachTowerEnergy; }
    void SetEachTowerCounter(int fEachTowerCounter)    {EachTowerCounter_ = fEachTowerCounter;}

    void SetTracksNonConeMuon03(int fTracksNonConeMuon03)    {TracksNonConeMuon03_ = fTracksNonConeMuon03;}
    void SetTracksNonConeElectron03(int fTracksNonConeElectron03)    {TracksNonConeElectron03_ = fTracksNonConeElectron03;}
    void SetTracksNonConepatMuon03(int fTracksNonConepatMuon03)    {TracksNonConepatMuon03_ = fTracksNonConepatMuon03;}
    void SetTracksNonConepatElectron03(int fTracksNonConepatElectron03)    {TracksNonConepatElectron03_ = fTracksNonConepatElectron03;}

    void SetTracksNonConeMuon04(int fTracksNonConeMuon04)    {TracksNonConeMuon04_ = fTracksNonConeMuon04;}
    void SetTracksNonConeElectron04(int fTracksNonConeElectron04)    {TracksNonConeElectron04_ = fTracksNonConeElectron04;}
    void SetTracksNonConepatMuon04(int fTracksNonConepatMuon04)    {TracksNonConepatMuon04_ = fTracksNonConepatMuon04;}
    void SetTracksNonConepatElectron04(int fTracksNonConepatElectron04)    {TracksNonConepatElectron04_ = fTracksNonConepatElectron04;}

    void SetTracksNonConeMuon05(int fTracksNonConeMuon05)    {TracksNonConeMuon05_ = fTracksNonConeMuon05;}
    void SetTracksNonConeElectron05(int fTracksNonConeElectron05)    {TracksNonConeElectron05_ = fTracksNonConeElectron05;}
    void SetTracksNonConepatMuon05(int fTracksNonConepatMuon05)    {TracksNonConepatMuon05_ = fTracksNonConepatMuon05;}
    void SetTracksNonConepatElectron05(int fTracksNonConepatElectron05)    {TracksNonConepatElectron05_ = fTracksNonConepatElectron05;}

    void SetLeadingElectronDeltaPhiTkClu(double fLeadingElectronDeltaPhiTkClu)    {LeadingElectronDeltaPhiTkClu_ = fLeadingElectronDeltaPhiTkClu;}
    void SetLeadingElectronDeltaEtaTkClu(double fLeadingElectronDeltaEtaTkClu)    {LeadingElectronDeltaEtaTkClu_ = fLeadingElectronDeltaEtaTkClu;}
    void SetLeadingElectronSigmaIeIe(double fLeadingElectronSigmaIeIe)    {LeadingElectronSigmaIeIe_ = fLeadingElectronSigmaIeIe;}
    void SetLeadingElectronDCot(double fLeadingElectronDCot)    {LeadingElectronDCot_ = fLeadingElectronDCot;}
    void SetLeadingElectronDist(double fLeadingElectronDist)    {LeadingElectronDist_ = fLeadingElectronDist;}
    void SetLeadingElectronInnerHits(double fLeadingElectronInnerHits)    {LeadingElectronInnerHits_ = fLeadingElectronInnerHits;}
    void SetLeadingElectronHE(double fLeadingElectronHE)    {LeadingElectronHE_ = fLeadingElectronHE;}

    void SetPatElectron1DeltaPhiTkClu(double fpatElectron1DeltaPhiTkClu)    {patElectron1DeltaPhiTkClu_ = fpatElectron1DeltaPhiTkClu;}
    void SetPatElectron1DeltaEtaTkClu(double fpatElectron1DeltaEtaTkClu)    {patElectron1DeltaEtaTkClu_ = fpatElectron1DeltaEtaTkClu;}
    void SetPatElectron1SigmaIeIe(double fpatElectron1SigmaIeIe)    {patElectron1SigmaIeIe_ = fpatElectron1SigmaIeIe;}
    void SetPatElectron1DCot(double fpatElectron1DCot)    {patElectron1DCot_ = fpatElectron1DCot;}
    void SetPatElectron1Dist(double fpatElectron1Dist)    {patElectron1Dist_ = fpatElectron1Dist;}
    void SetPatElectron1InnerHits(double fpatElectron1InnerHits)    {patElectron1InnerHits_ = fpatElectron1InnerHits;}
    void SetPatElectron1HE(double fpatElectron1HE)    {patElectron1HE_ = fpatElectron1HE;}

    int GetHLTPath(int idx)                    const { return hltTrigResults_[idx]; }
    double GetBosonElectronMass() const {return BosonElectronMass_;}
    double GetBosonElectronPt() const {return BosonElectronPt_;}
    double GetBosonElectronEta() const {return BosonElectronEta_;}
    double GetBosonElectronPhi() const {return BosonElectronPhi_;}
    double GetBosonMuonMass() const {return BosonMuonMass_;}
    double GetBosonMuonPt() const {return BosonMuonPt_;}
    double GetBosonMuonEta() const {return BosonMuonEta_;}
    double GetBosonMuonPhi() const {return BosonMuonPhi_;}

    double GetLeadingElectronPt() const {return LeadingElectronPt_;}
    double GetLeadingElectronEta() const {return LeadingElectronEta_;}
    double GetLeadingElectronPhi() const {return LeadingElectronPhi_;}
    const LorentzVector& GetLeadingElectronP4() const {return LeadingElectronP4_;}
    int GetLeadingElectronCharge() const {return LeadingElectronCharge_;}
    int GetElectronsN() const {return ElectronsN_;}
    double GetLeadingElectronTkDr03() const  {return LeadingElectronTkDr03_;}
    double GetLeadingElectronEcalDr03() const  {return LeadingElectronEcalDr03_;}
    double GetLeadingElectronHcalDr03() const  {return LeadingElectronHcalDr03_;}

    double GetLeadingElectronTkDr04() const  {return LeadingElectronTkDr04_;}
    double GetLeadingElectronEcalDr04() const  {return LeadingElectronEcalDr04_;}
    double GetLeadingElectronHcalDr04() const  {return LeadingElectronHcalDr04_;}

    double GetLeadingElectronrelIsoDr03() const {return LeadingElectronrelIsoDr03_;}
    double GetLeadingElectronrelIsoDr04() const {return LeadingElectronrelIsoDr04_;}

    double GetLeadingMuonPt() const {return LeadingMuonPt_;}
    double GetLeadingMuonEta() const {return LeadingMuonEta_;}
    double GetLeadingMuonPhi() const {return LeadingMuonPhi_;}
    const LorentzVector& GetLeadingMuonP4() const {return LeadingMuonP4_;}
    int GetLeadingMuonCharge() const {return LeadingMuonCharge_;}
    int GetMuonsN() const {return MuonsN_;}

    double GetLeadingMuonSumPtR03() const {return LeadingMuonSumPtR03_;}
    double GetLeadingMuonEmEtR03() const {return LeadingMuonEmEtR03_;}
    double GetLeadingMuonHadEtR03() const {return LeadingMuonHadEtR03_;}
    double GetLeadingMuonSumPtR05() const {return LeadingMuonSumPtR05_;}
    double GetLeadingMuonEmEtR05() const {return LeadingMuonEmEtR05_;}
    double GetLeadingMuonHadEtR05() const {return LeadingMuonHadEtR05_;}

    double GetLeadingMuonrelIsoDr03() const {return LeadingMuonrelIsoDr03_;}
    double GetLeadingMuonrelIsoDr05() const {return LeadingMuonrelIsoDr05_;}

    double GetLeadingMuonTrackerHits() const   {return LeadingMuonTrackerHits_;}
    double GetLeadingMuonPixelHits() const   {return LeadingMuonPixelHits_;}
    double GetLeadingMuonNormalizedChi2() const   {return LeadingMuonNormalizedChi2_;}
    double GetLeadingMuonMatchedStations() const   {return LeadingMuonMatchedStations_;}
    double GetLeadingMuonDxy() const    {return LeadingMuonDxy_;}
    bool GetLeadingMuonIsGlobal() const   {return LeadingMuonIsGlobal_;}
    bool GetLeadingMuonIsTracker() const   {return LeadingMuonIsTracker_;}

    double GetVertexMultiplicity(int i) const { return VertexMultiplicity_[i]; }
    double GetVertexChiNorm(int i) const { return VertexChiNorm_[i]; }
    double GetVertexNDOF(int i) const { return VertexNDOF_[i]; }
    double GetVz(int i) const { return Vz_[i]; }
    double GetVx(int i) const { return Vx_[i]; }
    double GetVy(int i) const { return Vy_[i]; }
    double GetTracksPt(int i,int j) const { return TracksPt_[i][j]; }
    double GetTrackEtaMin() const {return tracketamin_;}
    double GetTrackEtaMax() const {return tracketamax_;}
    double GetZDCdigifC(int i,int j) const { return ZDCdigifC_[i][j]; }

    double GetXiGenMinus() const {return xigenminus_;}
    double GetXiGenPlus() const {return xigenplus_;}

    double GetMxGenMinus() const { return mxgenminus_;}
    double GetMxGenPlus() const { return mxgenplus_;}
    double GetMx2GenMinus() const { return mx2genminus_;}
    double GetMx2GenPlus() const { return mx2genplus_;}
    double GetMxGenLeft() const { return mxgenleft_;}
    double GetMxGenRight() const { return mxgenright_;}
    double GetMx2GenLeft() const { return mx2genleft_;}
    double GetMx2GenRight() const { return mx2genright_;}
    double GetLrgGen() const { return lrggen_;}
    double GetEtaMaxGen() const { return etamaxgen_;}
    double GetEtaMinGen() const { return etamingen_;}
    double GetEpluspzGen() const { return epluspzgen_;}
    double GetEminuspzGen() const { return eminuspzgen_;}
    double GetEtExpoPlusGen() const { return etexpoplusgen_;}
    double GetEtExpoMinusGen() const { return etexpominusgen_;}
    double GetsumECastorMinusGen() const { return sumECastorMinusGen_;}
    double GetSumptGenLeft() const { return sumptgenleft_;}
    double GetSumptGenRight() const { return sumptgenright_;}

    double GetMxGenMinusCMS() const { return mxgenminusCMS_;}
    double GetMxGenPlusCMS() const { return mxgenplusCMS_;}
    double GetMx2GenMinusCMS() const { return mx2genminusCMS_;}
    double GetMx2GenPlusCMS() const { return mx2genplusCMS_;}
    double GetMxGenLeftCMS() const { return mxgenleftCMS_;}
    double GetMxGenRightCMS() const { return mxgenrightCMS_;}
    double GetMx2GenLeftCMS() const { return mx2genleftCMS_;}
    double GetMx2GenRightCMS() const { return mx2genrightCMS_;}
    double GetLrgGenCMS() const { return lrggenCMS_;}
    double GetEtaMaxGenCMS() const { return etamaxgenCMS_;}
    double GetEtaMinGenCMS() const { return etamingenCMS_;}
    double GetEpluspzGenCMS() const { return epluspzgenCMS_;}
    double GetEminuspzGenCMS() const { return eminuspzgenCMS_;}
    double GetEtExpoPlusGenCMS() const { return etexpoplusgenCMS_;}
    double GetEtExpoMinusGenCMS() const { return etexpominusgenCMS_;}
    double GetsumECastorMinusGenCMS() const { return sumECastorMinusGenCMS_;}
    double GetSumptGenLeftCMS() const { return sumptgenleftCMS_;}
    double GetSumptGenRightCMS() const { return sumptgenrightCMS_;}

    double GetBosonElectronMassPF() const {return BosonElectronMassPF_;}
    double GetBosonMuonMassPF() const {return BosonMuonMassPF_;}

    double GetSumEHFPlus()    const {return SumEHFPlus_;}
    double GetSumEHF_SPlus()    const {return SumEHF_SPlus_;}
    double GetSumEHF_LPlus()    const {return SumEHF_LPlus_;}
    double GetSumEtHFPlus()    const {return SumEtHFPlus_;}

    double GetSumEHFMinus()    const {return SumEHFMinus_;}
    double GetSumEHF_SMinus()    const {return SumEHF_SMinus_;}
    double GetSumEHF_LMinus()    const {return SumEHF_LMinus_;}
    double GetSumEtHFMinus()    const {return SumEtHFMinus_;}

    double GetSumEHEPlus()    const {return SumEHEPlus_;}
    double GetSumEtHEPlus()    const {return SumEtHEPlus_;}
    double GetSumEHEMinus()    const {return SumEHEMinus_;}
    double GetSumEtHEMinus()    const {return SumEtHEMinus_;}

    double GetSumEHBPlus()    const {return SumEHBPlus_;}
    double GetSumEtHBPlus()    const {return SumEtHBPlus_;}
    double GetSumEHBMinus()    const {return SumEHBMinus_;}
    double GetSumEtHBMinus()    const {return SumEtHBMinus_;}

    double GetSumEEEPlus()    const {return SumEEEPlus_;}
    double GetSumEtEEPlus()    const {return SumEtEEPlus_;}
    double GetSumEEEMinus()    const {return SumEEEMinus_;}
    double GetSumEtEEMinus()    const {return SumEtEEMinus_;}

    double GetSumEEBPlus()    const {return SumEEBPlus_;}
    double GetSumEtEBPlus()    const {return SumEtEBPlus_;}
    double GetSumEEBMinus()    const {return SumEEBMinus_;}
    double GetSumEtEBMinus()    const {return SumEtEBMinus_;}

    double GetEpluspzCalo()    const {return EPZCaloPlus_;}
    double GetEminuspzCalo()    const {return EPZCaloMinus_;}
    double GetEtExpoPlusCalo()    const {return EtEtaCaloPlus_;}
    double GetEtExpoPlusMinus()    const {return EtEtaCaloMinus_;}

    double GetEtaCaloMax()    const {return EtaCaloMax_;}
    double GetEtaCaloMin()    const {return EtaCaloMin_;}
    double GetLrgCalo() const { return lrgCalo_;}

    int GetMultiplicityHFPlus()    const {return MultiplicityHFPlus_;}
    int GetMultiplicityHEPlus()    const {return MultiplicityHEPlus_;}
    int GetMultiplicityEEPlus()    const {return MultiplicityEEPlus_;}
    int GetMultiplicityHFMinus()    const {return MultiplicityHFMinus_;}
    int GetMultiplicityHEMinus()    const {return MultiplicityHEMinus_;}
    int GetMultiplicityEEMinus()    const {return MultiplicityEEMinus_;}

    int GetVertex()    const {return Vertex_;}   
    double GetMxPFMinus() const { return mxpfminus_;}
    double GetMxPFPlus() const { return mxpfplus_;}
    double GetMx2PFMinus() const { return mx2pfminus_;}
    double GetMx2PFPlus() const { return mx2pfplus_;}
    double GetEtaMaxPF() const { return etamaxpf_;}
    double GetEtaMinPF() const { return etaminpf_;}
    double GetEpluspzPF() const { return epluspzpf_;}
    double GetEminuspzPF() const { return eminuspzpf_;}
    double GetEtExpoPlusPF() const { return etexpopluspf_;}
    double GetEtExpoMinusPF() const { return etexpominuspf_;}
    double GetLrgPF() const { return lrgPF_;}
    double GetMxPFLeft() const { return mxpfleft_;}
    double GetMxPFRight() const { return mxpfright_;}
    double GetMx2PFLeft() const { return mx2pfleft_;}
    double GetMx2PFRight() const { return mx2pfright_;}
    double GetSumptPFLeft() const { return sumptpfleft_;}
    double GetSumptPFRight() const { return sumptpfright_;}

    double GetMxPFNoWMinus() const { return mxpfnowminus_;}
    double GetMxPFNoWPlus() const { return mxpfnowplus_;}
    double GetMx2PFNoWMinus() const { return mx2pfnowminus_;}
    double GetMx2PFNoWPlus() const { return mx2pfnowplus_;}
    double GetEtaMaxPFNoW() const { return etamaxpfnow_;}
    double GetEtaMinPFNoW() const { return etaminpfnow_;}
    double GetEpluspzPFNoW() const { return epluspzpfnow_;}
    double GetEminuspzPFNoW() const { return eminuspzpfnow_;}
    double GetEtExpoPlusPFNoW() const { return etexpopluspfnow_;}
    double GetEtExpoMinusPFNoW() const { return etexpominuspfnow_;}
    double GetLrgPFNoW() const { return lrgPFnow_;}
    double GetMxPFNoWLeft() const { return mxpfnowleft_;}
    double GetMxPFNoWRight() const { return mxpfnowright_;}
    double GetMx2PFNoWLeft() const { return mx2pfnowleft_;}
    double GetMx2PFNoWRight() const { return mx2pfnowright_;}
    double GetSumptPFNoWLeft() const { return sumptpfnowleft_;}
    double GetSumptPFNoWRight() const { return sumptpfnowright_;}

    int GetPatNMuon() const {return patNMuon_;}
    double GetPatMuon1Pt() const {return patMuon1Pt_;}
    int GetPatMuon1Charge() const {return patMuon1Charge_;}
    double GetPatMuon1Phi() const {return patMuon1Phi_;}
    double GetPatMuon1Eta() const {return patMuon1Eta_;}
    double GetPatMuon1Et() const {return patMuon1Et_;}
    const LorentzVector& GetPatMuon1P4() const {return patMuon1P4_;}

    double GetPatMuon1SumPtR03() const {return patMuon1SumPtR03_;}
    double GetPatMuon1EmEtR03() const {return patMuon1EmEtR03_;}
    double GetPatMuon1HadEtR03() const {return patMuon1HadEtR03_;}    
    double GetPatMuon1SumPtR05() const {return patMuon1SumPtR05_;}
    double GetPatMuon1EmEtR05() const {return patMuon1EmEtR05_;}
    double GetPatMuon1HadEtR05() const {return patMuon1HadEtR05_;}    

    double GetPatMuon1relIsoDr03() const {return patMuon1relIsoDr03_;}
    double GetPatMuon1relIsoDr05() const {return patMuon1relIsoDr05_;}
    double GetPatMuon1relIso() const {return patMuon1relIso_;}

    double GetPatMuon1TrackerHits() const   {return patMuon1TrackerHits_;}
    double GetPatMuon1PixelHits() const   {return patMuon1PixelHits_;}
    double GetPatMuon1NormalizedChi2() const   {return patMuon1NormalizedChi2_;}
    double GetPatMuon1MatchedStations() const   {return patMuon1MatchedStations_;}
    double GetPatMuon1Dxy() const    {return patMuon1Dxy_;}
    bool GetPatMuon1IsGlobal() const   {return patMuon1IsGlobal_;}
    bool GetPatMuon1IsTracker() const   {return patMuon1IsTracker_;}

    double GetPatBosonMuonMass() const {return patBosonMuonMass_;}
    double GetPatBosonMuonPt() const {return patBosonMuonPt_;}
    double GetPatBosonMuonEta() const {return patBosonMuonEta_;}
    double GetPatBosonMuonPhi() const {return patBosonMuonPhi_;}

    int GetPatNElectron() const {return patNElectron_;}
    double GetPatElectron1Pt() const {return patElectron1Pt_;}
    int GetPatElectron1Charge() const {return patElectron1Charge_;}
    double GetPatElectron1Phi() const {return patElectron1Phi_;}
    double GetPatElectron1Eta() const {return patElectron1Eta_;}
    double GetPatElectron1Et() const {return patElectron1Et_;}
    const LorentzVector& GetPatElectron1P4() const {return patElectron1P4_;}

    double GetPatElectron1TkDr03() const  {return patElectron1TkDr03_;}    
    double GetPatElectron1EcalDr03() const  {return patElectron1EcalDr03_;}
    double GetPatElectron1HcalDr03() const  {return patElectron1HcalDr03_;}

    double GetPatElectron1TkDr04() const  {return patElectron1TkDr04_;}
    double GetPatElectron1EcalDr04() const  {return patElectron1EcalDr04_;}
    double GetPatElectron1HcalDr04() const  {return patElectron1HcalDr04_;}

    double GetPatElectron1relIsoDr03() const {return patElectron1relIsoDr03_;}
    double GetPatElectron1relIsoDr04() const {return patElectron1relIsoDr04_;}

    double GetPatBosonElectronMass() const {return patBosonElectronMass_;}
    double GetPatBosonElectronPt() const {return patBosonElectronPt_;}
    double GetPatBosonElectronEta() const {return patBosonElectronEta_;}
    double GetPatBosonElectronPhi() const {return patBosonElectronPhi_;}

    double GetEachTowerEta(int i) const { return EachTowerEta_[i]; }
    double GetEachTowerEnergy(int i) const { return EachTowerEnergy_[i]; }
    int GetEachTowerCounter() const {return EachTowerCounter_;}

    double GetCastorTowerEnergy(int i) const { return CastorTowerEnergy_[i]; }
    double GetCastorModule1Energy(int i) const { return CastorModule1Energy_[i]; }
    double GetCastorModule2Energy(int i) const { return CastorModule2Energy_[i]; }
    double GetCastorModule3Energy(int i) const { return CastorModule3Energy_[i]; }
    double GetCastorModule4Energy(int i) const { return CastorModule4Energy_[i]; }
    double GetCastorModule5Energy(int i) const { return CastorModule5Energy_[i]; }
    int GetCastorBadChannels(int i) const { return CastorBadChannels_[i]; }
    int GetCastorNumberBadChannels() const { return CastorNumberBadChannels_;}

    int GetTracksNonConeMuon03()    const {return TracksNonConeMuon03_;}
    int GetTracksNonConeElectron03()    const {return TracksNonConeElectron03_;}
    int GetTracksNonConepatMuon03()    const {return TracksNonConepatMuon03_;}
    int GetTracksNonConepatElectron03()    const {return TracksNonConepatElectron03_;}

    int GetTracksNonConeMuon04()    const {return TracksNonConeMuon04_;}
    int GetTracksNonConeElectron04()    const {return TracksNonConeElectron04_;}
    int GetTracksNonConepatMuon04()    const {return TracksNonConepatMuon04_;}
    int GetTracksNonConepatElectron04()    const {return TracksNonConepatElectron04_;}

    int GetTracksNonConeMuon05()    const {return TracksNonConeMuon05_;}
    int GetTracksNonConeElectron05()    const {return TracksNonConeElectron05_;}
    int GetTracksNonConepatMuon05()    const {return TracksNonConepatMuon05_;}
    int GetTracksNonConepatElectron05()    const {return TracksNonConepatElectron05_;}

    double GetLeadingElectronDeltaPhiTkClu() const {return LeadingElectronDeltaPhiTkClu_;}
    double GetLeadingElectronDeltaEtaTkClu() const {return LeadingElectronDeltaEtaTkClu_;}
    double GetLeadingElectronSigmaIeIe() const {return LeadingElectronSigmaIeIe_;}
    double GetLeadingElectronDCot() const {return LeadingElectronDCot_;}
    double GetLeadingElectronDist() const {return LeadingElectronDist_;}
    double GetLeadingElectronInnerHits() const {return LeadingElectronInnerHits_;}
    double GetLeadingElectronHE() const {return LeadingElectronHE_;}

    double GetPatElectron1DeltaPhiTkClu() const {return patElectron1DeltaPhiTkClu_;}
    double GetPatElectron1DeltaEtaTkClu() const {return patElectron1DeltaEtaTkClu_;}
    double GetPatElectron1SigmaIeIe() const {return patElectron1SigmaIeIe_;}
    double GetPatElectron1DCot() const {return patElectron1DCot_;}
    double GetPatElectron1Dist() const {return patElectron1Dist_;}
    double GetPatElectron1InnerHits() const {return patElectron1InnerHits_;}
    double GetPatElectron1HE() const {return patElectron1HE_;}

    double GetMETPt() const {return fmetPt_;}
    double GetMETPhi() const {return fmetPhi_;}
    double GetMETEt() const {return fmetEt_;}
    double GetMETSumEt() const {return fmetSumEt_;}
    double GetMETpx() const {return fmetpx_;}
    double GetMETpy() const {return fmetpy_;}
    const LorentzVector& GetMETP4() const {return fmetp4_;}

    double GetPatMETPt() const {return fpatmetPt_;}
    double GetPatMETPhi() const {return fpatmetPhi_;}
    double GetPatMETEt() const {return fpatmetEt_;}
    double GetPatMETSumEt() const {return fpatmetSumEt_;}
    double GetPatMETpx() const {return fpatmetpx_;}
    double GetPatMETpy() const {return fpatmetpy_;}
    const LorentzVector& GetPatMETP4() const {return fpatmetp4_;}

  private:
    friend class diffractiveWAnalysis::DiffractiveWAnalysis;

    void reset();

    int hltTrigResults_[20];
    double BosonElectronMass_;
    double BosonElectronPt_;
    double BosonElectronEta_;
    double BosonElectronPhi_;

    double BosonMuonMass_;
    double BosonMuonPt_;
    double BosonMuonEta_;
    double BosonMuonPhi_;

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
    bool LeadingMuonIsGlobal_;
    bool LeadingMuonIsTracker_;

    std::vector<double> VertexMultiplicity_;
    std::vector<double> VertexChiNorm_;
    std::vector<double> VertexNDOF_;
    std::vector<double> Vz_;
    std::vector<double> Vx_;
    std::vector<double> Vy_;
    std::vector<std::vector<double> > TracksPt_;
    std::vector<std::vector<double> > ZDCdigifC_;
    std::vector<double> EachTowerEta_;
    std::vector<double> EachTowerEnergy_;
    std::vector<double> CastorTowerEnergy_;
    std::vector<double> CastorModule1Energy_;
    std::vector<double> CastorModule2Energy_;
    std::vector<double> CastorModule3Energy_;
    std::vector<double> CastorModule4Energy_;
    std::vector<double> CastorModule5Energy_;
    std::vector<int> CastorBadChannels_;

    int CastorNumberBadChannels_;
    int EachTowerCounter_;

    double xigenminus_;
    double xigenplus_;
    double tracketamin_;
    double tracketamax_;

    double mxgenminus_;
    double mxgenplus_;
    double mx2genminus_;
    double mx2genplus_;
    double mxgenleft_;
    double mxgenright_;
    double mx2genleft_;
    double mx2genright_;
    double lrggen_;
    double etamaxgen_;
    double etamingen_;
    double epluspzgen_;
    double eminuspzgen_;
    double etexpoplusgen_;
    double etexpominusgen_;
    double sumECastorMinusGen_;
    double sumptgenleft_;
    double sumptgenright_;

    double mxgenminusCMS_;
    double mxgenplusCMS_;
    double mx2genminusCMS_;
    double mx2genplusCMS_;
    double mxgenleftCMS_;
    double mxgenrightCMS_;
    double mx2genleftCMS_;
    double mx2genrightCMS_;
    double lrggenCMS_;
    double etamaxgenCMS_;
    double etamingenCMS_;
    double epluspzgenCMS_;
    double eminuspzgenCMS_;
    double etexpoplusgenCMS_;
    double etexpominusgenCMS_;
    double sumECastorMinusGenCMS_;
    double sumptgenleftCMS_;
    double sumptgenrightCMS_;

    double BosonElectronMassPF_;
    double BosonMuonMassPF_;

    double SumEHFPlus_;
    double SumEHF_SPlus_;
    double SumEHF_LPlus_;
    double SumEtHFPlus_;
    double SumEHFMinus_;
    double SumEHF_SMinus_;
    double SumEHF_LMinus_;
    double SumEtHFMinus_;
    double SumEHEPlus_;
    double SumEtHEPlus_;
    double SumEHEMinus_;
    double SumEtHEMinus_;
    double SumEHBPlus_;
    double SumEtHBPlus_;
    double SumEHBMinus_;
    double SumEtHBMinus_;
    double SumEEEPlus_;
    double SumEtEEPlus_;
    double SumEEEMinus_;
    double SumEtEEMinus_;
    double SumEEBPlus_;
    double SumEtEBPlus_;
    double SumEEBMinus_;
    double SumEtEBMinus_;
    double EPZCaloPlus_;
    double EPZCaloMinus_;
    double EtEtaCaloPlus_;
    double EtEtaCaloMinus_;
    double EtaCaloMax_;
    double EtaCaloMin_;
    double lrgCalo_;
    int MultiplicityHFPlus_;
    int MultiplicityHEPlus_;
    int MultiplicityEEPlus_;
    int MultiplicityHFMinus_;
    int MultiplicityHEMinus_;
    int MultiplicityEEMinus_;

    int Vertex_;
    double mxpfminus_;
    double mxpfplus_;
    double mx2pfminus_;
    double mx2pfplus_;
    double etamaxpf_;
    double etaminpf_;
    double epluspzpf_; 
    double eminuspzpf_;
    double etexpopluspf_;
    double etexpominuspf_;
    double lrgPF_;
    double mxpfleft_;
    double mxpfright_;
    double mx2pfleft_;
    double mx2pfright_;
    double sumptpfleft_;
    double sumptpfright_;

    double mxpfnowminus_;
    double mxpfnowplus_;
    double mx2pfnowminus_;
    double mx2pfnowplus_;
    double etamaxpfnow_;
    double etaminpfnow_;
    double epluspzpfnow_;
    double eminuspzpfnow_;
    double etexpopluspfnow_;
    double etexpominuspfnow_;
    double lrgPFnow_;
    double mxpfnowleft_;
    double mxpfnowright_;
    double mx2pfnowleft_;
    double mx2pfnowright_;
    double sumptpfnowleft_;
    double sumptpfnowright_;

    int patNMuon_;

    double patMuon1Pt_;
    int patMuon1Charge_;
    double patMuon1Phi_;
    double patMuon1Eta_;
    double patMuon1Et_;
    LorentzVector patMuon1P4_;

    double patMuon1SumPtR03_;
    double patMuon1EmEtR03_;
    double patMuon1HadEtR03_;   
    double patMuon1SumPtR05_;
    double patMuon1EmEtR05_;
    double patMuon1HadEtR05_;   

    double patMuon1relIsoDr03_;
    double patMuon1relIsoDr05_;
    double patMuon1relIso_;

    double patMuon1TrackerHits_;
    double patMuon1PixelHits_;
    double patMuon1NormalizedChi2_;
    double patMuon1MatchedStations_;
    double patMuon1Dxy_;
    bool patMuon1IsGlobal_;
    bool patMuon1IsTracker_;

    double patBosonMuonMass_;
    double patBosonMuonPt_;
    double patBosonMuonEta_;
    double patBosonMuonPhi_;

    double patBosonElectronMass_;
    double patBosonElectronPt_;
    double patBosonElectronEta_;
    double patBosonElectronPhi_;

    int patNElectron_; 

    double patElectron1Pt_;
    int patElectron1Charge_;
    double patElectron1Phi_;
    double patElectron1Eta_;
    double patElectron1Et_;
    LorentzVector patElectron1P4_;

    double patElectron1TkDr03_;    
    double patElectron1EcalDr03_;
    double patElectron1HcalDr03_;

    double patElectron1TkDr04_;
    double patElectron1EcalDr04_;
    double patElectron1HcalDr04_;

    double patElectron1relIsoDr03_;
    double patElectron1relIsoDr04_;

    int TracksNonConeMuon03_;
    int TracksNonConeElectron03_;
    int TracksNonConepatMuon03_;
    int TracksNonConepatElectron03_;

    int TracksNonConeMuon04_;
    int TracksNonConeElectron04_;
    int TracksNonConepatMuon04_;
    int TracksNonConepatElectron04_;

    int TracksNonConeMuon05_;
    int TracksNonConeElectron05_;
    int TracksNonConepatMuon05_;
    int TracksNonConepatElectron05_;

    double LeadingElectronDeltaPhiTkClu_;
    double LeadingElectronDeltaEtaTkClu_;
    double LeadingElectronSigmaIeIe_;
    double LeadingElectronDCot_;
    double LeadingElectronDist_;
    double LeadingElectronInnerHits_;
    double LeadingElectronHE_;
    double patElectron1DeltaPhiTkClu_;
    double patElectron1DeltaEtaTkClu_;
    double patElectron1SigmaIeIe_;
    double patElectron1DCot_;
    double patElectron1Dist_;
    double patElectron1InnerHits_;
    double patElectron1HE_;

    double fmetPt_;
    double fmetPhi_;
    double fmetEt_;
    double fmetSumEt_;
    double fmetpx_;
    double fmetpy_;
    LorentzVector fmetp4_;

    double fpatmetPt_;
    double fpatmetPhi_;
    double fpatmetEt_;
    double fpatmetSumEt_;
    double fpatmetpx_;
    double fpatmetpy_;
    LorentzVector fpatmetp4_;

};

#endif    
