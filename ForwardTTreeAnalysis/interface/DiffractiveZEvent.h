#ifndef DiffractiveZEvent_h
#define DiffractiveZEvent_h

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/DetSetVector.h"

namespace diffractiveZAnalysis {
  class DiffractiveZAnalysis;
}

class DiffractiveZEvent {
  public:
    typedef diffractiveZAnalysis::DiffractiveZAnalysis analysis_type;
    static const char* name;

    typedef reco::Particle::LorentzVector LorentzVector;

    DiffractiveZEvent();
    ~DiffractiveZEvent();

    void SetHLTPath(int idx, int fHLTBit)         { hltTrigResults_[idx] = fHLTBit;}
    void SetDiElectronMass(double fDiElectronMass) { DiElectronMass_ = fDiElectronMass;}
    void SetDiElectronPt(double fDiElectronPt) { DiElectronPt_ = fDiElectronPt;}
    void SetDiElectronEta(double fDiElectronEta) { DiElectronEta_ = fDiElectronEta;}
    void SetDiElectronPhi(double fDiElectronPhi) { DiElectronPhi_ = fDiElectronPhi;}
    void SetDiMuonMass(double fDiMuonMass) { DiMuonMass_ = fDiMuonMass;}
    void SetDiMuonPt(double fDiMuonPt) { DiMuonPt_ = fDiMuonPt;}
    void SetDiMuonEta(double fDiMuonEta) { DiMuonEta_ = fDiMuonEta;}
    void SetDiMuonPhi(double fDiMuonPhi) { DiMuonPhi_ = fDiMuonPhi;}

    void SetLeadingElectronPt(double fLeadingElectronPt)    { LeadingElectronPt_     = fLeadingElectronPt;}
    void SetLeadingElectronEta(double fLeadingElectronEta)  { LeadingElectronEta_     = fLeadingElectronEta;}
    void SetLeadingElectronPhi(double fLeadingElectronPhi)  { LeadingElectronPhi_    = fLeadingElectronPhi;}
    void SetLeadingElectronP4(LorentzVector fLeadingElectronP4)    { LeadingElectronP4_     = fLeadingElectronP4;}
    void SetLeadingElectronCharge(int fLeadingElectronCharge)  { LeadingElectronCharge_     = fLeadingElectronCharge;}
    void SetSecondElectronPt(double fSecondElectronPt)    { SecondElectronPt_     = fSecondElectronPt;}
    void SetSecondElectronEta(double fSecondElectronEta)  { SecondElectronEta_     = fSecondElectronEta;}
    void SetSecondElectronPhi(double fSecondElectronPhi)  { SecondElectronPhi_    = fSecondElectronPhi;}
    void SetSecondElectronP4(LorentzVector fSecondElectronP4)    { SecondElectronP4_     = fSecondElectronP4;}
    void SetSecondElectronCharge(int fSecondElectronCharge)  { SecondElectronCharge_     = fSecondElectronCharge;}
    void SetElectronsN(int fElectronsN)  { ElectronsN_    = fElectronsN;}

    void SetLeadingElectronTkDr03(double fLeadingElectronTkDr03)    {LeadingElectronTkDr03_ = fLeadingElectronTkDr03;}
    void SetLeadingElectronEcalDr03(double fLeadingElectronEcalDr03)    {LeadingElectronEcalDr03_ = fLeadingElectronEcalDr03;}
    void SetLeadingElectronHcalDr03(double fLeadingElectronHcalDr03)    {LeadingElectronHcalDr03_ = fLeadingElectronHcalDr03;}
    void SetSecondElectronTkDr03(double fSecondElectronTkDr03)    {SecondElectronTkDr03_ = fSecondElectronTkDr03;}
    void SetSecondElectronEcalDr03(double fSecondElectronEcalDr03)    {SecondElectronEcalDr03_ = fSecondElectronEcalDr03;}
    void SetSecondElectronHcalDr03(double fSecondElectronHcalDr03)    {SecondElectronHcalDr03_ = fSecondElectronHcalDr03;}

    void SetLeadingElectronTkDr04(double fLeadingElectronTkDr04)    {LeadingElectronTkDr04_ = fLeadingElectronTkDr04;}
    void SetLeadingElectronEcalDr04(double fLeadingElectronEcalDr04)    {LeadingElectronEcalDr04_ = fLeadingElectronEcalDr04;}
    void SetLeadingElectronHcalDr04(double fLeadingElectronHcalDr04)    {LeadingElectronHcalDr04_ = fLeadingElectronHcalDr04;}
    void SetSecondElectronTkDr04(double fSecondElectronTkDr04)    {SecondElectronTkDr04_ = fSecondElectronTkDr04;}
    void SetSecondElectronEcalDr04(double fSecondElectronEcalDr04)    {SecondElectronEcalDr04_ = fSecondElectronEcalDr04;}
    void SetSecondElectronHcalDr04(double fSecondElectronHcalDr04)    {SecondElectronHcalDr04_ = fSecondElectronHcalDr04;}

    void SetLeadingElectronrelIsoDr03(double fLeadingElectronrelIsoDr03)    {LeadingElectronrelIsoDr03_ = fLeadingElectronrelIsoDr03;}
    void SetLeadingElectronrelIsoDr04(double fLeadingElectronrelIsoDr04)    {LeadingElectronrelIsoDr04_ = fLeadingElectronrelIsoDr04;}
    void SetSecondElectronrelIsoDr03(double fSecondElectronrelIsoDr03)    {SecondElectronrelIsoDr03_ = fSecondElectronrelIsoDr03;}
    void SetSecondElectronrelIsoDr04(double fSecondElectronrelIsoDr04)    {SecondElectronrelIsoDr04_ = fSecondElectronrelIsoDr04;}

    void SetLeadingMuonPt(double fLeadingMuonPt)    { LeadingMuonPt_     = fLeadingMuonPt;}
    void SetLeadingMuonEta(double fLeadingMuonEta)  { LeadingMuonEta_     = fLeadingMuonEta;}
    void SetLeadingMuonPhi(double fLeadingMuonPhi)  { LeadingMuonPhi_    = fLeadingMuonPhi;}
    void SetLeadingMuonP4(LorentzVector fLeadingMuonP4)    { LeadingMuonP4_     = fLeadingMuonP4;}
    void SetLeadingMuonCharge(int fLeadingMuonCharge)  { LeadingMuonCharge_     = fLeadingMuonCharge;}
    void SetSecondMuonPt(double fSecondMuonPt)    { SecondMuonPt_     = fSecondMuonPt;}
    void SetSecondMuonEta(double fSecondMuonEta)  { SecondMuonEta_     = fSecondMuonEta;}
    void SetSecondMuonPhi(double fSecondMuonPhi)  { SecondMuonPhi_    = fSecondMuonPhi;}
    void SetSecondMuonP4(LorentzVector fSecondMuonP4)    { SecondMuonP4_     = fSecondMuonP4;}
    void SetSecondMuonCharge(int fSecondMuonCharge)  { SecondMuonCharge_     = fSecondMuonCharge;}
    void SetMuonsN(int fMuonsN)  { MuonsN_    = fMuonsN;}

    void SetLeadingMuonSumPtR03(double fLeadingMuonSumPtR03)    {LeadingMuonSumPtR03_ = fLeadingMuonSumPtR03;}
    void SetLeadingMuonEmEtR03(double fLeadingMuonEmEtR03)    {LeadingMuonEmEtR03_ = fLeadingMuonEmEtR03;}
    void SetLeadingMuonHadEtR03(double fLeadingMuonHadEtR03)    {LeadingMuonHadEtR03_ = fLeadingMuonHadEtR03;}
    void SetLeadingMuonSumPtR05(double fLeadingMuonSumPtR05)    {LeadingMuonSumPtR05_ = fLeadingMuonSumPtR05;}
    void SetLeadingMuonEmEtR05(double fLeadingMuonEmEtR05)    {LeadingMuonEmEtR05_ = fLeadingMuonEmEtR05;}
    void SetLeadingMuonHadEtR05(double fLeadingMuonHadEtR05)    {LeadingMuonHadEtR05_ = fLeadingMuonHadEtR05;}

    void SetSecondMuonSumPtR03(double fSecondMuonSumPtR03)    {SecondMuonSumPtR03_ = fSecondMuonSumPtR03;}
    void SetSecondMuonEmEtR03(double fSecondMuonEmEtR03)    {SecondMuonEmEtR03_ = fSecondMuonEmEtR03;}
    void SetSecondMuonHadEtR03(double fSecondMuonHadEtR03)    {SecondMuonHadEtR03_ = fSecondMuonHadEtR03;}
    void SetSecondMuonSumPtR05(double fSecondMuonSumPtR05)    {SecondMuonSumPtR05_ = fSecondMuonSumPtR05;}
    void SetSecondMuonEmEtR05(double fSecondMuonEmEtR05)    {SecondMuonEmEtR05_ = fSecondMuonEmEtR05;}
    void SetSecondMuonHadEtR05(double fSecondMuonHadEtR05)    {SecondMuonHadEtR05_ = fSecondMuonHadEtR05;}

    void SetLeadingMuonrelIsoDr03(double fLeadingMuonrelIsoDr03)    {LeadingMuonrelIsoDr03_ = fLeadingMuonrelIsoDr03;}
    void SetSecondMuonrelIsoDr03(double fSecondMuonrelIsoDr03)    {SecondMuonrelIsoDr03_ = fSecondMuonrelIsoDr03;}
    void SetLeadingMuonrelIsoDr05(double fLeadingMuonrelIsoDr05)    {LeadingMuonrelIsoDr05_ = fLeadingMuonrelIsoDr05;}
    void SetSecondMuonrelIsoDr05(double fSecondMuonrelIsoDr05)    {SecondMuonrelIsoDr05_ = fSecondMuonrelIsoDr05;}

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

    void SetDiElectronMassPF(double fDiElectronMassPF) { DiElectronMassPF_ = fDiElectronMassPF;}
    void SetDiMuonMassPF(double fDiMuonMassPF) { DiMuonMassPF_ = fDiMuonMassPF;}

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

    void SetMxPFNoZMinus(double fmxpfnozminus)    { mxpfnozminus_     = fmxpfnozminus;} 
    void SetMxPFNoZPlus(double fmxpfnozplus)    { mxpfnozplus_     = fmxpfnozplus;}  
    void SetMx2PFNoZMinus(double fmx2pfnozminus)    { mx2pfnozminus_     = fmx2pfnozminus;}
    void SetMx2PFNoZPlus(double fmx2pfnozplus)    { mx2pfnozplus_     = fmx2pfnozplus;}
    void SetEtaMaxPFNoZ(double fetamaxpfnoz)    { etamaxpfnoz_     = fetamaxpfnoz;}
    void SetEtaMinPFNoZ(double fetaminpfnoz)    { etaminpfnoz_     = fetaminpfnoz;}
    void SetEpluspzPFNoZ(double fepluspzpfnoz)    { epluspzpfnoz_     = fepluspzpfnoz;}
    void SetEminuspzPFNoZ(double feminuspzpfnoz)    { eminuspzpfnoz_     = feminuspzpfnoz;}
    void SetEtExpoPlusPFNoZ(double fetexpopluspfnoz)    { etexpopluspfnoz_     = fetexpopluspfnoz;}
    void SetEtExpoMinusPFNoZ(double fetexpominuspfnoz)    { etexpominuspfnoz_     = fetexpominuspfnoz;}
    void SetLrgPFNoZ(double flrgPFnoz)    { lrgPFnoz_     = flrgPFnoz;}
    void SetMxPFNoZLeft(double fmxpfnozleft)    { mxpfnozleft_     = fmxpfnozleft;}
    void SetMxPFNoZRight(double fmxpfnozright)    { mxpfnozright_     = fmxpfnozright;}
    void SetMx2PFNoZLeft(double fmx2pfnozleft)    { mx2pfnozleft_     = fmx2pfnozleft;}
    void SetMx2PFNoZRight(double fmx2pfnozright)    { mx2pfnozright_     = fmx2pfnozright;}
    void SetSumptPFNoZLeft(double fsumptpfnozleft)    { sumptpfnozleft_     = fsumptpfnozleft;}
    void SetSumptPFNoZRight(double fsumptpfnozright)    { sumptpfnozright_     = fsumptpfnozright;}

    void SetPatNMuon(int fpatNMuon)    {patNMuon_ = fpatNMuon;}
    void SetPatMuon1Pt(double fpatMuon1Pt)    {patMuon1Pt_ = fpatMuon1Pt;}
    void SetPatMuon1Charge(int fpatMuon1Charge)    {patMuon1Charge_ = fpatMuon1Charge;}
    void SetPatMuon1Phi(double fpatMuon1Phi)    {patMuon1Phi_ = fpatMuon1Phi;}
    void SetPatMuon1Eta(double fpatMuon1Eta)    {patMuon1Eta_ = fpatMuon1Eta;}
    void SetPatMuon1Et(double fpatMuon1Et)    {patMuon1Et_ = fpatMuon1Et;}
    void SetPatMuon1P4(LorentzVector fpatMuon1P4)    {patMuon1P4_ = fpatMuon1P4;}

    void SetPatMuon2Pt(double fpatMuon2Pt)    {patMuon2Pt_ = fpatMuon2Pt;}
    void SetPatMuon2Charge(int fpatMuon2Charge)    {patMuon2Charge_ = fpatMuon2Charge;}
    void SetPatMuon2Phi(double fpatMuon2Phi)    {patMuon2Phi_ = fpatMuon2Phi;}
    void SetPatMuon2Eta(double fpatMuon2Eta)    {patMuon2Eta_ = fpatMuon2Eta;}
    void SetPatMuon2Et(double fpatMuon2Et)    {patMuon2Et_ = fpatMuon2Et;}
    void SetPatMuon2P4(LorentzVector fpatMuon2P4)    {patMuon2P4_ = fpatMuon2P4;}

    void SetPatMuon1SumPtR03(double fpatMuon1SumPtR03)    {patMuon1SumPtR03_ = fpatMuon1SumPtR03;}
    void SetPatMuon1EmEtR03(double fpatMuon1EmEtR03)    {patMuon1EmEtR03_ = fpatMuon1EmEtR03;}
    void SetPatMuon1HadEtR03(double fpatMuon1HadEtR03)    {patMuon1HadEtR03_ = fpatMuon1HadEtR03;}    
    void SetPatMuon1SumPtR05(double fpatMuon1SumPtR05)    {patMuon1SumPtR05_ = fpatMuon1SumPtR05;}
    void SetPatMuon1EmEtR05(double fpatMuon1EmEtR05)    {patMuon1EmEtR05_ = fpatMuon1EmEtR05;}
    void SetPatMuon1HadEtR05(double fpatMuon1HadEtR05)    {patMuon1HadEtR05_ = fpatMuon1HadEtR05;}    

    void SetPatMuon2SumPtR03(double fpatMuon2SumPtR03)    {patMuon2SumPtR03_ = fpatMuon2SumPtR03;}
    void SetPatMuon2EmEtR03(double fpatMuon2EmEtR03)    {patMuon2EmEtR03_ = fpatMuon2EmEtR03;}
    void SetPatMuon2HadEtR03(double fpatMuon2HadEtR03)    {patMuon2HadEtR03_ = fpatMuon2HadEtR03;}    
    void SetPatMuon2SumPtR05(double fpatMuon2SumPtR05)    {patMuon2SumPtR05_ = fpatMuon2SumPtR05;}
    void SetPatMuon2EmEtR05(double fpatMuon2EmEtR05)    {patMuon2EmEtR05_ = fpatMuon2EmEtR05;}
    void SetPatMuon2HadEtR05(double fpatMuon2HadEtR05)    {patMuon2HadEtR05_ = fpatMuon2HadEtR05;}  

    void SetPatMuon1relIsoDr03(double fpatMuon1relIsoDr03)    {patMuon1relIsoDr03_ = fpatMuon1relIsoDr03;}
    void SetPatMuon2relIsoDr03(double fpatMuon2relIsoDr03)    {patMuon2relIsoDr03_ = fpatMuon2relIsoDr03;}
    void SetPatMuon1relIsoDr05(double fpatMuon1relIsoDr05)    {patMuon1relIsoDr05_ = fpatMuon1relIsoDr05;}
    void SetPatMuon2relIsoDr05(double fpatMuon2relIsoDr05)    {patMuon2relIsoDr05_ = fpatMuon2relIsoDr05;}

    void SetPatMuon1relIso(double fpatMuon1relIso)    {patMuon1relIso_ = fpatMuon1relIso;}
    void SetPatMuon2relIso(double fpatMuon2relIso)    {patMuon2relIso_ = fpatMuon2relIso;}

    void SetPatDiMuonMass(double fpatDiMuonMass) { patDiMuonMass_ = fpatDiMuonMass;}
    void SetPatDiMuonPt(double fpatDiMuonPt) { patDiMuonPt_ = fpatDiMuonPt;}
    void SetPatDiMuonEta(double fpatDiMuonEta) { patDiMuonEta_ = fpatDiMuonEta;}
    void SetPatDiMuonPhi(double fpatDiMuonPhi) { patDiMuonPhi_ = fpatDiMuonPhi;}

    void SetPatNElectron(int fpatNElectron)    {patNElectron_ = fpatNElectron;}
    void SetPatElectron1Pt(double fpatElectron1Pt)    {patElectron1Pt_ = fpatElectron1Pt;}
    void SetPatElectron1Charge(int fpatElectron1Charge)    {patElectron1Charge_ = fpatElectron1Charge;}
    void SetPatElectron1Phi(double fpatElectron1Phi)    {patElectron1Phi_ = fpatElectron1Phi;}
    void SetPatElectron1Eta(double fpatElectron1Eta)    {patElectron1Eta_ = fpatElectron1Eta;}
    void SetPatElectron1Et(double fpatElectron1Et)    {patElectron1Et_ = fpatElectron1Et;}
    void SetPatElectron1P4(LorentzVector fpatElectron1P4)    {patElectron1P4_ = fpatElectron1P4;}

    void SetPatElectron2Pt(double fpatElectron2Pt)    {patElectron2Pt_ = fpatElectron2Pt;}
    void SetPatElectron2Charge(int fpatElectron2Charge)    {patElectron2Charge_ = fpatElectron2Charge;}
    void SetPatElectron2Phi(double fpatElectron2Phi)    {patElectron2Phi_ = fpatElectron2Phi;}
    void SetPatElectron2Eta(double fpatElectron2Eta)    {patElectron2Eta_ = fpatElectron2Eta;}
    void SetPatElectron2Et(double fpatElectron2Et)    {patElectron2Et_ = fpatElectron2Et;}
    void SetPatElectron2P4(LorentzVector fpatElectron2P4)    {patElectron2P4_ = fpatElectron2P4;}

    void SetPatElectron1TkDr03(double fpatElectron1TkDr03)    {patElectron1TkDr03_ = fpatElectron1TkDr03;}
    void SetPatElectron1EcalDr03(double fpatElectron1EcalDr03)    {patElectron1EcalDr03_ = fpatElectron1EcalDr03;}
    void SetPatElectron1HcalDr03(double fpatElectron1HcalDr03)    {patElectron1HcalDr03_ = fpatElectron1HcalDr03;}
    void SetPatElectron2TkDr03(double fpatElectron2TkDr03)    {patElectron2TkDr03_ = fpatElectron2TkDr03;}
    void SetPatElectron2EcalDr03(double fpatElectron2EcalDr03)    {patElectron2EcalDr03_ = fpatElectron2EcalDr03;}
    void SetPatElectron2HcalDr03(double fpatElectron2HcalDr03)    {patElectron2HcalDr03_ = fpatElectron2HcalDr03;}

    void SetPatElectron1TkDr04(double fpatElectron1TkDr04)    {patElectron1TkDr04_ = fpatElectron1TkDr04;}
    void SetPatElectron1EcalDr04(double fpatElectron1EcalDr04)    {patElectron1EcalDr04_ = fpatElectron1EcalDr04;}
    void SetPatElectron1HcalDr04(double fpatElectron1HcalDr04)    {patElectron1HcalDr04_ = fpatElectron1HcalDr04;}
    void SetPatElectron2TkDr04(double fpatElectron2TkDr04)    {patElectron2TkDr04_ = fpatElectron2TkDr04;}
    void SetPatElectron2EcalDr04(double fpatElectron2EcalDr04)    {patElectron2EcalDr04_ = fpatElectron2EcalDr04;}
    void SetPatElectron2HcalDr04(double fpatElectron2HcalDr04)    {patElectron2HcalDr04_ = fpatElectron2HcalDr04;}

    void SetPatElectron1relIsoDr03(double fpatElectron1relIsoDr03)    {patElectron1relIsoDr03_ = fpatElectron1relIsoDr03;}
    void SetPatElectron1relIsoDr04(double fpatElectron1relIsoDr04)    {patElectron1relIsoDr04_ = fpatElectron1relIsoDr04;}
    void SetPatElectron2relIsoDr03(double fpatElectron2relIsoDr03)    {patElectron2relIsoDr03_ = fpatElectron2relIsoDr03;}
    void SetPatElectron2relIsoDr04(double fpatElectron2relIsoDr04)    {patElectron2relIsoDr04_ = fpatElectron2relIsoDr04;}

    void SetPatDiElectronMass(double fpatDiElectronMass) { patDiElectronMass_ = fpatDiElectronMass;}
    void SetPatDiElectronPt(double fpatDiElectronPt) { patDiElectronPt_ = fpatDiElectronPt;}
    void SetPatDiElectronEta(double fpatDiElectronEta) { patDiElectronEta_ = fpatDiElectronEta;}
    void SetPatDiElectronPhi(double fpatDiElectronPhi) { patDiElectronPhi_ = fpatDiElectronPhi;}

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

    void SetSecondElectronDeltaPhiTkClu(double fSecondElectronDeltaPhiTkClu)    {SecondElectronDeltaPhiTkClu_ = fSecondElectronDeltaPhiTkClu;}
    void SetSecondElectronDeltaEtaTkClu(double fSecondElectronDeltaEtaTkClu)    {SecondElectronDeltaEtaTkClu_ = fSecondElectronDeltaEtaTkClu;}
    void SetSecondElectronSigmaIeIe(double fSecondElectronSigmaIeIe)    {SecondElectronSigmaIeIe_ = fSecondElectronSigmaIeIe;}
    void SetSecondElectronDCot(double fSecondElectronDCot)    {SecondElectronDCot_ = fSecondElectronDCot;}
    void SetSecondElectronDist(double fSecondElectronDist)    {SecondElectronDist_ = fSecondElectronDist;}
    void SetSecondElectronInnerHits(double fSecondElectronInnerHits)    {SecondElectronInnerHits_ = fSecondElectronInnerHits;}
    void SetSecondElectronHE(double fSecondElectronHE)    {SecondElectronHE_ = fSecondElectronHE;}

    void SetPatElectron1DeltaPhiTkClu(double fpatElectron1DeltaPhiTkClu)    {patElectron1DeltaPhiTkClu_ = fpatElectron1DeltaPhiTkClu;}
    void SetPatElectron1DeltaEtaTkClu(double fpatElectron1DeltaEtaTkClu)    {patElectron1DeltaEtaTkClu_ = fpatElectron1DeltaEtaTkClu;}
    void SetPatElectron1SigmaIeIe(double fpatElectron1SigmaIeIe)    {patElectron1SigmaIeIe_ = fpatElectron1SigmaIeIe;}
    void SetPatElectron1DCot(double fpatElectron1DCot)    {patElectron1DCot_ = fpatElectron1DCot;}
    void SetPatElectron1Dist(double fpatElectron1Dist)    {patElectron1Dist_ = fpatElectron1Dist;}
    void SetPatElectron1InnerHits(double fpatElectron1InnerHits)    {patElectron1InnerHits_ = fpatElectron1InnerHits;}
    void SetPatElectron1HE(double fpatElectron1HE)    {patElectron1HE_ = fpatElectron1HE;}

    void SetPatElectron2DeltaPhiTkClu(double fpatElectron2DeltaPhiTkClu)    {patElectron2DeltaPhiTkClu_ = fpatElectron2DeltaPhiTkClu;}
    void SetPatElectron2DeltaEtaTkClu(double fpatElectron2DeltaEtaTkClu)    {patElectron2DeltaEtaTkClu_ = fpatElectron2DeltaEtaTkClu;}
    void SetPatElectron2SigmaIeIe(double fpatElectron2SigmaIeIe)    {patElectron2SigmaIeIe_ = fpatElectron2SigmaIeIe;}
    void SetPatElectron2DCot(double fpatElectron2DCot)    {patElectron2DCot_ = fpatElectron2DCot;}
    void SetPatElectron2Dist(double fpatElectron2Dist)    {patElectron2Dist_ = fpatElectron2Dist;}
    void SetPatElectron2InnerHits(double fpatElectron2InnerHits)    {patElectron2InnerHits_ = fpatElectron2InnerHits;}
    void SetPatElectron2HE(double fpatElectron2HE)    {patElectron2HE_ = fpatElectron2HE;}

    int GetHLTPath(int idx)                    const { return hltTrigResults_[idx]; }
    double GetDiElectronMass() const {return DiElectronMass_;}
    double GetDiElectronPt() const {return DiElectronPt_;}
    double GetDiElectronEta() const {return DiElectronEta_;}
    double GetDiElectronPhi() const {return DiElectronPhi_;}
    double GetDiMuonMass() const {return DiMuonMass_;}
    double GetDiMuonPt() const {return DiMuonPt_;}
    double GetDiMuonEta() const {return DiMuonEta_;}
    double GetDiMuonPhi() const {return DiMuonPhi_;}

    double GetLeadingElectronPt() const {return LeadingElectronPt_;}
    double GetLeadingElectronEta() const {return LeadingElectronEta_;}
    double GetLeadingElectronPhi() const {return LeadingElectronPhi_;}
    const LorentzVector& GetLeadingElectronP4() const {return LeadingElectronP4_;}
    int GetLeadingElectronCharge() const {return LeadingElectronCharge_;}
    double GetSecondElectronPt() const {return SecondElectronPt_;}
    double GetSecondElectronEta() const {return SecondElectronEta_;}
    double GetSecondElectronPhi() const {return SecondElectronPhi_;}
    const LorentzVector& GetSecondElectronP4() const {return SecondElectronP4_;}
    int GetSecondElectronCharge() const {return SecondElectronCharge_;}
    int GetElectronsN() const {return ElectronsN_;}
    double GetLeadingElectronTkDr03() const  {return LeadingElectronTkDr03_;}
    double GetLeadingElectronEcalDr03() const  {return LeadingElectronEcalDr03_;}
    double GetLeadingElectronHcalDr03() const  {return LeadingElectronHcalDr03_;}
    double GetSecondElectronTkDr03() const  {return SecondElectronTkDr03_;}
    double GetSecondElectronEcalDr03() const  {return SecondElectronEcalDr03_;}
    double GetSecondElectronHcalDr03() const  {return SecondElectronHcalDr03_;}

    double GetLeadingElectronTkDr04() const  {return LeadingElectronTkDr04_;}
    double GetLeadingElectronEcalDr04() const  {return LeadingElectronEcalDr04_;}
    double GetLeadingElectronHcalDr04() const  {return LeadingElectronHcalDr04_;}
    double GetSecondElectronTkDr04() const  {return SecondElectronTkDr04_;}
    double GetSecondElectronEcalDr04() const  {return SecondElectronEcalDr04_;}
    double GetSecondElectronHcalDr04() const  {return SecondElectronHcalDr04_;}

    double GetLeadingElectronrelIsoDr03() const {return LeadingElectronrelIsoDr03_;}
    double GetLeadingElectronrelIsoDr04() const {return LeadingElectronrelIsoDr04_;}
    double GetSecondElectronrelIsoDr03() const {return SecondElectronrelIsoDr03_;}
    double GetSecondElectronrelIsoDr04() const {return SecondElectronrelIsoDr04_;}

    double GetLeadingMuonPt() const {return LeadingMuonPt_;}
    double GetLeadingMuonEta() const {return LeadingMuonEta_;}
    double GetLeadingMuonPhi() const {return LeadingMuonPhi_;}
    const LorentzVector& GetLeadingMuonP4() const {return LeadingMuonP4_;}
    int GetLeadingMuonCharge() const {return LeadingMuonCharge_;}
    double GetSecondMuonPt() const {return SecondMuonPt_;}
    double GetSecondMuonEta() const {return SecondMuonEta_;}
    double GetSecondMuonPhi() const {return SecondMuonPhi_;}
    const LorentzVector& GetSecondMuonP4() const {return SecondMuonP4_;}
    int GetSecondMuonCharge() const {return SecondMuonCharge_;}
    int GetMuonsN() const {return MuonsN_;}

    double GetLeadingMuonSumPtR03() const {return LeadingMuonSumPtR03_;}
    double GetLeadingMuonEmEtR03() const {return LeadingMuonEmEtR03_;}
    double GetLeadingMuonHadEtR03() const {return LeadingMuonHadEtR03_;}
    double GetLeadingMuonSumPtR05() const {return LeadingMuonSumPtR05_;}
    double GetLeadingMuonEmEtR05() const {return LeadingMuonEmEtR05_;}
    double GetLeadingMuonHadEtR05() const {return LeadingMuonHadEtR05_;}

    double GetSecondMuonSumPtR03() const {return SecondMuonSumPtR03_;}
    double GetSecondMuonEmEtR03() const {return SecondMuonEmEtR03_;}
    double GetSecondMuonHadEtR03() const {return SecondMuonHadEtR03_;}
    double GetSecondMuonSumPtR05() const {return SecondMuonSumPtR05_;}
    double GetSecondMuonEmEtR05() const {return SecondMuonEmEtR05_;}
    double GetSecondMuonHadEtR05() const {return SecondMuonHadEtR05_;}

    double GetLeadingMuonrelIsoDr03() const {return LeadingMuonrelIsoDr03_;}
    double GetSecondMuonrelIsoDr03() const {return SecondMuonrelIsoDr03_;}
    double GetLeadingMuonrelIsoDr05() const {return LeadingMuonrelIsoDr05_;}
    double GetSecondMuonrelIsoDr05() const {return SecondMuonrelIsoDr05_;}

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

    double GetDiElectronMassPF() const {return DiElectronMassPF_;}
    double GetDiMuonMassPF() const {return DiMuonMassPF_;}

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

    double GetMxPFNoZMinus() const { return mxpfnozminus_;}
    double GetMxPFNoZPlus() const { return mxpfnozplus_;}
    double GetMx2PFNoZMinus() const { return mx2pfnozminus_;}
    double GetMx2PFNoZPlus() const { return mx2pfnozplus_;}
    double GetEtaMaxPFNoZ() const { return etamaxpfnoz_;}
    double GetEtaMinPFNoZ() const { return etaminpfnoz_;}
    double GetEpluspzPFNoZ() const { return epluspzpfnoz_;}
    double GetEminuspzPFNoZ() const { return eminuspzpfnoz_;}
    double GetEtExpoPlusPFNoZ() const { return etexpopluspfnoz_;}
    double GetEtExpoMinusPFNoZ() const { return etexpominuspfnoz_;}
    double GetLrgPFNoZ() const { return lrgPFnoz_;}
    double GetMxPFNoZLeft() const { return mxpfnozleft_;}
    double GetMxPFNoZRight() const { return mxpfnozright_;}
    double GetMx2PFNoZLeft() const { return mx2pfnozleft_;}
    double GetMx2PFNoZRight() const { return mx2pfnozright_;}
    double GetSumptPFNoZLeft() const { return sumptpfnozleft_;}
    double GetSumptPFNoZRight() const { return sumptpfnozright_;}

    int GetPatNMuon() const {return patNMuon_;}
    double GetPatMuon1Pt() const {return patMuon1Pt_;}
    int GetPatMuon1Charge() const {return patMuon1Charge_;}
    double GetPatMuon1Phi() const {return patMuon1Phi_;}
    double GetPatMuon1Eta() const {return patMuon1Eta_;}
    double GetPatMuon1Et() const {return patMuon1Et_;}
    const LorentzVector& GetPatMuon1P4() const {return patMuon1P4_;}

    double GetPatMuon2Pt() const {return patMuon2Pt_;}
    int GetPatMuon2Charge() const {return patMuon2Charge_;}
    double GetPatMuon2Phi() const {return patMuon2Phi_;}
    double GetPatMuon2Eta() const {return patMuon2Eta_;}
    double GetPatMuon2Et() const {return patMuon2Et_;}
    const LorentzVector& GetPatMuon2P4() const {return patMuon2P4_;}

    double GetPatMuon1SumPtR03() const {return patMuon1SumPtR03_;}
    double GetPatMuon1EmEtR03() const {return patMuon1EmEtR03_;}
    double GetPatMuon1HadEtR03() const {return patMuon1HadEtR03_;}    
    double GetPatMuon1SumPtR05() const {return patMuon1SumPtR05_;}
    double GetPatMuon1EmEtR05() const {return patMuon1EmEtR05_;}
    double GetPatMuon1HadEtR05() const {return patMuon1HadEtR05_;}    

    double GetPatMuon2SumPtR03() const {return patMuon2SumPtR03_;}
    double GetPatMuon2EmEtR03() const {return patMuon2EmEtR03_;}
    double GetPatMuon2HadEtR03() const {return patMuon2HadEtR03_;}    
    double GetPatMuon2SumPtR05() const {return patMuon2SumPtR05_;}
    double GetPatMuon2EmEtR05() const {return patMuon2EmEtR05_;}
    double GetPatMuon2HadEtR05() const {return patMuon2HadEtR05_;}  

    double GetPatMuon1relIsoDr03() const {return patMuon1relIsoDr03_;}
    double GetPatMuon2relIsoDr03() const {return patMuon2relIsoDr03_;}
    double GetPatMuon1relIsoDr05() const {return patMuon1relIsoDr05_;}
    double GetPatMuon2relIsoDr05() const {return patMuon2relIsoDr05_;}

    double GetPatMuon1relIso() const {return patMuon1relIso_;}
    double GetPatMuon2relIso() const {return patMuon2relIso_;}

    double GetPatDiMuonMass() const {return patDiMuonMass_;}
    double GetPatDiMuonPt() const {return patDiMuonPt_;}
    double GetPatDiMuonEta() const {return patDiMuonEta_;}
    double GetPatDiMuonPhi() const {return patDiMuonPhi_;}

    int GetPatNElectron() const {return patNElectron_;}
    double GetPatElectron1Pt() const {return patElectron1Pt_;}
    int GetPatElectron1Charge() const {return patElectron1Charge_;}
    double GetPatElectron1Phi() const {return patElectron1Phi_;}
    double GetPatElectron1Eta() const {return patElectron1Eta_;}
    double GetPatElectron1Et() const {return patElectron1Et_;}
    const LorentzVector& GetPatElectron1P4() const {return patElectron1P4_;}

    double GetPatElectron2Pt() const {return patElectron2Pt_;}
    int GetPatElectron2Charge() const {return patElectron2Charge_;}
    double GetPatElectron2Phi() const {return patElectron2Phi_;}
    double GetPatElectron2Eta() const {return patElectron2Eta_;}
    double GetPatElectron2Et() const {return patElectron2Et_;}
    const LorentzVector& GetPatElectron2P4() const {return patElectron2P4_;}

    double GetPatElectron1TkDr03() const  {return patElectron1TkDr03_;}    
    double GetPatElectron1EcalDr03() const  {return patElectron1EcalDr03_;}
    double GetPatElectron1HcalDr03() const  {return patElectron1HcalDr03_;}
    double GetPatElectron2TkDr03() const  {return patElectron2TkDr03_;}
    double GetPatElectron2EcalDr03() const  {return patElectron2EcalDr03_;}
    double GetPatElectron2HcalDr03() const  {return patElectron2HcalDr03_;}

    double GetPatElectron1TkDr04() const  {return patElectron1TkDr04_;}
    double GetPatElectron1EcalDr04() const  {return patElectron1EcalDr04_;}
    double GetPatElectron1HcalDr04() const  {return patElectron1HcalDr04_;}
    double GetPatElectron2TkDr04() const  {return patElectron2TkDr04_;}
    double GetPatElectron2EcalDr04() const  {return patElectron2EcalDr04_;}
    double GetPatElectron2HcalDr04() const  {return patElectron2HcalDr04_;}

    double GetPatElectron1relIsoDr03() const {return patElectron1relIsoDr03_;}
    double GetPatElectron1relIsoDr04() const {return patElectron1relIsoDr04_;}
    double GetPatElectron2relIsoDr03() const {return patElectron2relIsoDr03_;}
    double GetPatElectron2relIsoDr04() const {return patElectron2relIsoDr04_;}

    double GetPatDiElectronMass() const {return patDiElectronMass_;}
    double GetPatDiElectronPt() const {return patDiElectronPt_;}
    double GetPatDiElectronEta() const {return patDiElectronEta_;}
    double GetPatDiElectronPhi() const {return patDiElectronPhi_;}

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

    double GetSecondElectronDeltaPhiTkClu() const {return SecondElectronDeltaPhiTkClu_;}
    double GetSecondElectronDeltaEtaTkClu() const {return SecondElectronDeltaEtaTkClu_;}
    double GetSecondElectronSigmaIeIe() const {return SecondElectronSigmaIeIe_;}
    double GetSecondElectronDCot() const {return SecondElectronDCot_;}
    double GetSecondElectronDist() const {return SecondElectronDist_;}
    double GetSecondElectronInnerHits() const {return SecondElectronInnerHits_;}
    double GetSecondElectronHE() const {return SecondElectronHE_;}

    double GetPatElectron1DeltaPhiTkClu() const {return patElectron1DeltaPhiTkClu_;}
    double GetPatElectron1DeltaEtaTkClu() const {return patElectron1DeltaEtaTkClu_;}
    double GetPatElectron1SigmaIeIe() const {return patElectron1SigmaIeIe_;}
    double GetPatElectron1DCot() const {return patElectron1DCot_;}
    double GetPatElectron1Dist() const {return patElectron1Dist_;}
    double GetPatElectron1InnerHits() const {return patElectron1InnerHits_;}
    double GetPatElectron1HE() const {return patElectron1HE_;}

    double GetPatElectron2DeltaPhiTkClu() const {return patElectron2DeltaPhiTkClu_;}
    double GetPatElectron2DeltaEtaTkClu() const {return patElectron2DeltaEtaTkClu_;}
    double GetPatElectron2SigmaIeIe() const {return patElectron2SigmaIeIe_;}
    double GetPatElectron2DCot() const {return patElectron2DCot_;}
    double GetPatElectron2Dist() const {return patElectron2Dist_;}
    double GetPatElectron2InnerHits() const {return patElectron2InnerHits_;}
    double GetPatElectron2HE() const {return patElectron2HE_;}

  private:
    friend class diffractiveZAnalysis::DiffractiveZAnalysis;

    void reset();

    int hltTrigResults_[20];
    double DiElectronMass_;
    double DiElectronPt_;
    double DiElectronEta_;
    double DiElectronPhi_;

    double DiMuonMass_;
    double DiMuonPt_;
    double DiMuonEta_;
    double DiMuonPhi_;

    double LeadingElectronPt_;
    double LeadingElectronEta_;
    double LeadingElectronPhi_;
    LorentzVector LeadingElectronP4_;
    int LeadingElectronCharge_;
    double SecondElectronPt_;
    double SecondElectronEta_;
    double SecondElectronPhi_;
    LorentzVector SecondElectronP4_;
    int SecondElectronCharge_;
    int ElectronsN_;

    double LeadingElectronTkDr03_;
    double LeadingElectronEcalDr03_;
    double LeadingElectronHcalDr03_;
    double SecondElectronTkDr03_;
    double SecondElectronEcalDr03_;
    double SecondElectronHcalDr03_;

    double LeadingElectronTkDr04_;
    double LeadingElectronEcalDr04_;
    double LeadingElectronHcalDr04_;
    double SecondElectronTkDr04_;
    double SecondElectronEcalDr04_;
    double SecondElectronHcalDr04_;

    double LeadingElectronrelIsoDr03_;
    double LeadingElectronrelIsoDr04_;
    double SecondElectronrelIsoDr03_;
    double SecondElectronrelIsoDr04_;

    double LeadingMuonPt_;
    double LeadingMuonEta_;
    double LeadingMuonPhi_;
    LorentzVector LeadingMuonP4_;
    int LeadingMuonCharge_;
    double SecondMuonPt_;
    double SecondMuonEta_;
    double SecondMuonPhi_;
    LorentzVector SecondMuonP4_;
    int SecondMuonCharge_;
    int MuonsN_;

    double LeadingMuonSumPtR03_;
    double LeadingMuonEmEtR03_;
    double LeadingMuonHadEtR03_;
    double LeadingMuonSumPtR05_;
    double LeadingMuonEmEtR05_;
    double LeadingMuonHadEtR05_;

    double SecondMuonSumPtR03_;
    double SecondMuonEmEtR03_;
    double SecondMuonHadEtR03_;
    double SecondMuonSumPtR05_;
    double SecondMuonEmEtR05_;
    double SecondMuonHadEtR05_;

    double LeadingMuonrelIsoDr03_;
    double SecondMuonrelIsoDr03_;
    double LeadingMuonrelIsoDr05_;
    double SecondMuonrelIsoDr05_;

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

    double DiElectronMassPF_;
    double DiMuonMassPF_;

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

    double mxpfnozminus_;
    double mxpfnozplus_;
    double mx2pfnozminus_;
    double mx2pfnozplus_;
    double etamaxpfnoz_;
    double etaminpfnoz_;
    double epluspzpfnoz_;
    double eminuspzpfnoz_;
    double etexpopluspfnoz_;
    double etexpominuspfnoz_;
    double lrgPFnoz_;
    double mxpfnozleft_;
    double mxpfnozright_;
    double mx2pfnozleft_;
    double mx2pfnozright_;
    double sumptpfnozleft_;
    double sumptpfnozright_;

    int patNMuon_;

    double patMuon1Pt_;
    int patMuon1Charge_;
    double patMuon1Phi_;
    double patMuon1Eta_;
    double patMuon1Et_;
    LorentzVector patMuon1P4_;

    double patMuon2Pt_;
    int patMuon2Charge_;
    double patMuon2Phi_ ;
    double patMuon2Eta_;
    double patMuon2Et_;
    LorentzVector patMuon2P4_;

    double patMuon1SumPtR03_;
    double patMuon1EmEtR03_;
    double patMuon1HadEtR03_;   
    double patMuon1SumPtR05_;
    double patMuon1EmEtR05_;
    double patMuon1HadEtR05_;   

    double patMuon2SumPtR03_;
    double patMuon2EmEtR03_;
    double patMuon2HadEtR03_;    
    double patMuon2SumPtR05_;
    double patMuon2EmEtR05_;
    double patMuon2HadEtR05_;  

    double patMuon1relIsoDr03_;
    double patMuon2relIsoDr03_;
    double patMuon1relIsoDr05_;
    double patMuon2relIsoDr05_;

    double patMuon1relIso_;
    double patMuon2relIso_;

    double patDiMuonMass_;
    double patDiMuonPt_;
    double patDiMuonEta_;
    double patDiMuonPhi_;

    double patDiElectronMass_;
    double patDiElectronPt_;
    double patDiElectronEta_;
    double patDiElectronPhi_;

    int patNElectron_; 

    double patElectron1Pt_;
    int patElectron1Charge_;
    double patElectron1Phi_;
    double patElectron1Eta_;
    double patElectron1Et_;
    LorentzVector patElectron1P4_;

    double patElectron2Pt_;
    int patElectron2Charge_;
    double patElectron2Phi_;
    double patElectron2Eta_;
    double patElectron2Et_;
    LorentzVector patElectron2P4_;

    double patElectron1TkDr03_;    
    double patElectron1EcalDr03_;
    double patElectron1HcalDr03_;
    double patElectron2TkDr03_;
    double patElectron2EcalDr03_;
    double patElectron2HcalDr03_;

    double patElectron1TkDr04_;
    double patElectron1EcalDr04_;
    double patElectron1HcalDr04_;
    double patElectron2TkDr04_;
    double patElectron2EcalDr04_;
    double patElectron2HcalDr04_;

    double patElectron1relIsoDr03_;
    double patElectron1relIsoDr04_;
    double patElectron2relIsoDr03_;
    double patElectron2relIsoDr04_;

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
    double SecondElectronDeltaPhiTkClu_;
    double SecondElectronDeltaEtaTkClu_;
    double SecondElectronSigmaIeIe_;
    double SecondElectronDCot_ ;
    double SecondElectronDist_;
    double SecondElectronInnerHits_;
    double SecondElectronHE_;
    double patElectron1DeltaPhiTkClu_;
    double patElectron1DeltaEtaTkClu_;
    double patElectron1SigmaIeIe_;
    double patElectron1DCot_;
    double patElectron1Dist_;
    double patElectron1InnerHits_;
    double patElectron1HE_;
    double patElectron2DeltaPhiTkClu_;
    double patElectron2DeltaEtaTkClu_;
    double patElectron2SigmaIeIe_;
    double patElectron2DCot_;
    double patElectron2Dist_;
    double patElectron2InnerHits_;
    double patElectron2HE_;    

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
