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

    void SetVertexMultiplicity(const std::vector<double>& fVertexMultiplicity) { VertexMultiplicity_ = fVertexMultiplicity; }
    void SetVertexChiNorm(const std::vector<double>& fVertexChiNorm) { VertexChiNorm_ = fVertexChiNorm; }
    void SetVertexNDOF(const std::vector<double>& fVertexNDOF) { VertexNDOF_ = fVertexNDOF; }
    void SetVz(const std::vector<double>& fVz) { Vz_ = fVz; }
    void SetVx(const std::vector<double>& fVx) { Vx_ = fVx; }
    void SetVy(const std::vector<double>& fVy) { Vy_ = fVy; }
    void SetTracksPt(const std::vector<std::vector<double> >& fTracksPt) { TracksPt_ = fTracksPt; }
    void SetZDCdigifC(const std::vector<std::vector<double> >& fZDCdigifC) { ZDCdigifC_ = fZDCdigifC; }

    void SetPrimaryGapMaxGen(double fPrimaryGapMaxGen)    { PrimaryGapMaxGen_     = fPrimaryGapMaxGen;}
    void SetSecondGapMaxGen(double fSecondGapMaxGen)    { SecondGapMaxGen_    = fSecondGapMaxGen;}
    void SetTracksPtGen(const std::vector<double>& fTracksPtGen)    { TracksPtGen_     = fTracksPtGen;}
    void SetEtaOfTracksPtGen(const std::vector<double>& fEtaOfTracksPtGen)    { EtaOfTracksPtGen_     = fEtaOfTracksPtGen;}
    void SetNTracksGen(int fNTracksGen)    { NTracksGen_     = fNTracksGen;}
    void SetMx2PlusGen(double fMx2PlusGen)    { Mx2PlusGen_     = fMx2PlusGen;}
    void SetMx2MinusGen(double fMx2MinusGen)    { Mx2MinusGen_     = fMx2MinusGen;}
    void SetMx2Gen(double fMx2Gen)    { Mx2Gen_     = fMx2Gen;}
    void SetMx2ZGen(double fMx2ZGen)    { Mx2ZGen_     = fMx2ZGen;}
    void SetNMx2PlusGen(int fNMx2PlusGen)    { NMx2PlusGen_     = fNMx2PlusGen;}
    void SetNMx2MinusGen(int fNMx2MinusGen)    { NMx2MinusGen_     = fNMx2MinusGen;}
    void SetEtaGaplimPlusGen(double fEtaGaplimPlusGen)    { EtaGaplimPlusGen_     = fEtaGaplimPlusGen;}
    void SetEtaGaplimMinusGen(double fEtaGaplimMinusGen)    { EtaGaplimMinusGen_     = fEtaGaplimMinusGen;}
    void SetNParticlesGen(int fNParticlesGen)    { NParticlesGen_     = fNParticlesGen;}
    void SetsumECastorMinusGen(double fsumECastorMinusGen)    { sumECastorMinusGen_     = fsumECastorMinusGen;}
    void SetsumECastorPlusGen(double fsumECastorPlusGen)    { sumECastorPlusGen_     = fsumECastorPlusGen;}
    void SetsumEZDCMinusGen(double fsumEZDCMinusGen)    { sumEZDCMinusGen_     = fsumEZDCMinusGen;}
    void SetsumEZDCPlusGen(double fsumEZDCPlusGen)    { sumEZDCPlusGen_     = fsumEZDCPlusGen;}
    void SetEtaOutcomingProtonGen(double fEtaOutcomingProtonGen)    { EtaOutcomingProtonGen_     = fEtaOutcomingProtonGen;}
    void SetxLGen(double fxLGen)   { xLGen_     = fxLGen;}
    void SetxLMostEnergeticGen(double fxLMostEnergeticGen)    { xLMostEnergeticGen_     = fxLMostEnergeticGen;}
    void SetxiZMinusGen(double fxiZMinusGen)    { xiZMinusGen_     = fxiZMinusGen;}
    void SetxiZPlusGen(double fxiZPlusGen)    { xiZPlusGen_     = fxiZPlusGen;}
    void SetEtaZGen(double fEtaZGen)    { EtaZGen_     = fEtaZGen;}
    void SetEnergyZGen(double fEnergyZGen)    { EnergyZGen_     = fEnergyZGen;}
    void SetpDissMassGen(double fpDissMassGen)    { pDissMassGen_     = fpDissMassGen;}
    void SetxLpDissMass(double fxLpDissMass)    { xLpDissMass_     = fxLpDissMass;}

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

    void SetEPZCaloPlus(double fEPZCaloPlus)    { EPZCaloPlus_ = fEPZCaloPlus;}
    void SetEPZCaloMinus(double fEPZCaloMinus)    { EPZCaloMinus_ = fEPZCaloMinus;}
    void SetXiCaloPlus(double fXiCaloPlus)    { XiCaloPlus_ = fXiCaloPlus;}
    void SetXiCaloMinus(double fXiCaloMinus)    { XiCaloMinus_ = fXiCaloMinus;}

    void SetEtaCaloMax(double fEtaCaloMax)    { EtaCaloMax_ = fEtaCaloMax;}
    void SetEtaCaloMin(double fEtaCaloMin)    { EtaCaloMin_ = fEtaCaloMin;}

    void SetMultiplicityHFPlus(int fMultiplicityHFPlus)    { MultiplicityHFPlus_ = fMultiplicityHFPlus;}
    void SetMultiplicityHEPlus(int fMultiplicityHEPlus)    { MultiplicityHEPlus_ = fMultiplicityHEPlus;}
    void SetMultiplicityEEPlus(int fMultiplicityEEPlus)    { MultiplicityEEPlus_ = fMultiplicityEEPlus;}
    void SetMultiplicityHFMinus(int fMultiplicityHFMinus)    { MultiplicityHFMinus_ = fMultiplicityHFMinus;}
    void SetMultiplicityHEMinus(int fMultiplicityHEMinus)    { MultiplicityHEMinus_ = fMultiplicityHEMinus;}
    void SetMultiplicityEEMinus(int fMultiplicityEEMinus)    { MultiplicityEEMinus_ = fMultiplicityEEMinus;}

    void SetVertex(int fVertex)    { Vertex_ = fVertex;}
    void SetXi_PF_minus(double fXi_PF_minus)    { Xi_PF_minus_ = fXi_PF_minus;}
    void SetXi_PF_plus(double fXi_PF_plus)    { Xi_PF_plus_ = fXi_PF_plus;}
    void SetXiMass(double fXi_mass)    { Xi_mass_ = fXi_mass;}
    void SetXiMassNoW(double fXi_massnow)    { Xi_massnow_ = fXi_massnow;}

    void SetEpz_PF_minus(double fEpz_PF_minus)    { Epz_PF_minus_= fEpz_PF_minus;}
    void SetEpz_PF_plus(double fEpz_PF_plus)    { Epz_PF_plus_ = fEpz_PF_plus;}
    void SetMultiplicityPF(int fMultiplicityPF)    { MultiplicityPF_ = fMultiplicityPF;}
    void SetSumEtaTimesEnergyPF(double fSumEtaTimesEnergyPF)    { SumEtaTimesEnergyPF_ = fSumEtaTimesEnergyPF;}
    void SetSumpxModulePF(double fSumpxModulePF)    { SumpxModulePF_ = fSumpxModulePF;}
    void SetSumpyModulePF(double fSumpyModulePF)    { SumpyModulePF_ = fSumpyModulePF;}
    void SetSumpzModulePF(double fSumpzModulePF)    { SumpzModulePF_ = fSumpzModulePF;}
    void SetSumpxPF(double fSumpxPF)    { SumpxPF_ = fSumpxPF;}
    void SetSumpyPF(double fSumpyPF)    { SumpyPF_ = fSumpyPF;}
    void SetSumpzPF(double fSumpzPF)    { SumpzPF_ = fSumpzPF;}
    void SetSumEnergyPF(double fSumEnergyPF)    { SumEnergyPF_ = fSumEnergyPF;}
    void SetMuEnergyPF(double fMuEnergyPF)    { MuEnergyPF_ = fMuEnergyPF;}
    void SetElectronEnergyPF(double fElectronEnergyPF)    { ElectronEnergyPF_ = fElectronEnergyPF;}
    void SetMaxGapPF(double fMaxGapPF)    { MaxGapPF_ = fMaxGapPF;}
    void SetSecondMaxGapPF(double fSecondMaxGapPF)    { SecondMaxGapPF_ = fSecondMaxGapPF;}
    void SetLimPlusGapPF(double fLimPlusGapPF)    { LimPlusGapPF_ = fLimPlusGapPF;}
    void SetLimMinusGapPF(double fLimMinusGapPF)    { LimMinusGapPF_ = fLimMinusGapPF;}
    void SetPTMaxGapMaxPF(double fPTMaxGapMaxPF)    { PTMaxGapMaxPF_ = fPTMaxGapMaxPF;}
    void SetPTMinGapMaxPF(double fPTMinGapMaxPF)    { PTMinGapMaxPF_ = fPTMinGapMaxPF;}
    void SetMultiplicityGapPlusPF(int fMultiplicityGapPlusPF)    { MultiplicityGapPlusPF_ = fMultiplicityGapPlusPF;}
    void SetMultiplicityGapMinusPF(int fMultiplicityGapMinusPF)    { MultiplicityGapMinusPF_ = fMultiplicityGapMinusPF;}

    void SetPatNMuon(int fpatNMuon)    {patNMuon_ = fpatNMuon;}
    void SetPatMuon1Pt(double fpatMuon1Pt)    {patMuon1Pt_ = fpatMuon1Pt;}
    void SetPatMuon1Charge(int fpatMuon1Charge)    {patMuon1Charge_ = fpatMuon1Charge;}
    void SetPatMuon1Phi(double fpatMuon1Phi)    {patMuon1Phi_ = fpatMuon1Phi;}
    void SetPatMuon1Eta(double fpatMuon1Eta)    {patMuon1Eta_ = fpatMuon1Eta;}
    void SetPatMuon1Et(double fpatMuon1Et)    {patMuon1Et_ = fpatMuon1Et;}

    void SetPatMuon1SumPtR03(double fpatMuon1SumPtR03)    {patMuon1SumPtR03_ = fpatMuon1SumPtR03;}
    void SetPatMuon1EmEtR03(double fpatMuon1EmEtR03)    {patMuon1EmEtR03_ = fpatMuon1EmEtR03;}
    void SetPatMuon1HadEtR03(double fpatMuon1HadEtR03)    {patMuon1HadEtR03_ = fpatMuon1HadEtR03;}    
    void SetPatMuon1SumPtR05(double fpatMuon1SumPtR05)    {patMuon1SumPtR05_ = fpatMuon1SumPtR05;}
    void SetPatMuon1EmEtR05(double fpatMuon1EmEtR05)    {patMuon1EmEtR05_ = fpatMuon1EmEtR05;}
    void SetPatMuon1HadEtR05(double fpatMuon1HadEtR05)    {patMuon1HadEtR05_ = fpatMuon1HadEtR05;}    

    void SetPatMuon1relIsoDr03(double fpatMuon1relIsoDr03)    {patMuon1relIsoDr03_ = fpatMuon1relIsoDr03;}
    void SetPatMuon1relIsoDr05(double fpatMuon1relIsoDr05)    {patMuon1relIsoDr05_ = fpatMuon1relIsoDr05;}
    void SetPatMuon1relIso(double fpatMuon1relIso)    {patMuon1relIso_ = fpatMuon1relIso;}

    void SetPatNElectron(int fpatNElectron)    {patNElectron_ = fpatNElectron;}
    void SetPatElectron1Pt(double fpatElectron1Pt)    {patElectron1Pt_ = fpatElectron1Pt;}
    void SetPatElectron1Charge(int fpatElectron1Charge)    {patElectron1Charge_ = fpatElectron1Charge;}
    void SetPatElectron1Phi(double fpatElectron1Phi)    {patElectron1Phi_ = fpatElectron1Phi;}
    void SetPatElectron1Eta(double fpatElectron1Eta)    {patElectron1Eta_ = fpatElectron1Eta;}
    void SetPatElectron1Et(double fpatElectron1Et)    {patElectron1Et_ = fpatElectron1Et;}

    void SetPatElectron1TkDr03(double fpatElectron1TkDr03)    {patElectron1TkDr03_ = fpatElectron1TkDr03;}
    void SetPatElectron1EcalDr03(double fpatElectron1EcalDr03)    {patElectron1EcalDr03_ = fpatElectron1EcalDr03;}
    void SetPatElectron1HcalDr03(double fpatElectron1HcalDr03)    {patElectron1HcalDr03_ = fpatElectron1HcalDr03;}

    void SetPatElectron1TkDr04(double fpatElectron1TkDr04)    {patElectron1TkDr04_ = fpatElectron1TkDr04;}
    void SetPatElectron1EcalDr04(double fpatElectron1EcalDr04)    {patElectron1EcalDr04_ = fpatElectron1EcalDr04;}
    void SetPatElectron1HcalDr04(double fpatElectron1HcalDr04)    {patElectron1HcalDr04_ = fpatElectron1HcalDr04;}

    void SetPatElectron1relIsoDr03(double fpatElectron1relIsoDr03)    {patElectron1relIsoDr03_ = fpatElectron1relIsoDr03;}
    void SetPatElectron1relIsoDr04(double fpatElectron1relIsoDr04)    {patElectron1relIsoDr04_ = fpatElectron1relIsoDr04;}

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

    void SetMETPt(double fmetPt)    {fmetPt_ = fmetPt;}
    void SetMETPhi(double fmetPhi)    {fmetPhi_ = fmetPhi;}
    void SetMETEt(double fmetEt)    {fmetEt_ = fmetEt;}
    void SetMETSumEt(double fmetSumEt)    {fmetSumEt_ = fmetSumEt;}
    void SetMETpx(double fmetpx)    {fmetpx_ = fmetpx;}
    void SetMETpy(double fmetpy)    {fmetpy_ = fmetpy;}

    void SetPatMETPt(double fpatmetPt)    {fpatmetPt_ = fpatmetPt;}
    void SetPatMETPhi(double fpatmetPhi)    {fpatmetPhi_ = fpatmetPhi;}
    void SetPatMETEt(double fpatmetEt)    {fpatmetEt_ = fpatmetEt;}
    void SetPatMETSumEt(double fpatmetSumEt)    {fpatmetSumEt_ = fpatmetSumEt;}
    void SetPatMETpx(double fpatmetpx)    {fpatmetpx_ = fpatmetpx;}
    void SetPatMETpy(double fpatmetpy)    {fpatmetpy_ = fpatmetpy;}

    int GetHLTPath(int idx)                    const { return hltTrigResults_[idx]; }

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

    double GetVertexMultiplicity(int i) const { return VertexMultiplicity_[i]; }
    double GetVertexChiNorm(int i) const { return VertexChiNorm_[i]; }
    double GetVertexNDOF(int i) const { return VertexNDOF_[i]; }
    double GetVz(int i) const { return Vz_[i]; }
    double GetVx(int i) const { return Vx_[i]; }
    double GetVy(int i) const { return Vy_[i]; }
    double GetTracksPt(int i,int j) const { return TracksPt_[i][j]; }
    double GetZDCdigifC(int i,int j) const { return ZDCdigifC_[i][j]; }
    double GetPrimaryGapMaxGen()    const {return PrimaryGapMaxGen_;}
    double GetSecondGapMaxGen()    const {return SecondGapMaxGen_;}
    double GetTracksPtGen(int i)    const {return TracksPtGen_[i];}
    double GetEtaOfTracksPtGen(int i)    const {return EtaOfTracksPtGen_[i];}
    int GetNTracksGen()    const {return NTracksGen_;}
    double GetMx2PlusGen()    const {return Mx2PlusGen_;}
    double GetMx2MinusGen()    const {return Mx2MinusGen_;}
    double GetMx2Gen()    const {return Mx2Gen_;}
    double GetMx2ZGen()    const {return Mx2ZGen_;}
    int GetNMx2PlusGen()    const {return NMx2PlusGen_;}
    int GetNMx2MinusGen()    const {return NMx2MinusGen_;}
    double GetEtaGaplimPlusGen()    const {return EtaGaplimPlusGen_;}
    double GetEtaGaplimMinusGen()    const {return EtaGaplimMinusGen_;}
    int GetNParticlesGen()    const {return NParticlesGen_;}
    double GetsumECastorMinusGen()    const {return sumECastorMinusGen_;}
    double GetsumECastorPlusGen()    const {return sumECastorPlusGen_;}
    double GetsumEZDCMinusGen()    const {return sumEZDCMinusGen_;}
    double GetsumEZDCPlusGen()    const {return sumEZDCPlusGen_;}
    double GetEtaOutcomingProtonGen()    const {return EtaOutcomingProtonGen_;}
    double GetxLGen()   const {return xLGen_;}
    double GetxLMostEnergeticGen()    const {return xLMostEnergeticGen_;}
    double GetxiZMinusGen()    const {return xiZMinusGen_;}
    double GetxiZPlusGen()    const {return xiZPlusGen_;}
    double GetEtaZGen()    const {return EtaZGen_;}
    double GetEnergyZGen()    const {return EnergyZGen_;}
    double GetpDissMassGen()    const {return pDissMassGen_;}
    double GetxLpDissMass()    const {return xLpDissMass_;}

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

    double GetEPZCaloPlus()    const {return EPZCaloPlus_;}
    double GetEPZCaloMinus()    const {return EPZCaloMinus_;}
    double GetXiCaloPlus()    const {return XiCaloPlus_;}
    double GetXiCaloMinus()    const {return XiCaloMinus_;}

    double GetEtaCaloMax()    const {return EtaCaloMax_;}
    double GetEtaCaloMin()    const {return EtaCaloMin_;}

    int GetMultiplicityHFPlus()    const {return MultiplicityHFPlus_;}
    int GetMultiplicityHEPlus()    const {return MultiplicityHEPlus_;}
    int GetMultiplicityEEPlus()    const {return MultiplicityEEPlus_;}
    int GetMultiplicityHFMinus()    const {return MultiplicityHFMinus_;}
    int GetMultiplicityHEMinus()    const {return MultiplicityHEMinus_;}
    int GetMultiplicityEEMinus()    const {return MultiplicityEEMinus_;}

    int GetVertex()    const {return Vertex_;}
    double GetXi_PF_minus()    const {return Xi_PF_minus_;}
    double GetXi_PF_plus()    const {return Xi_PF_plus_;}
    double GetXiMass()    const {return Xi_mass_;}
    double GetXiMassNoW()    const {return Xi_massnow_;}

    double GetEpz_PF_minus()    const {return Epz_PF_minus_;}
    double GetEpz_PF_plus()    const {return Epz_PF_plus_;}
    int GetMultiplicityPF()    const {return MultiplicityPF_;}
    double GetSumEtaTimesEnergyPF()    const {return SumEtaTimesEnergyPF_;}
    double GetSumpxModulePF()    const {return SumpxModulePF_;}
    double GetSumpyModulePF()    const {return SumpyModulePF_;}
    double GetSumpzModulePF()    const {return SumpzModulePF_;}
    double GetSumpxPF()    const {return SumpxPF_;}
    double GetSumpyPF()    const {return SumpyPF_;}
    double GetSumpzPF()    const {return SumpzPF_;}
    double GetSumEnergyPF()    const {return SumEnergyPF_;}
    double GetMuEnergyPF()    const {return MuEnergyPF_;}
    double GetElectronEnergyPF()    const {return ElectronEnergyPF_;}
    double GetMaxGapPF()    const {return MaxGapPF_;}
    double GetSecondMaxGapPF()    const {return SecondMaxGapPF_;}
    double GetLimPlusGapPF()    const {return LimPlusGapPF_;}
    double GetLimMinusGapPF()    const {return LimMinusGapPF_;}
    double GetPTMaxGapMaxPF()    const {return PTMaxGapMaxPF_;}
    double GetPTMinGapMaxPF()    const {return PTMinGapMaxPF_;}
    int GetMultiplicityGapPlusPF()    const {return MultiplicityGapPlusPF_;}
    int GetMultiplicityGapMinusPF()    const {return MultiplicityGapMinusPF_;}

    int GetPatNMuon() const {return patNMuon_;}
    double GetPatMuon1Pt() const {return patMuon1Pt_;}
    int GetPatMuon1Charge() const {return patMuon1Charge_;}
    double GetPatMuon1Phi() const {return patMuon1Phi_;}
    double GetPatMuon1Eta() const {return patMuon1Eta_;}
    double GetPatMuon1Et() const {return patMuon1Et_;}

    double GetPatMuon1SumPtR03() const {return patMuon1SumPtR03_;}
    double GetPatMuon1EmEtR03() const {return patMuon1EmEtR03_;}
    double GetPatMuon1HadEtR03() const {return patMuon1HadEtR03_;}    
    double GetPatMuon1SumPtR05() const {return patMuon1SumPtR05_;}
    double GetPatMuon1EmEtR05() const {return patMuon1EmEtR05_;}
    double GetPatMuon1HadEtR05() const {return patMuon1HadEtR05_;}    

    double GetPatMuon1relIsoDr03() const {return patMuon1relIsoDr03_;}
    double GetPatMuon1relIsoDr05() const {return patMuon1relIsoDr05_;}

    double GetPatMuon1relIso() const {return patMuon1relIso_;}

    int GetPatNElectron() const {return patNElectron_;}
    double GetPatElectron1Pt() const {return patElectron1Pt_;}
    int GetPatElectron1Charge() const {return patElectron1Charge_;}
    double GetPatElectron1Phi() const {return patElectron1Phi_;}
    double GetPatElectron1Eta() const {return patElectron1Eta_;}
    double GetPatElectron1Et() const {return patElectron1Et_;}

    double GetPatElectron1TkDr03() const  {return patElectron1TkDr03_;}    
    double GetPatElectron1EcalDr03() const  {return patElectron1EcalDr03_;}
    double GetPatElectron1HcalDr03() const  {return patElectron1HcalDr03_;}

    double GetPatElectron1TkDr04() const  {return patElectron1TkDr04_;}
    double GetPatElectron1EcalDr04() const  {return patElectron1EcalDr04_;}
    double GetPatElectron1HcalDr04() const  {return patElectron1HcalDr04_;}

    double GetPatElectron1relIsoDr03() const {return patElectron1relIsoDr03_;}
    double GetPatElectron1relIsoDr04() const {return patElectron1relIsoDr04_;}

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

    double GetPatMETPt() const {return fpatmetPt_;}
    double GetPatMETPhi() const {return fpatmetPhi_;}
    double GetPatMETEt() const {return fpatmetEt_;}
    double GetPatMETSumEt() const {return fpatmetSumEt_;}
    double GetPatMETpx() const {return fpatmetpx_;}
    double GetPatMETpy() const {return fpatmetpy_;}

  private:
    friend class diffractiveWAnalysis::DiffractiveWAnalysis;

    void reset();

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

    double PrimaryGapMaxGen_;
    double SecondGapMaxGen_;
    std::vector<double> TracksPtGen_;
    std::vector<double> EtaOfTracksPtGen_;
    int NTracksGen_;
    double Mx2PlusGen_;
    double Mx2MinusGen_;
    double Mx2Gen_;
    double Mx2ZGen_;
    int NMx2PlusGen_;
    int NMx2MinusGen_;
    double EtaGaplimPlusGen_;
    double EtaGaplimMinusGen_;
    int NParticlesGen_;
    double sumECastorMinusGen_;
    double sumECastorPlusGen_;
    double sumEZDCMinusGen_;
    double sumEZDCPlusGen_;
    double EtaOutcomingProtonGen_;
    double xLGen_;
    double xLMostEnergeticGen_;
    double xiZMinusGen_;
    double xiZPlusGen_;
    double EtaZGen_;
    double EnergyZGen_;
    double pDissMassGen_;
    double xLpDissMass_;

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
    double XiCaloPlus_;
    double XiCaloMinus_;
    double EtaCaloMax_;
    double EtaCaloMin_;
    int MultiplicityHFPlus_;
    int MultiplicityHEPlus_;
    int MultiplicityEEPlus_;
    int MultiplicityHFMinus_;
    int MultiplicityHEMinus_;
    int MultiplicityEEMinus_;

    int Vertex_;
    double Xi_PF_minus_;
    double Xi_PF_plus_;
    double Xi_mass_;
    double Xi_massnow_;
    double Epz_PF_minus_;
    double Epz_PF_plus_;
    int MultiplicityPF_;
    double SumEtaTimesEnergyPF_;
    double SumpxModulePF_;
    double SumpyModulePF_;
    double SumpzModulePF_;
    double SumpxPF_;
    double SumpyPF_;
    double SumpzPF_;
    double SumEnergyPF_;
    double MuEnergyPF_;
    double ElectronEnergyPF_;
    double MaxGapPF_;
    double SecondMaxGapPF_;
    double LimPlusGapPF_;
    double LimMinusGapPF_;
    double PTMaxGapMaxPF_;
    double PTMinGapMaxPF_;
    int MultiplicityGapPlusPF_;
    int MultiplicityGapMinusPF_;

    int patNMuon_;

    double patMuon1Pt_;
    int patMuon1Charge_;
    double patMuon1Phi_;
    double patMuon1Eta_;
    double patMuon1Et_;

    double patMuon1SumPtR03_;
    double patMuon1EmEtR03_;
    double patMuon1HadEtR03_;   
    double patMuon1SumPtR05_;
    double patMuon1EmEtR05_;
    double patMuon1HadEtR05_;   

    double patMuon1relIsoDr03_;
    double patMuon1relIsoDr05_;

    double patMuon1relIso_;

    int patNElectron_; 

    double patElectron1Pt_;
    int patElectron1Charge_;
    double patElectron1Phi_;
    double patElectron1Eta_;
    double patElectron1Et_;

    double patElectron1TkDr03_;    
    double patElectron1EcalDr03_;
    double patElectron1HcalDr03_;

    double patElectron1TkDr04_;
    double patElectron1EcalDr04_;
    double patElectron1HcalDr04_;

    double patElectron1relIsoDr03_;
    double patElectron1relIsoDr04_;

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

    double fpatmetPt_;
    double fpatmetPhi_;
    double fpatmetEt_;
    double fpatmetSumEt_;
    double fpatmetpx_;
    double fpatmetpy_;

};

#endif    
