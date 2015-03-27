#ifndef ZeroBiasEvent_h
#define ZeroBiasEvent_h

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/DetSetVector.h"

namespace zerobiasAnalysis {
  class ZeroBiasAnalysis;
}

class ZeroBiasEvent {
  public:
    typedef zerobiasAnalysis::ZeroBiasAnalysis analysis_type;
    static const char* name;

    typedef reco::Particle::LorentzVector LorentzVector;

    ZeroBiasEvent();
    ~ZeroBiasEvent();

    void SetNVertex(int fvtx)    {vtx_ = fvtx;}
    void SetMultiplicityTracks(int fntrk) {ntrk_ = fntrk;}
    void SetHLTPath(int idx, int fHLTBit)         { hltTrigResults_[idx] = fHLTBit;}
    void SetEHFPlus(const std::vector<double>& fEHFPlus) { EHFPlus_ = fEHFPlus; }
    void SetEHFMinus(const std::vector<double>& fEHFMinus) { EHFMinus_ = fEHFMinus; }
    void SetEHEPlus(const std::vector<double>& fEHEPlus) { EHEPlus_ = fEHEPlus; }
    void SetEHEMinus(const std::vector<double>& fEHEMinus) { EHEMinus_ = fEHEMinus; }
    void SetEHBPlus(const std::vector<double>& fEHBPlus) { EHBPlus_ = fEHBPlus; }
    void SetEHBMinus(const std::vector<double>& fEHBMinus) { EHBMinus_ = fEHBMinus; }
    void SetEEEPlus(const std::vector<double>& fEEEPlus) { EEEPlus_ = fEEEPlus; }
    void SetEEEMinus(const std::vector<double>& fEEEMinus) { EEEMinus_ = fEEEMinus; }
    void SetEEBPlus(const std::vector<double>& fEEBPlus) { EEBPlus_ = fEEBPlus; }
    void SetEEBMinus(const std::vector<double>& fEEBMinus) { EEBMinus_ = fEEBMinus; }
    void SetPFCharged(const std::vector<double>& fPFCharged) { PFCharged_ = fPFCharged; }
    void SetPFBarrelPlus(const std::vector<double>& fPFBarrelPlus) { PFBarrelPlus_ = fPFBarrelPlus; }
    void SetPFBarrelMinus(const std::vector<double>& fPFBarrelMinus) { PFBarrelMinus_ = fPFBarrelMinus; }
    void SetPFEndCapPlus(const std::vector<double>& fPFEndCapPlus) { PFEndCapPlus_ = fPFEndCapPlus; }
    void SetPFEndCapMinus(const std::vector<double>& fPFEndCapMinus) { PFEndCapMinus_ = fPFEndCapMinus; }
    void SetPFHFPlus(const std::vector<double>& fPFHFPlus) { PFHFPlus_ = fPFHFPlus; }
    void SetPFHFMinus(const std::vector<double>& fPFHFMinus) { PFHFMinus_ = fPFHFMinus; }
    void SetNHFPlus(int fNHFPlus) { NHFPlus_ = fNHFPlus; }
    void SetNHFMinus(int fNHFMinus) { NHFMinus_ = fNHFMinus; }
    void SetNHEPlus(int fNHEPlus) { NHEPlus_ = fNHEPlus; }
    void SetNHEMinus(int fNHEMinus) { NHEMinus_ = fNHEMinus; }
    void SetNHBPlus(int fNHBPlus) { NHBPlus_ = fNHBPlus; }
    void SetNHBMinus(int fNHBMinus) { NHBMinus_ = fNHBMinus; }
    void SetNEEPlus(int fNEEPlus) { NEEPlus_ = fNEEPlus; }
    void SetNEEMinus(int fNEEMinus) { NEEMinus_ = fNEEMinus; }
    void SetNEBPlus(int fNEBPlus) { NEBPlus_ = fNEBPlus; }
    void SetNEBMinus(int fNEBMinus) { NEBMinus_ = fNEBMinus; }
    void SetPFNCharged(int fPFNCharged) { PFNCharged_ = fPFNCharged; }
    void SetPFNBarrelPlus(int fPFNBarrelPlus) { PFNBarrelPlus_ = fPFNBarrelPlus; }
    void SetPFNBarrelMinus(int fPFNBarrelMinus) { PFNBarrelMinus_ = fPFNBarrelMinus; }
    void SetPFNEndCapPlus(int fPFNEndCapPlus) { PFNEndCapPlus_ = fPFNEndCapPlus; }
    void SetPFNEndCapMinus(int fPFNEndCapMinus) { PFNEndCapMinus_ = fPFNEndCapMinus; }
    void SetPFNHFPlus(int fPFNHFPlus) { PFNHFPlus_ = fPFNHFPlus; }
    void SetPFNHFMinus(int fPFNHFMinus) { PFNHFMinus_ = fPFNHFMinus; }
    void SetCastorTowerEnergy(const std::vector<double>& fCastorTowerEnergy) { CastorTowerEnergy_ = fCastorTowerEnergy; }
    void SetCastorModule1Energy(const std::vector<double>& fCastorModule1Energy) { CastorModule1Energy_ = fCastorModule1Energy; }
    void SetCastorModule2Energy(const std::vector<double>& fCastorModule2Energy) { CastorModule2Energy_ = fCastorModule2Energy; }
    void SetCastorModule3Energy(const std::vector<double>& fCastorModule3Energy) { CastorModule3Energy_ = fCastorModule3Energy; }
    void SetCastorModule4Energy(const std::vector<double>& fCastorModule4Energy) { CastorModule4Energy_ = fCastorModule4Energy; }
    void SetCastorModule5Energy(const std::vector<double>& fCastorModule5Energy) { CastorModule5Energy_ = fCastorModule5Energy; }
    void SetCastorBadChannels(const std::vector<int>& fCastorBadChannels) { CastorBadChannels_ = fCastorBadChannels; }
    void SetCastorNumberBadChannels(int fCastorNumberBadChannels) { CastorNumberBadChannels_ = fCastorNumberBadChannels;}

    int GetNVertex() const {return vtx_;}
    int GetMultiplicityTracks() const {return ntrk_;}
    int GetHLTPath(int idx) const { return hltTrigResults_[idx]; }
    double GetEHFPlus(int i) const { return EHFPlus_[i]; }
    double GetEHFMinus(int i) const { return EHFMinus_[i]; }
    double GetEHEPlus(int i) const { return EHEPlus_[i]; }
    double GetEHEMinus(int i) const { return EHEMinus_[i]; }
    double GetEHBPlus(int i) const { return EHBPlus_[i]; }
    double GetEHBMinus(int i) const { return EHBMinus_[i]; }
    double GetEEEPlus(int i) const { return EEEPlus_[i]; }
    double GetEEEMinus(int i) const { return EEEMinus_[i]; }
    double GetEEBPlus(int i) const { return EEBPlus_[i]; }
    double GetEEBMinus(int i) const { return EEBMinus_[i]; }
    double GetPFCharged(int i) const { return  PFCharged_[i]; }
    double GetPFBarrelPlus(int i) const { return  PFBarrelPlus_[i]; }
    double GetPFBarrelMinus(int i) const { return  PFBarrelMinus_[i]; }
    double GetPFEndCapPlus(int i) const { return  PFEndCapPlus_[i]; }
    double GetPFEndCapMinus(int i) const { return  PFEndCapMinus_[i]; }
    double GetPFHFPlus(int i) const { return  PFHFPlus_[i]; }
    double GetPFHFMinus(int i) const { return  PFHFMinus_[i]; }
    int GetNHFPlus() const { return NHFPlus_;}
    int GetNHFMinus() const { return NHFMinus_;}
    int GetNHEPlus() const { return NHEPlus_;}
    int GetNHEMinus() const { return NHEMinus_;}
    int GetNHBPlus() const { return NHBPlus_;}
    int GetNHBMinus() const { return NHBMinus_;}
    int GetNEEPlus() const { return NEEPlus_;}
    int GetNEEMinus() const { return NEEMinus_;}
    int GetNEBPlus() const { return NEBPlus_;}
    int GetNEBMinus() const { return NEBMinus_;}
    int GetPFNCharged() const { return PFNCharged_;}
    int GetPFNBarrelPlus() const { return PFNBarrelPlus_;}
    int GetPFNBarrelMinus() const { return PFNBarrelMinus_;}
    int GetPFNEndCapPlus() const { return PFNEndCapPlus_;}
    int GetPFNEndCapMinus() const { return PFNEndCapMinus_;}
    int GetPFNHFPlus() const { return PFNHFPlus_;}
    int GetPFNHFMinus() const { return PFNHFMinus_;}
    double GetCastorTowerEnergy(int i) const { return CastorTowerEnergy_[i]; }
    double GetCastorModule1Energy(int i) const { return CastorModule1Energy_[i]; }
    double GetCastorModule2Energy(int i) const { return CastorModule2Energy_[i]; }
    double GetCastorModule3Energy(int i) const { return CastorModule3Energy_[i]; }
    double GetCastorModule4Energy(int i) const { return CastorModule4Energy_[i]; }
    double GetCastorModule5Energy(int i) const { return CastorModule5Energy_[i]; }
    int GetCastorBadChannels(int i) const { return CastorBadChannels_[i]; }
    int GetCastorNumberBadChannels() const { return CastorNumberBadChannels_;}
    int GetLooseNoiseFilter() const {return LooseNoiseFilter_;}
    int GetTightNoiseFilter() const {return TightNoiseFilter_;}
    int GetBeamHaloLooseId() const {return BeamHaloLooseId_;}
    int GetBeamHaloTightId() const {return BeamHaloTightId_;}

  private:
    friend class zerobiasAnalysis::ZeroBiasAnalysis;

    void reset();

    int vtx_;
    int ntrk_;
    int hltTrigResults_[20];
    std::vector<double> EHFPlus_;
    std::vector<double> EHFMinus_;
    std::vector<double> EHEPlus_;
    std::vector<double> EHEMinus_;
    std::vector<double> EHBPlus_;
    std::vector<double> EHBMinus_;
    std::vector<double> EEEPlus_;
    std::vector<double> EEEMinus_;
    std::vector<double> EEBPlus_;
    std::vector<double> EEBMinus_;
    std::vector<double> PFCharged_;
    std::vector<double> PFBarrelPlus_;
    std::vector<double> PFBarrelMinus_;
    std::vector<double> PFEndCapPlus_;
    std::vector<double> PFEndCapMinus_;
    std::vector<double> PFHFPlus_;
    std::vector<double> PFHFMinus_;
    int NHFPlus_;
    int NHFMinus_;
    int NHEPlus_;
    int NHEMinus_;
    int NHBPlus_;
    int NHBMinus_;
    int NEEPlus_;
    int NEEMinus_;
    int NEBPlus_;
    int NEBMinus_;
    int PFNCharged_;
    int PFNBarrelPlus_;
    int PFNBarrelMinus_;
    int PFNEndCapPlus_;
    int PFNEndCapMinus_;
    int PFNHFPlus_;
    int PFNHFMinus_;
    std::vector<double> CastorTowerEnergy_;
    std::vector<double> CastorModule1Energy_;
    std::vector<double> CastorModule2Energy_;
    std::vector<double> CastorModule3Energy_;
    std::vector<double> CastorModule4Energy_;
    std::vector<double> CastorModule5Energy_;
    std::vector<int> CastorBadChannels_;
    int CastorNumberBadChannels_;
    int LooseNoiseFilter_;
    int TightNoiseFilter_;
    int BeamHaloLooseId_;
    int BeamHaloTightId_;

};

#endif    
