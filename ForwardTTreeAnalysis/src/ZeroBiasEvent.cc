#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ZeroBiasEvent.h"
#include <cstdio>

const char* ZeroBiasEvent::name = "ZeroBiasEvent";

ZeroBiasEvent::ZeroBiasEvent() {}

ZeroBiasEvent::~ZeroBiasEvent() {}

void ZeroBiasEvent::reset(){

  size_t len_hltTrigResults = sizeof(hltTrigResults_)/sizeof(int);
  for(size_t k = 0; k < len_hltTrigResults; ++k) hltTrigResults_[k] = 0;

  vtx_ = -999;
  ntrk_ = -999;
  EHFPlus_.clear();
  EHFMinus_.clear();
  EHEPlus_.clear();
  EHEMinus_.clear();
  EHBPlus_.clear();
  EHBMinus_.clear();
  EEEPlus_.clear();
  EEEMinus_.clear();
  EEBPlus_.clear();
  EEBMinus_.clear();
  PFCharged_.clear();
  PFBarrelPlus_.clear();
  PFBarrelMinus_.clear();
  PFEndCapPlus_.clear();
  PFEndCapMinus_.clear();
  PFHFPlus_.clear();
  PFHFMinus_.clear();
  CastorTowerEnergy_.clear();
  CastorModule1Energy_.clear();
  CastorModule2Energy_.clear();
  CastorModule3Energy_.clear();
  CastorModule4Energy_.clear();
  CastorModule5Energy_.clear();
  CastorBadChannels_.clear();
  CastorNumberBadChannels_ = -999;
  LooseNoiseFilter_ = -1;
  TightNoiseFilter_ = -1;
  BeamHaloLooseId_ = -1;
  BeamHaloTightId_ = -1;
  NHFPlus_=-999;
  NHFMinus_=-999;
  NHEPlus_=-999;
  NHEMinus_=-999;
  NHBPlus_=-999;
  NHBMinus_=-999;
  NEEPlus_=-999;
  NEEMinus_=-999;
  NEBPlus_=-999;
  NEBMinus_=-999;
  PFNCharged_=-999;
  PFNBarrelPlus_=-999;
  PFNBarrelMinus_=-999;
  PFNEndCapPlus_=-999;
  PFNEndCapMinus_=-999;
  PFNHFPlus_=-999;
  PFNHFMinus_=-999;

}
