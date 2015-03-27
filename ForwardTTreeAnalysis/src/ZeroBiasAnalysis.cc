#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ZeroBiasAnalysis.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ZeroBiasEvent.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"

#include "TDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "ForwardAnalysis/Utilities/interface/CastorEnergy.h"

#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include <stdio.h>
#include <math.h> 
#include <cmath>

using zerobiasAnalysis::ZeroBiasAnalysis;

const char* ZeroBiasAnalysis::name = "ZeroBiasAnalysis";

ZeroBiasAnalysis::ZeroBiasAnalysis(const edm::ParameterSet& pset):
  triggerResultsTag_(pset.getParameter<edm::InputTag>("TriggerResultsTag")),
  hltPathNames_(pset.getParameter<std::vector<std::string> >("hltPaths")),
  PVtxCollectionTag_(pset.getParameter<edm::InputTag>("PVtxCollectionTag")),
  trackTag_(pset.getParameter<edm::InputTag>("TrackTag")),
  caloTowerTag_(pset.getParameter<edm::InputTag>("CaloTowerTag")),
  energyThresholdHB_(pset.getParameter<double>("energyThresholdHB")),
  energyThresholdHE_(pset.getParameter<double>("energyThresholdHE")),
  energyThresholdHF_(pset.getParameter<double>("energyThresholdHF")),
  energyThresholdEB_(pset.getParameter<double>("energyThresholdEB")),
  energyThresholdEE_(pset.getParameter<double>("energyThresholdEE")),
  pfTag_(pset.getParameter<edm::InputTag>("pfTag")),
  pTPFThresholdCharged_(pset.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(pset.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(pset.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(pset.getParameter<double>("energyPFThresholdHF")),
  castorHitsTag_(pset.getParameter<edm::InputTag>("castorHitsTag")),
  RunA_(pset.getUntrackedParameter<Bool_t>("RunA", false)),
  RunB_(pset.getUntrackedParameter<Bool_t>("RunB", false)),
  fCGeVCastor_(pset.getParameter<double>("fCGeVCastor"))
{
}

void ZeroBiasAnalysis::setTFileService(){

  edm::Service<TFileService> fs;
  std::ostringstream oss;

  TFileDirectory triggerDir = fs->mkdir("TriggerInfo");
  hltTriggerNamesHisto_ = triggerDir.make<TH1F>("HLTTriggerNames","HLTTriggerNames",1,0,1);
  hltTriggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned k=0; k < hltPathNames_.size(); ++k){
    oss << "Using HLT reference trigger " << hltPathNames_[k] << std::endl;
    hltTriggerNamesHisto_->Fill(hltPathNames_[k].c_str(),1);
  }
  edm::LogVerbatim("Analysis") << oss.str();

  hltTriggerPassHisto_ = triggerDir.make<TH1F>("HLTTriggerPass","HLTTriggerPass",1,0,1);
  hltTriggerPassHisto_->SetBit(TH1::kCanRebin);

  TFileDirectory castorDir = fs->mkdir("CastorInfo");
  CastorChannelHisto_ = castorDir.make<TH1F>("CastorChannelWorking","Working Channel; Channel (id); # Times of working",240,0,240);
  for (int chan=1; chan <=224; chan++){
    char castor_channels[300];
    char castor_title[300];
    sprintf(castor_channels,"Castor_Channel_%d",chan);
    sprintf(castor_title,"Castor Channel %d Energy Distribution; Energy [GeV]; NEvents",chan);
    histo_castor_channels = castorDir.make<TH1F>(castor_channels,castor_title,1000,0,500);
    m_hVector_histo_castor_channels.push_back(histo_castor_channels);
  }

}

ZeroBiasAnalysis::~ZeroBiasAnalysis(){}

void ZeroBiasAnalysis::begin() {
  setTFileService();
}

void ZeroBiasAnalysis::begin(const edm::Run& run, const edm::EventSetup& setup) {}

void ZeroBiasAnalysis::end() {}

void ZeroBiasAnalysis::fill(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  eventData.reset();

  fillTriggerInfo(eventData,event,setup);
  fillVertexInfo(eventData,event,setup);
  fillTrackInfo(eventData,event,setup);
  fillDetectorVariables(eventData,event,setup);
  fillParticleFlow(eventData,event,setup);
  fillCastor(eventData,event,setup);
  fillCastorDebug(eventData,event,setup);
  fillNoiseInfo(eventData,event,setup);

}

// Fill Trigger
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ZeroBiasAnalysis::fillTriggerInfo(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(triggerResultsTag_, triggerResults);

  if( triggerResults.isValid() ){

    int nSize = triggerResults->size();
    const edm::TriggerNames& triggerNames = event.triggerNames(*triggerResults);

    size_t idxpath = 0;
    std::vector<std::string>::const_iterator hltpath = hltPathNames_.begin();
    std::vector<std::string>::const_iterator hltpaths_end = hltPathNames_.end();
    for(; hltpath != hltpaths_end; ++hltpath,++idxpath){
      std::string resolvedPathName;
      if( edm::is_glob( *hltpath ) ){
	std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(triggerNames.triggerNames(), *hltpath);
	if( matches.empty() ){
	  if (debug) edm::LogWarning("Configuration") << "Could not find trigger " << *hltpath << " in the path list.\n";
	}
	else if( matches.size() > 1)
	  throw cms::Exception("Configuration") << "HLT path type " << *hltpath << " not unique\n";
	else resolvedPathName = *(matches[0]);
      } else{
	resolvedPathName = *hltpath;
      }

      int idx_HLT = triggerNames.triggerIndex(resolvedPathName);

      if (idx_HLT >= 0 && idx_HLT < nSize){
	if (debug) std::cout << "Error accept_HLT?" << std::endl;
	int accept_HLT = ( triggerResults->wasrun(idx_HLT) && triggerResults->accept(idx_HLT) ) ? 1 : 0;
	if (debug) std::cout << "No... , accept_HLT: " << accept_HLT << std::endl;
	if (debug) std::cout << "Error eventData.SetHLTPath?" << std::endl;
	eventData.SetHLTPath(idxpath, accept_HLT);
	if (debug) std::cout << "No..." << std::endl;
	hltTriggerPassHisto_->Fill( (*hltpath).c_str(), 1 );
      }else{
	eventData.SetHLTPath(idxpath, -1);
	hltTriggerPassHisto_->Fill( (*hltpath).c_str(), -1 );
      }

    }

  }else{
    if (debug) std::cout << "\n No valid trigger result.\n" <<std::endl;
  }

}

// Fill Vertex
//////////////

void ZeroBiasAnalysis::fillVertexInfo(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Access vertex collection
  edm::Handle<edm::View<reco::Vertex> > vertexCollectionH;
  event.getByLabel(PVtxCollectionTag_,vertexCollectionH);
  const edm::View<reco::Vertex>& vtxColl = *(vertexCollectionH.product());

  // Find number of good vertices
  int nGoodVertices = 0;
  for(edm::View<reco::Vertex>::const_iterator vtx = vtxColl.begin(); vtx != vtxColl.end(); ++vtx){
    if(!vtx->isValid()) continue; // skip non-valid vertices
    if(vtx->isFake()) continue; // skip vertex from beam spot
    ++nGoodVertices;
  }

  eventData.SetNVertex(nGoodVertices);

}

// Fill Tracks
//////////////

void ZeroBiasAnalysis::fillTrackInfo(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Access collection
  edm::Handle<edm::View<reco::Track> > trackCollectionH;
  event.getByLabel(trackTag_,trackCollectionH);
  const edm::View<reco::Track>& trackColl = *(trackCollectionH.product());
  int nTracks = 0;
  edm::View<reco::Track>::const_iterator track = trackColl.begin();
  edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
  for(; track != tracks_end; ++track){
    ++nTracks;
  }
  eventData.SetMultiplicityTracks(nTracks);

}


//
// Fill Detector Variables
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ZeroBiasAnalysis::fillDetectorVariables(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  std::vector<double> VectorHFEnergyPlus;
  std::vector<double> VectorHFEnergyMinus;
  std::vector<double> VectorHEEnergyPlus;
  std::vector<double> VectorHEEnergyMinus;
  std::vector<double> VectorHBEnergyPlus;
  std::vector<double> VectorHBEnergyMinus;
  std::vector<double> VectorEEEnergyPlus;
  std::vector<double> VectorEEEnergyMinus;
  std::vector<double> VectorEBEnergyPlus;
  std::vector<double> VectorEBEnergyMinus;

  VectorHFEnergyPlus.clear();
  VectorHFEnergyMinus.clear();
  VectorHEEnergyPlus.clear();
  VectorHEEnergyMinus.clear();
  VectorHBEnergyPlus.clear();
  VectorHBEnergyMinus.clear();
  VectorEEEnergyPlus.clear();
  VectorEEEnergyMinus.clear();
  VectorEBEnergyPlus.clear();
  VectorEBEnergyMinus.clear();

  edm::Handle<CaloTowerCollection> towerCollectionH;
  event.getByLabel(caloTowerTag_,towerCollectionH);

  int towersize = towerCollectionH->size();
  int itTower;

  for(itTower=0; itTower < towersize; ++itTower){

    const CaloTower* calotower = &((*towerCollectionH)[itTower]);

    // Excluding HF Calorimeter Rings 12, 13, 29, 30, 40, 41
    if( ( (fabs(calotower->eta()) >= 2.866) && (fabs(calotower->eta()) < 3.152) ) || (fabs(calotower->eta()) >= 4.730) ) continue;

    bool hasHCAL = false;
    bool hasHF = false;
    bool hasHE = false;
    bool hasHB = false;
    bool hasHO = false;
    bool hasECAL = false;
    bool hasEE = false;
    bool hasEB = false;

    for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){
      DetId adetId = calotower->constituent(iconst);
      if(adetId.det()==DetId::Hcal){
	hasHCAL = true;
	HcalDetId hcalDetId(adetId);
	if(hcalDetId.subdet()==HcalForward) hasHF = true;
	else if(hcalDetId.subdet()==HcalEndcap) hasHE = true;
	else if(hcalDetId.subdet()==HcalBarrel) hasHB = true;
	else if(hcalDetId.subdet()==HcalOuter) hasHO = true;
      }
      else if(adetId.det()==DetId::Ecal){
	hasECAL = true;
	EcalSubdetector ecalSubDet = (EcalSubdetector)adetId.subdetId();
	if(ecalSubDet == EcalEndcap) hasEE = true;
	else if(ecalSubDet == EcalBarrel) hasEB = true;
      }
    }

    int zside = calotower->zside();
    double caloTowerEnergy = calotower->energy();
    double caloTowerEmEnergy = calotower->emEnergy();
    double caloTowerHadEnergy = calotower->hadEnergy();

    if( hasHF && !hasHE )
    {
      if( caloTowerEnergy > energyThresholdHF_ && fabs(calotower->eta())> 2.98 ) //// excluding HF ring1
      {
	if (debug) std::cout << "HF>> " << calotower->id() << " HAD energy " << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  VectorHFEnergyPlus.push_back(caloTowerEnergy);
	}
	else
	{
	  VectorHFEnergyMinus.push_back(caloTowerEnergy);
	}
      }
    }
    else if( hasHE && !hasHF && !hasHB )
    {
      if( caloTowerHadEnergy > energyThresholdHE_)
      {
	if (debug) std::cout << "HE>> " << calotower->id() << " HAD energy " << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;
	if(zside >= 0)
	{
	  VectorHEEnergyPlus.push_back(caloTowerEnergy);
	}
	else
	{
	  VectorHEEnergyMinus.push_back(caloTowerEnergy);
	}
      }
    }
    else if( hasHB && !hasHE )
    {
      if( caloTowerHadEnergy > energyThresholdHB_)
      {
	if (debug) std::cout << "HB>> " << calotower->id() << " HAD energy " << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;
	if(zside >= 0)
	{
	  VectorHBEnergyPlus.push_back(caloTowerEnergy);
	}
	else
	{
	  VectorHBEnergyMinus.push_back(caloTowerEnergy);
	}
      }
    }
    if( hasEE && !hasEB )
    {
      if( caloTowerEmEnergy >= energyThresholdEE_)
      {
	if (debug) std::cout << "EE>> " << calotower->id() << " HAD energy " << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;
	if(zside >= 0)
	{
	  VectorEEEnergyPlus.push_back(caloTowerEnergy);
	}
	else
	{
	  VectorEEEnergyMinus.push_back(caloTowerEnergy);
	}
      }
    }
    else if( hasEB && !hasEE )
    {
      if( caloTowerEmEnergy >= energyThresholdEB_)
      {
	if (debug) std::cout << "EB>> " << calotower->id() << " HAD energy " << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;
	if(zside >= 0)
	{
	  VectorEBEnergyPlus.push_back(caloTowerEnergy);
	}
	else
	{
	  VectorEBEnergyMinus.push_back(caloTowerEnergy);
	}
      }
    }
  } ////has to close calotower loop

  //Fill Variables
  eventData.SetEHFPlus(VectorHFEnergyPlus);
  eventData.SetEHFMinus(VectorHFEnergyMinus);
  eventData.SetEHEPlus(VectorHEEnergyPlus);
  eventData.SetEHEMinus(VectorHEEnergyMinus);
  eventData.SetEHBPlus(VectorHBEnergyPlus);
  eventData.SetEHBMinus(VectorHBEnergyMinus);
  eventData.SetEEEPlus(VectorEEEnergyPlus);
  eventData.SetEEEMinus(VectorEEEnergyMinus);
  eventData.SetEEBPlus(VectorEBEnergyPlus);
  eventData.SetEEBMinus(VectorEBEnergyMinus);

  eventData.SetNHFPlus(VectorHFEnergyPlus.size());
  eventData.SetNHFMinus(VectorHFEnergyMinus.size());
  eventData.SetNHEPlus(VectorHEEnergyPlus.size());
  eventData.SetNHEMinus(VectorHEEnergyMinus.size());
  eventData.SetNHBPlus(VectorHBEnergyPlus.size());
  eventData.SetNHBMinus(VectorHBEnergyMinus.size());
  eventData.SetNEEPlus(VectorEEEnergyPlus.size());
  eventData.SetNEEMinus(VectorEEEnergyMinus.size());
  eventData.SetNEBPlus(VectorEBEnergyPlus.size());
  eventData.SetNEBMinus(VectorEBEnergyMinus.size());

}

//
// Fill Particle Flow
//

void ZeroBiasAnalysis::fillParticleFlow(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  std::vector<double> VectorPFCharged;
  std::vector<double> VectorPFBarrelPlus;
  std::vector<double> VectorPFBarrelMinus;
  std::vector<double> VectorPFEndCapPlus;
  std::vector<double> VectorPFEndCapMinus;
  std::vector<double> VectorPFHFPlus;
  std::vector<double> VectorPFHFMinus;

  VectorPFCharged.clear();
  VectorPFBarrelPlus.clear();
  VectorPFBarrelMinus.clear();
  VectorPFEndCapPlus.clear();
  VectorPFEndCapMinus.clear();
  VectorPFHFPlus.clear();
  VectorPFHFMinus.clear();

  edm::Handle <reco::PFCandidateCollection> PFCandidates;
  event.getByLabel(pfTag_,PFCandidates);

  int pfsize = PFCandidates->size();
  int itPF;

  if(PFCandidates->size()>0){
    for(itPF=0; itPF < pfsize; ++itPF){
      const reco::PFCandidate* pfAll = &((*PFCandidates)[itPF]);

      // Excluding HF Calorimeter Rings 12, 13, 29, 30, 40, 41
      if( ( (fabs(pfAll->eta()) >= 2.866) && (fabs(pfAll->eta()) < 3.152) ) || (fabs(pfAll->eta()) >= 4.730) ) continue;

      if (fabs(pfAll->charge()) > 0 && pfAll->pt() > pTPFThresholdCharged_ ){
	VectorPFCharged.push_back(pfAll->energy());
      }

      if ( (fabs(pfAll->charge()) == 0 && ( pfAll->eta() <= 1.5 && pfAll->energy() > energyPFThresholdBar_) )){
	VectorPFBarrelPlus.push_back(pfAll->energy());
      }

      if ( (fabs(pfAll->charge()) == 0 && ( pfAll->eta() >= -1.5 && pfAll->energy() > energyPFThresholdBar_) )){
	VectorPFBarrelMinus.push_back(pfAll->energy());
      }

      if(pfAll->eta() > 1.5 && pfAll->eta() < 3 && pfAll->energy() > energyPFThresholdEnd_){
	VectorPFEndCapPlus.push_back(pfAll->energy());
      }

      if(pfAll->eta() < -1.5 && pfAll->eta() > -3 && pfAll->energy() > energyPFThresholdEnd_){
	VectorPFEndCapMinus.push_back(pfAll->energy());
      }

      if(pfAll->eta() > 3 && pfAll->energy() >energyPFThresholdHF_){
	VectorPFHFPlus.push_back(pfAll->energy());
      }

      if(pfAll->eta() < -3 && pfAll->energy() >energyPFThresholdHF_){
	VectorPFHFMinus.push_back(pfAll->energy());
      }

    }
  }

  eventData.SetPFCharged(VectorPFCharged);
  eventData.SetPFBarrelPlus(VectorPFBarrelPlus);
  eventData.SetPFBarrelMinus(VectorPFBarrelMinus);
  eventData.SetPFEndCapPlus(VectorPFEndCapPlus);
  eventData.SetPFEndCapMinus(VectorPFEndCapMinus);
  eventData.SetPFHFPlus(VectorPFHFPlus);
  eventData.SetPFHFMinus(VectorPFHFMinus);

  eventData.SetPFNCharged(VectorPFCharged.size());
  eventData.SetPFNBarrelPlus(VectorPFBarrelPlus.size());
  eventData.SetPFNBarrelMinus(VectorPFBarrelMinus.size());
  eventData.SetPFNEndCapPlus(VectorPFEndCapPlus.size());
  eventData.SetPFNEndCapMinus(VectorPFEndCapMinus.size());
  eventData.SetPFNHFPlus(VectorPFHFPlus.size());
  eventData.SetPFNHFMinus(VectorPFHFMinus.size());

}

//
// Fill Castor 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ZeroBiasAnalysis::fillCastor(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Phi: 16 modules, rh.id().sector(); 
  // W: 14 modules, rh.id().module(); 
  // Channel definition: 16*(rh.id().module()-1) + rh.id().sector(); 
  // For 2010, Castor uses only first five modules. 

  bool debug = false;
  bool debug_deep = false;
  bool idmod1 = false;
  bool idmod2 = false;
  bool idmod3 = false;
  bool idmod4 = false;
  bool idmod5 = false;

  std::vector<double> castor_tower;
  std::vector<double> castor_tower_module1;
  std::vector<double> castor_tower_module2;
  std::vector<double> castor_tower_module3;
  std::vector<double> castor_tower_module4;
  std::vector<double> castor_tower_module5;

  edm::Handle<CastorRecHitCollection> CastorRecHits;
  event.getByLabel(castorHitsTag_,CastorRecHits); 

  double sumCastorTower[16];
  double energyModule1[16];
  double energyModule2[16];
  double energyModule3[16];
  double energyModule4[16];
  double energyModule5[16];

  for(int isec = 0; isec < 16; isec++) {
    sumCastorTower[isec] = 0.;
    energyModule1[isec] = -999.;
    energyModule2[isec] = -999.;
    energyModule3[isec] = -999.;
    energyModule4[isec] = -999.;
    energyModule5[isec] = -999.; 
  }

  if( CastorRecHits.isValid() ){

    for (size_t i = 0; i < CastorRecHits->size(); ++i) {

      bool used_cha = false;
      const CastorRecHit & rh = (*CastorRecHits)[i];
      int cha = 16*(rh.id().module()-1) + rh.id().sector();    

      m_hVector_histo_castor_channels.at(cha-1)->Fill(rh.energy()*fCGeVCastor_);

      if(RunA_ && !RunB_){
	if(cha != 5 && cha != 6 && cha !=11 && cha !=12) used_cha = true;
	if (rh.id().module() > 5 ) continue;
      }
      if(!RunA_ && RunB_){
	if(cha != 13 && cha != 14 && (cha >=73 && cha <=80)) used_cha = true;
	if (rh.id().module() > 5 ) continue;
      }

      if((RunA_ && RunB_) || (!RunA_ && !RunB_)){
	used_cha = true;
      }

      if(used_cha == false) continue;

      if (debug_deep) std::cout << "Channel: " << cha << std::endl;
      if (debug_deep) std::cout << "Energy: " << rh.energy()*fCGeVCastor_ << " | Sector: " << rh.id().sector() << " | Module: " << rh.id().module() << " | Channel: " << cha << std::endl;

      for(int isec = 0; isec < 16; isec++) {
	if (rh.id().sector()== isec+1){
	  sumCastorTower[isec]+=rh.energy()*fCGeVCastor_; 
	  if (rh.id().module() == 1){
	    energyModule1[isec] = rh.energy()*fCGeVCastor_;
	    if (idmod1) std::cout << "Module " << rh.id().module() << ", Channel " << cha << ", isec " << isec << "." << std::endl;
	  }
	  if (rh.id().module() == 2){
	    energyModule2[isec] = rh.energy()*fCGeVCastor_;
	    if (idmod2) std::cout << "Module " << rh.id().module() << ", Channel " << cha << ", isec " << isec << "." << std::endl;
	  }
	  if (rh.id().module() == 3){
	    energyModule3[isec] = rh.energy()*fCGeVCastor_;
	    if (idmod3) std::cout << "Module " << rh.id().module() << ", Channel " << cha << ", isec " << isec << "." << std::endl;
	  }
	  if (rh.id().module() == 4){
	    energyModule4[isec] = rh.energy()*fCGeVCastor_;
	    if (idmod4) std::cout << "Module " << rh.id().module() << ", Channel " << cha << ", isec " << isec << "." << std::endl;
	  }
	  if (rh.id().module() == 5){
	    energyModule5[isec] = rh.energy()*fCGeVCastor_;
	    if (idmod5) std::cout << "Module " << rh.id().module() << ", Channel " << cha << ", isec " << isec << "." << std::endl;
	  }
	}
      }
    }

    for (int isec=0;isec<16;isec++){
      castor_tower.push_back(sumCastorTower[isec]);
      castor_tower_module1.push_back(energyModule1[isec]);
      castor_tower_module2.push_back(energyModule2[isec]);
      castor_tower_module3.push_back(energyModule3[isec]);
      castor_tower_module4.push_back(energyModule4[isec]);
      castor_tower_module5.push_back(energyModule5[isec]);
      if (debug) {
	std::cout << "Sector "<< isec+1 << ", Module 1, Energy [GeV]: " << energyModule1[isec] << std::endl;
	std::cout << "Sector "<< isec+1 << ", Module 2, Energy [GeV]: " << energyModule2[isec] << std::endl;
	std::cout << "Sector "<< isec+1 << ", Module 3, Energy [GeV]: " << energyModule3[isec] << std::endl;
	std::cout << "Sector "<< isec+1 << ", Module 4, Energy [GeV]: " << energyModule4[isec] << std::endl;
	std::cout << "Sector "<< isec+1 << ", Module 5, Energy [GeV]: " << energyModule5[isec] << std::endl;
	std::cout << "Sector "<< isec+1 << ", Total Energy [GeV]: " << sumCastorTower[isec] << std::endl;
      }
    }
    eventData.SetCastorTowerEnergy(castor_tower);
    eventData.SetCastorModule1Energy(castor_tower_module1);
    eventData.SetCastorModule2Energy(castor_tower_module2);
    eventData.SetCastorModule3Energy(castor_tower_module3);
    eventData.SetCastorModule4Energy(castor_tower_module4);
    eventData.SetCastorModule5Energy(castor_tower_module5);
  }else{
    if (debug) std::cout << "There is no Castor valid recHitSector "<< std::cout;
  }

}

//
// Fill Castor Check Channels 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ZeroBiasAnalysis::fillCastorDebug(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Phi: 16 modules, rh.id().sector();
  // W: 14 modules, rh.id().module();
  // Channel definition: 16*(rh.id().module()-1) + rh.id().sector();
  // For 2010, Castor uses only first five modules.

  bool debug = false;
  bool debug_deep = false;

  int NRecHits = 0;
  int NRecHitsPartial = 0;
  int BadChannels = 0;
  std::vector<int> Channels;
  std::vector<int> BChannels;

  edm::Handle<CastorRecHitCollection> CastorRecHits;
  event.getByLabel(castorHitsTag_,CastorRecHits);

  if( CastorRecHits.isValid() ){

    for (size_t i = 0; i < CastorRecHits->size(); ++i) {

      const CastorRecHit & rh = (*CastorRecHits)[i];
      int cha = 16*(rh.id().module()-1) + rh.id().sector();

      CastorChannelHisto_->Fill(cha);
      Channels.push_back(cha);

      ++NRecHits;
      if (rh.id().module() > 5 ) ++NRecHitsPartial;
      if (debug_deep){
	std::cout << "Channel: " << cha << std::endl;
      }
    }

    // Search Bad Channels
    const int size = (int) Channels.size();
    for (int i=1; i<=224; i++) {
      bool found=false;
      for (int j=0; j<size; j++){
	if (Channels[j]==i) {
	  if (debug) std::cout << "There is channel " << Channels[j] << std::endl;
	  found=true; 
	  break;
	}
      }
      if (!found) {
	++BadChannels;
	BChannels.push_back(i);
	if (debug) std::cout << "Channel " << i << " was not working." << std::endl;
      }
    }
    if (BadChannels < 1){
      BChannels.push_back(-999);
    }
    eventData.SetCastorNumberBadChannels(BadChannels); 
    eventData.SetCastorBadChannels(BChannels);  
  }else{
    if (debug) std::cout << "There is no Castor valid recHitSector "<< std::cout;
  }

}


//
// Fill Noise Info
//////////////////

void ZeroBiasAnalysis::fillNoiseInfo(ZeroBiasEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  edm::Handle<HcalNoiseSummary> noiseSummaryH;
  event.getByLabel("hcalnoise",noiseSummaryH);   

  bool passNoiseLoose = noiseSummaryH->passLooseNoiseFilter();
  bool passNoiseTight = noiseSummaryH->passTightNoiseFilter();

  edm::Handle<reco::BeamHaloSummary> beamHaloSummaryH;
  event.getByLabel("BeamHaloSummary",beamHaloSummaryH);

  bool beamHaloLooseId = beamHaloSummaryH->LooseId(); 
  bool beamHaloTightId = beamHaloSummaryH->TightId();

  eventData.LooseNoiseFilter_ = passNoiseLoose ? 1 : 0;
  eventData.TightNoiseFilter_ = passNoiseTight ? 1 : 0;

  eventData.BeamHaloLooseId_ = beamHaloLooseId ? 1 : 0;
  eventData.BeamHaloTightId_ = beamHaloTightId ? 1 : 0;
}


