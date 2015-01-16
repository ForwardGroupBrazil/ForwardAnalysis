#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveZAnalysis.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveZEvent.h"

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

#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "ForwardAnalysis/Utilities/interface/CastorEnergy.h"

#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include <stdio.h>
#include <math.h> 
#include <cmath>

using diffractiveZAnalysis::DiffractiveZAnalysis;

const char* DiffractiveZAnalysis::name = "DiffractiveZAnalysis";

DiffractiveZAnalysis::DiffractiveZAnalysis(const edm::ParameterSet& pset):
  triggerResultsTag_(pset.getParameter<edm::InputTag>("TriggerResultsTag")),
  hltPathNames_(pset.getParameter<std::vector<std::string> >("hltPaths")),
  electronTag_(pset.getParameter<edm::InputTag>("electronTag")),
  muonTag_(pset.getParameter<edm::InputTag>("muonTag")),
  pfTag_(pset.getParameter<edm::InputTag>("pfTag")),
  genTag_(pset.getParameter<edm::InputTag>("genTag")),
  PVtxCollectionTag_(pset.getParameter<edm::InputTag>("PVtxCollectionTag")),
  castorHitsTag_(pset.getParameter<edm::InputTag>("castorHitsTag")),
  zdcHitsTag_(pset.getParameter<edm::InputTag>("zdcHitsTag")),
  RunCastor_(pset.getUntrackedParameter<Bool_t>("RunCastor", false)),
  RunZDC_(pset.getUntrackedParameter<Bool_t>("RunZDC", false)),
  RunMC_(pset.getUntrackedParameter<Bool_t>("RunMC", false)),
  RunZPat_(pset.getUntrackedParameter<Bool_t>("RunZPat", false)),
  RunA_(pset.getUntrackedParameter<Bool_t>("RunA", false)),
  RunB_(pset.getUntrackedParameter<Bool_t>("RunB", false)),
  EachTower_(pset.getUntrackedParameter<Bool_t>("EachTower", false)),
  pTPFThresholdCharged_(pset.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(pset.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(pset.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(pset.getParameter<double>("energyPFThresholdHF")),
  energyThresholdHB_(pset.getParameter<double>("energyThresholdHB")),
  energyThresholdHE_(pset.getParameter<double>("energyThresholdHE")),
  energyThresholdHF_(pset.getParameter<double>("energyThresholdHF")),
  energyThresholdEB_(pset.getParameter<double>("energyThresholdEB")),
  energyThresholdEE_(pset.getParameter<double>("energyThresholdEE")),
  castorThreshold_(pset.getParameter<double>("castorThreshold")),
  fCGeVCastor_(pset.getParameter<double>("fCGeVCastor")),
  caloTowerTag_(pset.getParameter<edm::InputTag>("CaloTowerTag")),
  trackTag_(pset.getParameter<edm::InputTag>("TrackTag")),
  beamEnergy_(pset.getParameter<double>("beamEnergy"))  //one beam energy, IP = 2*beamEnergy_
{
}

void DiffractiveZAnalysis::setTFileService(){

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

DiffractiveZAnalysis::~DiffractiveZAnalysis(){}

void DiffractiveZAnalysis::begin() {
  setTFileService();
}

void DiffractiveZAnalysis::begin(const edm::Run& run, const edm::EventSetup& setup) {}

void DiffractiveZAnalysis::end() {}

void DiffractiveZAnalysis::fill(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  eventData.reset();

  fillTriggerInfo(eventData,event,setup);
  fillMuonsInfo(eventData,event,setup);
  fillElectronsInfo(eventData,event,setup);
  if (RunZPat_) fillZPat(eventData,event,setup);
  fillTracksInfo(eventData,event,setup);
  fillDetectorVariables(eventData,event,setup);
  fillVariables(eventData,event,setup);
  if (RunMC_) fillGenInfo(eventData,event,setup); 
  if (RunCastor_){
    fillCastor(eventData,event,setup);
    fillCastorDebug(eventData,event,setup);
  }
  if (RunZDC_) fillZDC(eventData,event,setup);
}

// Fill Trigger
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillTriggerInfo(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

// Fill Reco::Electron
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillElectronsInfo(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;
  ElectronVector.clear();

  edm::Handle<reco::GsfElectronCollection> electrons;
  event.getByLabel(electronTag_,electrons);

  int electronsize = electrons->size();
  int itElectron;
  if(electrons->size()>0){
    for(itElectron=0; itElectron < electronsize; ++itElectron){
      const reco::GsfElectron* electronAll = &((*electrons)[itElectron]);
      ElectronVector.push_back(electronAll);
    }
  }

  // Sorting Vector
  std::sort(ElectronVector.begin(), ElectronVector.end(), orderPT());

  if(ElectronVector.size()>1){

    double relIsoFirstElectronDr03 = (ElectronVector[0]->dr03TkSumPt()+ElectronVector[0]->dr03EcalRecHitSumEt()+ElectronVector[0]->dr03HcalTowerSumEt())/ElectronVector[0]->et();
    double relIsoFirstElectronDr04 = (ElectronVector[0]->dr04TkSumPt()+ElectronVector[0]->dr04EcalRecHitSumEt()+ElectronVector[0]->dr04HcalTowerSumEt())/ElectronVector[0]->et();
    double relIsoSecondElectronDr03 = (ElectronVector[1]->dr03TkSumPt()+ElectronVector[1]->dr03EcalRecHitSumEt()+ElectronVector[1]->dr03HcalTowerSumEt())/ElectronVector[1]->et();
    double relIsoSecondElectronDr04 = (ElectronVector[1]->dr04TkSumPt()+ElectronVector[1]->dr04EcalRecHitSumEt()+ElectronVector[1]->dr04HcalTowerSumEt())/ElectronVector[1]->et();
    double InnerHits1 = ElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    double InnerHits2 = ElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    // Dielectron Mass
    math::XYZTLorentzVector DielectronSystem(0.,0.,0.,0.);
    DielectronSystem += ElectronVector[0]->p4();
    DielectronSystem += ElectronVector[1]->p4();

    eventData.SetDiElectronMass(DielectronSystem.M());
    eventData.SetDiElectronPt(DielectronSystem.pt());
    eventData.SetDiElectronEta(DielectronSystem.eta());
    eventData.SetDiElectronPhi(DielectronSystem.phi());

    eventData.SetElectronsN(ElectronVector.size());
    eventData.SetLeadingElectronPt(ElectronVector[0]->pt());
    eventData.SetLeadingElectronEta(ElectronVector[0]->eta());
    eventData.SetLeadingElectronPhi(ElectronVector[0]->phi());
    eventData.SetLeadingElectronCharge(ElectronVector[0]->charge());
    eventData.SetLeadingElectronP4(ElectronVector[0]->p4());
    eventData.SetSecondElectronPt(ElectronVector[1]->pt());
    eventData.SetSecondElectronEta(ElectronVector[1]->eta());
    eventData.SetSecondElectronPhi(ElectronVector[1]->phi());
    eventData.SetSecondElectronCharge(ElectronVector[1]->charge());
    eventData.SetSecondElectronP4(ElectronVector[1]->p4());

    eventData.SetLeadingElectronTkDr03(ElectronVector[0]->dr03TkSumPt());
    eventData.SetLeadingElectronEcalDr03(ElectronVector[0]->dr03EcalRecHitSumEt());
    eventData.SetLeadingElectronHcalDr03(ElectronVector[0]->dr03HcalTowerSumEt());

    eventData.SetSecondElectronTkDr03(ElectronVector[1]->dr03TkSumPt());
    eventData.SetSecondElectronEcalDr03(ElectronVector[1]->dr03EcalRecHitSumEt());
    eventData.SetSecondElectronHcalDr03(ElectronVector[1]->dr03HcalTowerSumEt());

    eventData.SetLeadingElectronTkDr04(ElectronVector[0]->dr04TkSumPt());
    eventData.SetLeadingElectronEcalDr04(ElectronVector[0]->dr04EcalRecHitSumEt());
    eventData.SetLeadingElectronHcalDr04(ElectronVector[0]->dr04HcalTowerSumEt());

    eventData.SetSecondElectronTkDr04(ElectronVector[1]->dr04TkSumPt());
    eventData.SetSecondElectronEcalDr04(ElectronVector[1]->dr04EcalRecHitSumEt());
    eventData.SetSecondElectronHcalDr04(ElectronVector[1]->dr04HcalTowerSumEt());

    eventData.SetLeadingElectronrelIsoDr03(relIsoFirstElectronDr03);
    eventData.SetLeadingElectronrelIsoDr04(relIsoFirstElectronDr04);
    eventData.SetSecondElectronrelIsoDr03(relIsoSecondElectronDr03);
    eventData.SetSecondElectronrelIsoDr04(relIsoSecondElectronDr04);

    eventData.SetLeadingElectronDeltaPhiTkClu(ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetLeadingElectronDeltaEtaTkClu(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetLeadingElectronSigmaIeIe(ElectronVector[0]->sigmaIetaIeta());
    eventData.SetLeadingElectronDCot(ElectronVector[0]->convDcot());
    eventData.SetLeadingElectronDist(ElectronVector[0]->convDist());
    eventData.SetLeadingElectronInnerHits(InnerHits1);
    eventData.SetLeadingElectronHE(ElectronVector[0]->hadronicOverEm());

    eventData.SetSecondElectronDeltaPhiTkClu(ElectronVector[1]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetSecondElectronDeltaEtaTkClu(ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetSecondElectronSigmaIeIe(ElectronVector[1]->sigmaIetaIeta());
    eventData.SetSecondElectronDCot(ElectronVector[1]->convDcot());
    eventData.SetSecondElectronDist(ElectronVector[1]->convDist());
    eventData.SetSecondElectronInnerHits(InnerHits2);
    eventData.SetSecondElectronHE(ElectronVector[1]->hadronicOverEm());

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCount03 = 0;
    int goodTracksCount04 = 0;
    int goodTracksCount05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),ElectronVector[0]->eta(),ElectronVector[0]->phi()) > 0.3) && (deltaR(track->eta(),track->phi(),ElectronVector[1]->eta(),ElectronVector[1]->phi()) > 0.3))
      {
	goodTracksCount03++;
      }

      if ((deltaR(track->eta(),track->phi(),ElectronVector[0]->eta(),ElectronVector[0]->phi()) > 0.4) && (deltaR(track->eta(),track->phi(),ElectronVector[1]->eta(),ElectronVector[1]->phi()) > 0.4))
      {
	goodTracksCount04++;
      }

      if ((deltaR(track->eta(),track->phi(),ElectronVector[0]->eta(),ElectronVector[0]->phi()) > 0.5) && (deltaR(track->eta(),track->phi(),ElectronVector[1]->eta(),ElectronVector[1]->phi()) > 0.5))
      {
	goodTracksCount05++;
      }

    }

    eventData.SetTracksNonConeElectron03(goodTracksCount03);
    eventData.SetTracksNonConeElectron04(goodTracksCount04);
    eventData.SetTracksNonConeElectron05(goodTracksCount05);

    if (debug){
      std::cout << ">>> Reco Electron" << std::endl;
      std::cout << "electron1 -> dr03 TK: " << ElectronVector[0]->dr03TkSumPt() << "| dr03 Ecal: " << ElectronVector[0]->dr03EcalRecHitSumEt() << " | dr03 Hcal: " << ElectronVector[0]->dr03HcalTowerSumEt() << std::endl;
      std::cout << "electron1 -> dr04 TK: " << ElectronVector[0]->dr04TkSumPt() << "| dr04 Ecal: " << ElectronVector[0]->dr04EcalRecHitSumEt() << " | dr04 Hcal: " << ElectronVector[0]->dr04HcalTowerSumEt() <<  std::endl;
      std::cout << "electron2 -> dr03 TK: " << ElectronVector[1]->dr03TkSumPt() << "| dr03 Ecal: " << ElectronVector[1]->dr03EcalRecHitSumEt() << " | dr03 Hcal: " << ElectronVector[1]->dr03HcalTowerSumEt() << std::endl;
      std::cout << "electron2 -> dr04 TK: " << ElectronVector[1]->dr04TkSumPt() << "| dr04 Ecal: " << ElectronVector[1]->dr04EcalRecHitSumEt() << " | dr04 Hcal: " << ElectronVector[1]->dr04HcalTowerSumEt() <<  std::endl;
      std::cout << "Electron Isolation: " << relIsoFirstElectronDr03 << " | " << relIsoFirstElectronDr04 << std::endl;
      std::cout << "N Electrons: " << ElectronVector.size() << std::endl;
      std::cout << "Electron, pT 1: " << ElectronVector[0]->pt() << std::endl;
      std::cout << "Electron, pT 2: " << ElectronVector[1]->pt() << std::endl;
      std::cout << "Electron, eta 1: " << ElectronVector[0]->eta() << std::endl;
      std::cout << "Electron, eta 2: " << ElectronVector[1]->eta() << std::endl;
      std::cout << "Eta Z: " << DielectronSystem.eta() << std::endl;
      std::cout << "Phi Z: " << DielectronSystem.phi() << std::endl;
      std::cout << "pT Z: " << DielectronSystem.pt() << std::endl;
      std::cout << "energy Z: " << DielectronSystem.energy() << std::endl;
      std::cout << "pz Z: " << DielectronSystem.pz() << std::endl;
      std::cout << "DeltaPhiTkClu, electron1: " << ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx() << std::endl;
      std::cout << "DeltaEtaTkClu, electron1: " << ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx() << std::endl;
      std::cout << "SigmaIeIe, electron1: " << ElectronVector[0]->sigmaIetaIeta() << std::endl;
      std::cout << "Dcot, electron1: " << ElectronVector[0]->convDcot() << std::endl;
      std::cout << "Dist, electron1: " << ElectronVector[0]->convDist() << std::endl;
      std::cout << "Number Of Expected Inner Hits, electron1: " << ElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << std::endl;
      std::cout << "H/E, electron1: " << ElectronVector[0]->hadronicOverEm() << std::endl;
      std::cout << "DeltaPhiTkClu, electron2: " << ElectronVector[1]->deltaPhiSuperClusterTrackAtVtx() << std::endl;
      std::cout << "DeltaEtaTkClu, electron2: " << ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx() << std::endl;
      std::cout << "SigmaIeIe, electron2: " << ElectronVector[1]->sigmaIetaIeta() << std::endl;
      std::cout << "Dcot, electron2: " << ElectronVector[1]->convDcot() << std::endl;
      std::cout << "Dist, electron2: " << ElectronVector[1]->convDist() << std::endl;
      std::cout << "Number Of Expected Inner Hits, electron2: " << ElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << std::endl;
      std::cout << "H/E, electron2: " << ElectronVector[1]->hadronicOverEm() << std::endl;
      std::cout << "" << std::endl;
    }
  }
  else {
    eventData.SetDiElectronMass(-999.);
    eventData.SetDiElectronPt(-999.);
    eventData.SetDiElectronEta(-999.);
    eventData.SetDiElectronPhi(-999.);
    eventData.SetElectronsN(-1);
    eventData.SetLeadingElectronPt(-999.);
    eventData.SetLeadingElectronEta(-999.);
    eventData.SetLeadingElectronPhi(-999.);
    eventData.SetLeadingElectronCharge(-999);
    eventData.SetSecondElectronPt(-999.);
    eventData.SetSecondElectronEta(-999.);
    eventData.SetSecondElectronPhi(-999.);
    eventData.SetSecondElectronCharge(-999);

    eventData.SetLeadingElectronTkDr03(-999.);
    eventData.SetLeadingElectronEcalDr03(-999.);
    eventData.SetLeadingElectronHcalDr03(-999.);

    eventData.SetSecondElectronTkDr03(-999.);
    eventData.SetSecondElectronEcalDr03(-999.);
    eventData.SetSecondElectronHcalDr03(-999.);

    eventData.SetLeadingElectronTkDr04(-999.);
    eventData.SetLeadingElectronEcalDr04(-999.);
    eventData.SetLeadingElectronHcalDr04(-999.);

    eventData.SetSecondElectronTkDr04(-999.);
    eventData.SetSecondElectronEcalDr04(-999.);
    eventData.SetSecondElectronHcalDr04(-999.);

    eventData.SetLeadingElectronrelIsoDr03(-999.);
    eventData.SetLeadingElectronrelIsoDr04(-999.);
    eventData.SetSecondElectronrelIsoDr03(-999.);
    eventData.SetSecondElectronrelIsoDr04(-999.);

    eventData.SetLeadingElectronDeltaPhiTkClu(-999.);
    eventData.SetLeadingElectronDeltaEtaTkClu(-999.);
    eventData.SetLeadingElectronSigmaIeIe(-999.);
    eventData.SetLeadingElectronDCot(-999.);
    eventData.SetLeadingElectronDist(-999.);
    eventData.SetLeadingElectronInnerHits(-999.);
    eventData.SetLeadingElectronHE(-999.);
    eventData.SetSecondElectronDeltaPhiTkClu(-999.);
    eventData.SetSecondElectronDeltaEtaTkClu(-999.);
    eventData.SetSecondElectronSigmaIeIe(-999.);
    eventData.SetSecondElectronDCot(-999.);
    eventData.SetSecondElectronDist(-999.);
    eventData.SetSecondElectronInnerHits(-999.);
    eventData.SetSecondElectronHE(-999.);

  }
}

// Fill Reco::Muon
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillMuonsInfo(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;
  MuonVector.clear();

  edm::Handle<reco::MuonCollection> muons;
  event.getByLabel(muonTag_,muons);

  edm::Handle<reco::VertexCollection>  vertex;
  event.getByLabel(PVtxCollectionTag_, vertex);

  int muonsize = muons->size();
  int itMuon;
  if(muons->size()>0){
    for(itMuon=0; itMuon < muonsize; ++itMuon){
      const reco::Muon* muonAll = &((*muons)[itMuon]);
      MuonVector.push_back(muonAll);
    }
  }

  // Sorting Vector
  std::sort(MuonVector.begin(), MuonVector.end(), orderPT());

  if(MuonVector.size()>1){

    double muon1SumPtR03 = MuonVector[0]->isolationR03().sumPt;
    double muon1EmEtR03 = MuonVector[0]->isolationR03().emEt;
    double muon1HadEtR03 = MuonVector[0]->isolationR03().hadEt;
    double muon1SumPtR05 = MuonVector[0]->isolationR05().sumPt;
    double muon1EmEtR05 = MuonVector[0]->isolationR05().emEt;
    double muon1HadEtR05 = MuonVector[0]->isolationR05().hadEt;

    double muon2SumPtR03 = MuonVector[1]->isolationR03().sumPt;
    double muon2EmEtR03 = MuonVector[1]->isolationR03().emEt;
    double muon2HadEtR03 = MuonVector[1]->isolationR03().hadEt;
    double muon2SumPtR05 = MuonVector[1]->isolationR05().sumPt;
    double muon2EmEtR05 = MuonVector[1]->isolationR05().emEt;
    double muon2HadEtR05 = MuonVector[1]->isolationR05().hadEt;

    double relIsoFirstMuonDr03 = (muon1SumPtR03 + muon1EmEtR03 + muon1HadEtR03)/MuonVector[0]->pt();
    double relIsoSecondMuonDr03 = (muon2SumPtR03 + muon2EmEtR03 + muon2HadEtR03)/MuonVector[1]->pt();
    double relIsoFirstMuonDr05 = (muon1SumPtR05 + muon1EmEtR05 + muon1HadEtR05)/MuonVector[0]->pt();
    double relIsoSecondMuonDr05 = (muon2SumPtR05 + muon2EmEtR05 + muon2HadEtR05)/MuonVector[1]->pt();

    // Dimuon Mass
    math::XYZTLorentzVector DimuonSystem(0.,0.,0.,0.);
    DimuonSystem += MuonVector[0]->p4();
    DimuonSystem += MuonVector[1]->p4();
    eventData.SetDiMuonMass(DimuonSystem.M());
    eventData.SetDiMuonPt(DimuonSystem.pt());
    eventData.SetDiMuonEta(DimuonSystem.eta());
    eventData.SetDiMuonPhi(DimuonSystem.phi());
    eventData.SetMuonsN(MuonVector.size());
    eventData.SetLeadingMuonPt(MuonVector[0]->pt());
    eventData.SetLeadingMuonEta(MuonVector[0]->eta());
    eventData.SetLeadingMuonPhi(MuonVector[0]->phi());
    eventData.SetLeadingMuonCharge(MuonVector[0]->charge());
    eventData.SetLeadingMuonP4(MuonVector[0]->p4());
    eventData.SetSecondMuonPt(MuonVector[1]->pt());
    eventData.SetSecondMuonEta(MuonVector[1]->eta());
    eventData.SetSecondMuonPhi(MuonVector[1]->phi());
    eventData.SetSecondMuonCharge(MuonVector[1]->charge());
    eventData.SetSecondMuonP4(MuonVector[1]->p4());

    eventData.SetLeadingMuonSumPtR03(muon1SumPtR03);
    eventData.SetLeadingMuonEmEtR03(muon1EmEtR03);
    eventData.SetLeadingMuonHadEtR03(muon1HadEtR03);
    eventData.SetLeadingMuonSumPtR05(muon1SumPtR05);
    eventData.SetLeadingMuonEmEtR05(muon1EmEtR05);
    eventData.SetLeadingMuonHadEtR05(muon1HadEtR05);

    eventData.SetSecondMuonSumPtR03(muon2SumPtR03);
    eventData.SetSecondMuonEmEtR03(muon2EmEtR03);
    eventData.SetSecondMuonHadEtR03(muon2HadEtR03);
    eventData.SetSecondMuonSumPtR05(muon2SumPtR05);
    eventData.SetSecondMuonEmEtR05(muon2EmEtR05);
    eventData.SetSecondMuonHadEtR05(muon2HadEtR05);

    eventData.SetLeadingMuonrelIsoDr03(relIsoFirstMuonDr03);
    eventData.SetSecondMuonrelIsoDr03(relIsoSecondMuonDr03);
    eventData.SetLeadingMuonrelIsoDr05(relIsoFirstMuonDr05);
    eventData.SetSecondMuonrelIsoDr05(relIsoSecondMuonDr05);
    eventData.SetLeadingMuonTrackerHits(MuonVector[0]->track()->hitPattern().trackerLayersWithMeasurement());
    eventData.SetLeadingMuonPixelHits(MuonVector[0]->innerTrack()->hitPattern().numberOfValidPixelHits());
    eventData.SetLeadingMuonNormalizedChi2(MuonVector[0]->globalTrack()->normalizedChi2());
    eventData.SetLeadingMuonMatchedStations(MuonVector[0]->numberOfMatchedStations());
    eventData.SetLeadingMuonDxy(MuonVector[0]->innerTrack()->dxy(vertex->at(0).position()));
    eventData.SetLeadingMuonIsGlobal(MuonVector[0]->isGlobalMuon());
    eventData.SetLeadingMuonIsTracker(MuonVector[0]->isTrackerMuon());
    eventData.SetSecondMuonTrackerHits(MuonVector[1]->track()->hitPattern().trackerLayersWithMeasurement());
    eventData.SetSecondMuonPixelHits(MuonVector[1]->innerTrack()->hitPattern().numberOfValidPixelHits());
    eventData.SetSecondMuonNormalizedChi2(MuonVector[1]->globalTrack()->normalizedChi2());
    eventData.SetSecondMuonMatchedStations(MuonVector[1]->numberOfMatchedStations());
    eventData.SetSecondMuonDxy(MuonVector[1]->innerTrack()->dxy(vertex->at(0).position()));
    eventData.SetSecondMuonIsGlobal(MuonVector[1]->isGlobalMuon());
    eventData.SetSecondMuonIsTracker(MuonVector[1]->isTrackerMuon());

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCount03 = 0;
    int goodTracksCount04 = 0;
    int goodTracksCount05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),MuonVector[0]->eta(),MuonVector[0]->phi()) > 0.3) && (deltaR(track->eta(),track->phi(),MuonVector[1]->eta(),MuonVector[1]->phi()) > 0.3))
      {
	goodTracksCount03++;
      }

      if ((deltaR(track->eta(),track->phi(),MuonVector[0]->eta(),MuonVector[0]->phi()) > 0.4) && (deltaR(track->eta(),track->phi(),MuonVector[1]->eta(),MuonVector[1]->phi()) > 0.4))
      {
	goodTracksCount04++;
      }

      if ((deltaR(track->eta(),track->phi(),MuonVector[0]->eta(),MuonVector[0]->phi()) > 0.5) && (deltaR(track->eta(),track->phi(),MuonVector[1]->eta(),MuonVector[1]->phi()) > 0.5))
      {
	goodTracksCount05++;
      }

    }

    eventData.SetTracksNonConeMuon03(goodTracksCount03);
    eventData.SetTracksNonConeMuon04(goodTracksCount04);
    eventData.SetTracksNonConeMuon05(goodTracksCount05);

    if (debug){
      std::cout << "NMuons: " << MuonVector.size() << std::endl;
      std::cout << "Muon, pT 1: " << MuonVector[0]->pt() << std::endl;
      std::cout << "Muon, pT 2: " << MuonVector[1]->pt() << std::endl;
      std::cout << "Muon, eta 1: " << MuonVector[0]->eta() << std::endl;
      std::cout << "Muon, eta 2: " << MuonVector[1]->eta() << std::endl;
      std::cout << "Eta Z: " << DimuonSystem.eta() << std::endl;
      std::cout << "Phi Z: " << DimuonSystem.phi() << std::endl;
      std::cout << "pT Z: " << DimuonSystem.pt() << std::endl;
      std::cout << "energy Z: " << DimuonSystem.energy() << std::endl;
      std::cout << "pz Z: " << DimuonSystem.pz() << std::endl;
    }

  }
  else{
    eventData.SetDiMuonMass(-999.);
    eventData.SetDiMuonPt(-999.);
    eventData.SetDiMuonEta(-999.);
    eventData.SetDiMuonPhi(-999.);
    eventData.SetMuonsN(-1);
    eventData.SetLeadingMuonPt(-999.);
    eventData.SetLeadingMuonEta(-999.);
    eventData.SetLeadingMuonPhi(-999.);
    eventData.SetLeadingMuonCharge(-999);
    eventData.SetSecondMuonPt(-999.);
    eventData.SetSecondMuonEta(-999.);
    eventData.SetSecondMuonPhi(-999.);
    eventData.SetSecondMuonCharge(-999);

    eventData.SetLeadingMuonSumPtR03(-999.);
    eventData.SetLeadingMuonEmEtR03(-999.);
    eventData.SetLeadingMuonHadEtR03(-999.);
    eventData.SetLeadingMuonSumPtR05(-999.);
    eventData.SetLeadingMuonEmEtR05(-999.);
    eventData.SetLeadingMuonHadEtR05(-999.);

    eventData.SetSecondMuonSumPtR03(-999.);
    eventData.SetSecondMuonEmEtR03(-999.);
    eventData.SetSecondMuonHadEtR03(-999.);
    eventData.SetSecondMuonSumPtR05(-999.);
    eventData.SetSecondMuonEmEtR05(-999.);
    eventData.SetSecondMuonHadEtR05(-999.);

    eventData.SetLeadingMuonrelIsoDr03(-999.);
    eventData.SetSecondMuonrelIsoDr03(-999.);
    eventData.SetLeadingMuonrelIsoDr05(-999.);
    eventData.SetSecondMuonrelIsoDr05(-999.);
    eventData.SetLeadingMuonTrackerHits(-999.);
    eventData.SetLeadingMuonPixelHits(-999.);
    eventData.SetLeadingMuonNormalizedChi2(-999.);
    eventData.SetLeadingMuonMatchedStations(-999.);
    eventData.SetLeadingMuonDxy(-999.);
    eventData.SetLeadingMuonIsGlobal(false);
    eventData.SetLeadingMuonIsTracker(false);
    eventData.SetSecondMuonTrackerHits(-999.);
    eventData.SetSecondMuonPixelHits(-999.);
    eventData.SetSecondMuonNormalizedChi2(-999.);
    eventData.SetSecondMuonMatchedStations(-999.);
    eventData.SetSecondMuonDxy(-999.);
    eventData.SetSecondMuonIsGlobal(false);
    eventData.SetSecondMuonIsTracker(false);

  } 

}

// Fill Tracks Info
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillTracksInfo(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  std::vector<double> vertexNDOF;
  std::vector<double> vertexChiNorm;
  std::vector<double> vertexMultiplicity;
  double nhit=0;
  std::vector<double> V_x;
  std::vector<double> V_y;
  std::vector<double> V_z; 
  std::vector<double> pt;
  std::vector<double> tracks;
  std::vector<double> trackseta;
  std::vector<std::vector<double> > tracksPT;

  trackseta.clear();

  edm::Handle<reco::VertexCollection>  vertexCollectionHandle;
  event.getByLabel(PVtxCollectionTag_, vertexCollectionHandle);

  for(reco::VertexCollection::const_iterator vtx = vertexCollectionHandle->begin();vtx!=vertexCollectionHandle->end(); ++vtx)
  {
    reco::Vertex::trackRef_iterator it = vtx->tracks_begin();
    reco::Vertex::trackRef_iterator lastTrack = vtx->tracks_end();


    for(;it!=lastTrack;it++) {
      nhit+=(*it)->numberOfValidHits();
      pt.push_back((*it)->pt());
      trackseta.push_back((*it)->eta());
    }

    //Sorting the pt tracks, in order to take only the 31 most energetics
    const int  size = (int) pt.size();
    int *sorted = new int[size];
    double *v = new double[size];

    for (int i=0; i<size; i++) {
      v[i] = pt[i];
    }
    TMath::Sort(size, v, sorted, true);
    for (int i=0; i<size; i++) {
      tracks.push_back(pt[sorted[i]]);
      if (i>30) break;
    }

    tracksPT.push_back(tracks);
    tracks.clear();
    pt.clear();

    double ndof=vtx->ndof();
    double chiNorm=vtx->normalizedChi2();
    double NumbOfTracks=vtx->tracksSize();
    vertexNDOF.push_back(ndof);
    vertexChiNorm.push_back(chiNorm);
    vertexMultiplicity.push_back(NumbOfTracks);
    nhit=0;
    if ( ndof != 0 ) {
      V_x.push_back(vtx->x());
      V_y.push_back(vtx->y());
      V_z.push_back(vtx->z());
    } else {
      V_x.push_back(-999);
      V_y.push_back(-999);
      V_z.push_back(-999);
    }

  } // loop over vtx

  if(trackseta.size()>0){
    std::sort(trackseta.rbegin(), trackseta.rend());
    if(debug){
      for(unsigned int i=0;i<trackseta.size();i++){
	std::cout << "Tracks Eta: " << trackseta[i] << std::endl;  
      }
    }
    eventData.SetTrackEtaMax(trackseta[0]);
    eventData.SetTrackEtaMin(trackseta[trackseta.size()-1]);
  }else{
    eventData.SetTrackEtaMax(-999.);
    eventData.SetTrackEtaMin(-999.);
  }
  eventData.SetVertexMultiplicity(vertexMultiplicity);
  eventData.SetVertexChiNorm(vertexChiNorm);
  eventData.SetVertexNDOF(vertexNDOF);
  eventData.SetVz(V_z);
  eventData.SetVx(V_x);
  eventData.SetVy(V_y); 
  eventData.SetTracksPt(tracksPT);
}

// Fill Gen Level Information
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillGenInfo(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  double sumECastorGen = 0.;
  double sumECastorGenCMS = 0.;

  genVector.clear();
  genCMSVector.clear();
  genProtonVector.clear();

  // Fill Gen
  edm::Handle<reco::GenParticleCollection> genParticle;
  event.getByLabel(genTag_, genParticle);
  int gensize = genParticle->size();
  int itGen;
  if(genParticle->size()>0){
    for(itGen=0; itGen < gensize; ++itGen){
      const reco::GenParticle* genAll = &((*genParticle)[itGen]);
      if (genAll->status() != 1) continue;
      if (genAll->pdgId() == 2212 && genAll->pz()>0.7*beamEnergy_) genProtonVector.push_back(genAll);  
      if (fabs(genAll->pdgId()) == 2212) continue;
      if (genAll->eta()<-5.2 && genAll->eta()>-6.6) sumECastorGen+=genAll->energy();
      genVector.push_back(genAll);
      if( ( (fabs(genAll->eta()) >= 2.866) && (fabs(genAll->eta()) < 3.152) ) || (fabs(genAll->eta()) >= 4.730) ) continue;
      if ( ( fabs(genAll->eta()) <= 1.5 && genAll->energy() > energyPFThresholdBar_ ) ||
	  (fabs(genAll->eta()) > 1.5 && fabs(genAll->eta()) <= 3 && genAll->energy() > energyPFThresholdEnd_) ||
	  (fabs(genAll->eta()) > 3 && genAll->energy() >energyPFThresholdHF_) ) {
	if (genAll->eta()<-5.2 && genAll->eta()>-6.6) sumECastorGenCMS+=genAll->energy();
	genCMSVector.push_back(genAll);      
      }

    }
  }

  if(debug){
    std::cout << "\nGen Particles Info:" << std::endl;
    if(genVector.size()>0){
      for(unsigned int i=0;i<genVector.size();i++){
	std::cout << "GenParticle --> pdgId: " << genVector[i]->pdgId() << " | eta: " << genVector[i]->eta() << " | Energy [GeV]: " << genVector[i]->energy() << std::endl; 
      }
    }

    if(genCMSVector.size()>0){
      for(unsigned int i=0;i<genCMSVector.size();i++){
	std::cout << "GenParticle at CMS --> pdgId: " << genCMSVector[i]->pdgId() << " | eta: " << genCMSVector[i]->eta() << " | Energy [GeV]: " << genCMSVector[i]->energy() << std::endl;
      }
    }
  }

  double xi_p_gen_plus = -999.;
  double xi_p_gen_minus = -999.;

  // Gen Proton Info
  if(genProtonVector.size()>0){
    std::sort(genProtonVector.begin(), genProtonVector.end(), orderAbsolutPZ());
  }

  if(genProtonVector.size()==1){
    if(genProtonVector[0]->pz()>0){
      xi_p_gen_plus = (1. - (genProtonVector[0]->pz()/beamEnergy_));
      eventData.SetXiGenPlus(xi_p_gen_plus);
      eventData.SetXiGenMinus(-999.);
    }else{
      xi_p_gen_minus = (1. + (genProtonVector[0]->pz()/beamEnergy_));
      eventData.SetXiGenPlus(-999.);
      eventData.SetXiGenMinus(xi_p_gen_minus);
    }
  }

  if(genProtonVector.size()>1){
    if(genProtonVector[0]->pz()>0){
      xi_p_gen_plus = (1. - (genProtonVector[0]->pz()/beamEnergy_));
    }else{
      xi_p_gen_minus = (1. + (genProtonVector[0]->pz()/beamEnergy_));
    }
    if(genProtonVector[1]->pz()>0){
      xi_p_gen_plus = (1. - (genProtonVector[1]->pz()/beamEnergy_));
    }else{
      xi_p_gen_minus = (1. + (genProtonVector[1]->pz()/beamEnergy_));
    }
    eventData.SetXiGenPlus(xi_p_gen_plus);
    eventData.SetXiGenMinus(xi_p_gen_minus);
  }

  if(debug){
    std::cout << "\nProtons Gen Info:" << std::endl;
    std::cout << "Xi +: " << xi_p_gen_plus << std::endl;
    std::cout << "Xi -: " << xi_p_gen_minus << std::endl;
  }

  // Gen Info
  math::XYZTLorentzVector SPlus(0.,0.,0.,0.);
  math::XYZTLorentzVector SMinus(0.,0.,0.,0.);

  math::XYZTLorentzVector SPlusCMS(0.,0.,0.,0.);
  math::XYZTLorentzVector SMinusCMS(0.,0.,0.,0.);

  std::vector<std::pair<double, double> > GapVector;
  std::vector<std::pair<double, double> > GapCMSVector;
  GapVector.clear();
  GapCMSVector.clear();

  double E_pz_plus_gen = 0.;
  double E_pz_minus_gen = 0.;
  double et_expo_plus_gen = 0.;
  double et_expo_minus_gen = 0.;

  double E_pz_plus_gen_CMS = 0.;
  double E_pz_minus_gen_CMS = 0.;
  double et_expo_plus_gen_CMS = 0.;
  double et_expo_minus_gen_CMS = 0.;

  if(genVector.size()>0){
    std::sort(genVector.begin(), genVector.end(), orderETA());

    for(unsigned int i=0;i<genVector.size();i++){
      math::XYZTLorentzVector tmp(genVector[i]->px(), genVector[i]->py(), genVector[i]->pz(), genVector[i]->energy());
      E_pz_plus_gen += genVector[i]->energy()+genVector[i]->pz();
      et_expo_plus_gen += genVector[i]->et()*pow(2.71,genVector[i]->eta());
      E_pz_minus_gen += genVector[i]->energy()-genVector[i]->pz();
      et_expo_minus_gen += genVector[i]->et()*pow(2.71,-genVector[i]->eta());
      if(genVector[i]->eta()>0.){
	SPlus+=tmp;
      }else{
	SMinus+=tmp;
      }
    }

    if(debug){
      if(genVector.size()>0){
	std::cout << "\nGen Particles After Sort by Eta Info:" << std::endl;
	for(unsigned int i=0;i<genVector.size();i++){
	  std::cout << "Gen. Particles, eta: " << genVector[i]->eta() << std::endl;
	}
      }
      std::cout << "\nReconstructed Variables Gen Info:" << std::endl;
      std::cout << "Eta, Max: " << genVector[0]->eta() << std::endl;    
      std::cout << "Eta, Min: " << genVector[genVector.size()-1]->eta() << std::endl;
      std::cout << "PF Mass -: " << SMinus.M() << " [GeV]" << std::endl;
      std::cout << "PF Mass +: " << SPlus.M() << " [GeV]" << std::endl;
      std::cout << "PF Mass^2 -: " << SMinus.M2() << " [GeV]" << std::endl;
      std::cout << "PF Mass^2 +: " << SPlus.M2() << " [GeV]" << std::endl;
      std::cout << "E + pz: " << E_pz_plus_gen << " [GeV]" << std::endl;
      std::cout << "E - pz: " << E_pz_minus_gen << " [GeV]" << std::endl;
      std::cout << "Et*exp(+eta): " << et_expo_plus_gen << " [GeV]" << std::endl;
      std::cout << "Et*exp(-eta): " << et_expo_minus_gen << " [GeV]" << std::endl;
    }

    eventData.SetEtaMaxGen(genVector[0]->eta());
    eventData.SetEtaMinGen(genVector[genVector.size()-1]->eta());
    eventData.SetMxGenMinus(SMinus.M());
    eventData.SetMxGenPlus(SPlus.M());
    eventData.SetMx2GenMinus(SMinus.M2());
    eventData.SetMx2GenPlus(SPlus.M2());
    eventData.SetEpluspzGen(E_pz_plus_gen);
    eventData.SetEminuspzGen(E_pz_minus_gen);
    eventData.SetEtExpoPlusGen(et_expo_plus_gen);
    eventData.SetEtExpoMinusGen(et_expo_minus_gen);

  }else{

    eventData.SetEtaMaxGen(-999.);
    eventData.SetEtaMinGen(-999.);
    eventData.SetMxGenMinus(-999.);
    eventData.SetMxGenPlus(-999.);
    eventData.SetMx2GenMinus(-999.);
    eventData.SetMx2GenPlus(-999.);
    eventData.SetEpluspzGen(-999.);
    eventData.SetEminuspzGen(-999.);
    eventData.SetEtExpoPlusGen(-999.);
    eventData.SetEtExpoMinusGen(-999.);

  }

  if(genVector.size()>1){
    std::sort(genVector.begin(), genVector.end(), orderETA());
    for(unsigned int i=0;i<genVector.size()-1;i++){
      GapVector.push_back(std::make_pair(fabs(genVector[i+1]->eta()-genVector[i]->eta()),genVector[i]->eta()));
    }
  }

  // GEN Find LRG
  if(GapVector.size()>0){
    if(debug) std::cout << "\nLooking Gap" << std::endl;
    std::sort(GapVector.rbegin(), GapVector.rend());
    for(unsigned int i=0;i<GapVector.size()-1;i++){
      if (debug) std::cout << "GapSize: " << GapVector[i].first << " | eta edge: " << GapVector[i].second << std::endl;
    }
    eventData.SetLrgGen(GapVector[0].first);
  }else{
    eventData.SetLrgGen(-999.);
  } 

  math::XYZTLorentzVector systemX(0.,0.,0.,0.);
  math::XYZTLorentzVector systemY(0.,0.,0.,0.);
  //double xi_mx = 0;

  double sumpTX = 0.;
  double sumpTY = 0.;
  if(genVector.size()>0 && GapVector.size()>0){
    std::sort(genVector.begin(), genVector.end(), orderETA());
    std::sort(GapVector.rbegin(), GapVector.rend());
    for(unsigned int i=0;i<genVector.size();i++){
      math::XYZTLorentzVector tmp(genVector[i]->px(), genVector[i]->py(), genVector[i]->pz(), genVector[i]->energy());
      if(genVector[i]->eta()>GapVector[0].second){
	systemX += tmp;
	sumpTX += genVector[i]->pt();
      }else{
	systemY +=tmp;
	sumpTY += genVector[i]->pt();
      }
    }

    //if (systemX.M2() > systemY.M2() && systemX.M2() > 0.) xi_mx = systemX.M2()/(4*beamEnergy_*beamEnergy_);
    //else if (systemY.M2() > systemX.M2() && systemY.M2() > 0.) xi_mx = systemY.M2()/(4*beamEnergy_*beamEnergy_);

    if(debug){
      std::cout << "\nLooking Diffractive System" << std::endl;
      std::cout << "Mx: " << systemX.M() << " [GeV]" << std::endl;
      std::cout << "My: " << systemY.M() << " [GeV]" << std::endl;
    }
    eventData.SetMxGenLeft(systemY.M());
    eventData.SetMxGenRight(systemX.M());
    eventData.SetMx2GenLeft(systemY.M2());
    eventData.SetMx2GenRight(systemX.M2());
    eventData.SetSumptGenLeft(sumpTY);
    eventData.SetSumptGenRight(sumpTX);
  }else{
    eventData.SetMxGenLeft(-999.);
    eventData.SetMxGenRight(-999.);
    eventData.SetMx2GenLeft(-999.);
    eventData.SetMx2GenRight(-999.);
    eventData.SetSumptGenLeft(-999.);
    eventData.SetSumptGenRight(-999.);
  }

  // GEN CMS LRG
  if(genCMSVector.size()>0){
    std::sort(genCMSVector.begin(), genCMSVector.end(), orderETA());
    for(unsigned int i=0;i<genCMSVector.size();i++){
      math::XYZTLorentzVector tmp(genCMSVector[i]->px(), genCMSVector[i]->py(), genCMSVector[i]->pz(), genCMSVector[i]->energy());
      E_pz_plus_gen_CMS += genCMSVector[i]->energy()+genCMSVector[i]->pz();
      et_expo_plus_gen_CMS += genCMSVector[i]->et()*pow(2.71,genCMSVector[i]->eta());
      E_pz_minus_gen_CMS += genCMSVector[i]->energy()-genCMSVector[i]->pz();
      et_expo_minus_gen_CMS += genCMSVector[i]->et()*pow(2.71,-genCMSVector[i]->eta());
      if(genCMSVector[i]->eta()>0.){
	SPlusCMS+=tmp;
      }else{
	SMinusCMS+=tmp;
      }
    }

    if(debug){
      if(genVector.size()>0){
	std::cout << "\nGen CMS Particles After Sort by Eta Info:" << std::endl;
	for(unsigned int i=0;i<genCMSVector.size();i++){
	  std::cout << "Gen. Particles CMS, eta: " << genCMSVector[i]->eta() << std::endl;
	}
      }
      std::cout << "\nReconstructed Variables CMS Gen Info:" << std::endl;
      std::cout << "Eta, CMS Max: " << genCMSVector[0]->eta() << std::endl;
      std::cout << "Eta, CMS Min: " << genCMSVector[genCMSVector.size()-1]->eta() << std::endl;
      std::cout << "PF CMS Mass -: " << SMinusCMS.M() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass +: " << SPlusCMS.M() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass^2 -: " << SMinusCMS.M2() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass^2 +: " << SPlusCMS.M2() << " [GeV]" << std::endl;
      std::cout << "E + pz, CMS: " << E_pz_plus_gen_CMS << " [GeV]" << std::endl;
      std::cout << "E - pz, CMS: " << E_pz_minus_gen_CMS << " [GeV]" << std::endl;
      std::cout << "Et*exp(+eta), CMS: " << et_expo_plus_gen_CMS << " [GeV]" << std::endl;
      std::cout << "Et*exp(-eta), CMS: " << et_expo_minus_gen_CMS << " [GeV]" << std::endl;
    }
    eventData.SetEtaMaxGenCMS(genCMSVector[0]->eta());
    eventData.SetEtaMinGenCMS(genCMSVector[genCMSVector.size()-1]->eta());
    eventData.SetMxGenMinusCMS(SMinusCMS.M());
    eventData.SetMxGenPlusCMS(SPlusCMS.M());
    eventData.SetMx2GenMinusCMS(SMinusCMS.M2());
    eventData.SetMx2GenPlusCMS(SPlusCMS.M2());
    eventData.SetEpluspzGenCMS(E_pz_plus_gen_CMS);
    eventData.SetEminuspzGenCMS(E_pz_minus_gen_CMS);
    eventData.SetEtExpoPlusGenCMS(et_expo_plus_gen_CMS);
    eventData.SetEtExpoMinusGenCMS(et_expo_minus_gen_CMS);
  }else{
    eventData.SetEtaMaxGenCMS(-999.);
    eventData.SetEtaMinGenCMS(-999.);
    eventData.SetMxGenMinusCMS(-999.);
    eventData.SetMxGenPlusCMS(-999.);
    eventData.SetMx2GenMinusCMS(-999.);
    eventData.SetMx2GenPlusCMS(-999.);
    eventData.SetEpluspzGenCMS(-999.);
    eventData.SetEminuspzGenCMS(-999.);
    eventData.SetEtExpoPlusGenCMS(-999.);
    eventData.SetEtExpoMinusGenCMS(-999.);
  }

  if(genCMSVector.size()>1){
    std::sort(genCMSVector.begin(), genCMSVector.end(), orderETA());
    for(unsigned int i=0;i<genCMSVector.size()-1;i++){
      GapCMSVector.push_back(std::make_pair(fabs(genCMSVector[i+1]->eta()-genCMSVector[i]->eta()),genCMSVector[i]->eta()));
    }
  }

  // GEN CMS Find LRG
  if(GapCMSVector.size()>0){
    if (debug) std::cout << "\nLooking Gap CMS" << std::endl;
    std::sort(GapCMSVector.rbegin(), GapCMSVector.rend());
    for(unsigned int i=0;i<GapCMSVector.size()-1;i++){
      if (debug) std::cout << "GapSize CMS: " << GapCMSVector[i].first << " | eta edge: " << GapCMSVector[i].second << std::endl;
    }
    eventData.SetLrgGenCMS(GapCMSVector[0].first);
  }else{
    eventData.SetLrgGenCMS(-999.);
  }

  math::XYZTLorentzVector systemXCMS(0.,0.,0.,0.);
  math::XYZTLorentzVector systemYCMS(0.,0.,0.,0.);
  //double xi_mx_CMS = 0;


  double sumpTXCMS = 0.;
  double sumpTYCMS = 0.;
  if(genCMSVector.size()>0 && GapCMSVector.size()>0){
    std::sort(genCMSVector.begin(), genCMSVector.end(), orderETA());
    std::sort(GapCMSVector.rbegin(), GapCMSVector.rend());
    for(unsigned int i=0;i<genCMSVector.size();i++){
      math::XYZTLorentzVector tmp(genCMSVector[i]->px(), genCMSVector[i]->py(), genCMSVector[i]->pz(), genCMSVector[i]->energy());
      if(genCMSVector[i]->eta()>GapCMSVector[0].second){
	systemXCMS += tmp;
	sumpTXCMS += genCMSVector[i]->pt();
      }else{
	systemYCMS +=tmp;
	sumpTYCMS += genCMSVector[i]->pt();
      }
    }

    //if (systemXCMS.M2() > systemYCMS.M2() && systemXCMS.M2() > 0.) xi_mx_CMS = systemXCMS.M2()/(4*beamEnergy_*beamEnergy_);
    //else if (systemYCMS.M2() > systemXCMS.M2() && systemYCMS.M2() > 0.) xi_mx_CMS = systemYCMS.M2()/(4*beamEnergy_*beamEnergy_);

    if(debug){
      std::cout << "\nLooking CMS Diffractive System" << std::endl;
      std::cout << "Mx: " << systemXCMS.M() << " [GeV]" << std::endl;
      std::cout << "My: " << systemYCMS.M() << " [GeV]" << std::endl;
    }
    eventData.SetMxGenLeftCMS(systemYCMS.M());
    eventData.SetMxGenRightCMS(systemXCMS.M());
    eventData.SetMx2GenLeftCMS(systemYCMS.M2());
    eventData.SetMx2GenRightCMS(systemXCMS.M2());
    eventData.SetSumptGenLeftCMS(sumpTYCMS);
    eventData.SetSumptGenRightCMS(sumpTXCMS);
  }else{
    eventData.SetMxGenLeftCMS(-999.);
    eventData.SetMxGenRightCMS(-999.);
    eventData.SetMx2GenLeftCMS(-999.);
    eventData.SetMx2GenRightCMS(-999.);
    eventData.SetSumptGenLeftCMS(-999.);
    eventData.SetSumptGenRightCMS(-999.);
  }

  if(debug){
    std::cout << "\nCastor Energy" << std::endl;
    std::cout << "All particles: " << sumECastorGen << " [GeV]" << std::endl;
    std::cout << "CMS: " << sumECastorGenCMS << " [GeV]" << std::endl;
    std::cout << "" << std::endl;
  }

  eventData.SetsumECastorMinusGen(sumECastorGen);
  eventData.SetsumECastorMinusGenCMS(sumECastorGenCMS);

}

//
// Fill Detector Variables
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillDetectorVariables(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  std::vector<double> energy_tower;
  std::vector<double> eta_tower;

  energy_tower.clear();
  eta_tower.clear();
  towerVector.clear();

  double Epz_plus=0.;  
  double Epz_minus=0.;  
  double Et_eta_plus=0.;
  double Et_eta_minus=0.;

  int nTowersHF_plus = 0;
  int nTowersHF_minus = 0;
  int nTowersHE_plus = 0;
  int nTowersHE_minus = 0;
  int nTowersHB_plus = 0;
  int nTowersHB_minus = 0;
  int nTowersEE_plus = 0;
  int nTowersEE_minus = 0;
  int nTowersEB_plus = 0;
  int nTowersEB_minus = 0;

  //Sum(E)
  double sumEHF_plus = 0.;
  double sumEHF_minus = 0.;
  double sumEHF_L_plus = 0.;
  double sumEHF_L_minus = 0.;
  double sumEHF_S_plus = 0.;
  double sumEHF_S_minus = 0.;
  double sumEHE_plus = 0.;
  double sumEHE_minus = 0.;
  double sumEHB_plus = 0.;
  double sumEHB_minus = 0.;
  double sumEEE_plus = 0.;
  double sumEEE_minus = 0.;
  double sumEEB_plus = 0.;
  double sumEEB_minus = 0.;

  // Sum(ET)
  double sumETHF_plus = 0.;
  double sumETHF_minus = 0.;
  double sumETHE_plus = 0.;
  double sumETHE_minus = 0.;
  double sumETHB_plus = 0.;
  double sumETHB_minus = 0.;
  double sumETEB_plus = 0.;
  double sumETEB_minus = 0.;
  double sumETEE_plus = 0.;
  double sumETEE_minus = 0.;

  edm::Handle<CaloTowerCollection> towerCollectionH;
  event.getByLabel(caloTowerTag_,towerCollectionH);

  int towersize = towerCollectionH->size();
  int itTower;
  int counter_tower=0;

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
    double caloTowerPz = calotower->pz();
    double caloTowerEt = calotower->et();
    double caloTowerEmEt = calotower->emEt();
    double caloTowerHadEt = calotower->hadEt();
    double EHF_S = 0;
    double EHF_L = 0;

    bool CalAboveTh = false;

    if( hasHF && !hasHE )
    {
      if( caloTowerEnergy > energyThresholdHF_ && fabs(calotower->eta())> 2.98 )   //// excluding HF ring1
      {
	CalAboveTh = true;

	if (debug) std::cout << "HF>> " <<  calotower->id() << " HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl; 

	// Added Long and short fibers
	// emc=L-S
	// hac=2*S
	// Tot = L+S
	// S = hac/2
	// L = Tot - S

	EHF_S = caloTowerHadEnergy/2;
	EHF_L = caloTowerEnergy - caloTowerHadEnergy/2;

	if(zside >= 0)
	{
	  ++nTowersHF_plus;
	  sumEHF_plus += caloTowerEnergy;
	  sumEHF_S_plus += EHF_S;
	  sumEHF_L_plus += EHF_L;
	  sumETHF_plus += caloTowerEt;
	}
	else
	{
	  ++nTowersHF_minus;
	  sumEHF_minus += caloTowerEnergy;
	  sumEHF_S_minus += EHF_S;
	  sumEHF_L_minus += EHF_L;
	  sumETHF_minus += caloTowerEt;
	}
      }
    }
    else if( hasHE && !hasHF && !hasHB )
    {
      if( caloTowerHadEnergy > energyThresholdHE_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "HE>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  ++nTowersHE_plus;
	  sumEHE_plus += caloTowerHadEnergy;
	  sumETHE_plus += caloTowerHadEt;
	}
	else
	{
	  ++nTowersHE_minus;
	  sumEHE_minus += caloTowerHadEnergy;
	  sumETHE_minus += caloTowerHadEt;
	}
      }
    }
    else if( hasHB && !hasHE )
    {
      if( caloTowerHadEnergy > energyThresholdHB_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "HB>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  ++nTowersHB_plus;
	  sumEHB_plus += caloTowerHadEnergy;
	  sumETHB_plus += caloTowerHadEt;
	}
	else
	{
	  ++nTowersHB_minus;
	  sumEHB_minus += caloTowerHadEnergy;
	  sumETHB_minus += caloTowerHadEt;
	}
      }
    }

    if( hasEE && !hasEB )
    {
      if( caloTowerEmEnergy >= energyThresholdEE_)
      {

	CalAboveTh = true;

	if (debug) std::cout << "EE>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  ++nTowersEE_plus;
	  sumEEE_plus += caloTowerEmEnergy;
	  sumETEE_plus += caloTowerEmEt;
	}
	else
	{
	  ++nTowersEE_minus;
	  sumEEE_minus += caloTowerEmEnergy;
	  sumETEE_minus += caloTowerEmEt;
	}
      }
    }
    else if( hasEB && !hasEE )
    {
      if( caloTowerEmEnergy >= energyThresholdEB_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "EB>> " <<  calotower->id() << " HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl; 

	if(zside >= 0)
	{
	  ++nTowersEB_plus;
	  sumEEB_plus += caloTowerEmEnergy;
	  sumETEB_plus += caloTowerEmEt;
	}
	else
	{
	  ++nTowersEB_minus;
	  sumEEB_minus += caloTowerEmEnergy;
	  sumETEB_minus += caloTowerEmEt;
	}
      }
    }

    if(CalAboveTh)
    {
      ++counter_tower;
      energy_tower.push_back(calotower->energy());
      eta_tower.push_back(calotower->eta());
      towerVector.push_back(calotower);
      Et_eta_minus += caloTowerEt * pow(2.71,-calotower->eta());
      Et_eta_plus += caloTowerEt * pow(2.71,calotower->eta());
      Epz_plus  += caloTowerEnergy + caloTowerPz;
      Epz_minus += caloTowerEnergy - caloTowerPz;
    }

  }  ////has to close calotower loop

  // Towers, Calorimeters, ordered
  std::sort(towerVector.begin(), towerVector.end(), orderETA()); 
  if(towerVector.size()>0){
    eventData.SetEtaCaloMax(towerVector[0]->eta());
    eventData.SetEtaCaloMin(towerVector[towerVector.size()-1]->eta());
    if(debug){
      std::cout << "\nTowers Info: " << std::endl;
      for(unsigned int i=0;i<towerVector.size();i++){
	std::cout << "Tower Eta: " << towerVector[i]->eta() << " | Phi: " << towerVector[i]->phi() << std::endl;
      }
      std::cout << "Calorimeter, Eta Max: " << towerVector[0]->eta() << " | Eta Min: " << towerVector[towerVector.size()-1]->eta() << std::endl;
      std::cout << "Calorimeter, Et*exp(eta): " << Et_eta_plus << std::endl;
      std::cout << "Calorimeter, Et*exp(-eta): " << Et_eta_minus << std::endl;
      std::cout << "Calorimeter, E + pz: " << Epz_plus << std::endl;
      std::cout << "Calorimeter, E - pz: " << Epz_minus << std::endl;
      std::cout << "Calorimeter, Sum HF+: " << sumEHF_plus << " [GeV]" << std::endl;
      std::cout << "Calorimeter, Sum HF-: " << sumEHF_minus << " [GeV]"<< std::endl;
    }
  }else{
    eventData.SetEtaCaloMax(999.);
    eventData.SetEtaCaloMin(-999.);
  }

  // Find LRG CaloTower
  std::vector<std::pair<double, double> > GapCaloVector;
  GapCaloVector.clear();

  if(towerVector.size()>1){
    std::sort(towerVector.begin(), towerVector.end(), orderETA());
    for(unsigned int i=0;i<towerVector.size()-1;i++){
      GapCaloVector.push_back(std::make_pair(fabs(towerVector[i+1]->eta()-towerVector[i]->eta()),towerVector[i]->eta()));
    }
  }

  if(GapCaloVector.size()>0){
    std::sort(GapCaloVector.rbegin(), GapCaloVector.rend());
    eventData.SetLrgCalo(GapCaloVector[0].first);
    if (debug) std::cout << "\nLooking Calorimeter Gap CMS" << std::endl;
    for(unsigned int i=0;i<GapCaloVector.size()-1;i++){
      if (debug) std::cout << "Calorimeter GapSize CMS: " << GapCaloVector[i].first << " | eta edge: " << GapCaloVector[i].second << std::endl;
    }
  }else{
    eventData.SetLrgCalo(-999.);
  }

  //Fill Variables
  eventData.SetSumEHFPlus(sumEHF_plus);
  eventData.SetSumEHF_SPlus(sumEHF_S_plus);
  eventData.SetSumEHF_LPlus(sumEHF_L_plus);
  eventData.SetSumEtHFPlus(sumETHF_plus);

  eventData.SetSumEHFMinus(sumEHF_minus);
  eventData.SetSumEHF_SMinus(sumEHF_S_minus);
  eventData.SetSumEHF_LMinus(sumEHF_L_minus);
  eventData.SetSumEtHFMinus(sumETHF_minus);

  eventData.SetSumEHEPlus(sumEHE_plus);
  eventData.SetSumEtHEPlus(sumETHE_plus);
  eventData.SetSumEHEMinus(sumEHE_minus);
  eventData.SetSumEtHEMinus(sumETHE_minus);

  eventData.SetSumEHBPlus(sumEHB_plus);
  eventData.SetSumEtHBPlus(sumETHB_plus);
  eventData.SetSumEHBMinus(sumEHB_minus);
  eventData.SetSumEtHBMinus(sumETHB_minus);

  eventData.SetSumEEEPlus(sumEEE_plus);
  eventData.SetSumEtEEPlus(sumETEE_plus);
  eventData.SetSumEEEMinus(sumEEE_minus);
  eventData.SetSumEtEEMinus(sumETEE_minus);

  eventData.SetSumEEBPlus(sumEEB_plus);
  eventData.SetSumEtEBPlus(sumETEB_plus);
  eventData.SetSumEEBMinus(sumEEB_minus);
  eventData.SetSumEtEBMinus(sumETEB_minus);

  eventData.SetEpluspzCalo(Epz_plus);
  eventData.SetEminuspzCalo(Epz_minus);
  eventData.SetEtExpoPlusCalo(Et_eta_plus);
  eventData.SetEtExpoMinusCalo(Et_eta_minus);

  eventData.SetMultiplicityHFPlus(nTowersHF_plus);
  eventData.SetMultiplicityHEPlus(nTowersHE_plus);
  eventData.SetMultiplicityEEPlus(nTowersEE_plus);
  eventData.SetMultiplicityHFMinus(nTowersHF_minus);
  eventData.SetMultiplicityHEMinus(nTowersHE_minus);
  eventData.SetMultiplicityEEMinus(nTowersEE_minus);
  if(EachTower_){
    eventData.SetEachTowerCounter(counter_tower);
    eventData.SetEachTowerEta(eta_tower);
    eventData.SetEachTowerEnergy(energy_tower);
  }

}

//
// Fill Physics Variables
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillVariables(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug=false;
  bool debugDeep = false;

  PFMuonVector.clear();
  PFElectronVector.clear();
  PFVector.clear();
  PFNoZVector.clear();
  PFHFVector.clear();

  edm::Handle<reco::VertexCollection> Vertexes;
  event.getByLabel(PVtxCollectionTag_, Vertexes);

  eventData.SetVertex(Vertexes->size());

  edm::Handle <reco::PFCandidateCollection> PFCandidates;
  event.getByLabel(pfTag_,PFCandidates);
  int pfsize = PFCandidates->size();
  int itPF;

  double E_pz_plus = 0.;
  double E_pz_minus = 0.;
  double et_expo_plus = 0.;
  double et_expo_minus = 0.;

  math::XYZTLorentzVector SPlusCMS(0.,0.,0.,0.);
  math::XYZTLorentzVector SMinusCMS(0.,0.,0.,0.);

  if(PFCandidates->size()>0){
    for(itPF=0; itPF < pfsize; ++itPF){
      const reco::PFCandidate* pfAll = &((*PFCandidates)[itPF]);
      // Excluding HF Calorimeter Rings 12, 13, 29, 30, 40, 41
      if( ( (fabs(pfAll->eta()) >= 2.866) && (fabs(pfAll->eta()) < 3.152) ) || (fabs(pfAll->eta()) >= 4.730) ) continue;

      if ( (fabs(pfAll->charge()) > 0 && pfAll->pt() > pTPFThresholdCharged_ ) || (fabs(pfAll->charge()) == 0 && ( (fabs(pfAll->eta()) <= 1.5 && pfAll->energy() > energyPFThresholdBar_) || 
	      (fabs(pfAll->eta()) > 1.5 && fabs(pfAll->eta()) <= 3 && pfAll->energy() > energyPFThresholdEnd_) || (fabs(pfAll->eta()) > 3 && pfAll->energy() >energyPFThresholdHF_) ) ) ){

	PFVector.push_back(pfAll);
	if (pfAll->particleId()==reco::PFCandidate::e) PFElectronVector.push_back(pfAll);
	if (pfAll->particleId()==reco::PFCandidate::mu) PFMuonVector.push_back(pfAll);
	if (pfAll->particleId()==reco::PFCandidate::h_HF || pfAll->particleId()==reco::PFCandidate::egamma_HF) PFHFVector.push_back(pfAll);

      }

    }
  }

  // Computing PF Variables
  if(PFVector.size()>0){
    for(unsigned int i=0;i<PFVector.size();i++){
      math::XYZTLorentzVector tmp(PFVector[i]->px(), PFVector[i]->py(), PFVector[i]->pz(), PFVector[i]->energy());
      E_pz_plus += PFVector[i]->energy()+PFVector[i]->pz();
      et_expo_plus += PFVector[i]->et()*pow(2.71,PFVector[i]->eta());
      E_pz_minus += PFVector[i]->energy()-PFVector[i]->pz();
      et_expo_minus += PFVector[i]->et()*pow(2.71,-PFVector[i]->eta());
      if(PFVector[i]->eta()>0.){
	SPlusCMS+=tmp;
      }else{
	SMinusCMS+=tmp;
      }
    }
    eventData.SetMxPFMinus(SMinusCMS.M());
    eventData.SetMxPFPlus(SPlusCMS.M());
    eventData.SetMx2PFMinus(SMinusCMS.M2());
    eventData.SetMx2PFPlus(SPlusCMS.M2());
    eventData.SetEpluspzPF(E_pz_plus);
    eventData.SetEminuspzPF(E_pz_minus);
    eventData.SetEtExpoPlusPF(et_expo_plus);
    eventData.SetEtExpoMinusPF(et_expo_minus);
  }else{
    eventData.SetMxPFMinus(-999.);
    eventData.SetMxPFPlus(-999.);
    eventData.SetMx2PFMinus(-999.);
    eventData.SetMx2PFPlus(-999.);
    eventData.SetEpluspzPF(-999.);
    eventData.SetEminuspzPF(-999.);
    eventData.SetEtExpoPlusPF(-999.);
    eventData.SetEtExpoMinusPF(-999.);
  }

  if(PFVector.size()>0){
    std::sort(PFVector.begin(), PFVector.end(), orderETA());
    eventData.SetEtaMaxPF(PFVector[0]->eta());
    eventData.SetEtaMinPF(PFVector[PFVector.size()-1]->eta());
    if(debug){
      std::cout << "\nCMS Particles After Sort by Eta Info:" << std::endl;
      std::cout << "Reconstructed Variables CMS PF Info:" << std::endl;
      std::cout << "Eta, CMS Max: " << PFVector[0]->eta() << std::endl;
      std::cout << "Eta, CMS Min: " << PFVector[PFVector.size()-1]->eta() << std::endl;
      std::cout << "PF CMS Mass -: " << SMinusCMS.M() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass +: " << SPlusCMS.M() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass^2 -: " << SMinusCMS.M2() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass^2 +: " << SPlusCMS.M2() << " [GeV]" << std::endl;
      std::cout << "E + pz, CMS: " << E_pz_plus << " [GeV]" << std::endl;
      std::cout << "E - pz, CMS: " << E_pz_minus << " [GeV]" << std::endl;
      std::cout << "Et*exp(+eta), CMS: " << et_expo_plus << " [GeV]" << std::endl;
      std::cout << "Et*exp(-eta), CMS: " << et_expo_minus << " [GeV]" << std::endl;
      if(debugDeep){
	for(unsigned int i=0;i<PFVector.size();i++){
	  std::cout << "All PF Information --> pT: " << PFVector[i]->pt() << " [GeV] | eta: " << PFVector[i]->eta() << " | phi: " << PFVector[i]->phi() << " | id: " << PFVector[i]->particleId() << " | # Particles: " << PFVector.size() << std::endl;
	}
      }
    }
  }else{
    eventData.SetEtaMaxPF(-999.);
    eventData.SetEtaMinPF(-999.);
  }

  if(PFElectronVector.size()>0){
    std::sort(PFElectronVector.begin(), PFElectronVector.end(), orderPT());
    if(debugDeep){
      for(unsigned int i=0;i<PFElectronVector.size();i++){
	std::cout << "Electron PF Information --> pT: " << PFElectronVector[i]->pt() << " [GeV] | eta: " << PFElectronVector[i]->eta() << " | phi: " << PFElectronVector[i]->phi() << std::endl;
      }
    }
  }

  if(PFMuonVector.size()>0){
    std::sort(PFMuonVector.begin(), PFMuonVector.end(), orderPT());
    if(debugDeep){
      for(unsigned int i=0;i<PFMuonVector.size();i++){
	std::cout << "Muon PF Information --> pT: " << PFMuonVector[i]->pt() << " [GeV] | eta: " << PFMuonVector[i]->eta() << " | phi: " << PFMuonVector[i]->phi() << std::endl;
      }
    }
  }

  double PFHFEnergyPlus = 0.;
  double PFHFEnergyMinus = 0.;

  if(PFHFVector.size()>0){
    std::sort(PFHFVector.begin(), PFHFVector.end(), orderETA());
    for(unsigned int i=0;i<PFHFVector.size();i++){
      if(PFHFVector[i]->eta()>0){
	PFHFEnergyPlus += PFHFVector[i]->energy();
      }else{
	PFHFEnergyMinus += PFHFVector[i]->energy();
      }
      if (debugDeep) {
	std::cout << "PF HF Information --> pT: " << PFHFVector[i]->pt() << " [GeV] | eta: " << PFHFVector[i]->eta() << " | phi: " << PFHFVector[i]->phi() << " | id: " << PFVector[i]->particleId() << std::endl;
      }
    }
    if(debug){
      std::cout << "HF+ Energy from PF: " << PFHFEnergyPlus  << " [GeV]" << std::endl;
      std::cout << "HF- Energy from PF: " << PFHFEnergyMinus  << " [GeV]" << std::endl;
    }
  }

  bool ZMuon=false;
  bool ZElectron=false;
  bool ZPFMuon=false;
  bool ZPFElectron=false;

  // Z Candidates from Electron, Muon and PF Colections
  if(PFMuonVector.size()>1){
    if(InvariantMass(PFMuonVector[0],PFMuonVector[1])>60. && InvariantMass(PFMuonVector[0],PFMuonVector[1])<110.) ZPFMuon=true;
    if(debug) std::cout << "Invariant Z Mass, dimuon PF: " << InvariantMass(PFMuonVector[0],PFMuonVector[1]) << " [GeV]" << std::endl;
  }

  if(PFElectronVector.size()>1){
    if(InvariantMass(PFElectronVector[0],PFElectronVector[1])>60. && InvariantMass(PFElectronVector[0],PFElectronVector[1])<110.) ZPFElectron=true;
    if(debug) std::cout << "Invariant Z Mass, dielectron PF: " << InvariantMass(PFElectronVector[0],PFElectronVector[1]) << " [GeV]" << std::endl;
  }

  if(MuonVector.size()>1){
    if(InvariantMass(MuonVector[0],MuonVector[1])>60. && InvariantMass(MuonVector[0],MuonVector[1])<110.) ZMuon=true;
    if(debug) std::cout << "Invariant Z Mass, dimuon: " << InvariantMass(MuonVector[0],MuonVector[1]) << " [GeV]" << std::endl;
  }

  if(ElectronVector.size()>1){
    if(InvariantMass(ElectronVector[0],ElectronVector[1])>60. && InvariantMass(ElectronVector[0],ElectronVector[1])<110.) ZElectron=true;
    if(debug) std::cout << "Invariant Z Mass, dielectron: " << InvariantMass(ElectronVector[0],ElectronVector[1]) << " [GeV]" << std::endl;
  }

  // Fill PF Vector excluing Z candidate decay products
  if(PFVector.size()>1 && (ZPFMuon || ZMuon)){
    for(unsigned int i=0;i<PFVector.size();i++){
      if( ZPFMuon && (PFVector[i]->pt() == PFMuonVector[0]->pt() || PFVector[i]->pt() == PFMuonVector[1]->pt()) ) continue;
      if( ZMuon && (PFVector[i]->pt() == MuonVector[0]->pt() || PFVector[i]->pt() == MuonVector[1]->pt()) ) continue;
      PFNoZVector.push_back(PFVector[i]);
    }
  }
  else if(PFVector.size()>1 && (ZPFElectron || ZElectron)){
    for(unsigned int i=0;i<PFVector.size();i++){
      if( ZPFElectron && (PFVector[i]->pt() == PFElectronVector[0]->pt() || PFVector[i]->pt() == PFElectronVector[1]->pt()) ) continue;
      if( ZElectron && (PFVector[i]->pt() == ElectronVector[0]->pt() || PFVector[i]->pt() == ElectronVector[1]->pt()) ) continue;
      PFNoZVector.push_back(PFVector[i]);
    }
  }else{
    if(PFVector.size()>1){
      for(unsigned int i=0;i<PFVector.size();i++){
	PFNoZVector.push_back(PFVector[i]);
      }
    }
  }

  // PF Without Z products
  double E_pz_plus_noz = 0.;
  double E_pz_minus_noz = 0.;
  double et_expo_plus_noz = 0.;
  double et_expo_minus_noz = 0.;

  math::XYZTLorentzVector SPlusCMSnoz(0.,0.,0.,0.);
  math::XYZTLorentzVector SMinusCMSnoz(0.,0.,0.,0.);

  // Computing PF Variables
  if(PFNoZVector.size()>0){
    for(unsigned int i=0;i<PFNoZVector.size();i++){
      math::XYZTLorentzVector tmp(PFNoZVector[i]->px(), PFNoZVector[i]->py(), PFNoZVector[i]->pz(), PFNoZVector[i]->energy());
      E_pz_plus_noz += PFNoZVector[i]->energy()+PFNoZVector[i]->pz();
      et_expo_plus_noz += PFNoZVector[i]->et()*pow(2.71,PFNoZVector[i]->eta());
      E_pz_minus_noz += PFNoZVector[i]->energy()-PFNoZVector[i]->pz();
      et_expo_minus_noz += PFNoZVector[i]->et()*pow(2.71,-PFNoZVector[i]->eta());
      if(PFNoZVector[i]->eta()>0.){
	SPlusCMSnoz+=tmp;
      }else{
	SMinusCMSnoz+=tmp;
      }
    }
    eventData.SetMxPFNoZMinus(SMinusCMSnoz.M());
    eventData.SetMxPFNoZPlus(SPlusCMSnoz.M());
    eventData.SetMx2PFNoZMinus(SMinusCMSnoz.M2());
    eventData.SetMx2PFNoZPlus(SPlusCMSnoz.M2());
    eventData.SetEpluspzPFNoZ(E_pz_plus_noz);
    eventData.SetEminuspzPFNoZ(E_pz_minus_noz);
    eventData.SetEtExpoPlusPFNoZ(et_expo_plus_noz);
    eventData.SetEtExpoMinusPFNoZ(et_expo_minus_noz);
  }else{
    eventData.SetMxPFNoZMinus(-999.);
    eventData.SetMxPFNoZPlus(-999.);
    eventData.SetMx2PFNoZMinus(-999.);
    eventData.SetMx2PFNoZPlus(-999.);
    eventData.SetEpluspzPFNoZ(-999.);
    eventData.SetEminuspzPFNoZ(-999.);
    eventData.SetEtExpoPlusPFNoZ(-999.);
    eventData.SetEtExpoMinusPFNoZ(-999.);
  }

  if(PFNoZVector.size()>0){
    std::sort(PFNoZVector.begin(), PFNoZVector.end(), orderETA());
    eventData.SetEtaMaxPFNoZ(PFNoZVector[0]->eta());
    eventData.SetEtaMinPFNoZ(PFNoZVector[PFNoZVector.size()-1]->eta());
    if(debug){
      std::cout << "\nCMS Particles After Sort by Eta Info, excluding Z Boson products:" << std::endl;
      std::cout << "Reconstructed Variables CMS PF Info:" << std::endl;
      std::cout << "Eta, CMS Max: " << PFNoZVector[0]->eta() << std::endl;
      std::cout << "Eta, CMS Min: " << PFNoZVector[PFNoZVector.size()-1]->eta() << std::endl;
      std::cout << "PF CMS Mass -: " << SMinusCMSnoz.M() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass +: " << SPlusCMSnoz.M() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass^2 -: " << SMinusCMSnoz.M2() << " [GeV]" << std::endl;
      std::cout << "PF CMS Mass^2 +: " << SPlusCMSnoz.M2() << " [GeV]" << std::endl;
      std::cout << "E + pz, CMS: " << E_pz_plus_noz << " [GeV]" << std::endl;
      std::cout << "E - pz, CMS: " << E_pz_minus_noz << " [GeV]" << std::endl;
      std::cout << "Et*exp(+eta), CMS: " << et_expo_plus_noz << " [GeV]" << std::endl;
      std::cout << "Et*exp(-eta), CMS: " << et_expo_minus_noz << " [GeV]" << std::endl;
      if(debugDeep){
	if(PFNoZVector.size()>0){
	  for(unsigned int i=0;i<PFNoZVector.size();i++){
	    std::cout << "All PF Information, no Z products --> pT: " << PFNoZVector[i]->pt() << " [GeV] | eta: " << PFNoZVector[i]->eta() << " | phi: " << PFNoZVector[i]->phi() << " | id: " << PFNoZVector[i]->particleId() << " | # Particles: " << PFNoZVector.size() << std::endl;
	  }
	}
      }
    }
  }else{
    eventData.SetEtaMaxPFNoZ(-999.);
    eventData.SetEtaMinPFNoZ(-999.);
  }

  std::vector<std::pair<double, double> > GapPFVector;
  GapPFVector.clear();

  std::vector<std::pair<double, double> > GapPFVectornoz;
  GapPFVectornoz.clear();

  // Find LRG
  if(PFVector.size()>1){
    std::sort(PFVector.begin(), PFVector.end(), orderETA());
    for(unsigned int i=0;i<PFVector.size()-1;i++){
      GapPFVector.push_back(std::make_pair(fabs(PFVector[i+1]->eta()-PFVector[i]->eta()),PFVector[i]->eta()));
    }
  }

  if(GapPFVector.size()>0){
    std::sort(GapPFVector.rbegin(), GapPFVector.rend());
    eventData.SetLrgPF(GapPFVector[0].first);
    if (debugDeep) std::cout << "\nLooking Gap CMS" << std::endl;
    for(unsigned int i=0;i<GapPFVector.size()-1;i++){
      if (debugDeep) std::cout << "GapSize CMS: " << GapPFVector[i].first << " | eta edge: " << GapPFVector[i].second << std::endl;
    }
  }else{
    eventData.SetLrgPF(-999.);
  }

  // Find LRG No Z
  if(PFNoZVector.size()>1){
    std::sort(PFNoZVector.begin(), PFNoZVector.end(), orderETA());
    for(unsigned int i=0;i<PFNoZVector.size()-1;i++){
      GapPFVectornoz.push_back(std::make_pair(fabs(PFNoZVector[i+1]->eta()-PFNoZVector[i]->eta()),PFNoZVector[i]->eta()));
    }
  }

  if(GapPFVectornoz.size()>0){
    std::sort(GapPFVectornoz.rbegin(), GapPFVectornoz.rend());
    eventData.SetLrgPFNoZ(GapPFVectornoz[0].first);
    if (debugDeep) std::cout << "\nLooking Gap CMS, excluding Z decay products" << std::endl;
    for(unsigned int i=0;i<GapPFVectornoz.size()-1;i++){
      if (debugDeep) std::cout << "GapSize CMS: " << GapPFVectornoz[i].first << " | eta edge: " << GapPFVectornoz[i].second << std::endl;
    }
  }else{
    eventData.SetLrgPFNoZ(-999.);
  }

  // Calculate Sum pT and Mx left and right edge of the LRG
  math::XYZTLorentzVector systemXCMS(0.,0.,0.,0.);
  math::XYZTLorentzVector systemYCMS(0.,0.,0.,0.);
  double sumPTX = 0.;
  double sumPTY = 0.;

  if(PFVector.size()>0 && GapPFVector.size()>0){
    std::sort(PFVector.begin(), PFVector.end(), orderETA());
    std::sort(GapPFVector.rbegin(), GapPFVector.rend());
    for(unsigned int i=0;i<PFVector.size();i++){
      math::XYZTLorentzVector tmp(PFVector[i]->px(), PFVector[i]->py(), PFVector[i]->pz(), PFVector[i]->energy());
      if(PFVector[i]->eta()>GapPFVector[0].second){
	systemXCMS += tmp;
	sumPTX += PFVector[i]->pt();
      }else{
	systemYCMS +=tmp;
	sumPTY += PFVector[i]->pt();
      }
    }
    if(debug){
      std::cout << "\nLooking CMS Diffractive System" << std::endl;
      std::cout << "Mx: " << systemXCMS.M() << " [GeV]" << std::endl;
      std::cout << "My: " << systemYCMS.M() << " [GeV]" << std::endl;
      std::cout << "Sum pT X system: " << sumPTX << " [GeV]" << std::endl;
      std::cout << "Sum pT Y system: " << sumPTY << " [GeV]" << std::endl;
    }
    eventData.SetMxPFLeft(systemYCMS.M());
    eventData.SetMxPFRight(systemXCMS.M());
    eventData.SetMx2PFLeft(systemYCMS.M2());
    eventData.SetMx2PFRight(systemXCMS.M2());
    eventData.SetSumptPFLeft(sumPTY);
    eventData.SetSumptPFRight(sumPTX);
  }else{
    eventData.SetMxPFLeft(-999.);
    eventData.SetMxPFRight(-999.);
    eventData.SetMx2PFLeft(-999.);
    eventData.SetMx2PFRight(-999.);
    eventData.SetSumptPFLeft(-999.);
    eventData.SetSumptPFRight(-999.);
  }

  math::XYZTLorentzVector systemXCMSnoz(0.,0.,0.,0.);
  math::XYZTLorentzVector systemYCMSnoz(0.,0.,0.,0.);
  double sumPTXnoz = 0.;
  double sumPTYnoz = 0.;

  if(PFNoZVector.size()>0 && GapPFVectornoz.size()>0){
    std::sort(PFNoZVector.begin(), PFNoZVector.end(), orderETA());
    std::sort(GapPFVectornoz.rbegin(), GapPFVectornoz.rend());
    for(unsigned int i=0;i<PFNoZVector.size();i++){
      math::XYZTLorentzVector tmp(PFNoZVector[i]->px(), PFNoZVector[i]->py(), PFNoZVector[i]->pz(), PFNoZVector[i]->energy());
      if(PFNoZVector[i]->eta()>GapPFVectornoz[0].second){
	systemXCMSnoz += tmp;
	sumPTXnoz += PFNoZVector[i]->pt();
      }else{
	systemYCMSnoz +=tmp;
	sumPTYnoz += PFNoZVector[i]->pt();
      }
    }
    if(debug){
      std::cout << "\nLooking CMS Diffractive System, Excluding Z decay products" << std::endl;
      std::cout << "Mx: " << systemXCMSnoz.M() << " [GeV]" << std::endl;
      std::cout << "My: " << systemYCMSnoz.M() << " [GeV]" << std::endl;
      std::cout << "Sum pT X system: " << sumPTXnoz << " [GeV]" << std::endl;
      std::cout << "Sum pT Y system: " << sumPTYnoz << " [GeV]" << std::endl;
    }
    eventData.SetMxPFNoZLeft(systemYCMSnoz.M());
    eventData.SetMxPFNoZRight(systemXCMSnoz.M());
    eventData.SetMx2PFNoZLeft(systemYCMSnoz.M2());
    eventData.SetMx2PFNoZRight(systemXCMSnoz.M2());
    eventData.SetSumptPFNoZLeft(sumPTYnoz);
    eventData.SetSumptPFNoZRight(sumPTXnoz);
  }else{
    eventData.SetMxPFNoZLeft(-999.);
    eventData.SetMxPFNoZRight(-999.);
    eventData.SetMx2PFNoZLeft(-999.);
    eventData.SetMx2PFNoZRight(-999.);
    eventData.SetSumptPFNoZLeft(-999.);
    eventData.SetSumptPFNoZRight(-999.);
  }

}

//
// Fill Pat:Muon and Pat:Electron objects 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillZPat(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Declaring Variables
  bool debug = false;

  PatElectronVector.clear();
  PatMuonVector.clear();

  // Detector Objects and Candidates
  edm::Handle<std::vector<pat::Muon> > muons;
  event.getByLabel("patMuons", muons);

  edm::Handle<std::vector<pat::Electron> > electrons;
  event.getByLabel("patElectrons", electrons);

  edm::Handle<reco::VertexCollection>  vertex;
  event.getByLabel(PVtxCollectionTag_, vertex);

  int electronsize = electrons->size();
  int itElectron;
  if(electrons->size()>0){
    for(itElectron=0; itElectron < electronsize; ++itElectron){
      const pat::Electron* electronAll = &((*electrons)[itElectron]);
      PatElectronVector.push_back(electronAll);
    }
  }

  int muonsize = muons->size();
  int itMuon;
  if(muons->size()>0){
    for(itMuon=0; itMuon < muonsize; ++itMuon){
      const pat::Muon* muonAll = &((*muons)[itMuon]);
      PatMuonVector.push_back(muonAll);
    }
  }

  // Sorting Vector
  std::sort(PatElectronVector.begin(), PatElectronVector.end(), orderPT());
  std::sort(PatMuonVector.begin(), PatMuonVector.end(), orderPT());

  if(PatElectronVector.size()>0){
    for(unsigned int i=0;i<PatElectronVector.size();i++){
      if (debug) {
	std::cout << "PAT Electron --> pT: " << PatElectronVector[i]->pt() << " [GeV] | eta: " << PatElectronVector[i]->eta() << " | phi: " << PatElectronVector[i]->phi() << std::endl;
      }
    }
  }

  if(PatMuonVector.size()>0){
    for(unsigned int i=0;i<PatMuonVector.size();i++){
      if (debug) {
	std::cout << "PAT Muon --> pT: " << PatMuonVector[i]->pt() << " [GeV] | eta: " << PatMuonVector[i]->eta() << " | phi: " << PatMuonVector[i]->phi() << std::endl;
      }
    }
  }


  if(PatMuonVector.size()>1){
    double muon1SumPtR03 = PatMuonVector[0]->isolationR03().sumPt;
    double muon1EmEtR03 = PatMuonVector[0]->isolationR03().emEt;
    double muon1HadEtR03 = PatMuonVector[0]->isolationR03().hadEt;    
    double muon1SumPtR05 = PatMuonVector[0]->isolationR05().sumPt;
    double muon1EmEtR05 = PatMuonVector[0]->isolationR05().emEt;
    double muon1HadEtR05 = PatMuonVector[0]->isolationR05().hadEt;    
    double muon2SumPtR03 = PatMuonVector[1]->isolationR03().sumPt;
    double muon2EmEtR03 = PatMuonVector[1]->isolationR03().emEt;
    double muon2HadEtR03 = PatMuonVector[1]->isolationR03().hadEt;
    double muon2SumPtR05 = PatMuonVector[1]->isolationR05().sumPt;
    double muon2EmEtR05 = PatMuonVector[1]->isolationR05().emEt;
    double muon2HadEtR05 = PatMuonVector[1]->isolationR05().hadEt;
    double relIsoFirstMuonDr03 = (muon1SumPtR03 + muon1EmEtR03 + muon1HadEtR03)/PatMuonVector[0]->pt();
    double relIsoSecondMuonDr03 = (muon2SumPtR03 + muon2EmEtR03 + muon2HadEtR03)/PatMuonVector[1]->pt();
    double relIsoFirstMuonDr05 = (muon1SumPtR05 + muon1EmEtR05 + muon1HadEtR05)/PatMuonVector[0]->pt();
    double relIsoSecondMuonDr05 = (muon2SumPtR05 + muon2EmEtR05 + muon2HadEtR05)/PatMuonVector[1]->pt();

    double relIsoFirstMuon = (PatMuonVector[0]->trackIso()+PatMuonVector[0]->ecalIso()+PatMuonVector[0]->hcalIso())/PatMuonVector[0]->pt();
    double relIsoSecondMuon = (PatMuonVector[1]->trackIso()+PatMuonVector[1]->ecalIso()+PatMuonVector[1]->hcalIso())/PatMuonVector[1]->pt();

    // DiMuon Mass
    math::XYZTLorentzVector DipatMuonSystem(0.,0.,0.,0.);
    DipatMuonSystem += PatMuonVector[0]->p4();
    DipatMuonSystem += PatMuonVector[1]->p4();

    eventData.SetPatDiMuonMass(DipatMuonSystem.M());
    eventData.SetPatDiMuonEta(DipatMuonSystem.eta());
    eventData.SetPatDiMuonPhi(DipatMuonSystem.phi());
    eventData.SetPatDiMuonPt(DipatMuonSystem.pt());
    eventData.SetPatNMuon(PatMuonVector.size());
    eventData.SetPatMuon1Pt(PatMuonVector[0]->pt());
    eventData.SetPatMuon1Charge(PatMuonVector[0]->charge());
    eventData.SetPatMuon1Phi(PatMuonVector[0]->phi());
    eventData.SetPatMuon1Eta(PatMuonVector[0]->eta());
    eventData.SetPatMuon1Et(PatMuonVector[0]->et());
    eventData.SetPatMuon1P4(PatMuonVector[0]->p4());
    eventData.SetPatMuon2Pt(PatMuonVector[1]->pt());
    eventData.SetPatMuon2Charge(PatMuonVector[1]->charge());
    eventData.SetPatMuon2Phi(PatMuonVector[1]->phi());
    eventData.SetPatMuon2Eta(PatMuonVector[1]->eta());
    eventData.SetPatMuon2Et(PatMuonVector[1]->et());
    eventData.SetPatMuon2P4(PatMuonVector[1]->p4());
    eventData.SetPatMuon1SumPtR03(muon1SumPtR03);
    eventData.SetPatMuon1EmEtR03(muon1EmEtR03);
    eventData.SetPatMuon1HadEtR03(muon1HadEtR03);    
    eventData.SetPatMuon1SumPtR05(muon1SumPtR05);
    eventData.SetPatMuon1EmEtR05(muon1EmEtR05);
    eventData.SetPatMuon1HadEtR05(muon1HadEtR05);    
    eventData.SetPatMuon2SumPtR03(muon2SumPtR03);
    eventData.SetPatMuon2EmEtR03(muon2EmEtR03);
    eventData.SetPatMuon2HadEtR03(muon2HadEtR03);    
    eventData.SetPatMuon2SumPtR05(muon2SumPtR05);
    eventData.SetPatMuon2EmEtR05(muon2EmEtR05);
    eventData.SetPatMuon2HadEtR05(muon2HadEtR05);  
    eventData.SetPatMuon1relIsoDr03(relIsoFirstMuonDr03);
    eventData.SetPatMuon2relIsoDr03(relIsoSecondMuonDr03);
    eventData.SetPatMuon1relIsoDr05(relIsoFirstMuonDr05);
    eventData.SetPatMuon2relIsoDr05(relIsoSecondMuonDr05);
    eventData.SetPatMuon1relIso(relIsoFirstMuon);
    eventData.SetPatMuon2relIso(relIsoSecondMuon);
    eventData.SetPatMuon1TrackerHits(PatMuonVector[0]->track()->hitPattern().trackerLayersWithMeasurement());
    eventData.SetPatMuon1PixelHits(PatMuonVector[0]->innerTrack()->hitPattern().numberOfValidPixelHits());
    eventData.SetPatMuon1NormalizedChi2(PatMuonVector[0]->globalTrack()->normalizedChi2());
    eventData.SetPatMuon1MatchedStations(PatMuonVector[0]->numberOfMatchedStations());
    eventData.SetPatMuon1Dxy(PatMuonVector[0]->innerTrack()->dxy(vertex->at(0).position()));
    eventData.SetPatMuon1IsGlobal(PatMuonVector[0]->isGlobalMuon());
    eventData.SetPatMuon1IsTracker(PatMuonVector[0]->isTrackerMuon());
    eventData.SetPatMuon2TrackerHits(PatMuonVector[1]->track()->hitPattern().trackerLayersWithMeasurement());
    eventData.SetPatMuon2PixelHits(PatMuonVector[1]->innerTrack()->hitPattern().numberOfValidPixelHits());
    eventData.SetPatMuon2NormalizedChi2(PatMuonVector[1]->globalTrack()->normalizedChi2());
    eventData.SetPatMuon2MatchedStations(PatMuonVector[1]->numberOfMatchedStations());
    eventData.SetPatMuon2Dxy(PatMuonVector[1]->innerTrack()->dxy(vertex->at(0).position()));
    eventData.SetPatMuon2IsGlobal(PatMuonVector[1]->isGlobalMuon());
    eventData.SetPatMuon2IsTracker(PatMuonVector[1]->isTrackerMuon());

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCountm03= 0;
    int goodTracksCountm04= 0;
    int goodTracksCountm05= 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),PatMuonVector[0]->eta(),PatMuonVector[0]->phi()) > 0.3) && (deltaR(track->eta(),track->phi(),PatMuonVector[1]->eta(),PatMuonVector[1]->phi()) > 0.3))
      {
	goodTracksCountm03++;
      }

      if ((deltaR(track->eta(),track->phi(),PatMuonVector[0]->eta(),PatMuonVector[0]->phi()) > 0.4) && (deltaR(track->eta(),track->phi(),PatMuonVector[1]->eta(),PatMuonVector[1]->phi()) > 0.4))
      {
	goodTracksCountm04++;
      }

      if ((deltaR(track->eta(),track->phi(),PatMuonVector[0]->eta(),PatMuonVector[0]->phi()) > 0.5) && (deltaR(track->eta(),track->phi(),PatMuonVector[1]->eta(),PatMuonVector[1]->phi()) > 0.5))
      {
	goodTracksCountm05++;
      }

    }

    eventData.SetTracksNonConepatMuon03(goodTracksCountm03);
    eventData.SetTracksNonConepatMuon04(goodTracksCountm04);
    eventData.SetTracksNonConepatMuon05(goodTracksCountm05);

    if (debug){
      std::cout << ">>> Pat Muon" << std::endl;
      std::cout<<"Muon1 -> 0.3 Radion Rel Iso: "<<relIsoFirstMuonDr03<<" sumPt "<<muon1SumPtR03<<" emEt "<<muon1EmEtR03<<" hadEt "<<muon1HadEtR03<<std::endl;
      std::cout<<"Muon1 -> 0.5 Radion Rel Iso: "<<relIsoFirstMuonDr05<<" sumPt "<<muon1SumPtR05<<" emEt "<<muon1EmEtR05<<" hadEt "<<muon1HadEtR05<<std::endl;
      std::cout << "Muon1 -> trackIso(): " << PatMuonVector[0]->trackIso() << " | muon1 -> ecalIso(): " << PatMuonVector[0]->ecalIso() << " | muon1 -> hcalIso(): " << PatMuonVector[0]->hcalIso() << " | PatMuonVector[0]->Iso(): " << relIsoFirstMuon << std::endl; 
      std::cout<<"Muon2 -> 0.3 Radion Rel Iso: "<<relIsoSecondMuonDr03<<" sumPt "<<muon2SumPtR03<<" emEt "<<muon2EmEtR03<<" hadEt "<<muon2HadEtR03<<std::endl;
      std::cout<<"Muon2 -> 0.5 Radion Rel Iso: "<<relIsoSecondMuonDr05<<" sumPt "<<muon2SumPtR05<<" emEt "<<muon2EmEtR05<<" hadEt "<<muon2HadEtR05<<std::endl;
      std::cout << "Muon2 -> trackIso(): " << PatMuonVector[1]->trackIso() << " | muon2 -> ecalIso(): " << PatMuonVector[1]->ecalIso() << " | muon2 -> hcalIso(): " << PatMuonVector[1]->hcalIso() << " | PatMuonVector[1]->Iso(): " << relIsoSecondMuon << std::endl;  
      std::cout << "NSize: " << muons->size() << std::endl;
      std::cout << "Muon, pT 1: " << PatMuonVector[0]->pt() << std::endl;
      std::cout << "Muon, pT 2: " << PatMuonVector[1]->pt() << std::endl;
      std::cout << "Muon, eta 1: " << PatMuonVector[0]->eta() << std::endl;
      std::cout << "Muon, eta 2: " << PatMuonVector[1]->eta() << std::endl;
      std::cout << "Muon1, p4(): " << PatMuonVector[0]->p4() << std::endl;
      std::cout << "Muon2, p4(): " << PatMuonVector[1]->p4() << std::endl;
      std::cout << "DiMuon, M(): " << DipatMuonSystem.M() << std::endl;
      std::cout << "Eta Z: " << DipatMuonSystem.eta() << std::endl;
      std::cout << "Phi Z: " << DipatMuonSystem.phi() << std::endl;
      std::cout << "pT Z: " << DipatMuonSystem.pt() << std::endl;
      std::cout << "energy Z: " << DipatMuonSystem.energy() << std::endl;
      std::cout << "pz Z: " << DipatMuonSystem.pz() << std::endl;
      std::cout << "" << std::endl;
    }
  }
  else {
    eventData.SetPatMuon1Pt(-999.);
    eventData.SetPatMuon1Charge(-999);
    eventData.SetPatMuon1Phi(-999.);
    eventData.SetPatMuon1Eta(-999.);
    eventData.SetPatMuon1Et(-999.);
    eventData.SetPatMuon2Pt(-999.);
    eventData.SetPatMuon2Charge(-999);
    eventData.SetPatMuon2Phi(-999.);
    eventData.SetPatMuon2Eta(-999.);
    eventData.SetPatMuon2Et(-999.);
    eventData.SetPatMuon1SumPtR03(-999.);
    eventData.SetPatMuon1EmEtR03(-999.);
    eventData.SetPatMuon1HadEtR03(-999.);
    eventData.SetPatMuon1SumPtR05(-999.);
    eventData.SetPatMuon1EmEtR05(-999.);
    eventData.SetPatMuon1HadEtR05(-999.);
    eventData.SetPatMuon2SumPtR03(-999.);
    eventData.SetPatMuon2EmEtR03(-999.);
    eventData.SetPatMuon2HadEtR03(-999.);
    eventData.SetPatMuon2SumPtR05(-999.);
    eventData.SetPatMuon2EmEtR05(-999.);
    eventData.SetPatMuon2HadEtR05(-999.);
    eventData.SetPatMuon1relIsoDr03(-999.);
    eventData.SetPatMuon2relIsoDr03(-999.);
    eventData.SetPatMuon1relIsoDr05(-999.);
    eventData.SetPatMuon2relIsoDr05(-999.);
    eventData.SetPatMuon1relIso(-999.);
    eventData.SetPatMuon2relIso(-999.);
    eventData.SetPatMuon1TrackerHits(-999.);
    eventData.SetPatMuon1PixelHits(-999.);
    eventData.SetPatMuon1NormalizedChi2(-999.);
    eventData.SetPatMuon1MatchedStations(-999.);
    eventData.SetPatMuon1Dxy(-999.);
    eventData.SetPatMuon1IsGlobal(false);
    eventData.SetPatMuon1IsTracker(false);
    eventData.SetPatMuon2TrackerHits(-999.);
    eventData.SetPatMuon2PixelHits(-999.);
    eventData.SetPatMuon2NormalizedChi2(-999.);
    eventData.SetPatMuon2MatchedStations(-999.);
    eventData.SetPatMuon2Dxy(-999.);
    eventData.SetPatMuon2IsGlobal(false);
    eventData.SetPatMuon2IsTracker(false);
    eventData.SetPatDiMuonMass(-999.);
    eventData.SetPatDiMuonPt(-999.);
    eventData.SetPatDiMuonEta(-999.);
    eventData.SetPatDiMuonPhi(-999.);
  } 

  if(PatElectronVector.size()>1){
    double relIsoFirstElectronDr03 = (PatElectronVector[0]->dr03TkSumPt()+PatElectronVector[0]->dr03EcalRecHitSumEt()+PatElectronVector[0]->dr03HcalTowerSumEt())/PatElectronVector[0]->et();
    double relIsoFirstElectronDr04 = (PatElectronVector[0]->dr04TkSumPt()+PatElectronVector[0]->dr04EcalRecHitSumEt()+PatElectronVector[0]->dr04HcalTowerSumEt())/PatElectronVector[0]->et();
    double relIsoSecondElectronDr03 = (PatElectronVector[1]->dr03TkSumPt()+PatElectronVector[1]->dr03EcalRecHitSumEt()+PatElectronVector[1]->dr03HcalTowerSumEt())/PatElectronVector[1]->et();
    double relIsoSecondElectronDr04 = (PatElectronVector[1]->dr04TkSumPt()+PatElectronVector[1]->dr04EcalRecHitSumEt()+PatElectronVector[1]->dr04HcalTowerSumEt())/PatElectronVector[1]->et();
    double InnerHits1 = PatElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    double InnerHits2 = PatElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    // Dielectron Mass
    math::XYZTLorentzVector DipatElectronSystem(0.,0.,0.,0.);
    DipatElectronSystem += PatElectronVector[0]->p4();
    DipatElectronSystem += PatElectronVector[1]->p4();
    eventData.SetPatDiElectronMass(DipatElectronSystem.M());
    eventData.SetPatDiElectronEta(DipatElectronSystem.eta());
    eventData.SetPatDiElectronPhi(DipatElectronSystem.phi());
    eventData.SetPatDiElectronPt(DipatElectronSystem.pt());

    // Fill Electron Variables
    eventData.SetPatNElectron(PatElectronVector.size());
    eventData.SetPatElectron1Pt(PatElectronVector[0]->pt());
    eventData.SetPatElectron1Charge(PatElectronVector[0]->charge());
    eventData.SetPatElectron1Phi(PatElectronVector[0]->phi());
    eventData.SetPatElectron1Eta(PatElectronVector[0]->eta());
    eventData.SetPatElectron1Et(PatElectronVector[0]->et());
    eventData.SetPatElectron1P4(PatElectronVector[0]->p4());
    eventData.SetPatElectron2Pt(PatElectronVector[1]->pt());
    eventData.SetPatElectron2Charge(PatElectronVector[1]->charge());
    eventData.SetPatElectron2Phi(PatElectronVector[1]->phi());
    eventData.SetPatElectron2Eta(PatElectronVector[1]->eta());
    eventData.SetPatElectron2Et(PatElectronVector[1]->et());
    eventData.SetPatElectron2P4(PatElectronVector[1]->p4());
    eventData.SetPatElectron1TkDr03(PatElectronVector[0]->dr03TkSumPt());
    eventData.SetPatElectron1EcalDr03(PatElectronVector[0]->dr03EcalRecHitSumEt());
    eventData.SetPatElectron1HcalDr03(PatElectronVector[0]->dr03HcalTowerSumEt());
    eventData.SetPatElectron2TkDr03(PatElectronVector[1]->dr03TkSumPt());
    eventData.SetPatElectron2EcalDr03(PatElectronVector[1]->dr03EcalRecHitSumEt());
    eventData.SetPatElectron2HcalDr03(PatElectronVector[1]->dr03HcalTowerSumEt());
    eventData.SetPatElectron1TkDr04(PatElectronVector[0]->dr04TkSumPt());
    eventData.SetPatElectron1EcalDr04(PatElectronVector[0]->dr04EcalRecHitSumEt());
    eventData.SetPatElectron1HcalDr04(PatElectronVector[0]->dr04HcalTowerSumEt());
    eventData.SetPatElectron2TkDr04(PatElectronVector[1]->dr04TkSumPt());
    eventData.SetPatElectron2EcalDr04(PatElectronVector[1]->dr04EcalRecHitSumEt());
    eventData.SetPatElectron2HcalDr04(PatElectronVector[1]->dr04HcalTowerSumEt());
    eventData.SetPatElectron1relIsoDr03(relIsoFirstElectronDr03);
    eventData.SetPatElectron1relIsoDr04(relIsoFirstElectronDr04);
    eventData.SetPatElectron2relIsoDr03(relIsoSecondElectronDr03);
    eventData.SetPatElectron2relIsoDr04(relIsoSecondElectronDr04);
    eventData.SetPatElectron1DeltaPhiTkClu(PatElectronVector[0]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetPatElectron1DeltaEtaTkClu(PatElectronVector[0]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetPatElectron1SigmaIeIe(PatElectronVector[0]->sigmaIetaIeta());
    eventData.SetPatElectron1DCot(PatElectronVector[0]->convDcot());
    eventData.SetPatElectron1Dist(PatElectronVector[0]->convDist());
    eventData.SetPatElectron1InnerHits(InnerHits1);
    eventData.SetPatElectron1HE(PatElectronVector[0]->hadronicOverEm());
    eventData.SetPatElectron2DeltaPhiTkClu(PatElectronVector[1]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetPatElectron2DeltaEtaTkClu(PatElectronVector[1]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetPatElectron2SigmaIeIe(PatElectronVector[1]->sigmaIetaIeta());
    eventData.SetPatElectron2DCot(PatElectronVector[1]->convDcot());
    eventData.SetPatElectron2Dist(PatElectronVector[1]->convDist());
    eventData.SetPatElectron2InnerHits(InnerHits2);
    eventData.SetPatElectron2HE(PatElectronVector[1]->hadronicOverEm());

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCounte03 = 0;
    int goodTracksCounte04 = 0;
    int goodTracksCounte05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {

      if ((deltaR(track->eta(),track->phi(),PatElectronVector[0]->eta(),PatElectronVector[0]->phi()) > 0.3) && (deltaR(track->eta(),track->phi(),PatElectronVector[1]->eta(),PatElectronVector[1]->phi()) > 0.3))
      {
	goodTracksCounte03++;
      }

      if ((deltaR(track->eta(),track->phi(),PatElectronVector[0]->eta(),PatElectronVector[0]->phi()) > 0.4) && (deltaR(track->eta(),track->phi(),PatElectronVector[1]->eta(),PatElectronVector[1]->phi()) > 0.4))
      {
	goodTracksCounte04++;
      }

      if ((deltaR(track->eta(),track->phi(),PatElectronVector[0]->eta(),PatElectronVector[0]->phi()) > 0.5) && (deltaR(track->eta(),track->phi(),PatElectronVector[1]->eta(),PatElectronVector[1]->phi()) > 0.5))
      {
	goodTracksCounte05++;
      }

    }

    eventData.SetTracksNonConepatElectron03(goodTracksCounte03);
    eventData.SetTracksNonConepatElectron04(goodTracksCounte04);
    eventData.SetTracksNonConepatElectron05(goodTracksCounte05);

    if (debug) {
      std::cout << ">>> Pat Electron" << std::endl;
      std::cout << "electron1 -> dr03 TK: " << PatElectronVector[0]->dr03TkSumPt() << "| dr03 Ecal: " << PatElectronVector[0]->dr03EcalRecHitSumEt() << " | dr03 Hcal: " << PatElectronVector[0]->dr03HcalTowerSumEt() << std::endl;
      std::cout << "electron1 -> dr04 TK: " << PatElectronVector[0]->dr04TkSumPt() << "| dr04 Ecal: " << PatElectronVector[0]->dr04EcalRecHitSumEt() << " | dr04 Hcal: " << PatElectronVector[0]->dr04HcalTowerSumEt() <<  std::endl;
      std::cout << "electron2 -> dr03 TK: " << PatElectronVector[1]->dr03TkSumPt() << "| dr03 Ecal: " << PatElectronVector[1]->dr03EcalRecHitSumEt() << " | dr03 Hcal: " << PatElectronVector[1]->dr03HcalTowerSumEt() << std::endl;
      std::cout << "electron2 -> dr04 TK: " << PatElectronVector[1]->dr04TkSumPt() << "| dr04 Ecal: " << PatElectronVector[1]->dr04EcalRecHitSumEt() << " | dr04 Hcal: " << PatElectronVector[1]->dr04HcalTowerSumEt() <<  std::endl;
      std::cout << "Electron, pT 1: " << PatElectronVector[0]->pt() << std::endl;
      std::cout << "Electron, pT 2: " << PatElectronVector[1]->pt() << std::endl;
      std::cout << "Electron, eta 1: " << PatElectronVector[0]->eta() << std::endl;
      std::cout << "Electron, eta 2: " << PatElectronVector[1]->eta() << std::endl;
      std::cout << "Electron1, p4(): " << PatElectronVector[0]->p4() << std::endl;
      std::cout << "Electron2, p4(): " << PatElectronVector[1]->p4() << std::endl;
      std::cout << "DiElectron, M(): " << DipatElectronSystem.M() << std::endl;
      std::cout << "Eta Z: " << DipatElectronSystem.eta() << std::endl;
      std::cout << "Phi Z: " << DipatElectronSystem.phi() << std::endl;
      std::cout << "pT Z: " << DipatElectronSystem.pt() << std::endl;
      std::cout << "energy Z: " << DipatElectronSystem.energy() << std::endl;
      std::cout << "pz Z: " << DipatElectronSystem.pz() << std::endl;
      std::cout << "DeltaPhiTkClu, electron1: " << PatElectronVector[0]->deltaPhiSuperClusterTrackAtVtx() << std::endl;
      std::cout << "DeltaEtaTkClu, electron1: " << PatElectronVector[0]->deltaEtaSuperClusterTrackAtVtx() << std::endl;
      std::cout << "SigmaIeIe, electron1: " << PatElectronVector[0]->sigmaIetaIeta() << std::endl;
      std::cout << "Dcot, electron1: " << PatElectronVector[0]->convDcot() << std::endl;
      std::cout << "Dist, electron1: " << PatElectronVector[0]->convDist() << std::endl;
      std::cout << "Number Of Expected Inner Hits, electron1: " << PatElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << std::endl;
      std::cout << "H/E, electron1: " << PatElectronVector[0]->hadronicOverEm() << std::endl;
      std::cout << "DeltaPhiTkClu, electron2: " << PatElectronVector[1]->deltaPhiSuperClusterTrackAtVtx() << std::endl;
      std::cout << "DeltaEtaTkClu, electron2: " << PatElectronVector[1]->deltaEtaSuperClusterTrackAtVtx() << std::endl;
      std::cout << "SigmaIeIe, electron2: " << PatElectronVector[1]->sigmaIetaIeta() << std::endl;
      std::cout << "Dcot, electron2: " << PatElectronVector[1]->convDcot() << std::endl;
      std::cout << "Dist, electron2: " << PatElectronVector[1]->convDist() << std::endl;
      std::cout << "Number Of Expected Inner Hits, electron2: " << PatElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << std::endl;
      std::cout << "H/E, electron2: " << PatElectronVector[1]->hadronicOverEm() << std::endl;
      std::cout << "" << std::endl;
    }
  }
  else{
    eventData.SetPatElectron1Pt(-999.);
    eventData.SetPatElectron1Charge(-999);
    eventData.SetPatElectron1Phi(-999.);
    eventData.SetPatElectron1Eta(-999.);
    eventData.SetPatElectron1Et(-999.);
    eventData.SetPatElectron2Pt(-999.);
    eventData.SetPatElectron2Charge(-999);
    eventData.SetPatElectron2Phi(-999.);
    eventData.SetPatElectron2Eta(-999.);
    eventData.SetPatElectron2Et(-999.);
    eventData.SetPatElectron1TkDr03(-999.);
    eventData.SetPatElectron1EcalDr03(-999.);
    eventData.SetPatElectron1HcalDr03(-999.);
    eventData.SetPatElectron2TkDr03(-999.);
    eventData.SetPatElectron2EcalDr03(-999.);
    eventData.SetPatElectron2HcalDr03(-999.);
    eventData.SetPatElectron1TkDr04(-999.);
    eventData.SetPatElectron1EcalDr04(-999.);
    eventData.SetPatElectron1HcalDr04(-999.);
    eventData.SetPatElectron2TkDr04(-999.);
    eventData.SetPatElectron2EcalDr04(-999.);
    eventData.SetPatElectron2HcalDr04(-999.);
    eventData.SetPatElectron1relIsoDr03(-999.);
    eventData.SetPatElectron1relIsoDr04(-999.);
    eventData.SetPatElectron2relIsoDr03(-999.);
    eventData.SetPatElectron2relIsoDr04(-999.);
    eventData.SetPatDiElectronMass(-999.);
    eventData.SetPatDiElectronEta(-999.);
    eventData.SetPatDiElectronPhi(-999.);
    eventData.SetPatDiElectronPt(-999.);
    eventData.SetPatElectron1DeltaPhiTkClu(-999.);
    eventData.SetPatElectron1DeltaEtaTkClu(-999.);
    eventData.SetPatElectron1SigmaIeIe(-999.);
    eventData.SetPatElectron1DCot(-999.);
    eventData.SetPatElectron1Dist(-999.);
    eventData.SetPatElectron1InnerHits(-999.);
    eventData.SetPatElectron1HE(-999.);
    eventData.SetPatElectron2DeltaPhiTkClu(-999.);
    eventData.SetPatElectron2DeltaEtaTkClu(-999.);
    eventData.SetPatElectron2SigmaIeIe(-999.);
    eventData.SetPatElectron2DCot(-999.);
    eventData.SetPatElectron2Dist(-999.);
    eventData.SetPatElectron2InnerHits(-999.);
    eventData.SetPatElectron2HE(-999.);
  }

}

//
// Fill Castor 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillCastor(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Phi: 16 modules, rh.id().sector(); 
  // Z: 14 modules, rh.id().module(); 
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

void DiffractiveZAnalysis::fillCastorDebug(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Phi: 16 modules, rh.id().sector();
  // Z: 14 modules, rh.id().module();
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
// Fill ZDC
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveZAnalysis::fillZDC(DiffractiveZEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // ZDC have two sections: section 1 = EM, section 2 = HAD. EM has 5 modules. Had has 4 modules. 

  std::vector<std::vector<double> > digiAllPMT;
  std::vector<double> digiPMT;

  bool debug = false;
  bool debug_deep = false;

  double ZDCNSumEMEnergy = 0.;
  double ZDCNSumHADEnergy = 0.;
  double ZDCPSumEMEnergy = 0.;
  double ZDCPSumHADEnergy = 0.;

  double ZDCNSumEMTime = 0.;
  double ZDCNSumHADTime = 0.;
  double ZDCPSumEMTime = 0.;
  double ZDCPSumHADTime = 0.;

  int DigiDataADC[180];
  float DigiDatafC[180];

  edm::Handle <ZDCDigiCollection> zdc_digi_h;
  event.getByType(zdc_digi_h);
  edm::ESHandle<HcalDbService> conditions;
  const ZDCDigiCollection *zdc_digi = zdc_digi_h.failedToGet()? 0 : &*zdc_digi_h;

  edm::Handle <ZDCRecHitCollection> zdc_recHits_h;
  event.getByLabel(zdcHitsTag_, zdc_recHits_h);
  const ZDCRecHitCollection *zdc_recHits = zdc_recHits_h.failedToGet()? 0 : &*zdc_recHits_h;

  setup.get<HcalDbRecord>().get(conditions);

  if (zdc_recHits) {
    for (ZDCRecHitCollection::const_iterator zhit = zdc_recHits->begin(); zhit != zdc_recHits->end(); zhit++){		
      // Some Variables
      int ZDCSide      = (zhit->id()).zside();
      int ZDCSection   = (zhit->id()).section();
      //Float_t ZDCEnergy = zhit->energy();
      //Float_t ZDCRecHitTime = zhit->time();
      //int ZDCChannel   = (zhit->id()).channel();

      if (zhit->energy() >= 0.){
	if (ZDCSide == -1){
	  if (ZDCSection == 1 ){
	    ZDCNSumEMEnergy += zhit->energy();
	    ZDCNSumEMTime += zhit->time();
	  }
	  if (ZDCSection == 2 ){
	    ZDCNSumHADEnergy += zhit->energy();
	    ZDCNSumHADTime += zhit->time();
	  }
	}
	if (ZDCSide == 1){
	  if (ZDCSection == 1 ){
	    ZDCPSumEMEnergy += zhit->energy();
	    ZDCPSumEMTime += zhit->time();
	  }
	  if (ZDCSection == 2 ){
	    ZDCPSumHADEnergy += zhit->energy();
	    ZDCPSumHADTime += zhit->time();
	  }
	}
      }

    }

    if (debug){
      std::cout << "ZDC + | Total EM Energy: " << ZDCPSumEMEnergy << std::endl;
      std::cout << "ZDC + | Total HAD Energy: " << ZDCPSumHADEnergy << std::endl;
      std::cout << "ZDC + | EM <Time>: " << ZDCPSumEMTime/5. << std::endl;
      std::cout << "ZDC + | HAD <Time>: " << ZDCPSumHADTime/4. << std::endl;
      std::cout << "ZDC - | Total EM Energy: " << ZDCNSumEMEnergy << std::endl;
      std::cout << "ZDC - | Total HAD Energy: " << ZDCNSumHADEnergy << std::endl;
      std::cout << "ZDC - | EM <Time>: " << ZDCNSumEMTime/5. << std::endl;
      std::cout << "ZDC - | HAD <Time>: " << ZDCNSumHADTime/4. << std::endl;
    }

  }

  if (zdc_digi){
    for(int i=0; i<180; i++){DigiDatafC[i]=0;DigiDataADC[i]=0;}
    for (ZDCDigiCollection::const_iterator j=zdc_digi->begin();j!=zdc_digi->end();j++){
      const ZDCDataFrame digi = (const ZDCDataFrame)(*j);		
      int iSide      = digi.id().zside();
      int iSection   = digi.id().section();
      int iChannel   = digi.id().channel();
      int chid = (iSection-1)*5+(iSide+1)/2*9+(iChannel-1);

      const HcalQIEShape* qieshape=conditions->getHcalShape();
      const HcalQIECoder* qiecoder=conditions->getHcalCoder(digi.id());
      CaloSamples caldigi;
      HcalCoderDb coder(*qiecoder,*qieshape);

      coder.adc2fC(digi,caldigi);

      int fTS = digi.size();

      for (int i = 0; i < fTS; ++i) {
	DigiDatafC[i+chid*10] = caldigi[i];
	DigiDataADC[i+chid*10] = digi[i].adc();
	digiPMT.push_back(DigiDatafC[i+chid*10]);
	if (debug_deep){
	  std::cout << "Digi Size: " << fTS << std::endl;
	  std::cout << "DigiDataADC["<<i+chid*10<<"]: " << DigiDataADC[i+chid*10] << std::endl;
	  std::cout << "DigiDatafC["<<i+chid*10<<"]: " << DigiDatafC[i+chid*10] << std::endl;
	}
      }

      if (debug){
	std::cout << "iSide: " << iSide << std::endl;
	std::cout << "iSection: " << iSection << std::endl;
	std::cout << "iChannel: " << iChannel << std::endl;
	std::cout << "chid: " << chid << std::endl;
      }
      digiAllPMT.push_back(digiPMT);
    }
    eventData.SetZDCdigifC(digiAllPMT);
  }

}

template <class T, class W>
math::XYZTLorentzVector DiffractiveZAnalysis::DiSystem(T obj1, W obj2){
  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1->p4();
  DiObj += obj2->p4();
  return DiObj;
}

template <class T, class W>
double DiffractiveZAnalysis::InvariantMass(T lepton1, W lepton2){
  double mass_=-999.;
  math::XYZTLorentzVector DiSystem(0.,0.,0.,0.);
  DiSystem += lepton1->p4();
  DiSystem += lepton2->p4();
  mass_ = DiSystem.M();
  // Defense
  if (!std::isfinite(mass_)) {
    mass_ = -999.;
  }
  return mass_;
}
