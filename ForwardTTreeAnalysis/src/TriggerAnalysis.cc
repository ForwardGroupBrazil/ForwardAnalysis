#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/TriggerAnalysis.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/TriggerEvent.h"

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

#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include <stdio.h>
#include <math.h> 
#include <cmath>

using triggerAnalysis::TriggerAnalysis;

const char* TriggerAnalysis::name = "TriggerAnalysis";

TriggerAnalysis::TriggerAnalysis(const edm::ParameterSet& pset):
  triggerResultsTag_(pset.getParameter<edm::InputTag>("TriggerResultsTag")),
  hltPathNames_(pset.getParameter<std::vector<std::string> >("hltPaths")),
  electronTag_(pset.getParameter<edm::InputTag>("electronTag")),
  muonTag_(pset.getParameter<edm::InputTag>("muonTag")),
  metTag_(pset.getParameter<edm::InputTag>("metTag")),
  PVtxCollectionTag_(pset.getParameter<edm::InputTag>("PVtxCollectionTag")),
  trackTag_(pset.getParameter<edm::InputTag>("TrackTag"))
{
}

void TriggerAnalysis::setTFileService(){

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

}

TriggerAnalysis::~TriggerAnalysis(){}

void TriggerAnalysis::begin() {
  setTFileService();
}

void TriggerAnalysis::begin(const edm::Run& run, const edm::EventSetup& setup) {}

void TriggerAnalysis::end() {}

void TriggerAnalysis::fill(TriggerEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  eventData.reset();

  fillTriggerInfo(eventData,event,setup);
  fillMETInfo(eventData,event,setup);
  fillMuonsInfo(eventData,event,setup);
  fillElectronsInfo(eventData,event,setup);

}

// Fill Trigger
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TriggerAnalysis::fillTriggerInfo(TriggerEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

// Fill MET
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TriggerAnalysis::fillMETInfo(TriggerEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  // Fill pfMET
  edm::Handle<reco::PFMETCollection> met;
  event.getByLabel(metTag_,met);

  NeutrinoVector.clear();
  int neutrinosize = met->size();
  int itmet;

  if(met->size()>0){
    for(itmet=0; itmet < neutrinosize; ++itmet){
      const reco::PFMET* neutrinoAll = &((*met)[itmet]);
      NeutrinoVector.push_back(neutrinoAll);
    }

  }

  if(NeutrinoVector.size()>0){
    // Sorting Vector
    std::sort(NeutrinoVector.begin(), NeutrinoVector.end(), orderPT());
    for (unsigned int i=0;i<NeutrinoVector.size();i++){
      if (debug) std::cout << "reco::pfMET[" << i << "]\t\t---> pT [GeV]: " << NeutrinoVector[i]->pt() << " | eT [GeV]: " << NeutrinoVector[i]->et() << " | sum eT [GeV]: " << NeutrinoVector[i]->sumEt() << " | eta: " << NeutrinoVector[i]->eta() << " | phi: " << NeutrinoVector[i]->phi() << " | px [GeV]: " << NeutrinoVector[i]->px() << " | py [GeV]: " << NeutrinoVector[i]->py() << " | P4() [GeV]: " << NeutrinoVector[i]->p4()  << " | significance: " << NeutrinoVector[i]->significance() << std::endl;
    }
    eventData.SetMETPt(NeutrinoVector[0]->pt());
    eventData.SetMETPhi(NeutrinoVector[0]->phi());
    eventData.SetMETEt(NeutrinoVector[0]->et());
    eventData.SetMETSumEt(NeutrinoVector[0]->sumEt());
    eventData.SetMETpx(NeutrinoVector[0]->px());
    eventData.SetMETpy(NeutrinoVector[0]->py());
    eventData.SetMETP4(NeutrinoVector[0]->p4());
    eventData.SetMETsigma(NeutrinoVector[0]->significance());
  }else{
    eventData.SetMETPt(-999.);
    eventData.SetMETPhi(-999.);
    eventData.SetMETEt(-999.);
    eventData.SetMETSumEt(-999.);
    eventData.SetMETpx(-999.);
    eventData.SetMETpy(-999.);
    eventData.SetMETsigma(-999);
  }

}

// Fill Reco::Electron
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TriggerAnalysis::fillElectronsInfo(TriggerEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

  if(ElectronVector.size()>0){

    double relIsoFirstElectronDr03 = (ElectronVector[0]->dr03TkSumPt()+ElectronVector[0]->dr03EcalRecHitSumEt()+ElectronVector[0]->dr03HcalTowerSumEt())/ElectronVector[0]->et();
    double relIsoFirstElectronDr04 = (ElectronVector[0]->dr04TkSumPt()+ElectronVector[0]->dr04EcalRecHitSumEt()+ElectronVector[0]->dr04HcalTowerSumEt())/ElectronVector[0]->et();
    double InnerHits1 = ElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    eventData.SetElectronsN(ElectronVector.size());
    eventData.SetLeadingElectronPt(ElectronVector[0]->pt());
    eventData.SetLeadingElectronEta(ElectronVector[0]->eta());
    eventData.SetLeadingElectronPhi(ElectronVector[0]->phi());
    eventData.SetLeadingElectronCharge(ElectronVector[0]->charge());
    eventData.SetLeadingElectronP4(ElectronVector[0]->p4());

    eventData.SetLeadingElectronTkDr03(ElectronVector[0]->dr03TkSumPt());
    eventData.SetLeadingElectronEcalDr03(ElectronVector[0]->dr03EcalRecHitSumEt());
    eventData.SetLeadingElectronHcalDr03(ElectronVector[0]->dr03HcalTowerSumEt());

    eventData.SetLeadingElectronTkDr04(ElectronVector[0]->dr04TkSumPt());
    eventData.SetLeadingElectronEcalDr04(ElectronVector[0]->dr04EcalRecHitSumEt());
    eventData.SetLeadingElectronHcalDr04(ElectronVector[0]->dr04HcalTowerSumEt());

    eventData.SetLeadingElectronrelIsoDr03(relIsoFirstElectronDr03);
    eventData.SetLeadingElectronrelIsoDr04(relIsoFirstElectronDr04);

    eventData.SetLeadingElectronDeltaPhiTkClu(ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetLeadingElectronDeltaEtaTkClu(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetLeadingElectronSigmaIeIe(ElectronVector[0]->sigmaIetaIeta());
    eventData.SetLeadingElectronDCot(ElectronVector[0]->convDcot());
    eventData.SetLeadingElectronDist(ElectronVector[0]->convDist());
    eventData.SetLeadingElectronInnerHits(InnerHits1);
    eventData.SetLeadingElectronHE(ElectronVector[0]->hadronicOverEm());

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCountLeadingElectronR03 = 0;
    int goodTracksCountLeadingElectronR04 = 0;
    int goodTracksCountLeadingElectronR05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),ElectronVector[0]->eta(),ElectronVector[0]->phi()) > 0.3))
      {
	goodTracksCountLeadingElectronR03++;
      }

      if ((deltaR(track->eta(),track->phi(),ElectronVector[0]->eta(),ElectronVector[0]->phi()) > 0.4))
      {
	goodTracksCountLeadingElectronR04++;
      }

      if ((deltaR(track->eta(),track->phi(),ElectronVector[0]->eta(),ElectronVector[0]->phi()) > 0.5))
      {
	goodTracksCountLeadingElectronR05++;
      }

    }

    eventData.SetTracksNonConeLeadingElectronR03(goodTracksCountLeadingElectronR03);
    eventData.SetTracksNonConeLeadingElectronR04(goodTracksCountLeadingElectronR04);
    eventData.SetTracksNonConeLeadingElectronR05(goodTracksCountLeadingElectronR05);

    if (debug){
      std::cout << ">>> Reco Electron" << std::endl;
      std::cout << "electron1 -> dr03 TK: " << ElectronVector[0]->dr03TkSumPt() << "| dr03 Ecal: " << ElectronVector[0]->dr03EcalRecHitSumEt() << " | dr03 Hcal: " << ElectronVector[0]->dr03HcalTowerSumEt() << std::endl;
      std::cout << "electron1 -> dr04 TK: " << ElectronVector[0]->dr04TkSumPt() << "| dr04 Ecal: " << ElectronVector[0]->dr04EcalRecHitSumEt() << " | dr04 Hcal: " << ElectronVector[0]->dr04HcalTowerSumEt() <<  std::endl;
      std::cout << "Electron Isolation: " << relIsoFirstElectronDr03 << " | " << relIsoFirstElectronDr04 << std::endl;
      std::cout << "N Electrons: " << ElectronVector.size() << std::endl;
      std::cout << "Electron, pT 1: " << ElectronVector[0]->pt() << std::endl;
      std::cout << "Electron, eta 1: " << ElectronVector[0]->eta() << std::endl;
      std::cout << "DeltaPhiTkClu, electron1: " << ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx() << std::endl;
      std::cout << "DeltaEtaTkClu, electron1: " << ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx() << std::endl;
      std::cout << "SigmaIeIe, electron1: " << ElectronVector[0]->sigmaIetaIeta() << std::endl;
      std::cout << "Dcot, electron1: " << ElectronVector[0]->convDcot() << std::endl;
      std::cout << "Dist, electron1: " << ElectronVector[0]->convDist() << std::endl;
      std::cout << "Number Of Expected Inner Hits, electron1: " << ElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << std::endl;
      std::cout << "H/E, electron1: " << ElectronVector[0]->hadronicOverEm() << std::endl;
      std::cout << "" << std::endl;
    }

    double isoTk1 = ElectronVector[0]->dr03TkSumPt()/ElectronVector[0]->pt();
    double isoEcal1 = ElectronVector[0]->dr03EcalRecHitSumEt()/ElectronVector[0]->pt();
    double isoHcal1 = ElectronVector[0]->dr03HcalTowerSumEt()/ElectronVector[0]->pt();
    bool isoBarrel1WP95 = false;
    bool isoEndCap1WP95 = false;
    bool isolation1WP95 = false;
    bool eleEndCap1WP95 = false;
    bool eleBarrel1WP95 = false;
    bool candSel1WP95 = false;
    bool isoBarrel1WP80 = false;
    bool isoEndCap1WP80 = false;
    bool isolation1WP80 = false;
    bool eleEndCap1WP80 = false;
    bool eleBarrel1WP80 = false;
    bool candSel1WP80 = false;

    //Isolation Barrel
    if ((fabs(ElectronVector[0]->eta()) <= 1.4442) ){
      if (isoTk1<0.15 && isoEcal1<2.0 && isoHcal1<0.12) isoBarrel1WP95 = true;
      if (isoTk1<0.09 && isoEcal1<0.07 && isoHcal1<0.10) isoBarrel1WP80 = true;
    }

    // Isolation Endcap
    if ((fabs(ElectronVector[0]->eta()) >= 1.5660) && (fabs(ElectronVector[0]->eta()) <= 2.5)){
      if (isoTk1<0.08 && isoEcal1<0.06 && isoHcal1<0.05) isoEndCap1WP95 = true;
      if (isoTk1<0.04 && isoEcal1<0.05 && isoHcal1<0.025) isoEndCap1WP80 = true;
    }

    if ((isoEndCap1WP95 || isoBarrel1WP95)) isolation1WP95 = true;
    if ((isoEndCap1WP80 || isoBarrel1WP80)) isolation1WP80 = true;

    // Quality criteria Barrel
    if ((fabs(ElectronVector[0]->eta()) <= 1.4442) ){
      if (InnerHits1 <= 1 && fabs(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx()) < 0.007 && ElectronVector[0]->sigmaIetaIeta() < 0.01 && ElectronVector[0]->hadronicOverEm() < 0.15 ) eleBarrel1WP95 = true;
      if (InnerHits1 <= 0 && (fabs(ElectronVector[0]->convDcot()) >= 0.02 || fabs(ElectronVector[0]->convDist()) >= 0.02 ) && fabs(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx()) < 0.004 && fabs(ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx()) < 0.06 && ElectronVector[0]->sigmaIetaIeta() < 0.01 && ElectronVector[0]->hadronicOverEm() < 0.04 ) eleBarrel1WP80 = true;
    }

    // Quality criteria Endcap
    if ((fabs(ElectronVector[0]->eta()) >= 1.5660) && (fabs(ElectronVector[0]->eta()) <= 2.5)){
      if (InnerHits1 <= 1 && fabs(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx()) < 0.01 && ElectronVector[0]->sigmaIetaIeta() < 0.03 && ElectronVector[0]->hadronicOverEm() < 0.07) eleEndCap1WP95 = true;
      if (InnerHits1 <= 0 && (fabs(ElectronVector[0]->convDcot()) >= 0.02 || fabs(ElectronVector[0]->convDist()) >= 0.02 ) && fabs(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx()) < 0.007 && fabs(ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx()) < 0.03 && ElectronVector[0]->sigmaIetaIeta() < 0.03 && ElectronVector[0]->hadronicOverEm() < 0.025) eleEndCap1WP80 = true;
    }

    if ((eleEndCap1WP95 || eleBarrel1WP95)) candSel1WP95 = true;
    if ((eleEndCap1WP80 || eleBarrel1WP80)) candSel1WP80 = true;

    if(isolation1WP80 && candSel1WP80){
      eventData.SetLeadingElectronIsWP80(true);
    }else{
      eventData.SetLeadingElectronIsWP80(false);
    }

    if(isolation1WP95 && candSel1WP95){
      eventData.SetLeadingElectronIsWP95(true);
    }else{
      eventData.SetLeadingElectronIsWP95(false);
    }

  }
  else {
    eventData.SetElectronsN(-1);
    eventData.SetLeadingElectronPt(-999.);
    eventData.SetLeadingElectronEta(-999.);
    eventData.SetLeadingElectronPhi(-999.);
    eventData.SetLeadingElectronCharge(-999);

    eventData.SetLeadingElectronTkDr03(-999.);
    eventData.SetLeadingElectronEcalDr03(-999.);
    eventData.SetLeadingElectronHcalDr03(-999.);

    eventData.SetLeadingElectronTkDr04(-999.);
    eventData.SetLeadingElectronEcalDr04(-999.);
    eventData.SetLeadingElectronHcalDr04(-999.);

    eventData.SetLeadingElectronrelIsoDr03(-999.);
    eventData.SetLeadingElectronrelIsoDr04(-999.);

    eventData.SetLeadingElectronDeltaPhiTkClu(-999.);
    eventData.SetLeadingElectronDeltaEtaTkClu(-999.);
    eventData.SetLeadingElectronSigmaIeIe(-999.);
    eventData.SetLeadingElectronDCot(-999.);
    eventData.SetLeadingElectronDist(-999.);
    eventData.SetLeadingElectronInnerHits(-999.);
    eventData.SetLeadingElectronHE(-999.);

    eventData.SetLeadingElectronIsWP95(false);
    eventData.SetLeadingElectronIsWP80(false);

    eventData.SetTracksNonConeLeadingElectronR03(-999.);
    eventData.SetTracksNonConeLeadingElectronR04(-999.);
    eventData.SetTracksNonConeLeadingElectronR05(-999.);

  }

  // Second Lepton Info
  if(ElectronVector.size()>1){

    double relIsoSecondElectronDr03 = (ElectronVector[1]->dr03TkSumPt()+ElectronVector[1]->dr03EcalRecHitSumEt()+ElectronVector[1]->dr03HcalTowerSumEt())/ElectronVector[1]->et();
    double relIsoSecondElectronDr04 = (ElectronVector[1]->dr04TkSumPt()+ElectronVector[1]->dr04EcalRecHitSumEt()+ElectronVector[1]->dr04HcalTowerSumEt())/ElectronVector[1]->et();
    double InnerHits2 = ElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    eventData.SetSecondElectronPt(ElectronVector[1]->pt());
    eventData.SetSecondElectronEta(ElectronVector[1]->eta());
    eventData.SetSecondElectronPhi(ElectronVector[1]->phi());
    eventData.SetSecondElectronCharge(ElectronVector[1]->charge());
    eventData.SetSecondElectronP4(ElectronVector[1]->p4());
    eventData.SetSecondElectronTkDr03(ElectronVector[1]->dr03TkSumPt());
    eventData.SetSecondElectronEcalDr03(ElectronVector[1]->dr03EcalRecHitSumEt());
    eventData.SetSecondElectronHcalDr03(ElectronVector[1]->dr03HcalTowerSumEt());
    eventData.SetSecondElectronTkDr04(ElectronVector[1]->dr04TkSumPt());
    eventData.SetSecondElectronEcalDr04(ElectronVector[1]->dr04EcalRecHitSumEt());
    eventData.SetSecondElectronHcalDr04(ElectronVector[1]->dr04HcalTowerSumEt());
    eventData.SetSecondElectronrelIsoDr03(relIsoSecondElectronDr03);
    eventData.SetSecondElectronrelIsoDr04(relIsoSecondElectronDr04);
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

    int goodTracksCountSecondElectronR03 = 0;
    int goodTracksCountSecondElectronR04 = 0;
    int goodTracksCountSecondElectronR05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),ElectronVector[1]->eta(),ElectronVector[1]->phi()) > 0.3))
      {
	goodTracksCountSecondElectronR03++;
      }
      if ((deltaR(track->eta(),track->phi(),ElectronVector[1]->eta(),ElectronVector[1]->phi()) > 0.4))
      {
	goodTracksCountSecondElectronR04++;
      }
      if ((deltaR(track->eta(),track->phi(),ElectronVector[1]->eta(),ElectronVector[1]->phi()) > 0.5))
      {
	goodTracksCountSecondElectronR05++;
      }
    }
    eventData.SetTracksNonConeSecondElectronR03(goodTracksCountSecondElectronR03);
    eventData.SetTracksNonConeSecondElectronR04(goodTracksCountSecondElectronR04);
    eventData.SetTracksNonConeSecondElectronR05(goodTracksCountSecondElectronR05);

    if (debug){
      std::cout << ">>> Reco Electron" << std::endl;
      std::cout << "electron2 -> dr03 TK: " << ElectronVector[1]->dr03TkSumPt() << "| dr03 Ecal: " << ElectronVector[1]->dr03EcalRecHitSumEt() << " | dr03 Hcal: " << ElectronVector[1]->dr03HcalTowerSumEt() << std::endl;
      std::cout << "electron2 -> dr04 TK: " << ElectronVector[1]->dr04TkSumPt() << "| dr04 Ecal: " << ElectronVector[1]->dr04EcalRecHitSumEt() << " | dr04 Hcal: " << ElectronVector[1]->dr04HcalTowerSumEt() << std::endl;
      std::cout << "Electron Isolation: " << relIsoSecondElectronDr03 << " | " << relIsoSecondElectronDr04 << std::endl;
      std::cout << "N Electrons: " << ElectronVector.size() << std::endl;
      std::cout << "Electron, pT 2: " << ElectronVector[1]->pt() << std::endl;
      std::cout << "Electron, eta 2: " << ElectronVector[1]->eta() << std::endl;
      std::cout << "DeltaPhiTkClu, electron2: " << ElectronVector[1]->deltaPhiSuperClusterTrackAtVtx() << std::endl;
      std::cout << "DeltaEtaTkClu, electron2: " << ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx() << std::endl;
      std::cout << "SigmaIeIe, electron2: " << ElectronVector[1]->sigmaIetaIeta() << std::endl;
      std::cout << "Dcot, electron2: " << ElectronVector[1]->convDcot() << std::endl;
      std::cout << "Dist, electron2: " << ElectronVector[1]->convDist() << std::endl;
      std::cout << "Number Of Expected Inner Hits, electron2: " << ElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << std::endl;
      std::cout << "H/E, electron2: " << ElectronVector[1]->hadronicOverEm() << std::endl;
      std::cout << "" << std::endl;
    }
    double isoTk2 = ElectronVector[1]->dr03TkSumPt()/ElectronVector[1]->pt();
    double isoEcal2 = ElectronVector[1]->dr03EcalRecHitSumEt()/ElectronVector[1]->pt();
    double isoHcal1 = ElectronVector[1]->dr03HcalTowerSumEt()/ElectronVector[1]->pt();
    bool isoBarrel2WP95 = false;
    bool isoEndCap2WP95 = false;
    bool isolation2WP95 = false;
    bool eleEndCap2WP95 = false;
    bool eleBarrel2WP95 = false;
    bool candSel2WP95 = false;
    bool isoBarrel2WP80 = false;
    bool isoEndCap2WP80 = false;
    bool isolation2WP80 = false;
    bool eleEndCap2WP80 = false;
    bool eleBarrel2WP80 = false;
    bool candSel2WP80 = false;

    //Isolation Barrel
    if ((fabs(ElectronVector[1]->eta()) <= 1.4442) ){
      if (isoTk2<0.15 && isoEcal2<2.0 && isoHcal1<0.12) isoBarrel2WP95 = true;
      if (isoTk2<0.09 && isoEcal2<0.07 && isoHcal1<0.10) isoBarrel2WP80 = true;
    }

    // Isolation Endcap
    if ((fabs(ElectronVector[1]->eta()) >= 1.5660) && (fabs(ElectronVector[1]->eta()) <= 2.5)){
      if (isoTk2<0.08 && isoEcal2<0.06 && isoHcal1<0.05) isoEndCap2WP95 = true;
      if (isoTk2<0.04 && isoEcal2<0.05 && isoHcal1<0.025) isoEndCap2WP80 = true;
    }
    if ((isoEndCap2WP95 || isoBarrel2WP95)) isolation2WP95 = true;
    if ((isoEndCap2WP80 || isoBarrel2WP80)) isolation2WP80 = true;

    // Quality criteria Barrel
    if ((fabs(ElectronVector[1]->eta()) <= 1.4442) ){
      if (InnerHits2 <= 1 && fabs(ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx()) < 0.007 && ElectronVector[1]->sigmaIetaIeta() < 0.01 && ElectronVector[1]->hadronicOverEm() < 0.15 ) eleBarrel2WP95 = true;
      if (InnerHits2 <= 0 && (fabs(ElectronVector[1]->convDcot()) >= 0.02 || fabs(ElectronVector[1]->convDist()) >= 0.02 ) && fabs(ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx()) < 0.004 && fabs(ElectronVector[1]->deltaPhiSuperClusterTrackAtVtx()) < 0.06 && ElectronVector[1]->sigmaIetaIeta() < 0.01 && ElectronVector[1]->hadronicOverEm() < 0.04 ) eleBarrel2WP80 = true;
    }

    // Quality criteria Endcap
    if ((fabs(ElectronVector[1]->eta()) >= 1.5660) && (fabs(ElectronVector[1]->eta()) <= 2.5)){
      if (InnerHits2 <= 1 && fabs(ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx()) < 0.01 && ElectronVector[1]->sigmaIetaIeta() < 0.03 && ElectronVector[1]->hadronicOverEm() < 0.07) eleEndCap2WP95 = true;
      if (InnerHits2 <= 0 && (fabs(ElectronVector[1]->convDcot()) >= 0.02 || fabs(ElectronVector[1]->convDist()) >= 0.02 ) && fabs(ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx()) < 0.007 && fabs(ElectronVector[1]->deltaPhiSuperClusterTrackAtVtx()) < 0.03 && ElectronVector[1]->sigmaIetaIeta() < 0.03 && ElectronVector[1]->hadronicOverEm() < 0.025) eleEndCap2WP80 = true;
    }
    if ((eleEndCap2WP95 || eleBarrel2WP95)) candSel2WP95 = true;
    if ((eleEndCap2WP80 || eleBarrel2WP80)) candSel2WP80 = true;
    if(isolation2WP80 && candSel2WP80){
      eventData.SetSecondElectronIsWP80(true);
    }else{
      eventData.SetSecondElectronIsWP80(false);
    }
    if(isolation2WP95 && candSel2WP95){
      eventData.SetSecondElectronIsWP95(true);
    }else{
      eventData.SetSecondElectronIsWP95(false);
    }
  }
  else {
    eventData.SetSecondElectronPt(-999.);
    eventData.SetSecondElectronEta(-999.);
    eventData.SetSecondElectronPhi(-999.);
    eventData.SetSecondElectronCharge(-999);
    eventData.SetSecondElectronTkDr03(-999.);
    eventData.SetSecondElectronEcalDr03(-999.);
    eventData.SetSecondElectronHcalDr03(-999.);
    eventData.SetSecondElectronTkDr04(-999.);
    eventData.SetSecondElectronEcalDr04(-999.);
    eventData.SetSecondElectronHcalDr04(-999.);
    eventData.SetSecondElectronrelIsoDr03(-999.);
    eventData.SetSecondElectronrelIsoDr04(-999.);
    eventData.SetSecondElectronDeltaPhiTkClu(-999.);
    eventData.SetSecondElectronDeltaEtaTkClu(-999.);
    eventData.SetSecondElectronSigmaIeIe(-999.);
    eventData.SetSecondElectronDCot(-999.);
    eventData.SetSecondElectronDist(-999.);
    eventData.SetSecondElectronInnerHits(-999.);
    eventData.SetSecondElectronHE(-999.);
    eventData.SetSecondElectronIsWP95(false);
    eventData.SetSecondElectronIsWP80(false);
    eventData.SetTracksNonConeSecondElectronR03(-999);
    eventData.SetTracksNonConeSecondElectronR04(-999);
    eventData.SetTracksNonConeSecondElectronR05(-999);

  }

}

// Fill Reco::Muon
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TriggerAnalysis::fillMuonsInfo(TriggerEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

  if(MuonVector.size()>0){

    double muon1SumPtR03 = MuonVector[0]->isolationR03().sumPt;
    double muon1EmEtR03 = MuonVector[0]->isolationR03().emEt;
    double muon1HadEtR03 = MuonVector[0]->isolationR03().hadEt;
    double muon1SumPtR05 = MuonVector[0]->isolationR05().sumPt;
    double muon1EmEtR05 = MuonVector[0]->isolationR05().emEt;
    double muon1HadEtR05 = MuonVector[0]->isolationR05().hadEt;

    double relIsoFirstMuonDr03 = (muon1SumPtR03 + muon1EmEtR03 + muon1HadEtR03)/MuonVector[0]->pt();
    double relIsoFirstMuonDr05 = (muon1SumPtR05 + muon1EmEtR05 + muon1HadEtR05)/MuonVector[0]->pt();

    eventData.SetMuonsN(MuonVector.size());
    eventData.SetLeadingMuonPt(MuonVector[0]->pt());
    eventData.SetLeadingMuonEta(MuonVector[0]->eta());
    eventData.SetLeadingMuonPhi(MuonVector[0]->phi());
    eventData.SetLeadingMuonCharge(MuonVector[0]->charge());
    eventData.SetLeadingMuonP4(MuonVector[0]->p4());

    eventData.SetLeadingMuonSumPtR03(muon1SumPtR03);
    eventData.SetLeadingMuonEmEtR03(muon1EmEtR03);
    eventData.SetLeadingMuonHadEtR03(muon1HadEtR03);
    eventData.SetLeadingMuonSumPtR05(muon1SumPtR05);
    eventData.SetLeadingMuonEmEtR05(muon1EmEtR05);
    eventData.SetLeadingMuonHadEtR05(muon1HadEtR05);

    eventData.SetLeadingMuonrelIsoDr03(relIsoFirstMuonDr03);
    eventData.SetLeadingMuonrelIsoDr05(relIsoFirstMuonDr05);

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCountLeadingMuonR03 = 0;
    int goodTracksCountLeadingMuonR04 = 0;
    int goodTracksCountLeadingMuonR05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),MuonVector[0]->eta(),MuonVector[0]->phi()) > 0.3))
      {
	goodTracksCountLeadingMuonR03++;
      }

      if ((deltaR(track->eta(),track->phi(),MuonVector[0]->eta(),MuonVector[0]->phi()) > 0.4))
      {
	goodTracksCountLeadingMuonR04++;
      }

      if ((deltaR(track->eta(),track->phi(),MuonVector[0]->eta(),MuonVector[0]->phi()) > 0.5))
      {
	goodTracksCountLeadingMuonR05++;
      }

    }

    eventData.SetTracksNonConeLeadingMuonR03(goodTracksCountLeadingMuonR03);
    eventData.SetTracksNonConeLeadingMuonR04(goodTracksCountLeadingMuonR04);
    eventData.SetTracksNonConeLeadingMuonR05(goodTracksCountLeadingMuonR05);

    if (debug){
      std::cout << "NMuons: " << MuonVector.size() << std::endl;
      std::cout << "Muon, pT 1: " << MuonVector[0]->pt() << std::endl;
      std::cout << "Muon, eta 1: " << MuonVector[0]->eta() << std::endl;
      if(!MuonVector[0]->track().isNull()) std::cout << "Muon 1, Tracker Hits: " << MuonVector[0]->track()->hitPattern().trackerLayersWithMeasurement() << std::endl;
      if(!MuonVector[0]->globalTrack().isNull()) std::cout << "Muon 1, Chi2/ndof: " << MuonVector[0]->globalTrack()->normalizedChi2() << std::endl;
      std::cout << "Muon 1, Matched Stations: " << MuonVector[0]->numberOfMatchedStations() << std::endl;
      if(!MuonVector[0]->innerTrack().isNull()){
	std::cout << "Muon 1, dxy: " << MuonVector[0]->innerTrack()->dxy(vertex->at(0).position()) << std::endl;
	std::cout << "Muon 1, Number of Valid Pixel Hits: " << MuonVector[0]->innerTrack()->hitPattern().numberOfValidPixelHits() << std::endl;
      }
    }

    if(MuonVector[0]->isGlobalMuon() && MuonVector[0]->isTrackerMuon()){
      if(!MuonVector[0]->globalTrack().isNull() && MuonVector[0]->globalTrack()->hitPattern().numberOfValidMuonHits()>0){
	if(!MuonVector[0]->track().isNull() && MuonVector[0]->track()->hitPattern().trackerLayersWithMeasurement() > 10){
	  if(!MuonVector[0]->innerTrack().isNull() && MuonVector[0]->innerTrack()->hitPattern().numberOfValidPixelHits() > 0){
	    if(MuonVector[0]->numberOfMatchedStations() > 1 && MuonVector[0]->globalTrack()->normalizedChi2() < 10. && fabs(MuonVector[0]->innerTrack()->dxy(vertex->at(0).position())) < 0.2){
	      eventData.SetLeadingMuonIsGood(true);
	    }
	  }
	}
      }
    }else{
      eventData.SetLeadingMuonIsGood(false);
    }

    if(!MuonVector[0]->track().isNull()){
      eventData.SetLeadingMuonTrackerHits(MuonVector[0]->track()->hitPattern().trackerLayersWithMeasurement());
    }else{
      eventData.SetLeadingMuonTrackerHits(-999.);
    }
    if(!MuonVector[0]->innerTrack().isNull()){
      eventData.SetLeadingMuonPixelHits(MuonVector[0]->innerTrack()->hitPattern().numberOfValidPixelHits());
      eventData.SetLeadingMuonDxy(MuonVector[0]->innerTrack()->dxy(vertex->at(0).position()));
      eventData.SetLeadingMuonDz(MuonVector[0]->innerTrack()->dz(vertex->at(0).position()));
    }else{
      eventData.SetLeadingMuonPixelHits(-999.);
      eventData.SetLeadingMuonDxy(-999.);
      eventData.SetLeadingMuonDz(-999.);
    }
    if(!MuonVector[0]->globalTrack().isNull()){
      eventData.SetLeadingMuonNormalizedChi2(MuonVector[0]->globalTrack()->normalizedChi2());
    }else{
      eventData.SetLeadingMuonNormalizedChi2(-999.);
    }
    eventData.SetLeadingMuonMatchedStations(MuonVector[0]->numberOfMatchedStations());
    eventData.SetLeadingMuonIsGlobal(MuonVector[0]->isGlobalMuon());
    eventData.SetLeadingMuonIsTracker(MuonVector[0]->isTrackerMuon());

  }
  else{
    eventData.SetMuonsN(-1);
    eventData.SetLeadingMuonPt(-999.);
    eventData.SetLeadingMuonEta(-999.);
    eventData.SetLeadingMuonPhi(-999.);
    eventData.SetLeadingMuonCharge(-999);
    eventData.SetLeadingMuonSumPtR03(-999.);
    eventData.SetLeadingMuonEmEtR03(-999.);
    eventData.SetLeadingMuonHadEtR03(-999.);
    eventData.SetLeadingMuonSumPtR05(-999.);
    eventData.SetLeadingMuonEmEtR05(-999.);
    eventData.SetLeadingMuonHadEtR05(-999.);
    eventData.SetLeadingMuonrelIsoDr03(-999.);
    eventData.SetLeadingMuonrelIsoDr05(-999.);
    eventData.SetLeadingMuonTrackerHits(-999.);
    eventData.SetLeadingMuonPixelHits(-999.);
    eventData.SetLeadingMuonNormalizedChi2(-999.);
    eventData.SetLeadingMuonMatchedStations(-999.);
    eventData.SetLeadingMuonDxy(-999.);
    eventData.SetLeadingMuonDz(-999.);
    eventData.SetLeadingMuonIsGlobal(false);
    eventData.SetLeadingMuonIsTracker(false);
    eventData.SetLeadingMuonIsGood(false);
    eventData.SetTracksNonConeLeadingMuonR03(-999);
    eventData.SetTracksNonConeLeadingMuonR04(-999);
    eventData.SetTracksNonConeLeadingMuonR05(-999);
  } 

  // Second Muon
  if(MuonVector.size()>1){

    double muon2SumPtR03 = MuonVector[1]->isolationR03().sumPt;
    double muon2EmEtR03 = MuonVector[1]->isolationR03().emEt;
    double muon2HadEtR03 = MuonVector[1]->isolationR03().hadEt;
    double muon2SumPtR05 = MuonVector[1]->isolationR05().sumPt;
    double muon2EmEtR05 = MuonVector[1]->isolationR05().emEt;
    double muon2HadEtR05 = MuonVector[1]->isolationR05().hadEt;

    double relIsoSecondMuonDr03 = (muon2SumPtR03 + muon2EmEtR03 + muon2HadEtR03)/MuonVector[1]->pt();
    double relIsoSecondMuonDr05 = (muon2SumPtR05 + muon2EmEtR05 + muon2HadEtR05)/MuonVector[1]->pt();

    eventData.SetSecondMuonPt(MuonVector[1]->pt());
    eventData.SetSecondMuonEta(MuonVector[1]->eta());
    eventData.SetSecondMuonPhi(MuonVector[1]->phi());
    eventData.SetSecondMuonCharge(MuonVector[1]->charge());
    eventData.SetSecondMuonP4(MuonVector[1]->p4());
    eventData.SetSecondMuonSumPtR03(muon2SumPtR03);
    eventData.SetSecondMuonEmEtR03(muon2EmEtR03);
    eventData.SetSecondMuonHadEtR03(muon2HadEtR03);
    eventData.SetSecondMuonSumPtR05(muon2SumPtR05);
    eventData.SetSecondMuonEmEtR05(muon2EmEtR05);
    eventData.SetSecondMuonHadEtR05(muon2HadEtR05);
    eventData.SetSecondMuonrelIsoDr03(relIsoSecondMuonDr03);
    eventData.SetSecondMuonrelIsoDr05(relIsoSecondMuonDr05);

    edm::Handle<edm::View<reco::Track> > trackHandle;
    event.getByLabel(trackTag_,trackHandle);
    const edm::View<reco::Track>& trackColl = *(trackHandle.product());

    int goodTracksCountSecondMuonR03 = 0;
    int goodTracksCountSecondMuonR04 = 0;
    int goodTracksCountSecondMuonR05 = 0;

    // Tracks Outside Cone
    edm::View<reco::Track>::const_iterator track = trackColl.begin();
    edm::View<reco::Track>::const_iterator tracks_end = trackColl.end();
    for (; track != tracks_end; ++track)
    {
      if ((deltaR(track->eta(),track->phi(),MuonVector[1]->eta(),MuonVector[1]->phi()) > 0.3))
      {
	goodTracksCountSecondMuonR03++;
      }
      if ((deltaR(track->eta(),track->phi(),MuonVector[1]->eta(),MuonVector[1]->phi()) > 0.4))
      {
	goodTracksCountSecondMuonR04++;
      }
      if ((deltaR(track->eta(),track->phi(),MuonVector[1]->eta(),MuonVector[1]->phi()) > 0.5))
      {
	goodTracksCountSecondMuonR05++;
      }
    }
    eventData.SetTracksNonConeSecondMuonR03(goodTracksCountSecondMuonR03);
    eventData.SetTracksNonConeSecondMuonR04(goodTracksCountSecondMuonR04);
    eventData.SetTracksNonConeSecondMuonR05(goodTracksCountSecondMuonR05);

    if (debug){
      std::cout << "NMuons: " << MuonVector.size() << std::endl;
      std::cout << "Muon, pT 2: " << MuonVector[1]->pt() << std::endl;
      std::cout << "Muon, eta 2: " << MuonVector[1]->eta() << std::endl;
      if(!MuonVector[1]->track().isNull()) std::cout << "Muon 2, Tracker Hits: " << MuonVector[1]->track()->hitPattern().trackerLayersWithMeasurement() << std::endl;
      if(!MuonVector[1]->globalTrack().isNull()) std::cout << "Muon 2, Chi2/ndof: " << MuonVector[1]->globalTrack()->normalizedChi2() << std::endl;
      std::cout << "Muon 2, Matched Stations: " << MuonVector[1]->numberOfMatchedStations() << std::endl;
      if(!MuonVector[1]->innerTrack().isNull()){
	std::cout << "Muon 2, dxy: " << MuonVector[1]->innerTrack()->dxy(vertex->at(0).position()) << std::endl;
	std::cout << "Muon 2, Number of Valid Pixel Hits: " << MuonVector[1]->innerTrack()->hitPattern().numberOfValidPixelHits() << std::endl;
      }

    }
    if(MuonVector[1]->isGlobalMuon() && MuonVector[1]->isTrackerMuon()){
      if(!MuonVector[1]->globalTrack().isNull() && MuonVector[1]->globalTrack()->hitPattern().numberOfValidMuonHits()>0){
	if(!MuonVector[1]->track().isNull() && MuonVector[1]->track()->hitPattern().trackerLayersWithMeasurement() > 10){
	  if(!MuonVector[1]->innerTrack().isNull() && MuonVector[1]->innerTrack()->hitPattern().numberOfValidPixelHits() > 0){
	    if(MuonVector[1]->numberOfMatchedStations() > 1 && MuonVector[1]->globalTrack()->normalizedChi2() < 10. && fabs(MuonVector[1]->innerTrack()->dxy(vertex->at(0).position())) < 0.2){
	      eventData.SetSecondMuonIsGood(true);
	    }
	  }
	}
      }
    }else{
      eventData.SetSecondMuonIsGood(false);
    }
    if(!MuonVector[1]->track().isNull()){
      eventData.SetSecondMuonTrackerHits(MuonVector[1]->track()->hitPattern().trackerLayersWithMeasurement());
    }else{
      eventData.SetSecondMuonTrackerHits(-999.);
    }
    if(!MuonVector[1]->innerTrack().isNull()){
      eventData.SetSecondMuonPixelHits(MuonVector[1]->innerTrack()->hitPattern().numberOfValidPixelHits());
      eventData.SetSecondMuonDxy(MuonVector[1]->innerTrack()->dxy(vertex->at(0).position()));
      eventData.SetSecondMuonDz(MuonVector[1]->innerTrack()->dz(vertex->at(0).position()));
    }else{
      eventData.SetSecondMuonPixelHits(-999.);
      eventData.SetSecondMuonDxy(-999.);
      eventData.SetSecondMuonDz(-999.);
    }
    if(!MuonVector[1]->globalTrack().isNull()){
      eventData.SetSecondMuonNormalizedChi2(MuonVector[1]->globalTrack()->normalizedChi2());
    }else{
      eventData.SetSecondMuonNormalizedChi2(-999.);
    }
    eventData.SetSecondMuonMatchedStations(MuonVector[1]->numberOfMatchedStations());
    eventData.SetSecondMuonIsGlobal(MuonVector[1]->isGlobalMuon());
    eventData.SetSecondMuonIsTracker(MuonVector[1]->isTrackerMuon());
  }
  else{
    eventData.SetSecondMuonPt(-999.);
    eventData.SetSecondMuonEta(-999.);
    eventData.SetSecondMuonPhi(-999.);
    eventData.SetSecondMuonCharge(-999);
    eventData.SetSecondMuonSumPtR03(-999.);
    eventData.SetSecondMuonEmEtR03(-999.);
    eventData.SetSecondMuonHadEtR03(-999.);
    eventData.SetSecondMuonSumPtR05(-999.);
    eventData.SetSecondMuonEmEtR05(-999.);
    eventData.SetSecondMuonHadEtR05(-999.);
    eventData.SetSecondMuonrelIsoDr03(-999.);
    eventData.SetSecondMuonrelIsoDr05(-999.);
    eventData.SetSecondMuonTrackerHits(-999.);
    eventData.SetSecondMuonPixelHits(-999.);
    eventData.SetSecondMuonNormalizedChi2(-999.);
    eventData.SetSecondMuonMatchedStations(-999.);
    eventData.SetSecondMuonDxy(-999.);
    eventData.SetSecondMuonDz(-999.);
    eventData.SetSecondMuonIsGlobal(false);
    eventData.SetSecondMuonIsTracker(false);
    eventData.SetSecondMuonIsGood(false);
    eventData.SetTracksNonConeSecondMuonR03(-999);
    eventData.SetTracksNonConeSecondMuonR04(-999);
    eventData.SetTracksNonConeSecondMuonR05(-999);
  } 

}

// Templates
////////////

template <class T, class W>
math::XYZTLorentzVector TriggerAnalysis::DiSystem(T obj1, W obj2){
  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1->p4();
  DiObj += obj2->p4();
  return DiObj;
}

template <class T, class W>
double TriggerAnalysis::TransverseMass(T lepton, W met){

  //double tmass_=-999.;

  //double metT = sqrt(pow(met->px(),2) + pow(met->py(),2));
  //double lepT = sqrt(pow(lepton->px(),2) + pow(lepton->py(),2));
  //reco::Particle::LorentzVector TS = lepton->p4() + met->p4();
  //tmass_ = sqrt(pow(metT+lepT,2) - (TS.px()*TS.px()) - (TS.py()*TS.py()));

  //Defense
  //if (!std::isfinite(tmass_)) {
  //   tmass_ = -999.;
  //}

  double w_et = met->et() + lepton->pt();
  double w_px = met->px() + lepton->px();
  double w_py = met->py() + lepton->py();

  double tmass_ = w_et*w_et - w_px*w_px - w_py*w_py;
  tmass_ = (tmass_ > 0) ? sqrt(tmass_) : 0;

  return tmass_;

}
