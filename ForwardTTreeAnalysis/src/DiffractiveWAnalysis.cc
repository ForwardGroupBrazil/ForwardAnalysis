#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveWAnalysis.h"
#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/DiffractiveWEvent.h"

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

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

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

using diffractiveWAnalysis::DiffractiveWAnalysis;

const char* DiffractiveWAnalysis::name = "DiffractiveWAnalysis";

DiffractiveWAnalysis::DiffractiveWAnalysis(const edm::ParameterSet& pset):
  triggerResultsTag_(pset.getParameter<edm::InputTag>("TriggerResultsTag")),
  hltPathNames_(pset.getParameter<std::vector<std::string> >("hltPaths")),
  electronTag_(pset.getParameter<edm::InputTag>("electronTag")),
  muonTag_(pset.getParameter<edm::InputTag>("muonTag")),
  metTag_(pset.getParameter<edm::InputTag>("metTag")),
  patmetTag_(pset.getParameter<edm::InputTag>("patmetTag")),
  pfTag_(pset.getParameter<edm::InputTag>("pfTag")),
  genTag_(pset.getParameter<edm::InputTag>("genTag")),
  PVtxCollectionTag_(pset.getParameter<edm::InputTag>("PVtxCollectionTag")),
  castorHitsTag_(pset.getParameter<edm::InputTag>("castorHitsTag")),
  zdcHitsTag_(pset.getParameter<edm::InputTag>("zdcHitsTag")),
  RunCastor_(pset.getUntrackedParameter<Bool_t>("RunCastor", false)),
  RunZDC_(pset.getUntrackedParameter<Bool_t>("RunZDC", false)),
  RunMC_(pset.getUntrackedParameter<Bool_t>("RunMC", false)),
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
  trackTag_(pset.getParameter<edm::InputTag>("TrackTag"))
{
}

void DiffractiveWAnalysis::setTFileService(){

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

DiffractiveWAnalysis::~DiffractiveWAnalysis(){}

void DiffractiveWAnalysis::begin() {
  setTFileService();
}

void DiffractiveWAnalysis::begin(const edm::Run& run, const edm::EventSetup& setup) {}

void DiffractiveWAnalysis::end() {}

void DiffractiveWAnalysis::fill(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  eventData.reset();

  fillTriggerInfo(eventData,event,setup);
  fillElectronsInfo(eventData,event,setup);
  fillMuonsInfo(eventData,event,setup);
  fillMETInfo(eventData,event,setup);
  fillCollections(eventData,event,setup);
  fillTracksInfo(eventData,event,setup);
  fillDetectorVariables(eventData,event,setup);
  fillVariables(eventData,event,setup);
  if (RunMC_) fillGenInfo(eventData,event,setup); 
  if (RunCastor_){
    fillCastor(eventData,event,setup);
    fillCastorDebug(eventData,event,setup);
  }
  if (RunZDC_) fillZDC(eventData,event,setup);
  if (EachTower_) fillDetectorEnergyEtaInfo(eventData,event,setup);
}

// Fill Trigger
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillTriggerInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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


// F I L L   M U O N S   A N D   P A T :: M U O N
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillMuonsInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Fill reco::Muon
  edm::Handle<reco::MuonCollection> muons;
  event.getByLabel(muonTag_,muons);

  int muonsize = muons->size();
  int itMuon;
  MuonVector.clear();

  if(muons->size()>0){
    for(itMuon=0; itMuon < muonsize; ++itMuon){
      const reco::Muon* muonAll = &((*muons)[itMuon]);
      MuonVector.push_back(muonAll);
    }
  }

  // Fill pat::Muon
  edm::Handle<std::vector<pat::Muon> > patMuons;
  event.getByLabel("patMuons", patMuons);

  int patMuonsize = patMuons->size();
  int itpatMuon;
  PatMuonVector.clear();

  if(patMuons->size()>0){
    for(itpatMuon=0; itpatMuon < patMuonsize; ++itpatMuon){
      const pat::Muon* patMuonAll = &((*patMuons)[itpatMuon]);
      PatMuonVector.push_back(patMuonAll);    
    }
  }

}

// F I L L   E L E C T R O N    A N D   P A T :: E L E C T R O N
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillElectronsInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  // Fill reco::GsfElectron
  edm::Handle<reco::GsfElectronCollection> electrons;
  event.getByLabel(electronTag_,electrons);

  int electronsize = electrons->size();
  int itElectron;
  ElectronVector.clear();

  if(electrons->size()>0){
    for(itElectron=0; itElectron < electronsize; ++itElectron){
      const reco::GsfElectron* electronAll = &((*electrons)[itElectron]);
      ElectronVector.push_back(electronAll);
    }
  }

  // Fill pat::Electron
  edm::Handle<std::vector<pat::Electron> > patElectrons;
  event.getByLabel("patElectrons", patElectrons);

  int patElectronsize = patElectrons->size();
  int itpatElectron;
  PatElectronVector.clear();

  if(patElectrons->size()>0){
    for(itpatElectron=0; itpatElectron < patElectronsize; ++itpatElectron){
      const pat::Electron* patElectronAll = &((*patElectrons)[itpatElectron]);
      PatElectronVector.push_back(patElectronAll);
    }
  }

}

// F I L L   M E T
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillMETInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

  // Fill pat::MET
  // Apply Correction because Thresholds
  edm::Handle<std::vector<pat::MET> > patmet;
  event.getByLabel(patmetTag_,patmet);

  PatNeutrinoVector.clear();
  int patneutrinosize = patmet->size();
  int itpatmet;

  if(patmet->size()>0){
    for(itpatmet=0; itpatmet < patneutrinosize; ++itpatmet){
      const pat::MET* patneutrinoAll = &((*patmet)[itpatmet]);
      PatNeutrinoVector.push_back(patneutrinoAll);
    }
  }

}

// F I L L   C O L L E C T I O N S
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillCollections(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  if (debug) std::cout << "\n--BEGIN--" << std::endl;

  if(MuonVector.size()>0){

    // Sorting Vector by pT
    const int  MuonVectorSize = (int) MuonVector.size();
    int *sortMuonVector=   new int[MuonVectorSize];
    double *vmuon = new double[MuonVectorSize];

    for (int i=0; i<MuonVectorSize; i++) {
      vmuon[i] = MuonVector[i]->pt();
    }

    TMath::Sort(MuonVectorSize, vmuon, sortMuonVector, true);

    for (unsigned int i=0;i<MuonVector.size();i++){
      if (debug) std::cout << "ORDERED reco::Muon[" << sortMuonVector[i] << "]\t\t---> pT [GeV]: " << MuonVector[sortMuonVector[i]]->pt() << " | eT [GeV]: " << MuonVector[sortMuonVector[i]]->et() << " | eta: " << MuonVector[sortMuonVector[i]]->eta() << " | phi: " << MuonVector[sortMuonVector[i]]->phi() << " | px [GeV]: " << MuonVector[sortMuonVector[i]]->px() << " | py [GeV]: " << MuonVector[sortMuonVector[i]]->py() << std::endl;
    }

    double muon1SumPtR03 = MuonVector[sortMuonVector[0]]->isolationR03().sumPt;
    double muon1EmEtR03 = MuonVector[sortMuonVector[0]]->isolationR03().emEt;
    double muon1HadEtR03 = MuonVector[sortMuonVector[0]]->isolationR03().hadEt;
    double muon1SumPtR05 = MuonVector[sortMuonVector[0]]->isolationR05().sumPt;
    double muon1EmEtR05 = MuonVector[sortMuonVector[0]]->isolationR05().emEt;
    double muon1HadEtR05 = MuonVector[sortMuonVector[0]]->isolationR05().hadEt;

    double relIsoFirstMuonDr03 = (muon1SumPtR03 + muon1EmEtR03 + muon1HadEtR03)/MuonVector[sortMuonVector[0]]->pt();
    double relIsoFirstMuonDr05 = (muon1SumPtR05 + muon1EmEtR05 + muon1HadEtR05)/MuonVector[sortMuonVector[0]]->pt();

    eventData.SetMuonsN(MuonVectorSize);
    eventData.SetLeadingMuonPt(MuonVector[sortMuonVector[0]]->pt());
    eventData.SetLeadingMuonEta(MuonVector[sortMuonVector[0]]->eta());
    eventData.SetLeadingMuonPhi(MuonVector[sortMuonVector[0]]->phi());
    eventData.SetLeadingMuonCharge(MuonVector[sortMuonVector[0]]->charge());
    eventData.SetLeadingMuonP4(MuonVector[sortMuonVector[0]]->p4());
    eventData.SetLeadingMuonSumPtR03(muon1SumPtR03);
    eventData.SetLeadingMuonEmEtR03(muon1EmEtR03);
    eventData.SetLeadingMuonHadEtR03(muon1HadEtR03);
    eventData.SetLeadingMuonSumPtR05(muon1SumPtR05);
    eventData.SetLeadingMuonEmEtR05(muon1EmEtR05);
    eventData.SetLeadingMuonHadEtR05(muon1HadEtR05);
    eventData.SetLeadingMuonrelIsoDr03(relIsoFirstMuonDr03);
    eventData.SetLeadingMuonrelIsoDr05(relIsoFirstMuonDr05);

  }


  if(PatMuonVector.size()>0){

    // Sorting Vector by pT
    const int  PatMuonVectorSize = (int) PatMuonVector.size();
    int *sortPatMuonVector=   new int[PatMuonVectorSize];
    double *vpatmuon = new double[PatMuonVectorSize];

    for (int i=0; i<PatMuonVectorSize; i++) {
      vpatmuon[i] = PatMuonVector[i]->pt();
    }

    TMath::Sort(PatMuonVectorSize, vpatmuon, sortPatMuonVector, true);

    for (unsigned int i=0;i<PatMuonVector.size();i++){
      if (debug) std::cout << "ORDERED pat::Muon[" << sortPatMuonVector[i] << "]\t\t---> pT [GeV]: " << PatMuonVector[sortPatMuonVector[i]]->pt() << " | eT [GeV]: " << PatMuonVector[sortPatMuonVector[i]]->et() << " | eta: " << PatMuonVector[sortPatMuonVector[i]]->eta() << " | phi: " << PatMuonVector[sortPatMuonVector[i]]->phi() << " | px [GeV]: " << PatMuonVector[sortPatMuonVector[i]]->px() << " | py [GeV]: " << PatMuonVector[sortPatMuonVector[i]]->py() << std::endl;
    }

    eventData.SetPatNMuon(PatMuonVectorSize);
    eventData.SetPatMuon1Pt(PatMuonVector[sortPatMuonVector[0]]->pt());
    eventData.SetPatMuon1Charge(PatMuonVector[sortPatMuonVector[0]]->charge());
    eventData.SetPatMuon1Phi(PatMuonVector[sortPatMuonVector[0]]->phi());
    eventData.SetPatMuon1Eta(PatMuonVector[sortPatMuonVector[0]]->eta());
    eventData.SetPatMuon1Et(PatMuonVector[sortPatMuonVector[0]]->et());

    double Patmuon1SumPtR03 = PatMuonVector[sortPatMuonVector[0]]->isolationR03().sumPt;
    double Patmuon1EmEtR03 = PatMuonVector[sortPatMuonVector[0]]->isolationR03().emEt;
    double Patmuon1HadEtR03 = PatMuonVector[sortPatMuonVector[0]]->isolationR03().hadEt;
    double Patmuon1SumPtR05 = PatMuonVector[sortPatMuonVector[0]]->isolationR05().sumPt;
    double Patmuon1EmEtR05 = PatMuonVector[sortPatMuonVector[0]]->isolationR05().emEt;
    double Patmuon1HadEtR05 = PatMuonVector[sortPatMuonVector[0]]->isolationR05().hadEt; 

    eventData.SetPatMuon1SumPtR03(Patmuon1SumPtR03);
    eventData.SetPatMuon1EmEtR03(Patmuon1EmEtR03);
    eventData.SetPatMuon1HadEtR03(Patmuon1HadEtR03);
    eventData.SetPatMuon1SumPtR05(Patmuon1SumPtR05);
    eventData.SetPatMuon1EmEtR05(Patmuon1EmEtR05);
    eventData.SetPatMuon1HadEtR05(Patmuon1HadEtR05); 

    double PatrelIsoFirstMuonDr03 = (Patmuon1SumPtR03 + Patmuon1EmEtR03 + Patmuon1HadEtR03)/PatMuonVector[sortPatMuonVector[0]]->pt();
    double PatrelIsoFirstMuonDr05 = (Patmuon1SumPtR05 + Patmuon1EmEtR05 + Patmuon1HadEtR05)/PatMuonVector[sortPatMuonVector[0]]->pt();
    double PatrelIsoFirstMuon = (PatMuonVector[sortPatMuonVector[0]]->trackIso()+PatMuonVector[sortPatMuonVector[0]]->ecalIso()+PatMuonVector[sortPatMuonVector[0]]->hcalIso())/PatMuonVector[sortPatMuonVector[0]]->pt();

    eventData.SetPatMuon1relIsoDr03(PatrelIsoFirstMuonDr03);
    eventData.SetPatMuon1relIsoDr05(PatrelIsoFirstMuonDr05);
    eventData.SetPatMuon1relIso(PatrelIsoFirstMuon);

  }

  if(ElectronVector.size()>0){

    // Sorting Vector by pT
    const int  ElectronVectorSize = (int) ElectronVector.size();
    int *sortElectronVector=   new int[ElectronVectorSize];
    double *velectron = new double[ElectronVectorSize];

    for (int i=0; i<ElectronVectorSize; i++) {
      velectron[i] = ElectronVector[i]->pt();
    }

    TMath::Sort(ElectronVectorSize, velectron, sortElectronVector, true);

    for (unsigned int i=0;i<ElectronVector.size();i++){
      if (debug) std::cout << "ORDERED reco::Electron[" << sortElectronVector[i] << "]\t---> pT [GeV]: " << ElectronVector[sortElectronVector[i]]->pt() << " | eT [GeV]: " << ElectronVector[sortElectronVector[i]]->et() << " | eta: " << ElectronVector[sortElectronVector[i]]->eta() << " | phi: " << ElectronVector[sortElectronVector[i]]->phi() << " | px [GeV]: " << ElectronVector[sortElectronVector[i]]->px() << " | py [GeV]: " << ElectronVector[sortElectronVector[i]]->py() << std::endl;
    }

    // Fill Methods
    eventData.SetElectronsN(ElectronVectorSize);
    eventData.SetLeadingElectronPt(ElectronVector[sortElectronVector[0]]->pt());
    eventData.SetLeadingElectronEta(ElectronVector[sortElectronVector[0]]->eta());
    eventData.SetLeadingElectronPhi(ElectronVector[sortElectronVector[0]]->phi());
    eventData.SetLeadingElectronCharge(ElectronVector[sortElectronVector[0]]->charge());
    eventData.SetLeadingElectronP4(ElectronVector[sortElectronVector[0]]->p4());
    eventData.SetLeadingElectronTkDr03(ElectronVector[sortElectronVector[0]]->dr03TkSumPt());
    eventData.SetLeadingElectronEcalDr03(ElectronVector[sortElectronVector[0]]->dr03EcalRecHitSumEt());
    eventData.SetLeadingElectronHcalDr03(ElectronVector[sortElectronVector[0]]->dr03HcalTowerSumEt());
    eventData.SetLeadingElectronTkDr04(ElectronVector[sortElectronVector[0]]->dr04TkSumPt());
    eventData.SetLeadingElectronEcalDr04(ElectronVector[sortElectronVector[0]]->dr04EcalRecHitSumEt());
    eventData.SetLeadingElectronHcalDr04(ElectronVector[sortElectronVector[0]]->dr04HcalTowerSumEt());

    double relIsoFirstElectronDr03 = (ElectronVector[sortElectronVector[0]]->dr03TkSumPt()+ElectronVector[sortElectronVector[0]]->dr03EcalRecHitSumEt()+ElectronVector[sortElectronVector[0]]->dr03HcalTowerSumEt())/ElectronVector[sortElectronVector[0]]->et();
    double relIsoFirstElectronDr04 = (ElectronVector[sortElectronVector[0]]->dr04TkSumPt()+ElectronVector[sortElectronVector[0]]->dr04EcalRecHitSumEt()+ElectronVector[sortElectronVector[0]]->dr04HcalTowerSumEt())/ElectronVector[sortElectronVector[0]]->et();

    eventData.SetLeadingElectronrelIsoDr03(relIsoFirstElectronDr03);
    eventData.SetLeadingElectronrelIsoDr04(relIsoFirstElectronDr04);

    eventData.SetLeadingElectronDeltaPhiTkClu(ElectronVector[sortElectronVector[0]]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetLeadingElectronDeltaEtaTkClu(ElectronVector[sortElectronVector[0]]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetLeadingElectronSigmaIeIe(ElectronVector[sortElectronVector[0]]->sigmaIetaIeta());
    eventData.SetLeadingElectronDCot(ElectronVector[sortElectronVector[0]]->convDcot());
    eventData.SetLeadingElectronDist(ElectronVector[sortElectronVector[0]]->convDist());
    eventData.SetLeadingElectronInnerHits(ElectronVector[sortElectronVector[0]]->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
    eventData.SetLeadingElectronHE(ElectronVector[sortElectronVector[0]]->hadronicOverEm());

  }


  if(PatElectronVector.size()>0){

    // Sorting Vector by pT
    const int  PatElectronVectorSize = (int) PatElectronVector.size();
    int *sortPatElectronVector=   new int[PatElectronVectorSize];
    double *vpatelectron = new double[PatElectronVectorSize];

    for (int i=0; i<PatElectronVectorSize; i++) {
      vpatelectron[i] = PatElectronVector[i]->pt();
    }

    TMath::Sort(PatElectronVectorSize, vpatelectron, sortPatElectronVector, true);

    for (unsigned int i=0;i<PatElectronVector.size();i++){
      if (debug) std::cout << "ORDERED pat::Electron[" << sortPatElectronVector[i] << "]\t---> pT [GeV]: " << PatElectronVector[sortPatElectronVector[i]]->pt() << " | eT [GeV]: " << PatElectronVector[sortPatElectronVector[i]]->et() << " | eta: " << PatElectronVector[sortPatElectronVector[i]]->eta() << " | phi: " << PatElectronVector[sortPatElectronVector[i]]->phi() << " | px [GeV]: " << PatElectronVector[sortPatElectronVector[i]]->px() << " | py [GeV]: " << PatElectronVector[sortPatElectronVector[i]]->py() << std::endl;
    }


    // Fill Methods
    eventData.SetPatNElectron(PatElectronVectorSize);
    eventData.SetPatElectron1Pt(PatElectronVector[sortPatElectronVector[0]]->pt());
    eventData.SetPatElectron1Charge(PatElectronVector[sortPatElectronVector[0]]->charge());
    eventData.SetPatElectron1Phi(PatElectronVector[sortPatElectronVector[0]]->phi());
    eventData.SetPatElectron1Eta(PatElectronVector[sortPatElectronVector[0]]->eta());
    eventData.SetPatElectron1Et(PatElectronVector[sortPatElectronVector[0]]->et());
    eventData.SetPatElectron1TkDr03(PatElectronVector[sortPatElectronVector[0]]->dr03TkSumPt());
    eventData.SetPatElectron1EcalDr03(PatElectronVector[sortPatElectronVector[0]]->dr03EcalRecHitSumEt());
    eventData.SetPatElectron1HcalDr03(PatElectronVector[sortPatElectronVector[0]]->dr03HcalTowerSumEt());
    eventData.SetPatElectron1TkDr04(PatElectronVector[sortPatElectronVector[0]]->dr04TkSumPt());
    eventData.SetPatElectron1EcalDr04(PatElectronVector[sortPatElectronVector[0]]->dr04EcalRecHitSumEt());
    eventData.SetPatElectron1HcalDr04(PatElectronVector[sortPatElectronVector[0]]->dr04HcalTowerSumEt());

    double PatrelIsoFirstElectronDr03 = (PatElectronVector[sortPatElectronVector[0]]->dr03TkSumPt()+PatElectronVector[sortPatElectronVector[0]]->dr03EcalRecHitSumEt()+PatElectronVector[sortPatElectronVector[0]]->dr03HcalTowerSumEt())/PatElectronVector[sortPatElectronVector[0]]->et();
    double PatrelIsoFirstElectronDr04 = (PatElectronVector[sortPatElectronVector[0]]->dr04TkSumPt()+PatElectronVector[sortPatElectronVector[0]]->dr04EcalRecHitSumEt()+PatElectronVector[sortPatElectronVector[0]]->dr04HcalTowerSumEt())/PatElectronVector[sortPatElectronVector[0]]->et();

    eventData.SetPatElectron1relIsoDr03(PatrelIsoFirstElectronDr03);
    eventData.SetPatElectron1relIsoDr04(PatrelIsoFirstElectronDr04);

    eventData.SetPatElectron1DeltaPhiTkClu(PatElectronVector[sortPatElectronVector[0]]->deltaPhiSuperClusterTrackAtVtx());
    eventData.SetPatElectron1DeltaEtaTkClu(PatElectronVector[sortPatElectronVector[0]]->deltaEtaSuperClusterTrackAtVtx());
    eventData.SetPatElectron1SigmaIeIe(PatElectronVector[sortPatElectronVector[0]]->sigmaIetaIeta());
    eventData.SetPatElectron1DCot(PatElectronVector[sortPatElectronVector[0]]->convDcot());
    eventData.SetPatElectron1Dist(PatElectronVector[sortPatElectronVector[0]]->convDist());
    eventData.SetPatElectron1InnerHits(PatElectronVector[sortPatElectronVector[0]]->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
    eventData.SetPatElectron1HE(PatElectronVector[sortPatElectronVector[0]]->hadronicOverEm());

  }


  if(NeutrinoVector.size()>0){

    // Sorting Vector by pT
    const int  NeutrinoVectorSize = (int) NeutrinoVector.size();
    int *sortNeutrinoVector=   new int[NeutrinoVectorSize];
    double *vneutrino = new double[NeutrinoVectorSize];

    for (int i=0; i<NeutrinoVectorSize; i++) {
      vneutrino[i] = NeutrinoVector[i]->pt();
    }

    TMath::Sort(NeutrinoVectorSize, vneutrino, sortNeutrinoVector, true);

    for (unsigned int i=0;i<NeutrinoVector.size();i++){
      if (debug) std::cout << "ORDERED reco::pfMET[" << sortNeutrinoVector[i] << "]\t\t---> pT [GeV]: " << NeutrinoVector[sortNeutrinoVector[i]]->pt() << " | eT [GeV]: " << NeutrinoVector[sortNeutrinoVector[i]]->et() << " | sum eT [GeV]: " << NeutrinoVector[sortNeutrinoVector[i]]->sumEt() << " | eta: " << NeutrinoVector[sortNeutrinoVector[i]]->eta() << " | phi: " << NeutrinoVector[sortNeutrinoVector[i]]->phi() << " | px [GeV]: " << NeutrinoVector[sortNeutrinoVector[i]]->px() << " | py [GeV]: " << NeutrinoVector[sortNeutrinoVector[i]]->py() << std::endl;
    }

    eventData.SetMETPt(NeutrinoVector[sortNeutrinoVector[0]]->pt());
    eventData.SetMETPhi(NeutrinoVector[sortNeutrinoVector[0]]->phi());
    eventData.SetMETEt(NeutrinoVector[sortNeutrinoVector[0]]->et());
    eventData.SetMETSumEt(NeutrinoVector[sortNeutrinoVector[0]]->sumEt());
    eventData.SetMETpx(NeutrinoVector[sortNeutrinoVector[0]]->px());
    eventData.SetMETpy(NeutrinoVector[sortNeutrinoVector[0]]->py());

  }


  if(PatNeutrinoVector.size()>0){

    // Sorting Vector by pT
    const int  PatNeutrinoVectorSize = (int) PatNeutrinoVector.size();
    int *sortPatNeutrinoVector=   new int[PatNeutrinoVectorSize];
    double *vpatneutrino = new double[PatNeutrinoVectorSize];

    for (int i=0; i<PatNeutrinoVectorSize; i++) {
      vpatneutrino[i] = PatNeutrinoVector[i]->pt();
    }

    TMath::Sort(PatNeutrinoVectorSize, vpatneutrino, sortPatNeutrinoVector, true);

    for (unsigned int i=0;i<PatNeutrinoVector.size();i++){
      if (debug) std::cout << "ORDERED pat::MET[" << sortPatNeutrinoVector[i] << "]\t\t---> pT [GeV]: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->pt() << " | eT [GeV]: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->et() << " | sum eT [GeV]: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->sumEt() << " | eta: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->eta() << " | phi: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->phi() << " | px [GeV]: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->px() << " | py [GeV]: " << PatNeutrinoVector[sortPatNeutrinoVector[i]]->py() << std::endl;
    }

    eventData.SetPatMETPt(PatNeutrinoVector[sortPatNeutrinoVector[0]]->pt());
    eventData.SetPatMETPhi(PatNeutrinoVector[sortPatNeutrinoVector[0]]->phi());
    eventData.SetPatMETEt(PatNeutrinoVector[sortPatNeutrinoVector[0]]->et());
    eventData.SetPatMETSumEt(PatNeutrinoVector[sortPatNeutrinoVector[0]]->sumEt());
    eventData.SetPatMETpx(PatNeutrinoVector[sortPatNeutrinoVector[0]]->px());
    eventData.SetPatMETpy(PatNeutrinoVector[sortPatNeutrinoVector[0]]->py());

  }

  if (debug) std::cout << "--END--\n" << std::endl;


}


// Fill Tracks Info
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillTracksInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  std::vector<double> vertexNDOF;
  std::vector<double> vertexChiNorm;
  std::vector<double> vertexMultiplicity;
  double nhit=0;
  std::vector<double> V_x;
  std::vector<double> V_y;
  std::vector<double> V_z; 
  std::vector<double> pt;
  std::vector<double> tracks;
  std::vector<std::vector<double> > tracksPT;

  edm::Handle<reco::VertexCollection>  vertexCollectionHandle;
  event.getByLabel(PVtxCollectionTag_, vertexCollectionHandle);

  for(reco::VertexCollection::const_iterator vtx = vertexCollectionHandle->begin();vtx!=vertexCollectionHandle->end(); ++vtx)
  {
    reco::Vertex::trackRef_iterator it = vtx->tracks_begin();
    reco::Vertex::trackRef_iterator lastTrack = vtx->tracks_end();


    for(;it!=lastTrack;it++) {
      nhit+=(*it)->numberOfValidHits();
      pt.push_back((*it)->pt());
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

void DiffractiveWAnalysis::fillGenInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug=false;

  //Variable declaration
  int count=0;
  int Nstable_gen=0;
  math::XYZTLorentzVector part(0.,0.,0.,0.);
  math::XYZTLorentzVector partVis(0.,0.,0.,0.);
  math::XYZTLorentzVector partZ(0.,0.,0.,0.);
  double sumECastor_minus_gen=0;
  double sumECastor_plus_gen=0;
  double sumEZDC_minus_gen=0;
  double sumEZDC_plus_gen=0;
  double etaOutcomingProton=0;
  double energyOutcomingProton=0;
  double mostEnergeticXL=0;
  double mostEnergeticXLNum=0;
  std::vector<double> eta_gen_vec;
  double xi_Z_gen_minus=0;
  double xi_Z_gen_plus=0;
  double etaZ_gen=0;
  double energyZ_gen=0;
  double p_diss_mass=0;
  double p_diss=0;
  double xL_p_diss=0;

  double xL_etaGTP5=0;
  double xL_etaLTM5=0;
  int xL_LTM5Num=0;
  int xL_GTP5Num=0;

  std::vector<double> genpt;
  std::vector<double> tracks;
  std::vector<double> eta;
  std::vector<double> etaPT;
  std::vector<double> tracksPT;

  edm::Handle<reco::GenParticleCollection> genParticles;     
  event.getByLabel(genTag_,genParticles);  // standard PYTHIA collection

  for(size_t i = 0; i < genParticles->size(); ++ i) {

    const reco::GenParticle & p = (*genParticles)[i];

    int pdg = p.pdgId();
    int status = p.status();  	 
    double eta_gen = p.eta();
    double part_pt = p.pt();
    double ener_gen= p.energy();
    double px_gen  = p.px();
    double py_gen  = p.py();
    double pz_gen  = p.pz();
    double mass_gen= p.mass();
    int motherId=0;

    if (  p.mother() ) motherId =  p.mother()->pdgId();

    math::XYZTLorentzVector tmp( px_gen ,py_gen , pz_gen ,ener_gen );

    if (fabs(pdg)==11 || fabs(pdg)==22){ 	    
      if (motherId==23) {
	partZ+=tmp;
      }
    }

    if (count==2) {    /// only works for MC Pompyt dissociative
      p_diss_mass= mass_gen;
      p_diss= pz_gen;
      if ( pdg == 2212){
	etaOutcomingProton= eta_gen;
	energyOutcomingProton= ener_gen;
      }
    }

    if (( pdg == 23)){
      xi_Z_gen_minus=( ener_gen - pz_gen )/7000;
      xi_Z_gen_plus=( ener_gen + pz_gen )/7000;
      etaZ_gen=eta_gen;
      energyZ_gen= ener_gen;
    }

    if (status == 1) 
    {
      //// vector to find gaps (cut at 1 GeV in energy)
      if  ( ( fabs(eta_gen) <= 1.5  && ener_gen > energyPFThresholdBar_ )  ||
	  (fabs(eta_gen) > 1.5 && fabs(eta_gen) <= 3 && ener_gen > energyPFThresholdEnd_) ||
	  (fabs(eta_gen) > 3 && ener_gen >energyPFThresholdHF_)  ) {

	eta_gen_vec.push_back( eta_gen);
      }

      if (  count>2) {   
	part+=tmp;
      }

      if (  (fabs(eta_gen) < 4.7) && (part_pt > 0.10) ) {   // if particle has a chance to reach the detector ...
	partVis+=tmp;
      }

      // new xL_gen definition (after Sasha)
      if (count>=2 )
      {
	if (eta_gen > 4.7)  
	{
	  xL_etaGTP5 += pz_gen;
	  xL_GTP5Num++;
	}
	if (eta_gen < -4.7)  
	{
	  xL_etaLTM5 += pz_gen;
	  xL_LTM5Num++;
	}
      }

      if (count>=2 ){
	if (p_diss>0) {
	  if ( xL_p_diss < pz_gen ){
	    xL_p_diss= pz_gen;
	  }
	}
	if (p_diss<0) {
	  if ( xL_p_diss > pz_gen ){
	    xL_p_diss= pz_gen;
	  }
	}
      }

      if ( fabs(eta_gen)>5.2 && fabs(eta_gen)<6.6 ){
	if (debug) std::cout<<"Particle in Castor, having eta "<<eta_gen<<" and energy "<< ener_gen<<std::endl;
	if (eta_gen<0) sumECastor_minus_gen += ener_gen;
	if (eta_gen>0) sumECastor_plus_gen += ener_gen;
      }

      if ( fabs(eta_gen)>8.2  && ( pdg == 2112 || pdg == 22) ){
	if (debug) std::cout<<"Particle in ZDC, having eta "<<eta_gen<<" and energy "<< ener_gen<<std::endl;
	if (eta_gen<0) sumEZDC_minus_gen += ener_gen;
	if (eta_gen>0) sumEZDC_plus_gen += ener_gen;
      }      
      Nstable_gen++;

    }  // status =1
    count++;

  } // loop over particles

  //// Computing GAPs
  const int  size = (int) eta_gen_vec.size();
  int *sortedgen=   new int[size];
  double *vgen = new double[size];
  double eta_gap_limplus_gen = -10.0;
  double eta_gap_limminus_gen = -10.0;

  for (int i=0; i<size; i++) {
    vgen[i] = eta_gen_vec[i];
    if (debug) std::cout<<vgen[i]<<std::endl;
  }
  TMath::Sort(size, vgen, sortedgen, true);

  if (size > 1) {
    double *diff = new double[size-1];
    int *diffsorted = new int[size-1];
    for (int i=0; i<(size-1); i++) {
      diff[i] = fabs(eta_gen_vec[sortedgen[i+1]]-eta_gen_vec[sortedgen[i]]);
      if (debug) {
	std::cout<<" eta " << i << " size " << size << " diff "<< diff[i]<< std::endl;
	std::cout<<" GEN etas "  << " = " << eta_gen_vec[sortedgen[i+1]] << " - " <<  eta_gen_vec[sortedgen[i]] <<  " GAP diff "<< diff[i] << std::endl;
	std::cout<<" GEN etas "  << " = " << eta_gen_vec[sortedgen[i]] << std::endl;
      }
    }

    TMath::Sort(size-1, diff, diffsorted, true);

    //checking the max gap
    double max_eta_gap_gen=diff[diffsorted[0]];
    eta_gap_limminus_gen = eta_gen_vec[sortedgen[diffsorted[0]+1]] ;
    eta_gap_limplus_gen = eta_gen_vec[sortedgen[diffsorted[0]]] ;

    if (debug) std::cout << "GEN eta ranges " <<  eta_gap_limplus_gen  << " " <<  eta_gap_limminus_gen  << std::endl;
    eventData.SetPrimaryGapMaxGen(max_eta_gap_gen);

    if (size>2) {
      double max_second_eta_gap_gen=diff[diffsorted[1]];
      eventData.SetSecondGapMaxGen(max_second_eta_gap_gen);
      if (debug) std::cout<<" diff  " << diff[diffsorted[0]] << " sec " << diff[diffsorted[1]] << " diff size "<< diff[size-2] << std::endl;
    }

    delete [] diff;
    delete [] diffsorted;
  }

  delete [] sortedgen;
  delete [] vgen;

  math::XYZTLorentzVector dataMassG_plus(0.,0.,0.,0.);
  math::XYZTLorentzVector dataMassG_minus(0.,0.,0.,0.);
  int nplusG =0;
  int nminusG =0;
  int numseltracks =0;

  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    int status = p.status();  	 
    double eta_gen = p.eta();
    int charge = p.charge();
    double pt = p.pt();

    if( status == 1 && p.energy() > 1 ) {

      math::XYZTLorentzVector tmp( p.px(),p.py(),p.pz(),p.energy());

      if ( eta_gen >= eta_gap_limplus_gen ) {
	dataMassG_plus+=tmp;
	nplusG++;
      }
      else {
	dataMassG_minus+=tmp;
	nminusG++;
      }
    }

    if ( status == 1 ) {
      if ( charge && fabs(eta_gen)<2.6 &&  pt >= 0.1 ) {  // !!  condition for xsec per 3 charged prompt particles
	numseltracks++;
	genpt.push_back(pt);
	eta.push_back(eta_gen);
      }
    }
  } // end of genparticle loop

  double Mx2_gen=partVis.M2(); /// massaquadro visibile generata
  math::XYZTLorentzVector NOZ=partVis-partZ;
  double Mx2_NOZ_gen=NOZ.M2();
  if (debug) {
    std::cout << "Mx2_gen is "<< Mx2_gen<<" while eta of the outcoming proton is "<< etaOutcomingProton <<" and the energy "<< energyOutcomingProton << std::endl;
  }

  mostEnergeticXL = xL_etaGTP5/3500.;
  mostEnergeticXLNum = xL_GTP5Num ;
  if (fabs(xL_etaGTP5)<fabs(xL_etaLTM5)) 
  {
    mostEnergeticXL = xL_etaLTM5/3500.;
    mostEnergeticXLNum = xL_LTM5Num ;
  }

  if (debug) std::cout << "* XLgen " << mostEnergeticXL << " num " << mostEnergeticXLNum << " + " << xL_etaGTP5 << " - " << xL_etaLTM5 <<  std::endl;

  const int  size2 = (int) genpt.size();
  int *sorted = new int[size2];
  double *vv = new double[size2];
  for (int i=0; i<size2; i++) {
    vv[i] = genpt[i];
  }
  TMath::Sort(size2, vv, sorted, true);
  for (int i=0; i<size2; i++) {
    tracks.push_back(genpt[sorted[i]]);
    etaPT.push_back(eta[sorted[i]]);
    if (i>30) break;
  }  //  comes out size of 32!

  eventData.SetTracksPtGen(tracks);
  eventData.SetEtaOfTracksPtGen(etaPT);
  eventData.SetNTracksGen(tracks.size());

  genpt.clear();
  eta.clear();

  delete [] sorted;
  delete [] vv;

  eventData.SetMx2PlusGen(dataMassG_plus.M2());
  eventData.SetMx2MinusGen(dataMassG_minus.M2());
  eventData.SetMx2Gen(Mx2_NOZ_gen);
  eventData.SetMx2ZGen(Mx2_gen);
  eventData.SetNMx2PlusGen(nplusG);
  eventData.SetNMx2MinusGen(nminusG);
  eventData.SetEtaGaplimPlusGen(eta_gap_limplus_gen);
  eventData.SetEtaGaplimMinusGen(eta_gap_limminus_gen);
  eventData.SetNParticlesGen(Nstable_gen);
  eventData.SetsumECastorMinusGen(sumECastor_minus_gen);
  eventData.SetsumECastorPlusGen(sumECastor_plus_gen);
  eventData.SetsumEZDCMinusGen(sumEZDC_minus_gen);
  eventData.SetsumEZDCPlusGen(sumEZDC_plus_gen);
  eventData.SetEtaOutcomingProtonGen(etaOutcomingProton);
  eventData.SetxLGen(mostEnergeticXL);
  eventData.SetxLMostEnergeticGen(mostEnergeticXLNum);
  eventData.SetxiZMinusGen(xi_Z_gen_minus);
  eventData.SetxiZPlusGen(xi_Z_gen_plus);
  eventData.SetEtaZGen(etaZ_gen);
  eventData.SetEnergyZGen(energyZ_gen);
  eventData.SetpDissMassGen(p_diss_mass);
  eventData.SetxLpDissMass(xL_p_diss);

}

//
// Fill Detector Variables
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillDetectorVariables(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;

  double etaMax=-999;
  double etaMin=999;
  double Epz_plus=0;  
  double Epz_minus=0;  

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
  double xi_Calo_minus =0;
  double xi_Calo_plus =0;

  //bool hasHCAL;
  bool hasHF;
  bool hasHE;
  bool hasHB;
  //bool hasHO;
  //bool hasECAL;
  bool hasEE;
  bool hasEB;

  edm::Handle<CaloTowerCollection> towerCollectionH;
  event.getByLabel(caloTowerTag_,towerCollectionH);
  const CaloTowerCollection& towerCollection = *towerCollectionH;

  CaloTowerCollection::const_iterator calotower;
  calotower = towerCollection.begin();
  CaloTowerCollection::const_iterator calotowers_end = towerCollection.end();

  for(; calotower != calotowers_end; ++calotower) {

    if (fabs(calotower->eta())> 4.7) continue;   /// excluding ring12 and ring13 of HF

    //hasHCAL = false;
    hasHF = false;
    hasHE = false;
    hasHB = false;
    //hasHO = false;
    //hasECAL = false;
    hasEE = false;
    hasEB = false;     

    for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){
      DetId adetId = calotower->constituent(iconst);
      if(adetId.det()==DetId::Hcal){
	//hasHCAL = true;
	HcalDetId hcalDetId(adetId);
	if(hcalDetId.subdet()==HcalForward) hasHF = true;
	else if(hcalDetId.subdet()==HcalEndcap) hasHE = true;
	else if(hcalDetId.subdet()==HcalBarrel) hasHB = true;
	//else if(hcalDetId.subdet()==HcalOuter) hasHO = true;  
      } 
      else if(adetId.det()==DetId::Ecal){
	//hasECAL = true;
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
      if (calotower->eta() >= etaMax) etaMax=calotower->eta();
      if (calotower->eta() <= etaMin) etaMin=calotower->eta();
      xi_Calo_minus += caloTowerEt * pow(2.71,-calotower->eta()) / (7000);
      xi_Calo_plus += caloTowerEt * pow(2.71,calotower->eta()) / (7000);
      Epz_plus  += caloTowerEnergy + caloTowerPz;
      Epz_minus += caloTowerEnergy - caloTowerPz;
    }

  }  ////has to close calotower loop

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

  eventData.SetEPZCaloPlus(Epz_plus);
  eventData.SetEPZCaloMinus(Epz_minus);
  eventData.SetXiCaloPlus(xi_Calo_plus);
  eventData.SetXiCaloMinus(xi_Calo_minus);

  eventData.SetEtaCaloMax(etaMax);
  eventData.SetEtaCaloMin(etaMin);

  eventData.SetMultiplicityHFPlus(nTowersHF_plus);
  eventData.SetMultiplicityHEPlus(nTowersHE_plus);
  eventData.SetMultiplicityEEPlus(nTowersEE_plus);
  eventData.SetMultiplicityHFMinus(nTowersHF_minus);
  eventData.SetMultiplicityHEMinus(nTowersHE_minus);
  eventData.SetMultiplicityEEMinus(nTowersEE_minus); 

}

//
// Fill Physics Variables
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillVariables(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug=false;
  bool debugxi=false;
  bool debugOrder=false;

  std::vector<double> etas;
  double etaTimesEnergy=0.;
  double Epz_PF_plus=0.;
  double Epz_PF_minus=0.;
  double xi_PF_minus=0.;
  double xi_PF_plus=0.;
  double sumpx=0.;
  double sumpy=0.;
  double sumpz=0.;
  double sumpxModule=0.;
  double sumpyModule=0.;
  double sumpzModule=0.;
  double sumEnergyPF=0.;

  double MT_W_enu = 0.;
  double MT_W_munu = 0.;
  double MT_W_pfenu = 0.;
  double MT_W_pfmunu = 0.;

  int nPart_PF=0;

  std::vector<double> electronEnergy;
  std::vector<double> muEnergy;

  edm::Handle<reco::VertexCollection> Vertexes;
  event.getByLabel(PVtxCollectionTag_, Vertexes); 

  edm::Handle <reco::PFCandidateCollection> PFCandidates;
  event.getByLabel(pfTag_,PFCandidates);
  reco::PFCandidateCollection::const_iterator iter;

  eventData.SetVertex(Vertexes->size());

  PFMuonVector.clear();
  PFElectronVector.clear();

  // Fill All Electrons and Muons from Particle Flow Objects
  int NMuonsPF = 0;
  int NElectronsPF = 0;
  for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {

    const reco::PFCandidate *particle = &(*iter);

    //eta cut - excluding ring 12 13 HF  
    if (fabs(particle->eta())>4.7) continue;
    if (particle->particleId()==reco::PFCandidate::e) ++NElectronsPF;
    if (particle->particleId()==reco::PFCandidate::mu) ++NMuonsPF;

  }

  int pfsize = PFCandidates->size();
  for(int i=0; i < pfsize; ++i){
    const reco::PFCandidate* pfAll = &((*PFCandidates)[i]);

    if(NElectronsPF > 0){
      if (pfAll->particleId()==reco::PFCandidate::e){
	const reco::PFCandidate* pfAlle = &((*PFCandidates)[i]);
	PFElectronVector.push_back(pfAlle);
      }
    }

    if(NMuonsPF > 0){
      if (pfAll->particleId()==reco::PFCandidate::mu){
	const reco::PFCandidate* pfAllm = &((*PFCandidates)[i]);
	PFMuonVector.push_back(pfAllm);
      }
    }
  }


  if(PFElectronVector.size()>0){

    // Sorting Vector by pT
    const int  PFElectronVectorSize = (int) PFElectronVector.size();
    int *PFsortElectronVector=   new int[PFElectronVectorSize];
    double *PFvelectron = new double[PFElectronVectorSize];

    for (int i=0; i<PFElectronVectorSize; i++) {
      PFvelectron[i] = PFElectronVector[i]->pt();
    }

    TMath::Sort(PFElectronVectorSize, PFvelectron, PFsortElectronVector, true);

    for (unsigned int i=0;i<PFElectronVector.size();i++){
      if (debugOrder) std::cout << "ORDERED reco::Electron[" << PFsortElectronVector[i] << "]\t---> pT [GeV]: " << PFElectronVector[PFsortElectronVector[i]]->pt() << " | eT [GeV]: " << PFElectronVector[PFsortElectronVector[i]]->et() << " | eta: " << PFElectronVector[PFsortElectronVector[i]]->eta() << " | phi: " << PFElectronVector[PFsortElectronVector[i]]->phi() << " | px [GeV]: " << PFElectronVector[PFsortElectronVector[i]]->px() << " | py [GeV]: " << PFElectronVector[PFsortElectronVector[i]]->py() << std::endl;
    }

    MT_W_pfenu = TMath::Sqrt(2.*(PFElectronVector[0]->et()*NeutrinoVector[0]->et() - (PFElectronVector[0]->et()*TMath::Cos(PFElectronVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Cos(NeutrinoVector[0]->phi()) + PFElectronVector[0]->et()*TMath::Sin(PFElectronVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Sin(NeutrinoVector[0]->phi()) ) ) );

  }


  if(PFMuonVector.size()>0){

    // Sorting Vector by pT
    const int  PFMuonVectorSize = (int) PFMuonVector.size();
    int *PFsortMuonVector=   new int[PFMuonVectorSize];
    double *PFvMuon = new double[PFMuonVectorSize];

    for (int i=0; i<PFMuonVectorSize; i++) {
      PFvMuon[i] = PFMuonVector[i]->pt();
    }

    TMath::Sort(PFMuonVectorSize, PFvMuon, PFsortMuonVector, true);

    for (unsigned int i=0;i<PFMuonVector.size();i++){
      if (debugOrder) std::cout << "ORDERED reco::Muon[" << PFsortMuonVector[i] << "]\t---> pT [GeV]: " << PFMuonVector[PFsortMuonVector[i]]->pt() << " | eT [GeV]: " << PFMuonVector[PFsortMuonVector[i]]->et() << " | eta: " << PFMuonVector[PFsortMuonVector[i]]->eta() << " | phi: " << PFMuonVector[PFsortMuonVector[i]]->phi() << " | px [GeV]: " << PFMuonVector[PFsortMuonVector[i]]->px() << " | py [GeV]: " << PFMuonVector[PFsortMuonVector[i]]->py() << std::endl;
    }

    MT_W_pfmunu = TMath::Sqrt(2.*(PFMuonVector[0]->et()*NeutrinoVector[0]->et() - (PFMuonVector[0]->et()*TMath::Cos(PFMuonVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Cos(NeutrinoVector[0]->phi()) + PFMuonVector[0]->et()*TMath::Sin(PFMuonVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Sin(NeutrinoVector[0]->phi()) ) ) );

  }

  if(ElectronVector.size()>0){
    MT_W_enu = TMath::Sqrt(2.*(ElectronVector[0]->et()*NeutrinoVector[0]->et() - (ElectronVector[0]->et()*TMath::Cos(ElectronVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Cos(NeutrinoVector[0]->phi()) + ElectronVector[0]->et()*TMath::Sin(ElectronVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Sin(NeutrinoVector[0]->phi()) ) ) );
  }

  if(MuonVector.size()>0){
    MT_W_munu = TMath::Sqrt(2.*(MuonVector[0]->et()*NeutrinoVector[0]->et() - (MuonVector[0]->et()*TMath::Cos(MuonVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Cos(NeutrinoVector[0]->phi()) + MuonVector[0]->et()*TMath::Sin(MuonVector[0]->phi())*NeutrinoVector[0]->et()*TMath::Sin(NeutrinoVector[0]->phi()) ) ) );
  }

  bool w1 = false;
  bool w2 = false;
  bool w3 = false;
  bool w4 = false;

  if(PFMuonVector.size()>0){
    if (MT_W_pfmunu > 60. && MT_W_pfmunu < 110.){
      if(debugOrder)std::cout << " Mass Transverse W, muon PF: " << MT_W_pfmunu << std::endl;
      w1=true;
    }
  }

  if(PFElectronVector.size()>0){
    if (MT_W_pfenu > 60. && MT_W_pfenu < 110.){
      if(debugOrder)std::cout << " Mass Transverse W, electron PF: " << MT_W_pfenu << std::endl;
      w2=true;
    }
  }

  if(MuonVector.size()>0){
    if (MT_W_munu > 60. && MT_W_munu < 110.) {
      if(debugOrder)std::cout << " Mass Transverse W, muon: " << MT_W_munu << std::endl;
      w3=true;
    }
  }

  if(ElectronVector.size()>0){
    if (MT_W_enu > 60. && MT_W_enu < 110.){ 
      if(debugOrder)std::cout << " Mass Transverse W, electron: " << MT_W_enu << std::endl;
      w4=true;
    }
  }


  // Compute Gap Size Excluding W Candidates
  for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {

    const reco::PFCandidate *particle = &(*iter);
    double et=particle->et();
    double energy=particle->energy();
    double pt=particle->pt();
    double p=particle->p();
    double px=particle->px();
    double py=particle->py();
    double pz=particle->pz();
    double eta=particle->eta();
    double charge=particle->charge();
    double theta=particle->theta();

    // Fill 2D TTree (eta,energy);

    //eta cut - excluding ring 12 13 HF  
    if (fabs(eta)>4.7) continue;

    //int type=particle->particleId();

    TLorentzVector tmp(px,py,pz,energy);

    if  (  (fabs(charge) >0 && pt >  pTPFThresholdCharged_ ) ||
	(fabs(charge) == 0  && ( (fabs(eta) <= 1.5 && energy > energyPFThresholdBar_)  ||
				 (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > energyPFThresholdEnd_) ||
				 (fabs(eta) > 3 && energy >energyPFThresholdHF_) ) )   )
    {        

      nPart_PF++;

      Epz_PF_plus+=p+p*TMath::Cos(theta);
      Epz_PF_minus+=p-p*TMath::Cos(theta);
      xi_PF_minus += et * pow(2.71,-eta) / (7000);
      xi_PF_plus += et * pow(2.71,eta) / (7000);

      etaTimesEnergy+=eta*energy;
      sumpxModule +=fabs(px);
      sumpyModule +=fabs(py);
      sumpzModule +=fabs(pz);
      sumpx +=px;
      sumpy +=py;
      sumpz +=pz;
      sumEnergyPF +=energy;

      if(particle->particleId()==reco::PFCandidate::mu){
	if(PFMuonVector.size()>0){
	  if(w1 && (PFMuonVector[0]->pt()==pt)) continue;
	}
	if(MuonVector.size()>0){
	  if(w3 && (MuonVector[0]->pt()==pt)) continue;
	}
      }  

      if(particle->particleId()==reco::PFCandidate::e){ 
	if(PFElectronVector.size()>0){
	  if(w2 && (PFElectronVector[0]->pt()==pt)) continue;
	}
	if(ElectronVector.size()>0){
	  if(w4 && (ElectronVector[0]->pt()==pt)) continue;
	}
      }

      // Excluding W from LRG calculation
      if (debugOrder) std::cout << "Compute LRG, eta: " << eta << ", pT: "<< pt << " GeV" << std::endl;
      etas.push_back(eta);

    } 

  }

  eventData.SetXi_PF_minus(xi_PF_minus);
  eventData.SetXi_PF_plus(xi_PF_plus);
  eventData.SetEpz_PF_minus(Epz_PF_minus);
  eventData.SetEpz_PF_plus(Epz_PF_plus);
  eventData.SetMultiplicityPF(nPart_PF);
  eventData.SetSumEtaTimesEnergyPF(etaTimesEnergy);
  eventData.SetSumpxModulePF(sumpxModule);
  eventData.SetSumpyModulePF(sumpyModule);
  eventData.SetSumpzModulePF(sumpzModule);
  eventData.SetSumpxPF(sumpx);
  eventData.SetSumpyPF(sumpy);
  eventData.SetSumpzPF(sumpz);
  eventData.SetSumEnergyPF(sumEnergyPF);

  //// Computing GAPs
  //// adding two fake entries at +-4.9 in etas!!!
  etas.push_back(4.9);
  etas.push_back(-4.9);

  const int  size = (int) etas.size();
  int *sorted = new int[size];
  double *v = new double[size];
  double eta_gap_limplus = -10.0;
  double eta_gap_limminus = -10.0;

  for (int i=0; i<size; i++) {
    v[i] = etas[i];
    if (debug) std::cout<<v[i]<<std::endl;
  }
  TMath::Sort(size, v, sorted, true);

  if (size > 1) {
    double *diff = new double[size-1];
    int *diffsorted = new int[size-1];
    for (int i=0; i<(size-1); i++) {
      diff[i] = fabs(etas[sorted[i+1]]-etas[sorted[i]]);
    }

    TMath::Sort(size-1, diff, diffsorted, true);

    //checking the max gap
    double max_eta_gap=diff[diffsorted[0]];
    eta_gap_limminus = etas[sorted[diffsorted[0]+1]] ;
    eta_gap_limplus = etas[sorted[diffsorted[0]]] ;

    eventData.SetMaxGapPF(max_eta_gap);
    eventData.SetLimPlusGapPF(eta_gap_limplus);
    eventData.SetLimMinusGapPF(eta_gap_limminus);

    if (size>2) {
      double max_second_eta_gap=diff[diffsorted[1]];
      if (debug) std::cout<<" diff  " << diff[diffsorted[0]] << " sec " << diff[diffsorted[1]] << " diff size "<< diff[size-2] <<std::endl;
      eventData.SetSecondMaxGapPF(max_second_eta_gap);
    }

    else {
      eventData.SetSecondMaxGapPF(-999.);
    }

    delete [] diff;
    delete [] diffsorted;

  }

  else {

    eventData.SetMaxGapPF(-999.);
    eventData.SetSecondMaxGapPF(-999.);
    eventData.SetLimPlusGapPF(-999.);
    eventData.SetLimMinusGapPF(-999.);

  }

  delete [] sorted;
  delete [] v;

  //sorting electron energy
  const int  size3 = (int) electronEnergy.size();
  int *sorted3 = new int[size3];
  double *v3 = new double[size3];

  for (int i=0; i<size3; i++) {
    v3[i] = electronEnergy[i];
  }
  TMath::Sort(size3, v3, sorted3, true);
  for (int i=0; i<size3; i++) {
    electronEnergy[i] = v3[sorted3[i]];
  }

  //sorting muon energy
  const int  size4 = (int) muEnergy.size();
  int *sorted4 = new int[size4];
  double *v4 = new double[size4];

  for (int i=0; i<size4; i++) {
    v4[i] = muEnergy[i];
  }
  TMath::Sort(size4, v4, sorted4, true);
  for (int i=0; i<size4; i++) {
    muEnergy[i] = v4[sorted4[i]];
  }
  delete [] sorted3;
  delete [] v3;
  delete [] sorted4;
  delete [] v4;

  double MXsumPxPF = 0.;
  double MXsumPyPF = 0.;
  double MXsumPzPF = 0.;
  double MXsumEPF = 0.;
  double MYsumPxPF = 0.;
  double MYsumPyPF = 0.;
  double MYsumPzPF = 0.;
  double MYsumEPF = 0.;
  double xiMass = -999.;

  for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {

    const reco::PFCandidate *particle = &(*iter);
    double energy=particle->energy();
    double pt=particle->pt();
    double eta=particle->eta();
    double charge=particle->charge();

    //eta cut - excluding ring 12 13 HF  
    if (fabs(eta)>4.7) continue;

    if  (  (fabs(charge) >0 && pt >  pTPFThresholdCharged_ ) ||
	(fabs(charge) == 0  && ( (fabs(eta) <= 1.5 && energy > energyPFThresholdBar_)  ||
				 (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > energyPFThresholdEnd_) ||
				 (fabs(eta) > 3 && energy >energyPFThresholdHF_) ) )   )
    {        

      if ( particle->eta() >= eta_gap_limplus ){
	MXsumPxPF += particle->px();
	MXsumPyPF += particle->py();
	MXsumPzPF += particle->pz();
	MXsumEPF += particle->energy();
      }
      else {
	MYsumPxPF += particle->px();
	MYsumPyPF += particle->py();
	MYsumPzPF += particle->pz();
	MYsumEPF += particle->energy();
      }

    } 

  }

  TLorentzVector M_x(MXsumPxPF,MXsumPyPF,MXsumPzPF,MXsumEPF);
  TLorentzVector M_y(MYsumPxPF,MYsumPyPF,MYsumPzPF,MYsumEPF);

  double massX2 = pow(M_x.M(),2);
  double massY2 = pow(M_y.M(),2);

  if (massX2 > massY2 && massX2 > 0.) xiMass = massX2/(7000.*7000.);
  else if (massY2 > massX2 && massY2 > 0.) xiMass = massY2/(7000.*7000.);

  if (debugxi) {
    std::cout << "Xi Computation" << std::endl;
    std::cout << ">>>>> Xi, including Z: " << xiMass << std::endl;
  }

  eventData.SetXiMass(xiMass);

  //TLorentzVector dataMass_plus(0.,0.,0.,0.);
  //TLorentzVector dataMass_minus(0.,0.,0.,0.);
  int nplus =0;
  int nminus =0;

  double sumPTPFm = 0.;
  double sumPTPFp = 0.;
  double MXsumPxPFNoW = 0.;
  double MXsumPyPFNoW = 0.;
  double MXsumPzPFNoW = 0.;
  double MXsumEPFNoW = 0.;
  double MYsumPxPFNoW = 0.;
  double MYsumPyPFNoW = 0.;
  double MYsumPzPFNoW = 0.;
  double MYsumEPFNoW = 0.;
  double xiMassNoW = -999.;

  for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {
    const reco::PFCandidate *particle = &(*iter);
    double energy=particle->energy();
    double pt=particle->pt();
    double px=particle->px();
    double py=particle->py();
    double pz=particle->pz();
    double eta=particle->eta();
    double charge=particle->charge();

    //eta cut - excluding ring 12 13 HF  
    if (fabs(eta)>4.7) continue;

    TLorentzVector tmp(px,py,pz,energy); 

    if  (  (fabs(charge) >0 && pt >  pTPFThresholdCharged_ ) ||
	(fabs(charge) == 0  && ( (fabs(eta) <= 1.5 && energy > energyPFThresholdBar_)  ||
				 (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > energyPFThresholdEnd_) ||
				 (fabs(eta) > 3 && energy >energyPFThresholdHF_) ) )   )
    {

      if(particle->particleId()==reco::PFCandidate::mu){
	if(PFMuonVector.size()>0){
	  if(w1 && (PFMuonVector[0]->pt()==pt)) continue;
	}
	if(MuonVector.size()>0){
	  if(w3 && (MuonVector[0]->pt()==pt)) continue;
	}
      }

      if(particle->particleId()==reco::PFCandidate::e){
	if(PFElectronVector.size()>0){
	  if(w2 && (PFElectronVector[0]->pt()==pt)) continue;
	}
	if(ElectronVector.size()>0){
	  if(w4 && (ElectronVector[0]->pt()==pt)) continue;
	}
      }

      if ( eta >= eta_gap_limplus ){
	//dataMass_plus+=tmp;
	nplus++;
	sumPTPFp += pt;
	MXsumPxPFNoW += particle->px();
	MXsumPyPFNoW += particle->py();
	MXsumPzPFNoW += particle->pz();
	MXsumEPFNoW += particle->energy();
      }
      else {
	//dataMass_minus+=tmp;
	nminus++;
	sumPTPFm += pt;
	MYsumPxPFNoW += particle->px();
	MYsumPyPFNoW += particle->py();
	MYsumPzPFNoW += particle->pz();
	MYsumEPFNoW += particle->energy();
      }

      if (debugOrder) std::cout << "SUM, pT: Compute LRG, eta: " << eta << ", pT: "<< pt << " GeV" << std::endl;

    }
  }  // PF loop

  if (sumPTPFp > sumPTPFm){
    eventData.SetPTMaxGapMaxPF(sumPTPFp);
    eventData.SetPTMinGapMaxPF(sumPTPFm);
  }else{
    eventData.SetPTMaxGapMaxPF(sumPTPFm);
    eventData.SetPTMinGapMaxPF(sumPTPFp);
  }

  eventData.SetElectronEnergyPF(-999.); // First Electron, Fill Second Electron also. Eta, phi, pT and ISO from PF.
  eventData.SetMuEnergyPF(-999.); // First Muon, Fill Second Muon also. Eta, phi, pT and ISO from PF.
  eventData.SetMultiplicityGapPlusPF(nplus);
  eventData.SetMultiplicityGapMinusPF(nminus);

  TLorentzVector M_xNoz(MXsumPxPFNoW,MXsumPyPFNoW,MXsumPzPFNoW,MXsumEPFNoW);
  TLorentzVector M_yNoz(MYsumPxPFNoW,MYsumPyPFNoW,MYsumPzPFNoW,MYsumEPFNoW);

  double massX2NoW = pow(M_xNoz.M(),2);
  double massY2NoW = pow(M_yNoz.M(),2);

  if (massX2NoW > massY2NoW && massX2NoW > 0.) xiMassNoW = massX2NoW/(7000.*7000.);
  else if (massY2NoW > massX2NoW && massY2NoW > 0.) xiMassNoW = massY2NoW/(7000.*7000.);

  if (debugxi) {
    std::cout << "Xi Computation" << std::endl;
    std::cout << ">>>>> Xi, no Z: " << xiMassNoW << std::endl;
  }

  eventData.SetXiMassNoW(xiMassNoW);

}


//
// Fill Castor 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillCastor(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

void DiffractiveWAnalysis::fillCastorDebug(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

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

void DiffractiveWAnalysis::fillZDC(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  /*

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

*/

}

//
// Fill Tower Information Energy x Eta
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiffractiveWAnalysis::fillDetectorEnergyEtaInfo(DiffractiveWEvent& eventData, const edm::Event& event, const edm::EventSetup& setup){

  bool debug = false;
  bool debug_deep = false;

  //bool hasHCAL;
  bool hasHF;
  bool hasHE;
  bool hasHB;
  //bool hasHO;
  //bool hasECAL;
  bool hasEE;
  bool hasEB;

  std::vector<double> energy_tower;
  std::vector<double> eta_tower; 

  edm::Handle<CaloTowerCollection> towerCollectionH;
  event.getByLabel(caloTowerTag_,towerCollectionH);
  const CaloTowerCollection& towerCollection = *towerCollectionH;

  CaloTowerCollection::const_iterator calotower;
  calotower = towerCollection.begin();
  CaloTowerCollection::const_iterator calotowers_end = towerCollection.end();

  int counter_tower=0;

  for(; calotower != calotowers_end; ++calotower) {

    if (fabs(calotower->eta())> 4.7) continue;   /// excluding ring12 and ring13 of HF

    //hasHCAL = false;
    hasHF = false;
    hasHE = false;
    hasHB = false;
    //hasHO = false;
    //hasECAL = false;
    hasEE = false;
    hasEB = false;  

    for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){

      DetId adetId = calotower->constituent(iconst);
      if(adetId.det()==DetId::Hcal){
	//hasHCAL = true;
	if (debug_deep) std::cout << "HCAL is true." << std::endl;
	HcalDetId hcalDetId(adetId);
	if(hcalDetId.subdet()==HcalForward) {
	  hasHF = true;
	  if (debug_deep) std::cout << "HF is true." << std::endl;
	}
	else if(hcalDetId.subdet()==HcalEndcap) {
	  hasHE = true;
	  if (debug_deep) std::cout << "HE is true." << std::endl;
	}
	else if(hcalDetId.subdet()==HcalBarrel) {
	  hasHB = true;
	  if (debug_deep) std::cout << "HB is true." << std::endl;
	} 
	else if(hcalDetId.subdet()==HcalOuter) {
	  //hasHO = true;  
	  if (debug_deep) std::cout << "HO is true." << std::endl;
	}
      } 
      else if(adetId.det()==DetId::Ecal){
	//hasECAL = true;
	if (debug_deep) std::cout << "ECAL is true." << std::endl;
	EcalSubdetector ecalSubDet = (EcalSubdetector)adetId.subdetId();
	if(ecalSubDet == EcalEndcap) {
	  hasEE = true;
	  if (debug_deep) std::cout << "EE is true." << std::endl; 
	}
	else if(ecalSubDet == EcalBarrel) {
	  hasEB = true;
	  if (debug_deep) std::cout << "EB is true." << std::endl;
	}
      }
    }

    double caloTowerEnergy = calotower->energy();
    double caloTowerEta = calotower->eta();
    double caloTowerPhi = calotower->phi();
    double caloTowerEmEnergy = calotower->emEnergy();
    double caloTowerHadEnergy = calotower->hadEnergy();

    if( hasHF && !hasHE )
    {

      if (debug_deep) std::cout << "HF, no threshold." << std::endl;    

      if( caloTowerEnergy > energyThresholdHF_ && fabs(calotower->eta())> 2.98 )   //// excluding HF ring1
      {
	++counter_tower;
	energy_tower.push_back(caloTowerEnergy);
	eta_tower.push_back(caloTowerEta);

	if (debug) {
	  std::cout << "HF Energy for each CaloTower (GeV): " << caloTowerEnergy << " | Eta for each CaloTower: " << caloTowerEta << " | Phi for each CaloTower: " << caloTowerPhi << std::endl;
	}

      }
    }
    else if( hasHE && !hasHF && !hasHB )
    {

      if (debug_deep) std::cout << "HE, no threshold." << std::endl;
      if( caloTowerHadEnergy > energyThresholdHE_)
      {
	++counter_tower;
	energy_tower.push_back(caloTowerEnergy);
	eta_tower.push_back(caloTowerEta);

	if (debug) {
	  std::cout << "HE Energy for each CaloTower (GeV): " << caloTowerEnergy << " | Eta for each CaloTower: " << caloTowerEta << " | Phi for each CaloTower: " << caloTowerPhi << std::endl;
	}

      }
    }
    else if( hasHB && !hasHE )
    {
      if (debug_deep) std::cout << "HB, no threshold." << std::endl;
      if( caloTowerHadEnergy > energyThresholdHB_)
      {
	++counter_tower;
	energy_tower.push_back(caloTowerEnergy);
	eta_tower.push_back(caloTowerEta);

	if (debug) {
	  std::cout << "HB Energy for each CaloTower (GeV): " << caloTowerEnergy << " | Eta for each CaloTower: " << caloTowerEta << " | Phi for each CaloTower: " << caloTowerPhi << std::endl;
	}

      }
    }

    if( hasEE && !hasEB )
    {
      if (debug_deep) std::cout << "EE, no threshold." << std::endl;
      if( caloTowerEmEnergy >= energyThresholdEE_)
      {
	++counter_tower;
	energy_tower.push_back(caloTowerEnergy);
	eta_tower.push_back(caloTowerEta);
	if (debug) {
	  std::cout << "EB Energy for each CaloTower (GeV): " << caloTowerEnergy << " | Eta for each CaloTower: " << caloTowerEta << " | Phi for each CaloTower: " << caloTowerPhi << std::endl;
	}

      }
    }
    else if( hasEB && !hasEE )
    {
      if (debug_deep) std::cout << "EB, no threshold." << std::endl;
      if( caloTowerEmEnergy >= energyThresholdEB_)
      {
	++counter_tower;
	energy_tower.push_back(caloTowerEnergy);
	eta_tower.push_back(caloTowerEta);

	if (debug) {
	  std::cout << "EB Energy for each CaloTower (GeV): " << caloTowerEnergy << " | Eta for each CaloTower: " << caloTowerEta << " | Phi for each CaloTower: " << caloTowerPhi << std::endl;
	}
      }
    }
  }  ////has to close calotower loop

  if (debug) std::cout << "Active Towers: " << counter_tower << std::endl;
  eventData.SetEachTowerCounter(counter_tower);
  eventData.SetEachTowerEta(eta_tower);
  eventData.SetEachTowerEnergy(energy_tower);

}
