/*
   >>-------------------------<<
   Diffractive Z Filter
   >>-------------------------<<

Goal:
reduce NTuple size. Preselect events with at least 1 leptons and 1 MET.

Authors: D. Figueiredo, R. Arciadiacono and N. Cartiglia
 */

//--> System include files
#include <memory>

//--> Includes Default Files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

//--> RecoMuon and RecoElectron
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//--> PatMuon and PatElectron
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//--> MET Collection
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/PatCandidates/interface/MET.h"

using namespace edm;
using namespace std;
using namespace reco;


class diffractiveWFilter : public edm::EDFilter{
  public:
    explicit diffractiveWFilter(const edm::ParameterSet&);
    ~diffractiveWFilter();

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    int nLeptons_;
    edm::InputTag muonTag_;
    edm::InputTag electronTag_;
    edm::InputTag metTag_;
    edm::InputTag patmetTag_;

    std::vector<const reco::Muon*> MuonVector;
    std::vector<const reco::GsfElectron*> ElectronVector;
    std::vector<const pat::Muon*> PatMuonVector;
    std::vector<const pat::Electron*> PatElectronVector;
    std::vector<const reco::PFMET*> NeutrinoVector;
    std::vector<const pat::MET*> PatNeutrinoVector;
    std::vector<const reco::PFCandidate*> PFMuonVector;
    std::vector<const reco::PFCandidate*> PFElectronVector;

};


diffractiveWFilter::diffractiveWFilter(const edm::ParameterSet& iConfig):
  nLeptons_(iConfig.getUntrackedParameter<int>("nLeptons",1)),
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  metTag_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  patmetTag_(iConfig.getUntrackedParameter<edm::InputTag>("patmetTag"))
{

}

diffractiveWFilter::~diffractiveWFilter(){
}

void diffractiveWFilter::beginJob(){
}

bool diffractiveWFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool debug = true;

  // Fill reco::Muon
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(muonTag_,muons);

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
  iEvent.getByLabel("patMuons", patMuons);

  int patMuonsize = patMuons->size();
  int itpatMuon;
  PatMuonVector.clear();

  if(patMuons->size()>0){
    for(itpatMuon=0; itpatMuon < patMuonsize; ++itpatMuon){
      const pat::Muon* patMuonAll = &((*patMuons)[itpatMuon]);
      PatMuonVector.push_back(patMuonAll);
    }
  }

  // Fill reco::GsfElectron
  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel(electronTag_,electrons);

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
  iEvent.getByLabel("patElectrons", patElectrons);

  int patElectronsize = patElectrons->size();
  int itpatElectron;
  PatElectronVector.clear();

  if(patElectrons->size()>0){
    for(itpatElectron=0; itpatElectron < patElectronsize; ++itpatElectron){
      const pat::Electron* patElectronAll = &((*patElectrons)[itpatElectron]);
      PatElectronVector.push_back(patElectronAll);
    }
  }

  // Fill pfMET
  edm::Handle<reco::PFMETCollection> met;
  iEvent.getByLabel(metTag_,met);

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
  edm::Handle<std::vector<pat::MET> > patmet;
  iEvent.getByLabel(patmetTag_,patmet);

  PatNeutrinoVector.clear();
  int patneutrinosize = patmet->size();
  int itpatmet;

  if(patmet->size()>0){
    for(itpatmet=0; itpatmet < patneutrinosize; ++itpatmet){
      const pat::MET* patneutrinoAll = &((*patmet)[itpatmet]);
      PatNeutrinoVector.push_back(patneutrinoAll);
    }
  }


  // S O R T I N G   V E C T O R S
  ////////////////////////////////

  bool pTSel = false;
  bool NotSecondPt = true;

  if(MuonVector.size()>0){

    // Sorting Vector by pT
    const int MuonVectorSize = (int) MuonVector.size();
    int *sortMuonVector= new int[MuonVectorSize];
    double *vmuon = new double[MuonVectorSize];

    for (int i=0; i<MuonVectorSize; i++) {
      vmuon[i] = MuonVector[i]->pt();
    }

    TMath::Sort(MuonVectorSize, vmuon, sortMuonVector, true);

    if ( MuonVector[sortMuonVector[0]]->pt() > 20.) pTSel = true;
    if (MuonVector.size() > 1){
      if (MuonVector[sortMuonVector[1]]->pt() > 20.) NotSecondPt = false;
    }

  }

  if(PatMuonVector.size()>0){

    // Sorting Vector by pT
    const int PatMuonVectorSize = (int) PatMuonVector.size();
    int *sortPatMuonVector= new int[PatMuonVectorSize];
    double *vpatmuon = new double[PatMuonVectorSize];

    for (int i=0; i<PatMuonVectorSize; i++) {
      vpatmuon[i] = PatMuonVector[i]->pt();
    }

    TMath::Sort(PatMuonVectorSize, vpatmuon, sortPatMuonVector, true);

    if ( PatMuonVector[sortPatMuonVector[0]]->pt() > 20.) pTSel = true;
    if (PatMuonVector.size() > 1){
      if (PatMuonVector[sortPatMuonVector[1]]->pt() > 20.) NotSecondPt = false;
    }

  }


  if(ElectronVector.size()>0){

    // Sorting Vector by pT
    const int ElectronVectorSize = (int) ElectronVector.size();
    int *sortElectronVector= new int[ElectronVectorSize];
    double *velectron = new double[ElectronVectorSize];

    for (int i=0; i<ElectronVectorSize; i++) {
      velectron[i] = ElectronVector[i]->pt();
    }

    TMath::Sort(ElectronVectorSize, velectron, sortElectronVector, true);

    if ( ElectronVector[sortElectronVector[0]]->pt() > 20.) pTSel = true;
    if (ElectronVector.size() > 1){
      if (ElectronVector[sortElectronVector[1]]->pt() > 20.) NotSecondPt = false;
    }

  }


  if(PatElectronVector.size()>0){

    // Sorting Vector by pT
    const int PatElectronVectorSize = (int) PatElectronVector.size();
    int *sortPatElectronVector= new int[PatElectronVectorSize];
    double *vpatelectron = new double[PatElectronVectorSize];

    for (int i=0; i<PatElectronVectorSize; i++) {
      vpatelectron[i] = PatElectronVector[i]->pt();
    }

    TMath::Sort(PatElectronVectorSize, vpatelectron, sortPatElectronVector, true);

    if ( PatElectronVector[sortPatElectronVector[0]]->pt() > 20.) pTSel = true;
    if (PatElectronVector.size() > 1){
      if (PatElectronVector[sortPatElectronVector[1]]->pt() > 20.) NotSecondPt = false;
    }

  }

  bool NLepton = false;
  bool METpTSel = false;

  bool AllSelection = false;

  int ElectronSize = ElectronVector.size();
  int MuonSize = MuonVector.size();
  int PatElectronSize = PatElectronVector.size();
  int PatMuonSize = PatMuonVector.size();
  int NeutrinoSize = NeutrinoVector.size();
  int PatNeutrinoSize = PatNeutrinoVector.size();

  if (ElectronSize >= nLeptons_ || MuonSize >= nLeptons_ || PatElectronSize >= nLeptons_ || PatMuonSize >= nLeptons_) NLepton = true;

  if (NeutrinoSize > 0 || PatNeutrinoSize > 0){
    if (NeutrinoVector[0]->pt()>15. || PatNeutrinoVector[0]->pt()>15.) METpTSel = true;
  }

  AllSelection = NLepton & pTSel & NotSecondPt & METpTSel;

  if(debug){
    if (AllSelection) std::cout << "\n\n< Event Selected >\n\n" << std::endl;
  }

  return AllSelection;

}


void diffractiveWFilter::endJob(){
}


DEFINE_FWK_MODULE(diffractiveWFilter);
