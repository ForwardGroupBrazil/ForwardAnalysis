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

};


diffractiveWFilter::diffractiveWFilter(const edm::ParameterSet& iConfig):
  nLeptons_(iConfig.getUntrackedParameter<int>("nLeptons",1)),
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag"))
{

}

diffractiveWFilter::~diffractiveWFilter(){
}

void diffractiveWFilter::beginJob(){
}

bool diffractiveWFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool debug = false;
  bool Filter_diffractiveW = false;

  edm::Handle<std::vector<pat::Muon> > pmuons;
  iEvent.getByLabel("patMuons", pmuons);

  edm::Handle<std::vector<pat::Electron> > pelectrons;
  iEvent.getByLabel("patElectrons", pelectrons);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(muonTag_,muons);

  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel(electronTag_,electrons);

  edm::Handle<reco::PFMETCollection> met;
  iEvent.getByLabel("pfMet",met);

  int electronsize = electrons->size();
  int muonsize = muons->size();
  int pelectronsize = pelectrons->size();
  int pmuonsize = pmuons->size();
  int metsize = met->size();

  const reco::PFMET* neutrino = NULL;

  for(int itmet=0; itmet < metsize; ++itmet){
    neutrino = &((*met)[itmet]);
  }

  if ( (electronsize >= nLeptons_ || muonsize >= nLeptons_ || pelectronsize >= nLeptons_ || pmuonsize >= nLeptons_) && neutrino->et()>20.) {
    Filter_diffractiveW = true;
    if (debug) std::cout << "Event Selected." << std::endl;
  }

  return Filter_diffractiveW;

}


void diffractiveWFilter::endJob(){
}


DEFINE_FWK_MODULE(diffractiveWFilter);
