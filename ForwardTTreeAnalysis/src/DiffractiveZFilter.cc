/*
   >>-------------------------<<
   Diffractive Z Filter
   >>-------------------------<<

Goal:
reduce NTuple size. Preselect events with at least 2 leptons.

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

using namespace edm;
using namespace std;
using namespace reco;


class diffractiveZFilter : public edm::EDFilter{
  public:
    explicit diffractiveZFilter(const edm::ParameterSet&);
    ~diffractiveZFilter();

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    int nLeptons_;
    edm::InputTag muonTag_;
    edm::InputTag electronTag_;

};


diffractiveZFilter::diffractiveZFilter(const edm::ParameterSet& iConfig):
  nLeptons_(iConfig.getUntrackedParameter<int>("nLeptons",2)),
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag"))
{

}


diffractiveZFilter::~diffractiveZFilter(){
}

void diffractiveZFilter::beginJob(){
}

bool diffractiveZFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool debug = false;
  bool Filter_diffractiveZ = false;

  edm::Handle<std::vector<pat::Muon> > pmuons;
  iEvent.getByLabel("patMuons", pmuons);

  edm::Handle<std::vector<pat::Electron> > pelectrons;
  iEvent.getByLabel("patElectrons", pelectrons);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(muonTag_,muons);

  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel(electronTag_,electrons);

  int electronsize = electrons->size();
  int muonsize = muons->size();
  int pelectronsize = pelectrons->size();
  int pmuonsize = pmuons->size();

  if (electronsize >= nLeptons_ || muonsize >= nLeptons_ || pelectronsize >= nLeptons_ || pmuonsize >= nLeptons_) {
    Filter_diffractiveZ = true;
    if (debug) std::cout << "Event Selected." << std::endl;
  }

  return Filter_diffractiveZ;
}


void diffractiveZFilter::endJob(){
}


DEFINE_FWK_MODULE(diffractiveZFilter);
