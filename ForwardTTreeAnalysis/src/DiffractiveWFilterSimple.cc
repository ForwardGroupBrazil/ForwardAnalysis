/*
   >>-------------------------<<
   Diffractive W Filter
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

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


class diffractiveWFilterSimple : public edm::EDFilter{
  public:
    explicit diffractiveWFilterSimple(const edm::ParameterSet&);
    ~diffractiveWFilterSimple();

  private:
    void setTFileService();
    void PrintOrder();
    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

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

    TH1F *run_selected;
    TH1F *run_total;

    struct orderPT
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->pt() > vec2->pt());
	}
    };

    struct orderETA
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->eta() > vec2->eta());
	}
    };

    struct orderAbsolutPZ
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->pz()) > fabs(vec2->pz()));
	}
    };

};


diffractiveWFilterSimple::diffractiveWFilterSimple(const edm::ParameterSet& iConfig):
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  metTag_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  patmetTag_(iConfig.getUntrackedParameter<edm::InputTag>("patmetTag"))
{

}

diffractiveWFilterSimple::~diffractiveWFilterSimple(){
}

void diffractiveWFilterSimple::beginJob(){
  setTFileService();
}

bool diffractiveWFilterSimple::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool debug = false;

  // Fill reco::Muon
  //edm::Handle<reco::MuonCollection> muons;
  edm::Handle<edm::View<reco::Muon> > muons;
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

  // Fill reco::GsfElectron
  //edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<edm::View<reco::GsfElectron> >electrons;
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

  // Fill pfMET
  //edm::Handle<reco::PFMETCollection> met;
  edm::Handle<edm::View<reco::PFMET> > met;
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


  // S O R T I N G   V E C T O R S
  ////////////////////////////////

  std::sort(ElectronVector.begin(), ElectronVector.end(), orderPT());
  std::sort(MuonVector.begin(), MuonVector.end(), orderPT());
  std::sort(NeutrinoVector.begin(), NeutrinoVector.end(), orderPT());

  bool acceptmet = false;

  if (neutrinosize > 0){
    if (NeutrinoVector[0]->pt()>20.) acceptmet = true;
  }

  bool acceptlepton = false;

  if(muonsize > 0){
    if(MuonVector.size()==1){
      if(MuonVector[0]->pt()>10.) acceptlepton = true; // 20 GeV
    }
    if(MuonVector.size()>1){
      if(MuonVector[0]->pt()>10. && MuonVector[1]->pt()>10.) acceptlepton = false; // 20 GeV and 10 GeV
    }
  }

  if(electronsize >0){
    if(ElectronVector.size()==1){
      if(ElectronVector[0]->pt()>10.) acceptlepton = true; // 20 GeV
    }
    if(ElectronVector.size()>1){
      if(ElectronVector[0]->pt()>10. && ElectronVector[1]->pt()>5.) acceptlepton = false; // 20 GeV and 20 GeV
    }
  }

  bool acceptevent = false;
  acceptevent = acceptlepton && acceptmet;

  int runNumber = iEvent.id().run();  
  run_total->Fill(runNumber);
  if(acceptevent){
    run_selected->Fill(runNumber);
  }

  if(debug){
    PrintOrder();
    if(acceptevent) std::cout << "EVENT ACCEPTED" << std::endl;
  }

  return acceptevent;

}

void diffractiveWFilterSimple::endJob(){
}

void diffractiveWFilterSimple::setTFileService(){

  edm::Service<TFileService> fs;
  TFileDirectory runinfoDir = fs->mkdir("RunInfo");
  run_total = runinfoDir.make<TH1F>("RunTotal","; Run Id(); # Events per Run",1,0,1);
  run_total->SetBit(TH1::kCanRebin);
  run_selected = runinfoDir.make<TH1F>("RunSelected","Events per Run After Filter; Run Id(); # Events per Run",1,0,1);
  run_selected->SetBit(TH1::kCanRebin);

}

void diffractiveWFilterSimple::PrintOrder(){

  std::cout << "\n-- BEGIN EVENT --" << std::endl;

  if(ElectronVector.size()>0){
    for (unsigned int i=0;i<ElectronVector.size();i++){
      cout << "reco::Electron[" << i << "]\t---> pT [GeV]: " << ElectronVector[i]->pt() << " | eT [GeV]: " << ElectronVector[i]->et() << " | eta: " << ElectronVector[i]->eta() << " | phi: " << ElectronVector[i]->phi() << endl;
    }
  }

  if(MuonVector.size()>0){
    for (unsigned int i=0;i<MuonVector.size();i++){
      cout << "reco::Muon[" << i << "]\t---> pT [GeV]: " << MuonVector[i]->pt() << " | eT [GeV]: " << MuonVector[i]->et() << " | eta: " << MuonVector[i]->eta() << " | phi: " << MuonVector[i]->phi() << endl;
    }
  }

  if(PFElectronVector.size()>0){
    for (unsigned int i=0;i<PFElectronVector.size();i++){
      cout << "reco::PFElectron[" << i << "]\t---> pT [GeV]: " << PFElectronVector[i]->pt() << " | eT [GeV]: " << PFElectronVector[i]->et() << " | eta: " << PFElectronVector[i]->eta() << " | phi: " << PFElectronVector[i]->phi() << endl;
    }
  }

  if(PFMuonVector.size()>0){
    for (unsigned int i=0;i<PFMuonVector.size();i++){
      cout << "reco::PFMuon[" << i << "]\t---> pT [GeV]: " << PFMuonVector[i]->pt() << " | eT [GeV]: " << PFMuonVector[i]->et() << " | eta: " << PFMuonVector[i]->eta() << " | phi: " << PFMuonVector[i]->phi() << endl;
    }
  }

  if(PatElectronVector.size()>0){
    for (unsigned int i=0;i<PatElectronVector.size();i++){
      cout << "pat::Electron[" << i << "]\t---> pT [GeV]: " << PatElectronVector[i]->pt() << " | eT [GeV]: " << PatElectronVector[i]->et() << " | eta: " << PatElectronVector[i]->eta() << " | phi: " << PatElectronVector[i]->phi() << endl;
    }
  }

  if(PatMuonVector.size()>0){
    for (unsigned int i=0;i<PatMuonVector.size();i++){
      cout << "pat::Muon[" << i << "]\t---> pT [GeV]: " << PatMuonVector[i]->pt() << " | eT [GeV]: " << PatMuonVector[i]->et() << " | eta: " << PatMuonVector[i]->eta() << " | phi: " << PatMuonVector[i]->phi() << endl;
    }
  }

}

DEFINE_FWK_MODULE(diffractiveWFilterSimple);
