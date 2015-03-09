#ifndef DiffractiveZAnalysis_h
#define DiffractiveZAnalysis_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <math.h>
#include <cmath>

class DiffractiveZEvent;
class TH1F;
class TH2F;

namespace diffractiveZAnalysis {

  class DiffractiveZAnalysis {
    public:
      typedef DiffractiveZEvent event_type;
      static const char* name;

      DiffractiveZAnalysis() {} 
      DiffractiveZAnalysis(const edm::ParameterSet&);
      ~DiffractiveZAnalysis();

      void begin();
      void begin(const edm::Run&, const edm::EventSetup&);
      void fill(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void end();
    private:

      void setTFileService();
      void fillTriggerInfo(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillElectronsInfo(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillMuonsInfo(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillTracksInfo(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillGenInfo(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillZPat(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillDetectorVariables(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillVariables(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCastor(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCastorDebug(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void fillZDC(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);
      void VertexAssociation(DiffractiveZEvent&, const edm::Event&, const edm::EventSetup&);

      template <class T, class W>
	math::XYZTLorentzVector DiSystem(T obj1, W obj2);

      template <class T, class W>
	double InvariantMass(T lepton1, W lepton2);

      edm::InputTag triggerResultsTag_;
      std::vector<std::string> hltPathNames_;
      edm::InputTag electronTag_;
      edm::InputTag muonTag_;
      edm::InputTag pfTag_;
      edm::InputTag genTag_;
      edm::InputTag PVtxCollectionTag_;
      edm::InputTag castorHitsTag_;
      edm::InputTag zdcHitsTag_;
      bool RunCastor_;
      bool RunZDC_;
      bool RunMC_;
      bool RunZPat_;
      bool RunA_;
      bool RunB_;
      bool EachTower_;
      double pTPFThresholdCharged_;
      double energyPFThresholdBar_;
      double energyPFThresholdEnd_;
      double energyPFThresholdHF_;
      double energyThresholdHB_;
      double energyThresholdHE_;
      double energyThresholdHF_;
      double energyThresholdEB_;
      double energyThresholdEE_;
      double castorThreshold_;
      double fCGeVCastor_;
      edm::InputTag caloTowerTag_; 
      edm::InputTag trackTag_;
      double beamEnergy_;

      std::string selectionPathName_;

      std::vector<const reco::Muon*> MuonVector;
      std::vector<const reco::GsfElectron*> ElectronVector;
      std::vector<const pat::Muon*> PatMuonVector;
      std::vector<const pat::Electron*> PatElectronVector;
      std::vector<const reco::PFCandidate*> PFMuonVector;
      std::vector<const reco::PFCandidate*> PFElectronVector;
      std::vector<const reco::PFCandidate*> PFVector;
      std::vector<const reco::PFCandidate*> PFNoZVector;
      std::vector<const reco::PFCandidate*> PFHFVector;
      std::vector<const reco::GenParticle*> genVector;
      std::vector<const reco::GenParticle*> genCMSVector;
      std::vector<const reco::GenParticle*> genProtonVector;
      std::vector<const CaloTower*> towerVector;
      std::vector<const reco::Vertex*> VertexVector;

      TH1F *hltTriggerPassHisto_,*hltTriggerNamesHisto_;
      TH1F *CastorChannelHisto_;
      TH1F *histo_castor_channels;
      TH1F *hParticlesEta_, *hParticlesPhi_, *hParticlesPt_;
      std::vector<TH1F*> m_hVector_histo_castor_channels;

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

      struct orderVZ
      {
	template <class T, class W>
	  inline bool operator() (T vec1, W vec2)
	  {
	    return (fabs(vec1->z()) > fabs(vec2->z()));
	  }
      };

  };

} // namespace
#endif 
