#ifndef TriggerAnalysis_h
#define TriggerAnalysis_h

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

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <math.h>
#include <cmath>

class TriggerEvent;
class TH1F;
class TH2F;

namespace triggerAnalysis {

  class TriggerAnalysis {
    public:
      typedef TriggerEvent event_type;
      static const char* name;

      TriggerAnalysis() {} 
      TriggerAnalysis(const edm::ParameterSet&);
      ~TriggerAnalysis();

      void begin();
      void begin(const edm::Run&, const edm::EventSetup&);
      void fill(TriggerEvent&, const edm::Event&, const edm::EventSetup&);
      void end();
    private:

      void setTFileService();
      void fillVertexInfo(TriggerEvent&, const edm::Event&, const edm::EventSetup&);
      void fillTrackInfo(TriggerEvent&, const edm::Event&, const edm::EventSetup&);
      void fillTriggerInfo(TriggerEvent&, const edm::Event&, const edm::EventSetup&);
      void fillMETInfo(TriggerEvent&, const edm::Event&, const edm::EventSetup&);
      void fillElectronsInfo(TriggerEvent&, const edm::Event&, const edm::EventSetup&);
      void fillMuonsInfo(TriggerEvent&, const edm::Event&, const edm::EventSetup&);

      template <class T, class W>
	math::XYZTLorentzVector DiSystem(T obj1, W obj2);

      template <class T, class W>
	double TransverseMass(T lepton1, W lepton2);

      edm::InputTag triggerResultsTag_;
      std::vector<std::string> hltPathNames_;
      edm::InputTag electronTag_;
      edm::InputTag muonTag_;
      edm::InputTag metTag_;
      edm::InputTag PVtxCollectionTag_;
      edm::InputTag trackTag_;

      std::string selectionPathName_;

      std::vector<const reco::PFMET*> NeutrinoVector;
      std::vector<const reco::Muon*> MuonVector;
      std::vector<const reco::GsfElectron*> ElectronVector;

      TH1F *hltTriggerPassHisto_,*hltTriggerNamesHisto_;

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
