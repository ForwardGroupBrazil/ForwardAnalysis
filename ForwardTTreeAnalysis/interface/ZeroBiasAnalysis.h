#ifndef ZeroBiasAnalysis_h
#define ZeroBiasAnalysis_h

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

class ZeroBiasEvent;
class TH1F;
class TH2F;

namespace zerobiasAnalysis {

  class ZeroBiasAnalysis {
    public:
      typedef ZeroBiasEvent event_type;
      static const char* name;

      ZeroBiasAnalysis() {} 
      ZeroBiasAnalysis(const edm::ParameterSet&);
      ~ZeroBiasAnalysis();

      void begin();
      void begin(const edm::Run&, const edm::EventSetup&);
      void fill(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void end();
    private:

      void setTFileService();
      void fillVertexInfo(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillTrackInfo(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillTriggerInfo(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillDetectorVariables(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillParticleFlow(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCastor(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCastorDebug(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);
      void fillNoiseInfo(ZeroBiasEvent&, const edm::Event&, const edm::EventSetup&);

      edm::InputTag triggerResultsTag_;
      std::vector<std::string> hltPathNames_;
      edm::InputTag PVtxCollectionTag_;
      edm::InputTag trackTag_;
      edm::InputTag caloTowerTag_; 
      double energyThresholdHB_;
      double energyThresholdHE_;
      double energyThresholdHF_;
      double energyThresholdEB_;
      double energyThresholdEE_;
      edm::InputTag pfTag_;
      double pTPFThresholdCharged_;
      double energyPFThresholdBar_;
      double energyPFThresholdEnd_;
      double energyPFThresholdHF_;
      edm::InputTag castorHitsTag_;
      bool RunA_;
      bool RunB_;
      double fCGeVCastor_;

      std::string selectionPathName_;
      TH1F *hltTriggerPassHisto_,*hltTriggerNamesHisto_;
      TH1F *CastorChannelHisto_;
      TH1F *histo_castor_channels;
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
