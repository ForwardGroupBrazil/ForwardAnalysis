#ifndef DiffractiveWAnalysis_h
#define DiffractiveWAnalysis_h

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

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include <vector>
#include <string>
#include <map>

class DiffractiveWEvent;
class TH1F;
class TH2F;

namespace diffractiveWAnalysis {

  class DiffractiveWAnalysis {
    public:
      typedef DiffractiveWEvent event_type;
      static const char* name;

      DiffractiveWAnalysis() {} 
      DiffractiveWAnalysis(const edm::ParameterSet&);
      ~DiffractiveWAnalysis();

      void begin();
      void begin(const edm::Run&, const edm::EventSetup&);
      void fill(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void end();
    private:

      void setTFileService();
      void fillTriggerInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillTracksInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillGenInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillDetectorVariables(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillVariables(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCastor(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCastorDebug(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillZDC(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillDetectorEnergyEtaInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillElectronsInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillMuonsInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillMETInfo(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      void fillCollections(DiffractiveWEvent&, const edm::Event&, const edm::EventSetup&);
      
      template <class T, class W>
      double WBoson(T lepton, W met);

      edm::InputTag triggerResultsTag_;
      std::vector<std::string> hltPathNames_;
      edm::InputTag electronTag_;
      edm::InputTag muonTag_;
      edm::InputTag metTag_;
      edm::InputTag patmetTag_;
      edm::InputTag pfTag_;
      edm::InputTag genTag_;
      edm::InputTag PVtxCollectionTag_;
      edm::InputTag castorHitsTag_;
      edm::InputTag zdcHitsTag_;
      bool RunCastor_;
      bool RunZDC_;
      bool RunMC_;
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

      std::string selectionPathName_;

      std::vector<const reco::Muon*> MuonVector;
      std::vector<const reco::GsfElectron*> ElectronVector;
      std::vector<const pat::Muon*> PatMuonVector;
      std::vector<const pat::Electron*> PatElectronVector;
      std::vector<const reco::PFMET*> NeutrinoVector;
      std::vector<const pat::MET*> PatNeutrinoVector;
      std::vector<const reco::PFCandidate*> PFMuonVector;
      std::vector<const reco::PFCandidate*> PFElectronVector;

      TH1F *hltTriggerPassHisto_,*hltTriggerNamesHisto_;
      TH1F *CastorChannelHisto_;
      TH1F *histo_castor_channels;
      std::vector<TH1F*> m_hVector_histo_castor_channels;
      int indexE, indexM, indexpE, indexpM;
  };

} // namespace
#endif 
