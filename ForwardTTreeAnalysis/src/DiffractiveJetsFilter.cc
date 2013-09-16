//
// Author: dmf@cern.ch
//


// system include files
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//JetCollection
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//deltaR
#include "DataFormats/Math/interface/deltaR.h"

//ROOT
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>

using namespace edm;
using namespace std;
using namespace math;
using namespace reco;

class DiffractiveJetsFilter : public edm::EDFilter
{
  public:
    explicit DiffractiveJetsFilter(const edm::ParameterSet&);
    ~DiffractiveJetsFilter();

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    string calAlgoFilter;
    double PtCut1;
    double PtCut2;

};

DiffractiveJetsFilter::DiffractiveJetsFilter(const edm::ParameterSet& iConfig) :
  calAlgoFilter( iConfig.getUntrackedParameter<string>("calAlgoFilter")),
  PtCut1( iConfig.getParameter<double>( "PtCut1" )),
  PtCut2( iConfig.getParameter<double>( "PtCut2" ))
{
}


DiffractiveJetsFilter::~DiffractiveJetsFilter()
{
}

  bool
DiffractiveJetsFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool FilterResult = false;

  edm::Handle<edm::View<reco::Jet> > jetsis5;
  Handle<reco::Jet> jetsis5;
  iEvent.getByLabel(calAlgoFilter,jetsis5);

  int jet5size = jetsis5->size();
  int ijet5;

  double jet1Pt = 0;
  double jet2Pt = 0;
  const reco::Jet* jet1=NULL;
  const reco::Jet* jet2=NULL;

  if (jet5size >= 2){

    for(ijet5=0; ijet5 < jet5size; ++ijet5){
      const reco::Jet* jetAll = &((*jetsis5)[ijet5]);

      if (jetAll==NULL) continue;
      if (jet1==NULL) {jet1=jetAll; continue;}
      if (jetAll->pt()>jet1->pt()) {
	jet2=jet1;
	jet1=jetAll;
	continue;
      }

      if (jet2==NULL) {jet2=jetAll; continue;}
      if (jetAll->pt()>jet2->pt()) jet2 = jetAll;
    }

    jet1Pt = jet1->pt();
    jet2Pt = jet2->pt();

    if (jet1Pt >= PtCut1 && jet2Pt >= PtCut2){
      FilterResult = true;
    }

    else {
      FilterResult = false;
    }

  }
  return FilterResult;
}

void DiffractiveJetsFilter::beginJob(){
}

void DiffractiveJetsFilter::endJob(){
}

DEFINE_FWK_MODULE(DiffractiveJetsFilter);
