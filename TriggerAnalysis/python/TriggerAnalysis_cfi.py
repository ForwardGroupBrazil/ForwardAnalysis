import FWCore.ParameterSet.Config as cms

from ForwardAnalysis.TriggerAnalysis.pfThresholds_cfi import pfThresholds

TriggerAnalysis = cms.PSet(
    hltPaths = cms.vstring(''),
    TriggerResultsTag = cms.InputTag("TriggerResults::HLT"),
    electronTag = cms.InputTag("gsfElectrons"),
    muonTag = cms.InputTag("muons"),
    metTag = cms.InputTag("pfMet"),
    PVtxCollectionTag=cms.InputTag('goodOfflinePrimaryVertices'),
    TrackTag = cms.InputTag("analysisTracks")
    )
