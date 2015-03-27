import FWCore.ParameterSet.Config as cms

from ForwardAnalysis.ZeroBiasStudies.pfThresholds_cfi import pfThresholds

ZeroBiasAnalysis = cms.PSet(
    TriggerResultsTag = cms.InputTag("TriggerResults::HLT"),
    hltPaths = cms.vstring(''),
    PVtxCollectionTag=cms.InputTag('goodOfflinePrimaryVertices'),
    TrackTag = cms.InputTag("analysisTracks"),
    CaloTowerTag = cms.InputTag("towerMaker"),
    energyThresholdHB = cms.double(0.),
    energyThresholdHE = cms.double(0.),
    energyThresholdHF = cms.double(0.),
    energyThresholdEB = cms.double(0.),
    energyThresholdEE = cms.double(0.),
    pfTag = cms.InputTag("particleFlow"),
    pTPFThresholdCharged = cms.double(0.),
    energyPFThresholdBar = cms.double(0.),
    energyPFThresholdEnd = cms.double(0.),
    energyPFThresholdHF = cms.double(0.),
    castorHitsTag = cms.InputTag("castorRecHitCorrector"),
    RunA = cms.untracked.bool(True), # When RunA == RunB == True/False (All Channels). RunA=True and RunB=False (Good Channels Run A)/ RunA=False and RunB=True (Good Channels RunB)
    RunB = cms.untracked.bool(True),
    fCGeVCastor = cms.double(0.015) # for MC we need use 0.9375 (0.015/0.016)
    )
