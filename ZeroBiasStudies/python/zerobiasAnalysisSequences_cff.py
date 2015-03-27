import FWCore.ParameterSet.Config as cms

from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff import *
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
es_prefer_l1GtTriggerMaskAlgoTrig = cms.ESPrefer("L1GtTriggerMaskAlgoTrigTrivialProducer","l1GtTriggerMaskAlgoTrig")
es_prefer_l1GtTriggerMaskTechTrig = cms.ESPrefer("L1GtTriggerMaskTechTrigTrivialProducer","l1GtTriggerMaskTechTrig")

from ForwardAnalysis.AnalysisSequences.primaryVertexFilter_cfi import *
from ForwardAnalysis.AnalysisSequences.filterScraping_cfi import *
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *
from ForwardAnalysis.ZeroBiasStudies.zerobiasHLTPaths_cfi import *

from ForwardAnalysis.Utilities.analysisTracks_cfi import *
from ForwardAnalysis.Utilities.selectTracksAssociatedToPV_cfi import *
selectTracksAssociatedToPV.src = "analysisTracks"
selectTracksAssociatedToPV.vertexTag = "goodOfflinePrimaryVertices"
selectTracksAssociatedToPV.maxDistanceFromVertex = 0.4

from ForwardAnalysis.Utilities.tracksOutsideJets_cfi import *
tracksOutsideJets.src = "selectTracksAssociatedToPV" 
tracksOutsideJets.JetTag = "ak5PFJets"
tracksOutsideJets.JetConeSize = 0.5

from ForwardAnalysis.AnalysisSequences.tracksTransverseRegion_cfi import *
tracksTransverseRegion.src = "selectTracksAssociatedToPV"
tracksTransverseRegion.JetTag = "ak5PFJets"

from ForwardAnalysis.Utilities.trackMultiplicity_cfi import *
trackMultiplicityTransverseRegion = trackMultiplicity.clone( src = "tracksTransverseRegion" ) 

offlineSelection = cms.Sequence(primaryVertexFilter)
noiseRejection = cms.Sequence(filterScraping + HBHENoiseFilter)
eventSelection = cms.Sequence(noiseRejection + zerobiasHLTFilter)

tracks = cms.Sequence(analysisTracks*
                      selectTracksAssociatedToPV*
                      tracksOutsideJets*
                      tracksTransverseRegion)

analysisSequences = cms.Sequence(tracks)
