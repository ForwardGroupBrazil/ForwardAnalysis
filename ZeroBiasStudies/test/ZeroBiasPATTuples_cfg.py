#flake8: noqa

'''
>>---------------------<<
Zero Bias NTuple Producer
>>---------------------<<

Goal:
Produce your ZeroBias ntuple.

Usage:
    cmsRun ZeroBiasPATTuples_cfg.py Run=A

Example:
    cmsRun ZeroBiasPATTuples_cfg.py 

Optional arguments:
    Run = A, B, MC_none, MC_PU, MC_FlatWeight or MC_FlatWeight_and_PU

Authors: D. Figueiredo, R. Arciadiacono and N. Cartiglia
'''

import FWCore.ParameterSet.Config as cms
import os, sys
import atexit

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Run','Full',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: data or MC.")
options.register('condition','Full',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Channels CASTOR conditions.")
options.parseArguments()

process = cms.Process("Analysis")

class config: pass
config.verbose = True
config.writeEdmOutput = False
config.outputTTreeFile = 'ZeroBiasPATTuple.root'
config.runPATSequences = True
config.comEnergy = 7000.0
config.trackAnalyzerName = 'trackHistoAnalyzer'
config.trackTagName = 'analysisTracks'
config.NumberOfEvents = -1

#
# Define Options to Run
#
######################################################################################

if options.Run == "A":
  print("")
  print("#####")
  print("Run A")
  print("#####")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  l1list = 'L1_ZeroBias','L1_BptxMinus_NotBptxPlus'
  triggerlist = 'HLT_ZeroBias','HLT_L1_BPTX_PlusOnly','HLT_L1_BPTX_MinusOnly','HLT_L1_BPTX','HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_L1_BptxXOR_BscMinBiasOR','HLT_L1Tech_BSC_minBias_OR','HLT_L1Tech_BSC_minBias','HLT_L1Tech_BSC_halo','HLT_L1Tech_BSC_halo_forPhysicsBackground','HLT_L1Tech_BSC_HighMultiplicity'
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False

elif options.Run == "B":
  print("")
  print("#####")
  print("Run B")
  print("#####")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False
  l1list = 'L1_ZeroBias','L1_BptxMinus_NotBptxPlus'
  triggerlist = 'HLT_ZeroBias','HLT_L1_BPTX_PlusOnly','HLT_L1_BPTX_MinusOnly','HLT_L1_BPTX','HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_L1_BptxXOR_BscMinBiasOR','HLT_L1Tech_BSC_minBias_OR','HLT_L1Tech_BSC_minBias','HLT_L1Tech_BSC_halo','HLT_L1Tech_BSC_halo_forPhysicsBackground','HLT_L1Tech_BSC_HighMultiplicity'

elif options.Run == "Full":
  print("")
  print("############")
  print("RunA or RunB")
  print("############")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False
  l1list = 'L1_ZeroBias','L1_BptxMinus_NotBptxPlus'
  triggerlist = 'HLT_ZeroBias','HLT_L1_BPTX_PlusOnly','HLT_L1_BPTX_MinusOnly','HLT_L1_BPTX','HLT_L1_BscMinBiasOR_BptxPlusORMinus','HLT_L1_BptxXOR_BscMinBiasOR','HLT_L1Tech_BSC_minBias_OR','HLT_L1Tech_BSC_minBias','HLT_L1Tech_BSC_halo','HLT_L1Tech_BSC_halo_forPhysicsBackground','HLT_L1Tech_BSC_HighMultiplicity'

elif options.Run == "MC_FlatWeight_and_PU":
  print("")
  print("#####################")
  print("MC Flat Weight and PU")
  print("#####################")
  print("")
  config.globalTagNameMC = 'START42_V17D::All'
  config.TriggerOn = False
  triggerlist = 'no_trigger','no_trigger'
  l1list = 'no_trigger','no_trigger'
  config.runOnMC = True
  config.runPUMC = True
  config.runGen = True

elif options.Run == "MC_FlatWeight":
  print("")
  print("##############")
  print("MC Flat Weight")
  print("##############")
  print("")
  config.globalTagNameMC = 'START42_V17D::All'
  config.TriggerOn = False
  triggerlist = 'no_trigger','no_trigger'
  l1list = 'no_trigger','no_trigger'
  config.runOnMC = True
  config.runPUMC = False
  config.runGen = True

elif options.Run == "MC_PU":
  print("")
  print("#####")
  print("MC PU")
  print("#####")
  print("")
  config.globalTagNameMC = 'START42_V17D::All'
  config.TriggerOn = False
  triggerlist = 'no_trigger','no_trigger'
  l1list = 'no_trigger','no_trigger'
  config.runOnMC = True
  config.runPUMC = True
  config.runGen = False

elif options.Run == "MC_none":
  print("")
  print("#######")
  print("MC None")
  print("#######")
  print("")
  config.globalTagNameMC = 'START42_V17D::All'
  config.TriggerOn = False
  triggerlist = 'no_trigger','no_trigger'
  l1list = 'no_trigger','no_trigger'
  config.runOnMC = True
  config.runPUMC = False
  config.runGen = False

else:
  print("")
  print("")
  raise RuntimeError, "Unknown option. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
  print("")


print("")
print(">>> Input Options:")
print("Run with MC: %s" % config.runOnMC)
print("Run MC with Pile Up: %s" % config.runPUMC)
print("Run MC with Flat Weight: %s" % config.runGen)
print("Run with Trigger: %s" % config.TriggerOn)
if not config.runOnMC:
   print("Data Global Tag: " + config.globalTagNameData)
else:
   print("MC Global Tag: " + config.globalTagNameMC)
print("")

#
# Define Triggers and Input Files
#
######################################################################################

if config.runOnMC:
    config.l1Paths = (l1list)
    config.hltPaths =(triggerlist)
    config.inputFileName = 'root://eoscms.cern.ch//store/user/dmf/SamplesDebugCMSSW428/WMuNuPythia6.root'

else:
    config.l1Paths = (l1list)
    config.hltPaths = (triggerlist)
    config.inputFileName = '/storage/dmf/uerj-1/TestSamples/MinBias.root'

#
# CMSSW Main Code
#
######################################################################################

process = cms.Process("Analysis")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(False),
	SkipEvent = cms.untracked.vstring('ProductNotFound')
	)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 'file:%s' % config.inputFileName )
    #duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

#
# Output
#
######################################################################################

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(config.outputTTreeFile))


# Detector Conditions and Scales
#
######################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#
# Global Tag Input
#
######################################################################################

if config.runOnMC: process.GlobalTag.globaltag = config.globalTagNameMC
else: process.GlobalTag.globaltag = config.globalTagNameData

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

#
# PAT Sequences
#
######################################################################################

if config.runPATSequences:
    from ForwardAnalysis.Skimming.addPATSequences import addPATSequences
    addPATSequences(process,config.runOnMC)

    if config.runOnMC:
        process.patTrigger.addL1Algos = cms.bool( False )
        process.patJets.addTagInfos   = cms.bool( False )
    else:
        process.patTrigger.addL1Algos = cms.bool( True )
        process.patJets.addTagInfos   = cms.bool( True )

from ForwardAnalysis.Utilities.addCastorRecHitCorrector import addCastorRecHitCorrector
addCastorRecHitCorrector(process)

#
# CASTOR
#
######################################################################################

from ForwardAnalysis.Utilities.addCastorRecHitCorrector import addCastorRecHitCorrector
addCastorRecHitCorrector(process)

#
# Remove PAT MCMatching for Data
#
######################################################################################

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.coreTools import *

#if not config.runOnMC:
#    removeMCMatching(process, ['All'])

#
# PAT Muons and Electrons WorkFlow
#
######################################################################################

from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import *

from PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *

if not config.runOnMC:

    process.makePatElectrons = cms.Sequence(
        patElectrons
    )

    process.makePatMuons = cms.Sequence(
         patMuons
    )

else:
    process.makePatElectrons = cms.Sequence(
         electronMatch *
         patElectrons
    )

    process.makePatMuons = cms.Sequence(
         muonMatch*
         patMuons
    )


#
# Open Common Modules
#
######################################################################################

process.load("ForwardAnalysis.ZeroBiasStudies.zerobiasAnalysisSequences_cff")
process.tracksTransverseRegion.JetTag = "selectedPatJetsPFlow"

#
# Import PSET for each Module
#
######################################################################################

from ForwardAnalysis.ZeroBiasStudies.ZeroBiasAnalysis_cfi import ZeroBiasAnalysis

#
# Define Analyzers
#
######################################################################################

process.zerobiasHLTFilter.HLTPaths = config.hltPaths

process.zerobiasTTree = cms.EDAnalyzer("EventInfoZeroBiasAnalysisTTree",
        EventInfo = cms.PSet(
                    RunOnData = cms.untracked.bool(not config.runOnMC),
                    RunWithMCPU = cms.untracked.bool(config.runPUMC),
                    RunWithWeightGen = cms.untracked.bool(config.runGen)
        ),
        ZeroBiasAnalysis = ZeroBiasAnalysis
        )

process.zerobiasTTree.ZeroBiasAnalysis.hltPaths = config.hltPaths

if config.runOnMC:
     process.zerobiasTTree.ZeroBiasAnalysis.fCGeVCastor = 0.9375
     process.zerobiasTTree.ZeroBiasAnalysis.castorHitsTag = "castorreco"
else:
     process.zerobiasTTree.ZeroBiasAnalysis.fCGeVCastor = 0.015
     process.zerobiasTTree.ZeroBiasAnalysis.castorHitsTag = "castorRecHitCorrector"

if options.condition=="A":
     print("")
     print(">>>> RunA Castor Conditions")
     print("")
     process.zerobiasTTree.ZeroBiasAnalysis.RunA = True
     process.zerobiasTTree.ZeroBiasAnalysis.RunB = False

elif options.condition=="B":
     print("")
     print(">>>> RunB Castor Conditions")
     print("")
     process.zerobiasTTree.ZeroBiasAnalysis.RunA = False
     process.zerobiasTTree.ZeroBiasAnalysis.RunB = True
else:
     print("")
     print(">>>> Full Castor")
     print("")
     process.zerobiasTTree.ZeroBiasAnalysis.RunA = False
     process.zerobiasTTree.ZeroBiasAnalysis.RunB = False

#
# Define MC Access
#
######################################################################################

process.pat_Producer = cms.Path(process.makePatElectrons + process.makePatMuons)
process.castor_step = cms.Path(process.castorSequence)

if config.TriggerOn:
       print(">> With Trigger.")
       process.analysis_zerobiasAnalysisPATTriggerInfoTTree_step = cms.Path(
       process.analysisSequences + process.eventSelection + process.zerobiasTTree)

else:
       print(">> No Trigger.")
       process.analysis_zerobiasAnalysisPATTriggerInfoTTree_step = cms.Path(
       process.analysisSequences + process.noiseRejection + process.zerobiasTTree)
