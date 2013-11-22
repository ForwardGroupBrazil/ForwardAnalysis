#flake8: noqa

'''
>>--------------------------------------<<
Exclusive Dijets Zero Bias NTuple Producer
>>--------------------------------------<<

Goal:
Produce your ZeroBias Exclusive Dijets ntuple.

Usage:
    cmsRun ZeroBiasPATTuplesExclusiveDijets_cfg.py Run=A

Example:
    cmsRun ZeroBiasPATTuplesExclusiveDijets_cfg.py 

Optional arguments:
    Run = A, B, MC_none, MC_PU, MC_FlatWeight or MC_FlatWeight_and_PU

Authors: D. Figueiredo and E. Melo
'''

import FWCore.ParameterSet.Config as cms
import os, sys
import atexit

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Run','Full',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: data or MC.")
options.parseArguments()

process = cms.Process("Analysis")

class config: pass
config.verbose = True
config.writeEdmOutput = False
config.outputTTreeFile = 'ZeroBiasExclusiveDijets.root'
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
  print("##############")
  print("Run A && Run B")
  print("##############")
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
    config.inputFileName = '/afs/cern.ch/work/d/dmf/public/TestSamples/DyToMuMuPU2010/DyToMuMu.root'

else:
    config.l1Paths = (l1list)
    config.hltPaths = (triggerlist)
    config.inputFileName = '/afs/cern.ch/work/d/dmf/public/TestSamples/ZeroBias2010/minimumBias2010.root'

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
# Open Common Modules
#
######################################################################################

process.load("ForwardAnalysis.ExclusiveDijetsAnalysis.exclusiveDijetsAnalysisSequences_cff")
#process.pfCandidateNoiseThresholds.src = "pfNoPileUpPFlow"
process.pfCandidateNoiseThresholds.src = "particleFlow"
process.tracksTransverseRegion.JetTag = "selectedPatJetsPFlow"

#
# Import PSET for each Module
#
######################################################################################

from ForwardAnalysis.ForwardTTreeAnalysis.DiffractiveAnalysis_cfi import DiffractiveAnalysis
from ForwardAnalysis.ForwardTTreeAnalysis.ExclusiveDijetsAnalysis_cfi import ExclusiveDijetsAnalysis
from ForwardAnalysis.ForwardTTreeAnalysis.PATTriggerInfo_cfi import PATTriggerInfo
from ForwardAnalysis.ForwardTTreeAnalysis.DijetsTriggerAnalysis_cfi import DijetsTriggerAnalysis
from ForwardAnalysis.ForwardTTreeAnalysis.PFCandInfo_cfi import PFCandInfo

#PATTriggerInfo.L1AlgoBitName =  config.l1Paths
PATTriggerInfo.HLTAlgoBitName = config.hltPaths
PATTriggerInfo.runALLTriggerPath = True

PATTriggerInfo.HLTAlgoBitName = config.hltPaths
PATTriggerInfo.runALLTriggerPath = True

#
# Define Filter
#
######################################################################################

process.DiffractiveJetsFilter = cms.EDFilter("DiffractiveJetsFilter",
        calAlgoFilter = cms.untracked.string("selectedPatJetsPFlow"),
        PtCut1 = cms.double(15.0),
        PtCut2 = cms.double(15.0)
        )

#
# Define Analyzers
#
######################################################################################

process.exclusiveDijetsAnalysisTTree = cms.EDAnalyzer("EventInfoPFCandInfoDiffractiveExclusiveDijetsAnalysisTTree",
        EventInfo = cms.PSet(
                    RunOnData = cms.untracked.bool(not config.runOnMC),
                    RunWithMCPU = cms.untracked.bool(config.runPUMC),
                    RunWithWeightGen = cms.untracked.bool(config.runGen)
        ),
        DiffractiveAnalysis = DiffractiveAnalysis,
        ExclusiveDijetsAnalysis = ExclusiveDijetsAnalysis,
        PFCandInfo = PFCandInfo
        )

process.exclusiveDijetsHLTFilter.HLTPaths = config.hltPaths

process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.hltPath = ''
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.trackTag = 'analysisTracks'
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.vertexTag = "goodOfflinePrimaryVertices"
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.particleFlowTag = "pfCandidateNoiseThresholds"
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.jetTag = "selectedPatJetsPFlow"
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.energyThresholdHF = 7.0
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.accessCastorInfo = True
process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.accessZDCInfo = False

process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.hltPaths = config.hltPaths
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.TrackTag = 'analysisTracks'
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.VertexTag = "goodOfflinePrimaryVertices"
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.ParticleFlowTag = "pfCandidateNoiseThresholds"
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.JetTag = "selectedPatJetsPFlow"
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.JetNonCorrTag = "ak5PFJets"
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.EnergyThresholdHF = 7.0
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.PFlowThresholds.Transition.hadronHF.energy = 7.0
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.PFlowThresholds.Transition.emHF.energy = 7.0
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.PFlowThresholds.Forward.hadronHF.energy = 7.0
process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.PFlowThresholds.Forward.emHF.energy = 7.0

if options.Run=="A":
     print("")
     print(">>>> RunA Castor Conditions")
     print("")
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunA = True
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunB = False

elif options.Run=="B":
     print("")
     print(">>>> RunB Castor Conditions")
     print("")
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunA = False
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunB = True

elif options.Run=="Full":
     print("")
     print(">>>> Full Castor")
     print("")
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunA = False
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunB = False

if config.runOnMC:
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.fCGeVCastor = 0.9375
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.accessMCInfo = True
     process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.AccessMCInfo = True
     process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.CastorRecHitTag = "castorreco"
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.castorRecHitTag ="castorreco"
     process.gen_step = cms.Path(process.genChargedParticles+
                              process.genProtonDissociative*process.edmNtupleMxGen+
                              process.genStableParticles*
                              process.etaMaxGen+process.etaMinGen*
                              process.edmNtupleEtaMaxGen+process.edmNtupleEtaMinGen)

else:
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.accessMCInfo = False
     process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.AccessMCInfo = False
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.fCGeVCastor = 0.015
     process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.CastorRecHitTag = "castorRecHitCorrector"
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.castorRecHitTag ="castorRecHitCorrector"

process.exclusiveDijetsAnalysisTTreeBefore = process.exclusiveDijetsAnalysisTTree.clone()
process.exclusiveDijetsAnalysisTTreeAfter = process.exclusiveDijetsAnalysisTTree.clone()

#
# Run Path.
# If TriggerOn = True (Run with trigger)
#
########################################

process.castor_step = cms.Path(process.castorSequence)

if config.TriggerOn:
       print(">> With Trigger.")
       process.analysis_diffractiveExclusiveDijetsAnalysisPATTriggerInfoTTree_step = cms.Path(
       process.analysisSequences + process.eventSelectionOnlyHLT + process.exclusiveDijetsAnalysisTTreeBefore + process.eventSelection +  process.exclusiveDijetsAnalysisTTreeAfter)

else:
       print(">> No Trigger.")
       process.analysis_diffractiveDiffractiveZAnalysisPATTriggerInfoTTree_step = cms.Path(
       process.analysisSequences + process.exclusiveDijetsAnalysisTTreeBefore + process.eventSelection + process.exclusiveDijetsAnalysisTTreeAfter)

