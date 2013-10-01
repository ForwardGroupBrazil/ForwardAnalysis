#flake8: noqa

'''
>>------------------<<
Dijets NTuple Producer
>>------------------<<

Goal:
Produce your diffractive dijets ntuple. 

Usage:
    cmsRun DijetsScaleMultiple.py

Example:
    cmsRun DijetsScaleMultiple.py Run=data_MuonA condition=B

Optional arguments:
    Run = data_exclusiveB, data_inclusiveB, MC_none, MC_PU, MC_FlatWeight or MC_FlatWeight_and_PU 
    condition = A or B for Castor remove channels depending of Run. If any, no conditions.

Authors: D. Figueiredo, E. Melo, A. Vilela and S. Fonseca
'''

import FWCore.ParameterSet.Config as cms
import os, sys
import atexit

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Run','data_Multijet',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: data or MC.")
options.register('condition','Full',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Channels CASTOR conditions.")
options.parseArguments()

process = cms.Process("Analysis")

class config: pass
config.verbose = True
config.writeEdmOutput = False
config.outputTTreeFile = 'DijetsPATTuple.root'
config.runPATSequences = True
config.comEnergy = 7000.0
config.trackAnalyzerName = 'trackHistoAnalyzer'
config.trackTagName = 'analysisTracks'
config.NumberOfEvents = -1

#
# Define Options to Run
#
########################

if options.Run == "data_Multijet":
  print("")
  print("#############")
  print("Data Multijet")
  print("#############")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  triggerlist = 'HLT_Jet15U','HLT_Jet15U_v*','HLT_Jet30U','HLT_Jet30U_v*','HLT_Jet50U','HLT_Jet50U','HLT_ExclDiJet30U_HFOR_v*','HLT_ExclDiJet30U_HFAND_v*','HLT_DiJetAve15U','HLT_DiJetAve15U_v*','HLT_DiJetAve30U','HLT_DiJetAve30U_v*','HLT_DiJetAve50U','HLT_DiJetAve50U_v*' 
  l1list = 'L1_ZeroBias','L1_SingleEG5'
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False

if options.Run == "data_Jet":
  print("")
  print("########")
  print("Data Jet")
  print("########")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  triggerlist = 'HLT_Jet15U','HLT_Jet15U_v*','HLT_Jet30U','HLT_Jet30U_v*','HLT_Jet50U','HLT_Jet50U','HLT_ExclDiJet30U_HFOR_v*','HLT_ExclDiJet30U_HFAND_v*','HLT_DiJetAve15U','HLT_DiJetAve15U_v*','HLT_DiJetAve30U','HLT_DiJetAve30U_v*','HLT_DiJetAve50U','HLT_DiJetAve50U_v*'
  l1list = 'L1_ZeroBias','L1_SingleEG5'
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False

elif options.Run == "data_JetMet":
  print("")
  print("#####################")
  print("Data Inclusive JetMet ")
  print("#####################")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  triggerlist = 'HLT_Jet15U','HLT_Jet30U','HLT_Jet50U'
  l1list = 'L1_ZeroBias','L1_SingleEG5'
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False

elif options.Run == "data_JetMetTau":
  print("")
  print("########################")
  print("Data Inclusive JetMetTau")
  print("########################")
  print("")
  config.globalTagNameData = 'GR_R_42_V23::All'
  config.TriggerOn = True
  triggerlist = 'HLT_Jet15U','HLT_Jet30U','HLT_Jet50U'
  l1list = 'L1_ZeroBias','L1_SingleEG5'
  config.runOnMC = False
  config.runPUMC = False
  config.runGen = False

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
#################################

if config.runOnMC:
    config.l1Paths = (l1list)
    config.hltPaths =(triggerlist)
    config.inputFileName = '/afs/cern.ch/work/d/dmf/public/TestSamples/PYTHIA6_QCD_15to3000_private_SL_RECO2010/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_cff_py_RAW2DIGI_L1RECO.root'

else:
    config.l1Paths = (l1list)
    config.hltPaths = (triggerlist)
    config.inputFileName = '/afs/cern.ch/work/d/dmf/public/TestSamples/MultiJets2010/Multijet2010.root'

#
# CMSSW Main Code
#
#################

process = cms.Process("Analysis")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(False),
	SkipEvent = cms.untracked.vstring('ProductNotFound')
	)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(config.NumberOfEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 'file:%s' % config.inputFileName )
    #duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

#
# Output
#
########

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(config.outputTTreeFile))


#
# Detector Conditions and Scales
#
################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#
# Global Tag Input
#
##################

if config.runOnMC: process.GlobalTag.globaltag = config.globalTagNameMC
else: process.GlobalTag.globaltag = config.globalTagNameData

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

#
# PAT Sequences
#
###############

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
# Open Analysis Module
#
######################

process.load("ForwardAnalysis.ExclusiveDijetsAnalysis.exclusiveDijetsAnalysisSequences_cff")
process.pfCandidateNoiseThresholds.src = "pfNoPileUpPFlow"
process.tracksTransverseRegion.JetTag = "selectedPatJetsPFlow"


from ForwardAnalysis.ForwardTTreeAnalysis.DiffractiveAnalysis_cfi import DiffractiveAnalysis
from ForwardAnalysis.ForwardTTreeAnalysis.ExclusiveDijetsAnalysis_cfi import ExclusiveDijetsAnalysis
from ForwardAnalysis.ForwardTTreeAnalysis.PATTriggerInfo_cfi import PATTriggerInfo
from ForwardAnalysis.ForwardTTreeAnalysis.DijetsTriggerAnalysis_cfi import DijetsTriggerAnalysis  
from ForwardAnalysis.ForwardTTreeAnalysis.PFCandInfo_cfi import PFCandInfo

#PATTriggerInfo.L1AlgoBitName =  config.l1Paths 
PATTriggerInfo.HLTAlgoBitName = config.hltPaths 
PATTriggerInfo.runALLTriggerPath = True

#
# Module with EventInfo + DiffractiveAnalysis + ExclusiveDijets 
#
###############################################################

process.DiffractiveJetsFilter = cms.EDFilter("DiffractiveJetsFilter",
        calAlgoFilter = cms.untracked.string("selectedPatJetsPFlow"),
        PtCut1 = cms.double(15.0),
        PtCut2 = cms.double(15.0)
        )

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

if options.condition=="A":
     print("")
     print(">>>> RunA Castor Conditions")
     print("")
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunA = True
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunB = False

elif options.condition=="B":
     print("")
     print(">>>> RunB Castor Conditions")
     print("")
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunA = False
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunB = True
else:
     print("")
     print(">>>> Full Castor")
     print("")
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunA = False
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.RunB = False

if config.runOnMC:
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.fCGeVCastor = 0.9375
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.accessMCInfo = True
     process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.AccessMCInfo = True
     process.gen_step = cms.Path(process.genChargedParticles+
                              process.genProtonDissociative*process.edmNtupleMxGen+
                              process.genStableParticles*
                              process.etaMaxGen+process.etaMinGen*
                              process.edmNtupleEtaMaxGen+process.edmNtupleEtaMinGen)


else:
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.accessMCInfo = False
     process.exclusiveDijetsAnalysisTTree.ExclusiveDijetsAnalysis.AccessMCInfo = False
     process.exclusiveDijetsAnalysisTTree.DiffractiveAnalysis.fCGeVCastor = 0.015

process.exclusiveDijetsAnalysisTTreePFShiftedUp = process.exclusiveDijetsAnalysisTTree.clone()
process.exclusiveDijetsAnalysisTTreePFShiftedUp.DiffractiveAnalysis.particleFlowTag = "pfCandidateNoiseThresholdsShiftedUp"
process.exclusiveDijetsAnalysisTTreePFShiftedUp.DiffractiveAnalysis.edmNtupleEtaMaxTag = "edmNtupleEtaMaxShiftedUp"
process.exclusiveDijetsAnalysisTTreePFShiftedUp.DiffractiveAnalysis.edmNtupleEtaMinTag = "edmNtupleEtaMinShiftedUp"
process.exclusiveDijetsAnalysisTTreePFShiftedUp.ExclusiveDijetsAnalysis.ParticleFlowTag = "pfCandidateNoiseThresholdsShiftedUp"

process.exclusiveDijetsAnalysisTTreePFShiftedDown = process.exclusiveDijetsAnalysisTTree.clone() 	 
process.exclusiveDijetsAnalysisTTreePFShiftedDown.DiffractiveAnalysis.particleFlowTag = "pfCandidateNoiseThresholdsShiftedDown" 	 
process.exclusiveDijetsAnalysisTTreePFShiftedDown.DiffractiveAnalysis.edmNtupleEtaMaxTag = "edmNtupleEtaMaxShiftedDown" 	 
process.exclusiveDijetsAnalysisTTreePFShiftedDown.DiffractiveAnalysis.edmNtupleEtaMinTag = "edmNtupleEtaMinShiftedDown" 	 
process.exclusiveDijetsAnalysisTTreePFShiftedDown.ExclusiveDijetsAnalysis.ParticleFlowTag = "pfCandidateNoiseThresholdsShiftedDown"

#
# Run Path. 
# If TriggerOn = True (Run with trigger)
#
########################################

process.castor_step = cms.Path(process.castorSequence)

if config.TriggerOn:
    process.analysis_diffractiveExclusiveDijetsAnalysisPATTriggerInfoShiftedTTree_step = cms.Path(
    process.analysisSequencesShiftedUp + process.analysisSequencesShiftedDown + process.analysisSequences + process.eventSelectionHLT + process.DiffractiveJetsFilter +
    process.exclusiveDijetsAnalysisTTreePFShiftedUp + process.exclusiveDijetsAnalysisTTreePFShiftedDown + process.exclusiveDijetsAnalysisTTree)

else: 
    process.analysis_diffractiveExclusiveDijetsAnalysisPATTriggerInfoShiftedTTree_step = cms.Path(
    process.analysisSequencesShiftedUp + process.analysisSequencesShiftedDown + process.analysisSequences + process.eventSelection + process.DiffractiveJetsFilter + process.exclusiveDijetsAnalysisTTreePFShiftedUp + process.exclusiveDijetsAnalysisTTreePFShiftedDown + process.exclusiveDijetsAnalysisTTree)

