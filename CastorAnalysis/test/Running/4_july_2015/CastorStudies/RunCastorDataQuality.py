#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

#
# Electron Run2010B
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/Electron2010B/ElectronB.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_Electron2010B_Reco.root\" \"trigger_all_electron\" 0 25.0 25.0 1 \"no_weight\" 1.0 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_ZElectron_71mb.root\" ")


#
# Muon Run2010B
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/Muon_P2/MuonP2.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_muon2010B_p2_HLT_Mu9_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_weight\" 1.0 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP2_Z.root\" ")


#
# Muon Run2010B
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/Muon_P3/MuonP3.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_muon2010B_p3_HLT_Mu15_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_weight\" 1.0 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP3_Z.root\" ")


#
# Pomwig Plus DyToEE
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoEE_Plus/SingleZtoEEPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigPlus_electron_Reco.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"mc_lumi_pu_weight\" 0.001315134506364 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_ZElectron_71mb.root\" ")


#
# Pomwig Minus DyToEE
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoEE_Minus/SingleZtoEEMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigMinus_electron_Reco.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"mc_lumi_pu_weight\" 0.001329954557739 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_ZElectron_71mb.root\" ")

#
# Pomwig Plus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Plus/SingleZtoMuMuPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigPlus_muonP2_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 0.000209534870222 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP2_Z.root\"")


#
# Pomwig Plus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Plus/SingleZtoMuMuPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigPlus_muonP3_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 0.000986788687006 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP3_Z.root\" ")

#
# Pomwig Minus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Minus/SingleZtoMuMuMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigMinus_muonP2_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 0.000261107336324 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP2_Z.root\" ")

#
# Pomwig Minus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Minus/SingleZtoMuMuMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigMinus_muonP3_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 0.001229665331152 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP3_Z.root\" ")


#
# DyToEE
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/ZtoEE/ZtoEE.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_DyToEE_Reco_pu.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"mc_lumi_pu_weight\" 1. \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_ZElectron_71mb.root\" ")


#
# DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/ZtoMuMu/ZtoMuMu.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_DyToMuMu_Reco_pu_P2.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP2_Z.root\" ")


#
# DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/ZtoMuMu/ZtoMuMu.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_DyToMuMu_Reco_pu_P3.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP3_Z.root\" ")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
