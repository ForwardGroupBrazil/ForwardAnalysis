#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''


#
# Pomwig Plus DyToEE
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoEE_Plus/SingleZtoEEPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigPlus_electron_Reco.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"mc_lumi_pu_weight\" 1. \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_ZElectron_71mb.root\" ")


#
# Pomwig Minus DyToEE
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoEE_Minus/SingleZtoEEMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigMinus_electron_Reco.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"mc_lumi_pu_weight\" 1. \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_ZElectron_71mb.root\" ")

#
# Pomwig Plus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Plus/SingleZtoMuMuPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigPlus_muonP2_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP2_Z.root\"")


#
# Pomwig Plus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Plus/SingleZtoMuMuPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigPlus_muonP3_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP3_Z.root\" ")

#
# Pomwig Minus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Minus/SingleZtoMuMuMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigMinus_muonP2_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP2_Z.root\" ")

#
# Pomwig Minus DyToMuMu
#

os.system("./CastorAnalysis \"/storage/dmf/uerj-1/Samples/DiffractiveZ/11_june_2015/PomwigZtoMuMu_Minus/SingleZtoMuMuMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_castor_PomwigMinus_muonP3_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"puMCZ.root\" \"pu_71mbMuonP3_Z.root\" ")


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


#
# Pomwig Plus WtoENu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoE_Plus/SingleWtoEPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_electron_Reco.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\"")

#
# Pomwig Minus WtoENu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoE_Minus/SingleWtoEMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_electron_Reco.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\"")


#
# Pomwig Plus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Plus/SingleWtoMuPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_muonP2_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\"")


#
# Pomwig Plus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Plus/SingleWtoMuPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_muonP3_Reco.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.012941528859207 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP3_W.root\"")

#
# Pomwig Minus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Minus/SingleWtoMuMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_muon_Reco_P2.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\"")


#
# WToENu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoENu/WtoE.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoENu_Reco_pu.root\" \"no_trigger_nocorrection\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\"")

#
# WToMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoMuNu/WtoMu.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoMuNu_Reco_pu_P2.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\"")

#
# WToMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoMuNu/WtoMu.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoMuNu_Reco_pu_P3.root\" \"no_trigger_nocorrection\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 1. \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP3_W.root\"")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
