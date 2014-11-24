#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

#
# EG Run2010A
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010A/ElectronWA.root\" \"diffractiveWAnalysisTTree\" \"histo_EG2010A_Reco.root\" \"trigger_all_electron\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Electron Run2010B
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Electron2010B/ElectronWB.root\" \"diffractiveWAnalysisTTree\" \"histo_Electron2010B_Reco.root\" \"trigger_all_electron\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Muon Run2010A
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P1/MuonWP1.root\" \"diffractiveWAnalysisTTree\" \"histo_muon2010A_p1_HLT_Mu9_Reco.root\" \"trigger\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Muon Run2010B
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P2/MuonWP2.root\" \"diffractiveWAnalysisTTree\" \"histo_muon2010B_p2_HLT_Mu9_Reco.root\" \"trigger\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Muon Run2010B
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/Muon_P3/MuonWP3.root\" \"diffractiveWAnalysisTTree\" \"histo_muon2010B_p3_HLT_Mu15_Reco.root\" \"trigger\" 1 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Pomwig Plus WtoENu
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Plus/pomwigWtoEPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_electron_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Pomwig Minus WtoENu
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoE_Minus/pomwigWtoEMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_electron_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Pomwig Plus WtoMuNu
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Plus/pomwigWtoMuPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_muon_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# Pomwig Minus WtoMuNu
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/PomwigWtoMu_Minus/pomwigWtoMuMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_muon_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# WToENu
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoENu/pythiaWtoE.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoENu_Reco.root\" \"no_trigger_correction\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

#
# WToMuNu
#

os.system("./DiffractiveW \"/storage/dmf/uerj-1/Samples/DiffractiveW/4_november_2014/WtoMuNu/pythiaWtoMu.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoMuNu_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"puData.root\" ")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
