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

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/Electron2010A/ElectronA.root\" \"diffractiveZAnalysisTTree\" \"histo_EG2010A_Reco.root\" \"trigger_all_electron\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"no_mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Electron Run2010B
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/Electron2010B/ElectronB.root\" \"diffractiveZAnalysisTTree\" \"histo_Electron2010B_Reco.root\" \"trigger_all_electron\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"no_mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Muon Run2010A
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/Muon_P1/Muon_P1.root\" \"diffractiveZAnalysisTTree\" \"histo_muon2010A_p1_HLT_Mu9_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Muon Run2010B
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/Muon_P2/Muon_P2.root\" \"diffractiveZAnalysisTTree\" \"histo_muon2010B_p2_HLT_Mu9_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Muon Run2010B
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/Muon_P3/Muon_P3.root\" \"diffractiveZAnalysisTTree\" \"histo_muon2010B_p3_HLT_Mu15_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Pompyt Plus Electron
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/PompytZPlus/PompytZPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_PompytPlus_electron_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Pompyt Minus Electron
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/PompytZMinus/PompytZMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_PompytMinus_electron_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Pompyt Plus Muon
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/PompytZPlus/PompytZPlus.root\" \"diffractiveZAnalysisTTree\" \"histo_PompytPlus_muon_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# Pompyt Minus Muon
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/PompytZMinus/PompytZMinus.root\" \"diffractiveZAnalysisTTree\" \"histo_PompytMinus_muon_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# DyToEE
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/DyToEE/DyToEE.root\" \"diffractiveZAnalysisTTree\" \"histo_DyToEE_Reco.root\" \"no_trigger_correction\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoElectron\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

#
# DyToMuMu
#

os.system("./DiffractiveZ \"/storage1/dmf/Samples/DiffractiveZ/10_january_2014/DyToMuMu/DyToMuMu.root\" \"diffractiveZAnalysisTTree\" \"histo_DyToMuMu_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_weight\" 1.0 \"RecoMuon\" -2.0 1.0 \"channelcorrection.root\" \"hf\"")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
