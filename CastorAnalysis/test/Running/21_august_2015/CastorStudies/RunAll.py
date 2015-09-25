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

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/Electron2010B/ElectronB.root\" \"diffractiveWAnalysisTTree\" \"histo_Electron2010B_Reco.root\" \"trigger_all_electron\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\" \"no_fit_castor\" \"ToFit.root\"")

#
# Muon Run2010B
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/Muon_P2/MuonP2.root\" \"diffractiveWAnalysisTTree\" \"histo_muon2010B_p2_HLT_Mu9_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\" \"no_fit_castor\" \"ToFit.root\"")

#
# Muon Run2010B
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/Muon_P3/MuonP3.root\" \"diffractiveWAnalysisTTree\" \"histo_muon2010B_p3_HLT_Mu15_Reco.root\" \"trigger_all_muon\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"no_weight\" 1.0 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP3_W.root\" \"no_fit_castor\" \"ToFit.root\"")

#
# Pomwig Plus WtoENu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoE_Plus/SingleWtoEPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_electron_Reco.root\" \"no_trigger_correction\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.013734768233469 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\" \"fit_castor\" \"ToFit.root\"")

#
# Pomwig Minus WtoENu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoE_Minus/SingleWtoEMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_electron_Reco.root\" \"no_trigger_correction\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.014464369283855 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\" \"fit_castor\" \"ToFit.root\"")


#
# Pomwig Plus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Plus/SingleWtoMuPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_muonP2_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.002749051757767 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\" \"fit_castor\" \"ToFit.root\"")


#
# Pomwig Plus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Plus/SingleWtoMuPlus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigPlus_muonP3_Reco.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.012941528859207 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP3_W.root\" \"fit_castor\" \"ToFit.root\"")

#
# Pomwig Minus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Minus/SingleWtoMuMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_muon_Reco_P2.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.002395667944446 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\" \"fit_castor\" \"ToFit.root\"")

#
# Pomwig Minus WtoMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/PomwigWtoMu_Minus/SingleWtoMuMinus.root\" \"diffractiveWAnalysisTTree\" \"histo_PomwigMinus_muon_Reco_P3.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.011277927289848 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP3_W.root\" \"fit_castor\" \"ToFit.root\"")


#
# WToENu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoENu/WtoE.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoENu_Reco_pu.root\" \"no_trigger_correction\" 0 25.0 25.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.128010102029489 \"RecoElectron\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_WElectron_71mb.root\" \"fit_castor\" \"ToFit.root\"")

#
# WToMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoMuNu/WtoMu.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoMuNu_Reco_pu_P2.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.024618345583322 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP2_W.root\" \"fit_castor\" \"ToFit.root\"")

#
# WToMuNu
#

os.system("./DiffractiveW \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoMuNu/WtoMu.root\" \"diffractiveWAnalysisTTree\" \"histo_WtoMuNu_Reco_pu_P3.root\" \"no_trigger_correction\" 0 20.0 20.0 1 \"no_multiple_pileup\" \"mc_lumi_pu_weight\" 0.115894154750768 \"RecoMuon\" 1.5 -1.0 \"channelcorrector.root\" \"pf\" \"puMC.root\" \"pu_71mbMuonP3_W.root\" \"fit_castor\" \"ToFit.root\"")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
