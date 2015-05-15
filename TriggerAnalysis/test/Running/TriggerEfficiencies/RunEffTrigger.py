#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

os.system("./TriggerEfficiency \"/afs/cern.ch/user/d/dmf/private/uerjpc02/Samples/TriggerAnalysisMuon/30_march_2015/Muon_P2/MuonP2.root\" \"histo_trigger_eff_muonp2.root\" \"triggerAnalysisTTree\" \"Muon\" 8 3") # Ref: HLT_Ele10_SW_L1R(8), Trigger: HLT_Mu9(3)

os.system("./TriggerEfficiency \"/afs/cern.ch/user/d/dmf/private/uerjpc02/Samples/TriggerAnalysisMuon/30_march_2015/Muon_P3/MuonP3.root\" \"histo_trigger_eff_muonp3.root\" \"triggerAnalysisTTree\" \"Muon\" 12 8") # Ref: HLT_Ele10_LW_L1R_v*(12), Trigger: HLT_Mu15(8)

os.system("./TriggerEfficiency \"/afs/cern.ch/user/d/dmf/private/uerjpc02/Samples/TriggerAnalysisElectron/30_march_2015/Electron2010B/ElectronB.root\" \"histo_trigger_eff_electron.root\" \"triggerAnalysisTTree\" \"Electron\" 16 9") # Ref: HLT_Mu15_v*(16), Trigger: HLT_Ele17_SW_CaloEleId_L1R(9)

#----------------------------------------------------------->>>

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
