#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''


os.system("./TriggerEfficiency \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"histo_trigger_eff.root\" \"TriggerAnalysisTTree\" \"Muon\" 3 4")

#----------------------------------------------------------->>>

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
