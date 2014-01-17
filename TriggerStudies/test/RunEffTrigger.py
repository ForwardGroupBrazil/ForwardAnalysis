#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''


os.system("./TriggerEffCom \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"histo_effTriggerMultijetsRunB_RefOr30U_And30U_castor.root\" \"exclusiveDijetsAnalysisTTree\" 3 4 1 1 1 1.0")

#----------------------------------------------------------->>>

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
