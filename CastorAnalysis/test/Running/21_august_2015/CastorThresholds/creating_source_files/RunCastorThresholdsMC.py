#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''


# Unpaired

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/9_december_2013/MinBiasMC/MinBiasMC.root\" \"histo_Castor_MinBias_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpairedmc\" 1 1")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
