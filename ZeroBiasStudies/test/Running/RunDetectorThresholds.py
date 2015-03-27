#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''


# Unpaired

os.system("./DetectorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasB_v3/ZeroBiasB_v3.root\" \"histo_Castor_threshold_p3_unpaired.root\" \"zerobiasTTree\" \"unpaired\" 141956 149291")

os.system("./DetectorThreshold \"/storage1/dmf/Samples/DiffractiveZ/9_december_2013/MinBiasMC/MinBiasMC.root\" \"histo_Castor_MinBias_unpaired.root\" \"zerobiasTTree\" \"unpairedmc\" 1 1")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
