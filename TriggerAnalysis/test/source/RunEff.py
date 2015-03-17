#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

os.system("./EffMacroCom \"/storage1/dmf/Samples/ExclusiveDijets/20_january_2014/ZeroBiasB/ZeroBiasB.root\" \"histo_effCutsMinBias2010RunB_castor.root\" \"exclusiveDijetsAnalysisTTreeBefore\" 1 0 1 1 1 1.0")

#----------------------------------------------------------->>>

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
