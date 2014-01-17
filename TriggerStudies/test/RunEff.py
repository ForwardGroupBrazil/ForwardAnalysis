#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

os.system("./EffMacroCom \"ZeroBiasB_Example.root\" \"histo_effCutsMinBias2010RunB_castor.root\" \"exclusiveDijetsAnalysisTTreeBefore\" 1 0 1 1 1 1.0")

#----------------------------------------------------------->>>

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
