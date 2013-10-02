#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

#
# Herwig 
#

os.system("./MC \"/storage1/eliza/PATTuples_Mar2013/mc_HERWIG/mc_total/mc_herwig_lowpu_0.root\" \"puHistoHerwig.root\" \"exclusiveDijetsAnalysisTTree/ProcessedTree\" 60 60 1 0 0 0 1 0 1 25")

#
# Pythia 
#

os.system("./MC \"/storage1/eliza/PATTuples_Mar2013/mc_PYTHIA/LowPU/mc_total/mc_pythia_lowpu_0.root\" \"puHistoPythia6.root\" \"exclusiveDijetsAnalysisTTree/ProcessedTree\" 60 60 1 0 0 0 1 0 1 25")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
