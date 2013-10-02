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

os.system("./MC_PU_Distribution \"/storage1/eliza/PATTuples_HF7_Nov2012/data_MCHERWIG/data/herwigppPF7_total_0.root\" \"puHistoHerwig.root\" \"exclusiveDijetsAnalysisTTree/ProcessedTree\" 60 60 1 0 0 0 1 0 1 25")

#
# Pythia 
#

os.system("./MC_PU_Distribution \"/storage1/eliza/PATTuples_HF7_Nov2012/data_MCPYTHIAPF7/pythia6PF7_total_0.root\" \"puHistoPythia6.root\" \"exclusiveDijetsAnalysisTTree/ProcessedTree\" 60 60 1 0 0 0 1 0 1 25")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
