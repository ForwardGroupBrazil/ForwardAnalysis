#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

os.system("./mcPU_Distributions \"/storage/dmf/uerj-1/Samples/DiffractiveW/25_january_2015/WtoMuNu/pythiaWtoMu.root\" \"puMC.root\" \"diffractiveWAnalysisTTree\" 25")

#----------------------------------------------------------->>>

print '\n'
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
