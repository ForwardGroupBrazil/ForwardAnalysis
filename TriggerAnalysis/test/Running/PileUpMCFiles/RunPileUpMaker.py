#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

os.system("./mcPU_Distributions \"/afs/cern.ch/user/d/dmf/private/uerj-1/Samples/DiffractiveZ/10_march_2015/ZtoEE/ZtoEE.root\" \"puMCZ.root\" \"diffractiveZAnalysisTTree\" 25")

#----------------------------------------------------------->>>

print '\n'
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
