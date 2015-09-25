#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

os.system("./CastorBackground \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoMuNu/WtoMu.root\" \"diffractiveWAnalysisTTree\" \"/storage1/dmf/Samples/ZeroBias/30_march_2015/ZeroBiasA/ZeroBiasA.root\" \"zerobiasTTree\" \"histo_mc.root\" \"mc\" 1.5 1.")

os.system("./CastorBackground \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/WtoMuNu/WtoMu.root\" \"diffractiveWAnalysisTTree\" \"/storage1/dmf/Samples/ZeroBias/30_march_2015/ZeroBiasA/ZeroBiasA.root\" \"zerobiasTTree\" \"histo_mixture.root\" \"mixture\" 1.5 1.")

os.system("./CastorBackground \"/storage1/dmf/Samples/DiffractiveW/11_june_2015/Muon_P3/MuonP3.root\" \"diffractiveWAnalysisTTree\" \"/storage1/dmf/Samples/ZeroBias/30_march_2015/ZeroBiasA/ZeroBiasA.root\" \"zerobiasTTree\" \"histo_data.root\" \"data\" 1.5 1.")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
