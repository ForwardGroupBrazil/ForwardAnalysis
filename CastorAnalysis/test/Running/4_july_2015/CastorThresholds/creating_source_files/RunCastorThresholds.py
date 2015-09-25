#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''


# Unpaired

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA_v3/ZeroBiasA_v3.root\" \"histo_Castor_threshold_p1_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpaired\" 136035 138559")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA_v3/ZeroBiasA_v3.root\" \"histo_Castor_threshold_p2_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpaired\" 138560 141955")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasB_v3/ZeroBiasB_v3.root\" \"histo_Castor_threshold_p3_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpaired\" 141956 149291")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
