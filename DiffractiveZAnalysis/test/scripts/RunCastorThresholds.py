#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

# Collisions

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA/ZeroBiasA.root\" \"histo_Castor_threshold_p1_collisions.root\" \"diffractiveZAnalysisTTreeBefore\" \"collisions\" 136035 138559")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA/ZeroBiasA.root\" \"histo_Castor_threshold_p2_collisions.root\" \"diffractiveZAnalysisTTreeBefore\" \"collisions\" 138560 141955")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasB/ZeroBiasB.root\" \"histo_Castor_threshold_p3_collisions.root\" \"diffractiveZAnalysisTTreeBefore\" \"collisions\" 141956 149291")

# Unpaired

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA/ZeroBiasA.root\" \"histo_Castor_threshold_p1_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpaired\" 136035 138559")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA/ZeroBiasA.root\" \"histo_Castor_threshold_p2_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpaired\" 138560 141955")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasB/ZeroBiasB.root\" \"histo_Castor_threshold_p3_unpaired.root\" \"diffractiveZAnalysisTTreeBefore\" \"unpaired\" 141956 149291")

# No Collisions 

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA/ZeroBiasA.root\" \"histo_Castor_threshold_p1_nocollisions.root\" \"diffractiveZAnalysisTTreeBefore\" \"no_collisions\" 136035 138559")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasA/ZeroBiasA.root\" \"histo_Castor_threshold_p2_nocollisions.root\" \"diffractiveZAnalysisTTreeBefore\" \"no_collisions\" 138560 141955")

os.system("./CastorThreshold \"/storage1/dmf/Samples/DiffractiveZ/5_October_2013/ZeroBias/ZeroBiasB/ZeroBiasB.root\" \"histo_Castor_threshold_p3_nocollisions.root\" \"diffractiveZAnalysisTTreeBefore\" \"no_collisions\" 141956 149291")

#---------->

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
