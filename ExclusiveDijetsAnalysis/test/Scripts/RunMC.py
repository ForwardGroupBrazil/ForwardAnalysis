#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

#
# Pythia PU
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Pythia/Pythia6QCD.root\" \"histo_pythia6pu_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 167651.4 \"mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Herwig PU
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Herwig/HerwigQCD.root\" \"histo_herwigpu_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 123120.5 \"mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")


#----------------------------------------------------------->>>

#
# Pythia PU0
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Pythia/Pythia6QCD.root\" \"histo_pythia6pu_pT60_60_pu0.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 1139791.5 \"mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Herwig PU0
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Herwig/HerwigQCD.root\" \"histo_herwigpu_pT60_60_pu0.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 826050.4 \"mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")


#----------------------------------------------------------->>>

#
# ExHume 
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/ExHume/ExHume.root\" \"histo_exhume.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 0.0368  \"no_mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Pomwig 
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Pomwig/PomwigPDE.root\" \"histo_pomwig.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 9.23  \"no_mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Pompyt Plus 
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/PompytPlus/PompytPlus.root\" \"histo_pompytplus.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 105.4 \"no_mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Pompyt Minus 
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/PompytMinus/PompytMinus.root\" \"histo_pompytminus.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 100.2 \"no_mc_event_weight\" 1 0 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
