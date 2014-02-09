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

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_PYTHIA/LowPU/mc_total/mc_pythia_lowpu_0.root\" \"histo_pythia6pu_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 55781.0683578044 \"mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Herwig PU
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_HERWIG/mc_total/mc_herwig_lowpu_0.root\" \"histo_herwigpu_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 60060.4549155978 \"mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")


#----------------------------------------------------------->>>

#
# Pythia PU0
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_PYTHIA/LowPU/mc_total/mc_pythia_lowpu_0.root\" \"histo_pythia6pu_pT60_60_pu0.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 365676.168663604 \"mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Herwig PU0
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_HERWIG/mc_total/mc_herwig_lowpu_0.root\" \"histo_herwigpu_pT60_60_pu0.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 402676.822340036 \"mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")


#----------------------------------------------------------->>>

#
# ExHume 
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_EXHUME/mc_exhume_0.root\" \"histo_exhume.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 0.023802870674982 \"no_mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Pomwig 
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_POMWIG/mc_pomwig_dpe_0.root\" \"histo_pomwig.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 9.58245938037161 \"no_mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Pompyt Plus 
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_POMPYT/mc_pompyt_plus_0.root\" \"histo_pompytplus.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 109.678367396884 \"no_mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")

#
# Pompyt Minus 
#

os.system("./ExclusiveDijet \"/storage1/eliza/PATTuples_Mar2013/mc_POMPYT/mc_pompyt_minus_0.root\" \"histo_pompytminus.root\" \"exclusiveDijetsAnalysisTTree\" \"no_trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoHerwig.root\" \"effMinBias2010RunB_castor.root\" \"histo_effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"mc_lumi_weight\" 104.729539408482 \"no_mc_event_weight\" 1 0 60 60 \"no_castor\" \"channelcorrection.root\" 1.0 \"preselection\"")


print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
