#!/usr/bin/env python
import os

print ''
print '@@@@@@@@@@@@'
print 'Run Analysis'
print '@@@@@@@@@@@@'
print ''

#
# MultiJet
#

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"plushisto_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_triggercorr_effcorr_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTreePFShiftedUp\" \"trigger\" \"no_multiple_pileup\" \"plus\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"sigmaPlusHLTDijet50_And30U_pT60_castor.root\" \"cut_correction\" \"trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"minushisto_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_triggercorr_effcorr_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTreePFShiftedDown\" \"trigger\" \"no_multiple_pileup\" \"minus\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"sigmaMinusHLTDijet50_And30U_pT60_castor.root\" \"cut_correction\" \"trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"nonehisto_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_triggercorr_effcorr_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"cut_correction\" \"trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

######--------------------------------------------------------------------------------------->

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"nonehisto_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"plus_jets_histo_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"trigger\" \"no_multiple_pileup\" \"plus\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"sigmaMinusHLTDijet50_And30U_pT60_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"minus_jets_histo_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"trigger\" \"no_multiple_pileup\" \"minus\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"plus_pf_histo_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTreePFShiftedUp\" \"trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"sigmaMinusHLTDijet50_And30U_pT60_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"minus_pf_histo_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTreePFShiftedDown\" \"trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"no_trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"nonehisto_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_effcorr_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"cut_correction\" \"no_trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

os.system("./ExclusiveDijet \"/storage1/eliza/Samples/ExclusiveDijets/Multijets/Multijets.root\" \"nonehisto_multijet_HLTExclDiJet30U_AND_HFPreSel_GoodVertex_triggerRefDijet50And30_triggercorr_68mb_24bin_pT60_60.root\" \"exclusiveDijetsAnalysisTTree\" \"trigger\" \"no_multiple_pileup\" \"none\" \"no_pileup_correction\" \"pu_MultiJets_completejson_excldijets_march2013_68mb_bin24.root\" \"puHistoPythia6.root\" \"effMinBias2010RunB_castor.root\" \"effTriggerMultijetsRunB_RefDijetAve50U_And30U_castor.root\" \"no_cut_correction\" \"trigger_correction\" \"no_mc_lumi_weight\" 1.0 \"no_mc_event_weight\" 1 7 60 60 \"castor_no_correction\" \"channelcorrection.root\" 1.0 \"preselection\"")

print ''
print '@@@@@@@@@@@@@@@@@@'
print 'Finishing Analysis'
print '@@@@@@@@@@@@@@@@@@'
print ''
