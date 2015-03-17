import FWCore.ParameterSet.Config as cms

import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
triggerHLTFilter = copy.deepcopy(hltHighLevel)
triggerHLTFilter.throw = False
triggerHLTFilter.HLTPaths =['']
