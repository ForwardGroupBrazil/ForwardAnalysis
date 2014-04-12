import FWCore.ParameterSet.Config as cms

import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
diffractiveWHLTFilter = copy.deepcopy(hltHighLevel)
diffractiveWHLTFilter.throw = False
diffractiveWHLTFilter.HLTPaths =['']
