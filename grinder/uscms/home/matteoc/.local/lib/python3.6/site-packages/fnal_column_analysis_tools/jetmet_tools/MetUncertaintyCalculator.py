#if (jetPtJESUp > 20) {
#    MetXCorrjesUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
#    MetYCorrjesUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
#}
#if (jetPtJESDown > 20) {
#    MetXCorrjesDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
#    MetYCorrjesDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
#}
#if (jetPtJERUp > 20) {
#    MetXCorrjerUp += -1 * (thisJetJERUp.Px() - thisJet.Px());
#    MetYCorrjerUp += -1 * (thisJetJERUp.Py() - thisJet.Py());
#}
#if (jetPtJERDown > 20) {
#    MetXCorrjerDown += -1 * (thisJetJERDown.Px() - thisJet.Px());
#    MetYCorrjerDown += -1 * (thisJetJERDown.Py() - thisJet.Py());
#}

from ..analysis_objects.JaggedCandidateArray import JaggedCandidateArray

from fnal_column_analysis_tools.util import numpy as np
from fnal_column_analysis_tools.util import awkward

# returns corrected MET.x, corrected MET.y
def calculateType1MetXY(pfmet_in, correctedAK4Jets):
    pass
