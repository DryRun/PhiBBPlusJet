import os
import sys
import ROOT
ROOT.gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsRootUtils.so"))
ROOT.gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsAnalysisTools.so"))
ROOT.gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
from DAZSLE.PhiBBPlusJet.bacon_event_selector import *

def MakeSRSelector(jet_type, n2_ddt_cut=0., tau21_ddt_cut=None, jet_systematic="nominal", tag=None):
	selector_name = "EventSelector_SR"
	if jet_systematic != "nominal":
		selector_name += "_" + jet_systematic
	if tag:
		selector_name += "_" + tag
	event_selector = BaconEventSelector(selector_name)

	if jet_type == "AK8":
		event_selector.add_cut("Min_AK8Puppijet0_pt", {"Min_AK8Puppijet0_pt":450., "systematic":jet_systematic})
		event_selector.add_cut("Min_AK8Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("AK8Puppijet0_isTightVJet")
	elif jet_type == "CA15":
		event_selector.add_cut("Min_CA15Puppijet0_pt", {"Min_CA15Puppijet0_pt":450., "systematic":jet_systematic})
		event_selector.add_cut("Min_CA15Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("CA15Puppijet0_isTightVJet")

	event_selector.add_cut("Max_neleLoose", 0)
	event_selector.add_cut("Max_nmuLoose", 0)
	event_selector.add_cut("Max_ntau", 0)
	event_selector.add_cut("Max_pfmet", {"Max_pfmet":140., "systematic":jet_systematic})

	# AK8 or CA15 cuts
	if jet_type == "AK8":
		if tau21_ddt_cut != None:
			event_selector.add_cut("Max_AK8Puppijet0_tau21DDT", tau21_ddt_cut)
		else:
			event_selector.add_cut("Max_AK8Puppijet0_N2DDT", n2_ddt_cut)
	elif jet_type == "CA15":
		if tau21_ddt_cut != None:
			event_selector.add_cut("Max_CA15Puppijet0_tau21DDT", tau21_ddt_cut)
		else:
			event_selector.add_cut("Max_CA15Puppijet0_N2DDT", n2_ddt_cut)

	#if jet_type == "AK8":
	#	event_selector.add_cut("Min_AK8Puppijet0_rho", -6.0)
	#	event_selector.add_cut("Max_AK8Puppijet0_rho", -2.1)
	#elif jet_type == "CA15":
	#	event_selector.add_cut("Min_CA15Puppijet0_rho", -4.7)
	#	event_selector.add_cut("Max_CA15Puppijet0_rho", -1.0)

	return event_selector

def MakeMuCRSelector(jet_type, n2_ddt_cut=0., jet_systematic="nominal"):
	selector_name = "EventSelector_muCR"
	if jet_systematic != "nominal":
		selector_name += "_" + jet_systematic
	event_selector = BaconEventSelector(selector_name)

	if jet_type == "AK8":
		event_selector.add_cut("Min_AK8Puppijet0_pt", {"Min_AK8Puppijet0_pt":400., "systematic":jet_systematic})
		event_selector.add_cut("Min_AK8Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("AK8Puppijet0_isTightVJet")
	elif jet_type == "CA15":
		event_selector.add_cut("Min_CA15Puppijet0_pt", {"Min_CA15Puppijet0_pt":400., "systematic":jet_systematic})
		event_selector.add_cut("Min_CA15Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("CA15Puppijet0_isTightVJet")

	# General event cuts, independent of jet type
	event_selector.add_cut("Max_neleLoose", 0)
	event_selector.add_cut("Max_ntau", 0)
	event_selector.add_cut("Min_nmuLoose", 1)
	event_selector.add_cut("Max_nmuLoose", 1)
	event_selector.add_cut("Min_vmuoLoose0_pt", 55.)
	event_selector.add_cut("Max_vmuoLoose0_abseta", 2.1)
	event_selector.add_cut("Min_dphi_mu_jet", {"Min_dphi_mu_jet":2. * ROOT.TMath.Pi() / 3., "jet_type":jet_type})
	event_selector.add_cut("Min_nAK4PuppijetsMPt50dR08_0", 1)

	# AK8 or CA15 cuts
	if jet_type == "AK8":
		event_selector.add_cut("Max_AK8Puppijet0_N2DDT", n2_ddt_cut)
	elif jet_type == "CA15":
		event_selector.add_cut("Max_CA15Puppijet0_N2DDT", n2_ddt_cut)
	return event_selector

def MakePreselectionSelector(jet_type, n2_ddt_cut=0., tau21_ddt_cut=None, jet_systematic="nominal", tag=None):
	selector_name = "EventSelector_Preselection"
	if jet_systematic != "nominal":
		selector_name += "_" + jet_systematic
	if tag:
		selector_name += "_" + tag
	event_selector = BaconEventSelector(selector_name)

	if jet_type == "AK8":
		#event_selector.add_cut("Min_AK8Puppijet0_pt", {"Min_AK8Puppijet0_pt":450., "systematic":jet_systematic})
		#event_selector.add_cut("Min_AK8Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("AK8Puppijet0_isTightVJet")
	elif jet_type == "CA15":
		#event_selector.add_cut("Min_CA15Puppijet0_pt", {"Min_CA15Puppijet0_pt":450., "systematic":jet_systematic})
		#event_selector.add_cut("Min_CA15Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("CA15Puppijet0_isTightVJet")

	event_selector.add_cut("Max_neleLoose", 0)
	event_selector.add_cut("Max_nmuLoose", 0)
	event_selector.add_cut("Max_ntau", 0)
	event_selector.add_cut("Max_pfmet", {"Max_pfmet":140., "systematic":jet_systematic})

	return event_selector

def MakeN2CRSelector(jet_type, n2_ddt_cut=0., tau21_ddt_cut=None, jet_systematic="nominal", tag=None):
	selector_name = "EventSelector_N2CR"
	if jet_systematic != "nominal":
		selector_name += "_" + jet_systematic
	if tag:
		selector_name += "_" + tag
	event_selector = BaconEventSelector(selector_name)

	if jet_type == "AK8":
		event_selector.add_cut("Min_AK8Puppijet0_pt", {"Min_AK8Puppijet0_pt":450., "systematic":jet_systematic})
		event_selector.add_cut("Min_AK8Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("AK8Puppijet0_isTightVJet")
	elif jet_type == "CA15":
		event_selector.add_cut("Min_CA15Puppijet0_pt", {"Min_CA15Puppijet0_pt":450., "systematic":jet_systematic})
		event_selector.add_cut("Min_CA15Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("CA15Puppijet0_isTightVJet")

	event_selector.add_cut("Max_neleLoose", 0)
	event_selector.add_cut("Max_nmuLoose", 0)
	event_selector.add_cut("Max_ntau", 0)
	event_selector.add_cut("Max_pfmet", {"Max_pfmet":140., "systematic":jet_systematic})

	# AK8 or CA15 cuts
	if jet_type == "AK8":
		if tau21_ddt_cut != None:
			event_selector.add_cut("Min_AK8Puppijet0_tau21DDT", tau21_ddt_cut)
		else:
			event_selector.add_cut("Min_AK8Puppijet0_N2DDT", n2_ddt_cut)
	elif jet_type == "CA15":
		if tau21_ddt_cut != None:
			event_selector.add_cut("Min_CA15Puppijet0_tau21DDT", tau21_ddt_cut)
		else:
			event_selector.add_cut("Min_CA15Puppijet0_N2DDT", n2_ddt_cut)

	#if jet_type == "AK8":
	#	event_selector.add_cut("Min_AK8Puppijet0_rho", -6.0)
	#	event_selector.add_cut("Max_AK8Puppijet0_rho", -2.1)
	#elif jet_type == "CA15":
	#	event_selector.add_cut("Min_CA15Puppijet0_rho", -4.7)
	#	event_selector.add_cut("Max_CA15Puppijet0_rho", -1.0)

	return event_selector


def MakeN2SRSelector(jet_type, dbtag_cut=0., jet_systematic="nominal", tag=None):
	selector_name = "EventSelector_N2SR"
	if jet_systematic != "nominal":
		selector_name += "_" + jet_systematic
	if tag:
		selector_name += "_" + tag
	event_selector = BaconEventSelector(selector_name)

	if jet_type == "AK8":
		event_selector.add_cut("Min_AK8Puppijet0_pt", {"Min_AK8Puppijet0_pt":450., "systematic":jet_systematic})
		event_selector.add_cut("Min_AK8Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("AK8Puppijet0_isTightVJet")
	elif jet_type == "CA15":
		event_selector.add_cut("Min_CA15Puppijet0_pt", {"Min_CA15Puppijet0_pt":450., "systematic":jet_systematic})
		event_selector.add_cut("Min_CA15Puppijet0_msd_puppi", 40.)
		event_selector.add_cut("CA15Puppijet0_isTightVJet")

	event_selector.add_cut("Max_neleLoose", 0)
	event_selector.add_cut("Max_nmuLoose", 0)
	event_selector.add_cut("Max_ntau", 0)
	event_selector.add_cut("Max_pfmet", {"Max_pfmet":140., "systematic":jet_systematic})

	# AK8 or CA15 cuts
	if jet_type == "AK8":
		event_selector.add_cut("Min_AK8Puppijet0_doublecsv", dbtag_cut)
	elif jet_type == "CA15":
		event_selector.add_cut("Min_CA15Puppijet0_doublesub", dbtag_cut)

	#if jet_type == "AK8":
	#	event_selector.add_cut("Min_AK8Puppijet0_rho", -6.0)
	#	event_selector.add_cut("Max_AK8Puppijet0_rho", -2.1)
	#elif jet_type == "CA15":
	#	event_selector.add_cut("Min_CA15Puppijet0_rho", -4.7)
	#	event_selector.add_cut("Max_CA15Puppijet0_rho", -1.0)

	return event_selector
