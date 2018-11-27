import os
import sys
import pickle
sys.path.append(".")
from scale_factors import scale_factors
from DAZSLE.PhiBBPlusJet.analysis_configuration import background_names as backgrounds 
from DAZSLE.PhiBBPlusJet.analysis_configuration import simulated_signal_names as signals

from ROOT import *

# When run as a script, this file computes all the normalization systematics from the input histograms, and saves them to a pickle file.
# When loaded as a module, this file provides a dict of normalization systematics, and a list of shape systematics
# systematics_norm = {syst_name : {process : (pass unc., fail unc.)}}

def PrintNormalizationSystematics(jet_type, wp_name):
	systematics_norm = GetNormalizationSystematics(jet_type, wp_name)
	for syst_name in sorted(systematics_norm.keys()):
		print syst_name
		for process in sorted(systematics_norm[syst_name].keys()):
			print "\t" + process
			print "\t\t",
			print systematics_norm[syst_name][process]

def GetSystNormFile(jet_type, wp_name):
	return os.path.expandvars("$CMSSW_BASE/data/DAZSLE/PhiBBPlusJet/data/systematics/syst_{}_{}_norm.pkl".format(jet_type, wp_name))

def GetNormalizationSystematics(jet_type, wp_name):
	with open(GetSystNormFile(jet_type, wp_name), 'r') as f:
		systematics_norm = pickle.load(f)
	return systematics_norm
	
def GetShapeSystematics():
	# Shape systematics: dict store the 1/sigma values for the datacard. User is responsible for providing the shape-varied histograms in the appropriate format.
	systematics_shape = {}
	systematics_shape["smear"] = {}
	systematics_shape["scale"] = {}
	systematics_shape["scalept"] = {}
	for process in backgrounds + signals:
		if process == "qcd":
			continue
		systematics_shape["smear"][process] = 0.5
		systematics_shape["scale"][process] = 0.1
		systematics_shape["scalept"][process] = 0.4


def ComputeNormalizationSystematics(jet_type, wp_name):
	histogram_file = TFile("/uscms/home/dryu/PhiBB2017/data/histograms/optimization/histograms_SR_{}_{}.root".format(jet_type, wp_name), "READ")
	systematics_norm = {}

	systematics_norm["lumi"] = {}
	systematics_norm["muveto"] = {}
	systematics_norm["eleveto"] = {}
	for process in backgrounds:
		if not process in ["tqq", "qcd"]:
			systematics_norm["lumi"][process] = 1.025
		if not process in ["qcd"]:
			systematics_norm["muveto"][process] = 1.005
			systematics_norm["eleveto"][process] = 1.005

	# Tagging systematics
	systematics_norm["veff"] = {}
	for process in ["wqq", "zqq", "hqq125", "tthqq125", "vbfhqq125", "whqq125", "zhqq125"] + signals:
		systematics_norm["veff"][process] = {}
		for box in ["passn2_passdbtag", "passn2_faildbtag", "failn2_passdbtag", "failn2_faildbtag"]:
			if "passn2" in box:
				sf_unc = 1.0 + scale_factors[jet_type]["V_SF_ERR"] / scale_factors[jet_type]["V_SF"]
			else:
				# Calculation details:
				#p = s P
				#dp = ds P
				#p + f = P + F
				#f = P+F-sP = P(1-s)+F
				#df = -P ds = -dp
				#Combine wants 1+dp/p and 1+df/f
				#1+dp/p = 1 + ds/s
				#1+df/f = 1-dp/f = 1-(dp/p)(p/f) = 1-(ds/s)(p/f)
				pass_hname =  "{}_{}_pt_vs_msd".format(process.replace("failn2", "passn2"), box)
				pass_hist = histogram_file.Get(pass_hname)
				if not pass_hist:
					print "[systematics::ComputeNormalizationSystematics] ERROR : Couldn't find histogram {} in file {}".format(pass_hname, histogram_file.GetPath())
					sys.exit(1)
				pass_normalization = pass_hist.Integral() * scale_factors[jet_type]["V_SF"]

				fail_hname =  "{}_{}_pt_vs_msd".format(process.replace("failn2", "failn2"), box)
				fail_hist = histogram_file.Get(fail_hname)
				if not fail_hist:
					print "[systematics::ComputeNormalizationSystematics] ERROR : Couldn't find histogram {} in file {}".format(fail_hname, histogram_file.GetPath())
					sys.exit(1)
				fail_normalization = fail_hist.Integral() + pass_hist.Integral() * (1 - scale_factors[jet_type]["V_SF"])

				if fail_normalization > 0:
					sf_unc = 1.0 - (scale_factors[jet_type]["V_SF_ERR"] / scale_factors[jet_type]["V_SF"]) * (pass_normalization / fail_normalization)
				else:
					sf_unc = 1.0
			systematics_norm["veff"][process][box] = sf_unc

	systematics_norm["bbeff"] = {}
	for process in ["wqq", "zqq", "hqq125", "tthqq125", "vbfhqq125", "whqq125", "zhqq125"] + signals:
		systematics_norm["bbeff"][process] = {}
		for box in ["passn2_passdbtag", "passn2_faildbtag", "failn2_passdbtag", "failn2_faildbtag"]:
			if "passdbtag" in box:
				sf_unc = 1.0 + scale_factors[jet_type]["BB_SF_ERR"] / scale_factors[jet_type]["BB_SF"]
			else:
				pass_hname =  "{}_{}_pt_vs_msd".format(process.replace("faildbtag", "passdbtag"), box)
				pass_hist = histogram_file.Get(pass_hname)
				if not pass_hist:
					print "[systematics::ComputeNormalizationSystematics] ERROR : Couldn't find histogram {} in file {}".format(pass_hname, histogram_file.GetPath())
					sys.exit(1)
				pass_normalization = pass_hist.Integral() * scale_factors[jet_type]["BB_SF"]

				fail_hname =  "{}_{}_pt_vs_msd".format(process.replace("faildbtag", "faildbtag"), box)
				fail_hist = histogram_file.Get(fail_hname)
				if not fail_hist:
					print "[systematics::ComputeNormalizationSystematics] ERROR : Couldn't find histogram {} in file {}".format(fail_hname, histogram_file.GetPath())
					sys.exit(1)
				fail_normalization = fail_hist.Integral() + pass_hist.Integral() * (1 - scale_factors[jet_type]["BB_SF"])

				if fail_normalization > 0:
					sf_unc = 1.0 - (scale_factors[jet_type]["BB_SF_ERR"] / scale_factors[jet_type]["BB_SF"]) * (pass_normalization / fail_normalization)
				else:
					sf_unc = 1.0
			systematics_norm["bbeff"][process][box] = sf_unc

	# W and Z systetmatics
	systematics_norm["znormQ"] = {"zqq":1.1, "wqq":1.1}

	systematics_norm["znormEW"] = {"wqq": [1.15, 1.15, 1.25, 1.35, 1.35, 1.35], "zqq":[1.15, 1.15, 1.25, 1.35, 1.35, 1.35]}

	systematics_norm["wznormEW"] = {"wqq":[1.05, 1.05, 1.05, 1.15, 1.15, 1.15], "zqq":[1.05, 1.05, 1.05, 1.15, 1.15, 1.15]}

	# H systematics
	systematics_norm["hpt"] = {}
	for process in ["hqq125", "tthqq125", "vbfhqq125", "whqq125", "zhqq125"]:
		systematics_norm["hpt"][process] = 1.3

	# Systematics determined from variation in the event processing
	for systematic in ["JER", "JES", "PU", "Trigger"]:
		systematics_norm[systematic] = {}
		for process in backgrounds + signals:
			if process == "qcd":
				continue
			systematics_norm[systematic][process] = {}
			for box in ["passn2_passdbtag", "passn2_faildbtag", "failn2_passdbtag", "failn2_faildbtag"]:
				hname = "{}_{}_pt_vs_msd".format(process, box)
				if not histogram_file.Get(hname):
					print "Couldn't get histogram {} from file {}. Setting systematic uncertainty to zero.".format(hname, histogram_file.GetPath())
					systematics_norm[systematic][process][box] = 1.0
					continue
				nominal = histogram_file.Get(hname).Integral()

				hname = "{}_{}_pt_vs_msd_{}Up".format(process, box, systematic)
				if not histogram_file.Get(hname):
					print "Couldn't get histogram {} from file {}.".format(hname, histogram_file.GetPath())
					sys.exit(1)
				up = histogram_file.Get(hname).Integral()

				hname = "{}_{}_pt_vs_msd_{}Down".format(process, box, systematic)
				if not histogram_file.Get(hname):
					print "Couldn't get histogram {} from file {}.".format(hname, histogram_file.GetPath())
					sys.exit(1)
				down = histogram_file.Get(hname).Integral()

				if nominal > 0:
					systematics_norm[systematic][process][box] = 1.0 + max(abs(up - nominal) / nominal, abs(down - nominal) / nominal)
				else:
					systematics_norm[systematic][process][box] = 1.0

	with open(GetSystNormFile(jet_type, wp_name), 'w') as f:
		pickle.dump(systematics_norm, f)

	PrintNormalizationSystematics(jet_type, wp_name)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser("Systematics helper")
	parser.add_argument("--compute", action="store_true", help="Compute systematics")
	parser.add_argument("--jet_type", type=str, help="Jet type")
	parser.add_argument("--wp", type=str, default="dbtag0.9_n2wpdbpass0.15_n2wpdbfail0.05", help="Working point name")
	parser.add_argument("--printt", action="store_true", help="Print systematics to terminal")
	parser.add_argument("--print_table", action="store_true", help="Print systematics as LaTeX table")
	args = parser.parse_args()

	if args.compute:
		ComputeNormalizationSystematics(args.jet_type, args.wp)
	if args.printt:
		PrintNormalizationSystematics(args.jet_type, args.wp)
