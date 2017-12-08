import os
import sys
from array import array
import ROOT
import math
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libDAZSLEPhiBBPlusJet.so"))
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so"))
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
import DAZSLE.PhiBBPlusJet.style as style
seaborn = Root.SeabornInterface()
seaborn.Initialize()

formatted_samples = {}
formatted_labels = {
	"SR":{
		"Inclusive":"Inclusive",
		"Min_AK8Puppijet0_pt":"$p_{\\mathrm{T}}>\SI{450}{\\giga\\electronvolt}$",
		"Min_CA15Puppijet0_pt":"$p_{\\mathrm{T}}>\SI{450}{\\giga\\electronvolt}$",
		"Min_AK8Puppijet0_msd_puppi":"$m_{SD}>\SI{40}{\\giga\\electronvolt}$",
		"Min_CA15Puppijet0_msd_puppi":"$m_{SD}>\SI{40}{\\giga\\electronvolt}$",
		"AK8Puppijet0_isTightVJet":"Tight VJet",
		"CA15Puppijet0_isTightVJet":"Tight VJet",
		"Max_neleLoose":"$n_{e}==0$",
		"Max_nmuLoose":"$n_{\\mu}==0$",
		"Max_ntau":"$n_{\\tau}==0$",
		"Max_pfmet":"PF MET$<\\SI{140}{\\giga\\electronvolt}$",
		"Max_AK8Puppijet0_N2DDT":"$N_{2}^{\\mathrm{DDT}}<0$",
		"Max_CA15Puppijet0_N2DDT":"$N_{2}^{\\mathrm{DDT}}<0$",
	},
	"muCR":{
		"Inclusive":"Inclusive",
		"Min_AK8Puppijet0_pt":"$p_{\\mathrm{T}}>\\SI{400}{\\giga\\electronvolt}$",
		"Min_CA15Puppijet0_pt":"$p_{\\mathrm{T}}>\\SI{400}{\\giga\\electronvolt}$",
		"Min_AK8Puppijet0_msd_puppi":"$m_{\\mathrm{SD}}>\\SI{40}{\\giga\\electronvolt}$",
		"Min_CA15Puppijet0_msd_puppi":"$m_{\\mathrm{SD}}>\\SI{40}{\\giga\\electronvolt}$",
		"AK8Puppijet0_isTightVJet":"Tight VJet",
		"CA15Puppijet0_isTightVJet":"Tight VJet",
		"Max_neleLoose":"$n_{e}==0$",
		"Max_ntau":"$n_{\\tau}==0$",
		"Min_nmuLoose":"$n_{\\mu}\\geq1$",
		"Max_nmuLoose":"$n_{\\mu}\\leq1$",
		"Min_vmuoLoose0_pt":"$p_{\\mathrm{T}}(\\mu)>\\SI{55}{\\giga\\electronvolt}$",
		"Max_vmuoLoose0_abseta":"$|\\eta(\\mu)|<2.1$",
		"Min_dphi_mu_jet":"$\\Delta\\phi(\\mu,j)<\\frac{2\\pi}{3}$",
		"Min_nAK4PuppijetsMPt50dR08_0":"Min\_nAK4PuppijetsMPt50dR08\_0",
		"Max_AK8Puppijet0_N2DDT":"$N_{2}^{\\mathrm{DDT}}<0$",
		"Max_CA15Puppijet0_N2DDT":"$N_{2}^{\\mathrm{DDT}}<0$",
	}
}

def get_rhocut_hist(hist, rho_min, rho_max):
	rhocut_hist = hist.Clone()
	rhocut_hist.SetName(hist.GetName() + "_rhocut")
	for i in range(1,rhocut_hist.GetNbinsX()+1):
		for j in range(1,rhocut_hist.GetNbinsY()+1):
			massVal = rhocut_hist.GetXaxis().GetBinCenter(i)
			ptVal = rhocut_hist.GetYaxis().GetBinLowEdge(j) + rhocut_hist.GetYaxis().GetBinWidth(j) * 0.3
			rhoVal = math.log(massVal * massVal / ptVal / ptVal)
			if rhoVal < rho_min or rhoVal > rho_max:
				rhocut_hist.SetBinContent(i, j, 0.)
	return rhocut_hist


# Convert a cutflow histogram to a latex table row
def cutflow_histogram_to_row(name, cutflow_hist, pass_hist, fail_hist):
	if name in formatted_samples:
		row_string = "\t\t" + formatted_samples[name].replace("_","\_")
	else:
		row_string = "\t\t" + name.replace("_","\_")
	last_cutflow_bin = -1
	for bin in xrange(1, cutflow_hist.GetNbinsX() + 1):
		if cutflow_hist.GetXaxis().GetBinLabel(bin) != "":
			row_string += "\t&\t" + str(int(cutflow_hist.GetBinContent(bin)))
			last_cutflow_bin = cutflow_hist.GetBinContent(bin)

	# Rho and dbtag cuts
	if last_cutflow_bin != pass_hist.Integral(0, pass_hist.GetNbinsX()+1, 0, pass_hist.GetNbinsY()+1) + fail_hist.Integral(0, fail_hist.GetNbinsX()+1, 0, fail_hist.GetNbinsY()+1):
		print "WARNING : {}: last cutflow bin = {}, but pass+fail integral = {}".format(name, last_cutflow_bin, pass_hist.Integral(0, pass_hist.GetNbinsX()+1, 0, pass_hist.GetNbinsY()+1) + fail_hist.Integral(0, fail_hist.GetNbinsX()+1, 0, fail_hist.GetNbinsY()+1))
	pass_hist_rhocut = get_rhocut_hist(pass_hist, params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	fail_hist_rhocut = get_rhocut_hist(fail_hist, params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	row_string += "\t&\t" + str(int(pass_hist_rhocut.Integral() + fail_hist_rhocut.Integral()))
	row_string += "\t&\t" + str(int(pass_hist_rhocut.Integral()))

	row_string += "\t\\\\"
	return row_string

def cutflow_histogram_to_row_percent(name, cutflow_hist, pass_hist, fail_hist):
	if name in formatted_samples:
		row_string = "\t\t" + formatted_samples[name].replace("_","\_")
	else:
		row_string = "\t\t" + name.replace("_","\_")
	inclusive_counts = cutflow_hist.GetBinContent(1)
	if inclusive_counts == 0:
		print "WARNING : inclusive_counts = 0 for {}".format(name)
		return "NULL"
	last_cutflow_bin = -1
	for bin in xrange(2, cutflow_hist.GetNbinsX() + 1):
		if cutflow_hist.GetXaxis().GetBinLabel(bin) != "":
			eff = 1. * cutflow_hist.GetBinContent(bin) / inclusive_counts
			deff = math.sqrt(eff * (1. - eff) / inclusive_counts)
			row_string += "\t&\t$" + str(round(100. * eff, 1)) + " \\pm " + str(round(100. * deff, 1)) + "$"
			last_cutflow_bin = cutflow_hist.GetBinContent(bin)

	# Rho and dbtag cuts
	if last_cutflow_bin != pass_hist.Integral(0, pass_hist.GetNbinsX()+1, 0, pass_hist.GetNbinsY()+1) + fail_hist.Integral(0, fail_hist.GetNbinsX()+1, 0, fail_hist.GetNbinsY()+1):
		print "WARNING : {}: last cutflow bin = {}, but pass+fail integral = {}".format(name, last_cutflow_bin, pass_hist.Integral(0, pass_hist.GetNbinsX()+1, 0, pass_hist.GetNbinsY()+1) + fail_hist.Integral(0, fail_hist.GetNbinsX()+1, 0, fail_hist.GetNbinsY()+1))
	pass_hist_rhocut = get_rhocut_hist(pass_hist, params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	fail_hist_rhocut = get_rhocut_hist(fail_hist, params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	eff = float(pass_hist_rhocut.Integral() + fail_hist_rhocut.Integral()) / inclusive_counts
	deff = math.sqrt(eff * (1. - eff) / inclusive_counts)
	row_string += "\t&\t$" + str(round(100. * eff, 1)) + " \\pm " + str(round(100. * deff, 1)) + "$"
	eff = float(pass_hist_rhocut.Integral()) / inclusive_counts
	deff = math.sqrt(eff * (1. - eff) / inclusive_counts)
	row_string += "\t&\t$" + str(round(100. * eff, 1)) + " \\pm " + str(round(100. * deff, 1)) + "$"

	row_string += "\t\\\\"
	return row_string

def get_cutflow_headers(hist, region, jet_type):
	header_string = "\t\t Sample"
	columns = 1
	for bin in xrange(1, hist.GetNbinsX() + 1):
		if hist.GetXaxis().GetBinLabel(bin) != "":
			raw_label = hist.GetXaxis().GetBinLabel(bin)
			if raw_label in formatted_labels[region]:
				header_string += "\t&\t" + formatted_labels[region][raw_label]
			else:
				header_string += "\t&\t" + raw_label
			columns += 1
	header_string += "\t&\t${}<\\rho<{}$".format(params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	columns += 1
	header_string += "\t&\t2$times$b tag"
	columns += 1
	header_string += "\\\\"
	#print header_string
	return header_string, columns

def get_cutflow_headers_percent(hist, region, jet_type):
	header_string = "\t\t Sample"
	columns = 1
	for bin in xrange(2, hist.GetNbinsX() + 1):
		if hist.GetXaxis().GetBinLabel(bin) != "":
			raw_label = hist.GetXaxis().GetBinLabel(bin)
			if raw_label in formatted_labels[region]:
				header_string += "\t&\t" + formatted_labels[region][raw_label]
			else:
				header_string += "\t&\t" + raw_label
			columns += 1
	header_string += "\t&\t${}<\\rho<{}$".format(params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	columns += 1
	header_string += "\t&\t2$times$b tag"
	columns += 1
	header_string += "\\\\"
	#print header_string
	return header_string, columns

if __name__ == "__main__":
	jet_types = ["AK8", "CA15"]
	regions = ["SR", "muCR"]
	supersamples = {
		"SR":["data_obs", "Sbb50", "Sbb100", "Sbb125", "Sbb200", "Sbb300", "Sbb350", "qcd", "tqq", "zqq", "wqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125",],
		"muCR":["data_singlemu", "tqq", "qcd", "wlnu", "zll", "Sbb125"]
	}
	for jet_type in jet_types:
		for region in regions:
			# Make latex table
			table_file_path = os.path.expandvars("$HOME/DAZSLE/data/EventSelection/cutflows/{}_{}.tex".format(region, jet_type))
			print "Opening table file {}".format(table_file_path)
			table_file = open(table_file_path, 'w')
			table_file_percent_path = os.path.expandvars("$HOME/DAZSLE/data/EventSelection/cutflows/{}_{}_percent.tex".format(region, jet_type))
			print "Opening percent table file {}".format(table_file_percent_path)
			table_file_percent = open(table_file_percent_path, 'w')
			cutflow_histograms = {}
			pass_histograms = {}
			fail_histograms = {}
			row_strings = {}
			row_strings_percent = {}
			header_string = None
			header_string_percent = None
			columns = 0
			for supersample in supersamples[region]:
				for sample in config.samples[supersample]:
					fname = os.path.expandvars("$HOME/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root".format(sample, jet_type))
					f = TFile(fname, "READ")
					hname = "CutFlowCounter_EventSelector_{}".format(region)
					if not f.Get(hname):
						print "ERROR : Unable to get histogram {} from file {}".format(hname, f.GetPath())
						sys.exit(1)
					if not supersample in cutflow_histograms:
						cutflow_histograms[supersample] = f.Get(hname).Clone()
						cutflow_histograms[supersample].SetName("h_CutFlowCounter_{}_{}_{}".format(jet_type, region, supersample))
						cutflow_histograms[supersample].SetDirectory(0)
					else:
						cutflow_histograms[supersample].Add(f.Get(hname))

					hname = "h_{}_{}_pass_unweighted".format(region, jet_type)
					if not f.Get(hname):
						print "ERROR : Unable to get histogram {} from file {}".format(hname, f.GetPath())
						sys.exit(1)
					if not supersample in pass_histograms:
						pass_histograms[supersample] = f.Get(hname).Clone()
						pass_histograms[supersample].SetName(hname + "_" + supersample)
						pass_histograms[supersample].SetDirectory(0)
					else:
						pass_histograms[supersample].Add(f.Get(hname))

					hname = "h_{}_{}_fail_unweighted".format(region, jet_type)
					if not f.Get(hname):
						print "ERROR : Unable to get histogram {} from file {}".format(hname, f.GetPath())
						sys.exit(1)
					if not supersample in fail_histograms:
						fail_histograms[supersample] = f.Get(hname).Clone()
						fail_histograms[supersample].SetName(hname + "_" + supersample)
						fail_histograms[supersample].SetDirectory(0)
					else:
						fail_histograms[supersample].Add(f.Get(hname))
					f.Close()
				cutflow_histograms[supersample].SetDirectory(0)
				pass_histograms[supersample].SetDirectory(0)
				fail_histograms[supersample].SetDirectory(0)

				print "Converting histogram to string for {} {} {} {}".format(jet_type, region, supersample, sample)
				row_strings[supersample] = cutflow_histogram_to_row(supersample, cutflow_histograms[supersample], pass_histograms[supersample], fail_histograms[supersample])
				row_strings_percent[supersample] = cutflow_histogram_to_row_percent(supersample, cutflow_histograms[supersample], pass_histograms[supersample], fail_histograms[supersample])
				if not header_string:
					header_string, columns_percent = get_cutflow_headers(cutflow_histograms[supersample], region, jet_type)
				if not header_string_percent:
					header_string_percent, columns_percent = get_cutflow_headers_percent(cutflow_histograms[supersample], region, jet_type)
			table_file.write("\\begin{table}\n")
			table_file.write("\t\\begin{tabular}{|c|" + "c|".join(["" for x in xrange(columns)]) + "}\n")
			table_file.write("\t\t\\hline\n")
			table_file.write(header_string + "\n")
			table_file.write("\t\t\\hline\n")
			for supersample in supersamples[region]:
				table_file.write(row_strings[supersample] + "\n")
				table_file.write("\t\t\\hline\n")
			table_file.write("\t\\end{tabular}\n")
			table_file.write("\\end{table}\n")
			table_file.close()

			table_file_percent.write("\\begin{table}\n")
			table_file_percent.write("\t\\begin{tabular}{|c|" + "c|".join(["" for x in xrange(columns_percent)]) + "}\n")
			table_file_percent.write("\t\t\\hline\n")
			table_file_percent.write(header_string_percent + "\n")
			table_file_percent.write("\t\t\\hline\n")
			for supersample in supersamples[region]:
				table_file_percent.write(row_strings_percent[supersample] + "\n")
				table_file_percent.write("\t\t\\hline\n")
			table_file_percent.write("\t\\end{tabular}\n")
			table_file_percent.write("\\end{table}\n")
			table_file_percent.close()

			# Make pretty plot
			#cutflow_canvas = 
