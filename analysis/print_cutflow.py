import os
import sys
from array import array
import ROOT
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libDAZSLEPhiBBPlusJet.so"))
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so"))
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
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

# Convert a cutflow histogram to a latex table row
def cutflow_histogram_to_row(name, hist):
	if name in formatted_samples:
		row_string = "\t\t" + formatted_samples[name].replace("_","\_")
	else:
		row_string = "\t\t" + name.replace("_","\_")
	for bin in xrange(1, hist.GetNbinsX() + 1):
		if hist.GetXaxis().GetBinLabel(bin) != "":
			row_string += "\t&\t" + str(int(hist.GetBinContent(bin)))
	row_string += "\t\\\\"
	return row_string

def get_cutflow_headers(hist, region):
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
	header_string += "\\\\"
	print header_string
	return header_string, columns

if __name__ == "__main__":
	jet_types = ["AK8", "CA15"]
	regions = ["SR", "muCR"]
	supersamples = {
		"SR":["data_obs", "Sbb100", "Sbb200", "Sbb300", "qcd", "tqq", "zqq", "wqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125",],
		"muCR":["data_singlemu", "tqq", "qcd", "wlnu", "zll", "Sbb125"]
	}
	for jet_type in jet_types:
		for region in regions:
			# Make latex table
			table_file = open(os.path.expandvars("$HOME/DAZSLE/data/EventSelection/cutflows/{}_{}.tex".format(region, jet_type)), 'w')
			cutflow_histograms = {}
			row_strings = {}
			header_string = None
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
					f.Close()
				cutflow_histograms[supersample].SetDirectory(0)
				row_strings[supersample] = cutflow_histogram_to_row(supersample, cutflow_histograms[supersample])
				if not header_string:
					header_string, columns = get_cutflow_headers(cutflow_histograms[supersample], region)
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

			# Make pretty plot
			#cutflow_canvas = 
