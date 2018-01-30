# Background cutflow: THStack
# Plus a few signal cutflows

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

def make_background_stack(backgrounds, cutflow_hists):
	stack = THStack("background_cutflows", "background_cutflows")
	backgrounds.sort(key=lambda x: cutflow_hists[x].GetBinContent(1))
	for background in backgrounds:
		cutflow_hists[background].SetFillColor(style.background_colors[background])
		cutflow_hists[background].SetLineColor(1)
		cutflow_hists[background].SetLineStyle(1)
		cutflow_hists[background].SetLineWidth(2)
		stack.Add(cutflow_hists[background])
	return stack

def make_background_eff_stack(backgrounds, cutflow_hists):
	total_inclusive_counts = 0
	for background in backgrounds:
		total_inclusive_counts += cutflow_hists[background].GetBinContent(1)

	cutflow_eff_hists = {}
	for background in backgrounds: 
		cutflow_eff_hists[background] = cutflow_hists[background].Clone()
		cutflow_eff_hists[background].Scale(1. / total_inclusive_counts)

	stack = THStack("background_eff_stack", "background_eff_stack")
	for background in backgrounds:
		stack.Add(cutflow_eff_hists[background])

	# Total background efficiency histogram with statistical uncertainties
	total_background_hist = None
	for background in backgrounds:
		if not total_background_hist:
			total_background_hist = cutflow_hists[background].Clone()
		else:
			total_background_hist.Add(cutflow_hists[background])
	inclusive_counts_hist = total_background_hist.Clone()
	inclusive_counts_hist.Reset()
	for bin in xrange(1, inclusive_counts_hist.GetNbinsX() + 1):
		inclusive_counts_hist.SetBinContent(bin, total_background_hist.GetBinContent(1))
		inclusive_counts_hist.SetBinError(bin, total_background_hist.GetBinError(1))
	total_background_eff_hist = total_background_hist.Clone()
	total_background_eff_hist.Divide(total_background_hist, inclusive_counts_hist, 1., 1., "B")

	return stack, total_background_eff_hist

def convert_cutflow_to_efficiency(cutflow_hist):
	inclusive_hist = cutflow_hist.Clone()
	inclusive_hist.Reset()
	for bin in xrange(1, inclusive_hist.GetNbinsX() + 1):
		inclusive_hist.SetBinContent(bin, cutflow_hist.GetBinContent(1))
		inclusive_hist.SetBinError(bin, cutflow_hist.GetBinError(1))
	cutflow_eff_hist = cutflow_hist.Clone()
	cutflow_eff_hist.Reset()
	cutflow_eff_hist.Divide(cutflow_hist, inclusive_hist, 1., 1., "B")
	return cutflow_eff_hist

# Chop out unused bins, normalize to inclusive, add rho cut and 2xbtag cut
def process_cutflow_hist(cutflow_hist, pass_hist, fail_hist):
	nbins = 0
	cuts = []
	event_counts = {}
	devent_counts = {}
	effs = {}
	deffs = {}
	inclusive_counts = 0
	for bin in xrange(1, cutflow_hist.GetNbinsX() + 1):
		this_cut = cutflow_hist.GetXaxis().GetBinLabel(bin)
		if this_cut == "Inclusive":
			inclusive_counts = cutflow_hist.GetBinContent(bin)
		if not this_cut in cuts:
			cuts.append(this_cut)
			event_counts[this_cut] = cutflow_hist.GetBinContent(bin)
			devent_counts[this_cut] = cutflow_hist.GetBinError(bin)
			effs[this_cut] = 1. * cutflow_hist.GetBinContent(bin) / inclusive_counts
			deffs[this_cut] = math.sqrt(effs[this_cut] * (1. - effs[this_cut]) / inclusive_counts)

	# Rho cut
	pass_hist_rhocut = get_rhocut_hist(pass_hist, params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	pass_hist_rhocut.RebinX(pass_hist_rhocut.GetNbinsX()).RebinY(pass_hist_rhocut.GetNbinsY())
	fail_hist_rhocut = get_rhocut_hist(fail_hist, params[jet_type]["RHO"][0], params[jet_type]["RHO"][1])
	fail_hist_rhocut.RebinX(fail_hist_rhocut.GetNbinsX()).RebinY(fail_hist_rhocut.GetNbinsY())
	total_hist_rhocut = pass_hist_rhocut.Clone()
	total_hist_rhocut.Add(fail_hist_rhocut)
	event_counts["#rho window"] = total_hist_rhocut.GetBinContent(1,1)
	devent_counts["#rho window"] = total_hist_rhocut.GetBinError(1,1)
	cuts.append("#rho window")

	# Double b-tag
	event_counts["Double b tag"] = pass_hist_rhocut.GetBinContent(1,1)
	devent_counts["Double b tag"] = pass_hist_rhocut.GetBinError(1,1)
	cuts.append("Double b tag")

	# Efficiencies
	for cut in cuts:
		effs[cut] = 1. * event_counts[cut] / inclusive_counts
		deffs[cut] = math.sqrt(effs[cut] * (1. - effs[cut]) / inclusive_counts)

	# Make new histograms
	new_cutflow_hist = TH1D(cutflow_hist.GetName() + "_processed", cutflow_hist.GetName() + "_processed", len(cuts), 0.5, len(cuts)+0.5)
	new_cutflow_eff_hist = TH1D(cutflow_hist.GetName() + "_eff_processed", cutflow_hist.GetName() + "_eff_processed", len(cuts), 0.5, len(cuts)+0.5)
	for binm1, cut in enumerate(cuts):
		bin = binm1 + 1
		new_cutflow_hist.GetXaxis().SetBinLabel(bin, cut)
		new_cutflow_hist.SetBinContent(bin, event_counts[cut])
		new_cutflow_hist.SetBinError(bin, devent_counts[cut])

		new_cutflow_eff_hist.GetXaxis().SetBinLabel(bin, cut)
		new_cutflow_eff_hist.SetBinContent(bin, effs[cut])
		new_cutflow_eff_hist.SetBinError(bin, deffs[cut])

	return new_cutflow_hist, new_cutflow_eff_hist



if __name__ == "__main__":
	luminosity = 35900 # in pb^-1
	from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections

	jet_types = ["AK8", "CA15"]
	regions = ["SR"]#, "muCR"]
	supersamples = {
		"SR":["Sbb50", "Sbb100", "Sbb125", "Sbb200", "Sbb300", "Sbb350", "qcd", "tqq", "zqq", "wqq", "vvqq", "stqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125",],
		#"muCR":["tqq", "qcd", "wlnu", "zll", "Sbb125"]
	}
	for jet_type in jet_types:
		for region in regions:
			cutflow_histograms_raw = {}
			cutflow_histograms = {}
			cutflow_eff_histograms = {}
			pass_histograms = {}
			fail_histograms = {}
			for supersample in supersamples[region]:
				for sample in config.samples[supersample]:
					fname = os.path.expandvars("$HOME/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root".format(sample, jet_type))
					f = TFile(fname, "READ")
					n_input_events = f.Get("h_input_nevents").Integral()

					hname = "CutFlowCounter_EventSelector_{}".format(region)
					if not f.Get(hname):
						print "ERROR : Unable to get histogram {} from file {}".format(hname, f.GetPath())
						sys.exit(1)
					this_cutflow_histogram = f.Get(hname)
					this_cutflow_histogram.Scale(luminosity * cross_sections[sample] / n_input_events)
					if not supersample in cutflow_histograms_raw:
						cutflow_histograms_raw[supersample] = this_cutflow_histogram
						cutflow_histograms_raw[supersample].SetName("h_CutFlowCounter_{}_{}_{}".format(jet_type, region, supersample))
						cutflow_histograms_raw[supersample].SetDirectory(0)
					else:
						cutflow_histograms_raw[supersample].Add(this_cutflow_histogram)

					hname = "h_{}_{}_pass_unweighted".format(region, jet_type)
					if not f.Get(hname):
						print "ERROR : Unable to get histogram {} from file {}".format(hname, f.GetPath())
						sys.exit(1)
					this_pass_histogram = f.Get(hname)
					this_pass_histogram.Scale(luminosity * cross_sections[sample] / n_input_events)
					if not supersample in pass_histograms:
						pass_histograms[supersample] = this_pass_histogram.Clone()
						pass_histograms[supersample].SetName(hname + "_" + supersample)
						pass_histograms[supersample].SetDirectory(0)
					else:
						pass_histograms[supersample].Add(this_pass_histogram)

					hname = "h_{}_{}_fail_unweighted".format(region, jet_type)
					if not f.Get(hname):
						print "ERROR : Unable to get histogram {} from file {}".format(hname, f.GetPath())
						sys.exit(1)
					this_fail_histogram = f.Get(hname)
					this_fail_histogram.Scale(luminosity * cross_sections[sample] / n_input_events)
					if not supersample in fail_histograms:
						fail_histograms[supersample] = this_fail_histogram.Clone()
						fail_histograms[supersample].SetName(hname + "_" + supersample)
						fail_histograms[supersample].SetDirectory(0)
					else:
						fail_histograms[supersample].Add(this_fail_histogram)

					f.Close()
				cutflow_histograms_raw[supersample].SetDirectory(0)
				pass_histograms[supersample].SetDirectory(0)
				fail_histograms[supersample].SetDirectory(0)

				cutflow_histograms[supersample], cutflow_eff_histograms[supersample] = process_cutflow_hist(cutflow_histograms_raw[supersample], pass_histograms[supersample], fail_histograms[supersample])

			# Add up Higgs histograms
			higgs_subprocesses = ["hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]
			for higgs_subprocess in higgs_subprocesses:
				if not "hbb" in cutflow_histograms:
					cutflow_histograms_raw["hbb"] = cutflow_histograms_raw[higgs_subprocess].Clone()
					pass_histograms["hbb"] = pass_histograms[higgs_subprocess].Clone()
					fail_histograms["hbb"] = fail_histograms[higgs_subprocess].Clone()
				else:
					cutflow_histograms_raw["hbb"].Add(cutflow_histograms_raw[higgs_subprocess])
					pass_histograms["hbb"].Add(pass_histograms[higgs_subprocess])
					fail_histograms["hbb"].Add(fail_histograms[higgs_subprocess])
			cutflow_histograms["hbb"], cutflow_eff_histograms["hbb"] = process_cutflow_hist(cutflow_histograms_raw["hbb"], pass_histograms[supersample], fail_histograms[supersample])

			# Sort and style background histograms (do before making the stack, so the stacked stuff is in the right order)
			backgrounds = ["qcd","tqq","wqq","zqq","hbb","stqq","vvqq"]
			backgrounds.sort(key=lambda x: cutflow_histograms[x].GetBinContent(1))
			for background in backgrounds:
				cutflow_histograms[background].SetFillColor(style.background_colors[background])
				cutflow_histograms[background].SetLineColor(1)
				cutflow_histograms[background].SetLineStyle(1)
				cutflow_histograms[background].SetLineWidth(2)

			# Make background efficiency stack and total efficiency histogram
			background_eff_stack, total_background_eff_hist = make_background_eff_stack(backgrounds, cutflow_histograms)

			# Signal efficiency histograms
			if jet_type == "AK8":
				signals = ["Sbb50", "Sbb100", "Sbb125"]
			else:
				signals = ["Sbb125", "Sbb200", "Sbb300"]
			signal_eff_histograms = {}
			for i, signal in enumerate(signals):
				signal_eff_histograms[signal] = convert_cutflow_to_efficiency(cutflow_histograms[signal])

			# Draw
			c = TCanvas("c_cutflow_eff_{}_{}".format(region, jet_type))
			l = TLegend(0.6, 0.6, 0.85, 0.85)
			l.SetFillColor(0)
			l.SetBorderSize(0)

			background_eff_stack.SetMinimum(-0.05)
			background_eff_stack.SetMaximum(1.1)
			background_eff_stack.Draw()

			total_background_eff_hist.SetMarkerStyle(20)
			total_background_eff_hist.SetMarkerSize(0)
			total_background_eff_hist.SetLineColor(seaborn.GetColorRoot("deep", 2))
			#total_background_eff_hist.SetFillColor(seaborn.GetColorRoot("pastel", 2))
			#total_background_eff_hist.SetFillStyle(3001)
			total_background_eff_hist.Draw("hist same")
			l.AddEntry(total_background_eff_hist, "Total background", "lf")
			for background in backgrounds:
				l.AddEntry(cutflow_histograms[background], style.background_legends[background], "lf")

			for signal in signals:
				signal_eff_histograms[signal].SetMarkerStyle(20 + i)
				signal_eff_histograms[signal].SetMarkerColor(seaborn.GetColorRoot("default", i))
				signal_eff_histograms[signal].SetLineColor(seaborn.GetColorRoot("default", i))
				signal_eff_histograms[signal].SetLineWidth(2)
				signal_eff_histograms[signal].Draw("hist same")
				l.AddEntry(signal_eff_histograms[signal], signal, "l")
			l.Draw()
			c.SaveAs("/uscms/home/dryu/DAZSLE/data/EventSelection/figures/{}.pdf".format(c.GetName()))

