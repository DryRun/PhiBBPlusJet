# Get signal acceptance and efficiencies
import os
import sys
import time
from ROOT import *
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
gROOT.SetBatch(True)


if __name__ == "__main__":
	timestamp = str(time.time()).replace(".","_")
	# Get the signal efficiencies from cutflow histograms and print
	jet_types = ["AK8", "CA15"]
	backgrounds = ["QCD"]

	#this_background_efficiencies = {}
	this_background_efficiencies_opt = {}

	#for jet_type in jet_types:
	#	this_background_efficiencies[jet_type] = {}
	#	for background in backgrounds:			
	#		histogram_file = TFile("~/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root".format(background, #jet_type), "READ")
	#		cutflow_histogram = histogram_file.Get("CutFlowCounter_EventSelector_SR_weighted")
	#		inclusive = cutflow_histogram.GetBinContent(1)
	#		pass_histogram = histogram_file.Get("h_SR_{}_pass".format(jet_type))
	#		if not pass_histogram:
	#			print "ERROR : Couldn't get histogram {} from file {}".format("h_SR_{}_pass".format(jet_type), #histogram_file.GetPath())
	#			sys.exit(1)
	#		final = pass_histogram.Integral()
	#		this_background_efficiencies[jet_type][background] = float(final)/float(inclusive)

	# tau21 optimization
	tau21_values = [0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7]
	dcsv_values = [0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]
	for tau21_value in tau21_values:
		this_background_efficiencies_opt[tau21_value] = {}
		for dcsv_value in dcsv_values:
			this_background_efficiencies_opt[tau21_value][dcsv_value] = {}
			for jet_type in jet_types:
				this_background_efficiencies_opt[tau21_value][dcsv_value][jet_type] = {}
				for background in backgrounds:
					histogram_file = TFile("~/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root".format(background, jet_type), "READ")
					cutflow_histogram = histogram_file.Get("CutFlowCounter_EventSelector_SR_tau21ddt{}_weighted".format(tau21_value))
					inclusive = cutflow_histogram.GetBinContent(1)
					pass_histogram = histogram_file.Get("h_SR_tau21ddt{}_{}_pass_dcsv{}".format(tau21_value, jet_type, dcsv_value))
					if not pass_histogram:
						print "ERROR : Couldn't get histogram {} from file {}".format("h_SR_tau21ddt{}_{}_pass_dcsv{}".format(tau21_value, jet_type, dcsv_value), histogram_file.GetPath())
						sys.exit(1)
					final = pass_histogram.Integral()
					this_background_efficiencies_opt[tau21_value][dcsv_value][jet_type][background] = float(final)/float(inclusive)

	#print this_background_efficiencies
	#print this_background_efficiencies_opt

	for jet_type in jet_types:
		for background in backgrounds:
			h = TH2D("h_{}".format(background), "h_{}".format(background), len(tau21_values), -0.5, len(tau21_values) - 0.5, len(dcsv_values), -0.5, len(dcsv_values) - 0.5)
			for xval, tau21_value in enumerate(tau21_values):
				xbin = xval + 1
				h.GetXaxis().SetBinLabel(xbin, str(tau21_value))
				for yval, dcsv_value in enumerate(dcsv_values):
					ybin = yval + 1
					h.GetYaxis().SetBinLabel(ybin, str(dcsv_value))
					h.SetBinContent(xbin, ybin, this_background_efficiencies_opt[tau21_value][dcsv_value][jet_type][background])
			c = TCanvas("c_opt_background_eff_{}_{}".format(background, jet_type), "c_opt_background_eff_{}_{}".format(background, jet_type), 800, 600)
			c.SetRightMargin(0.2)
			h.GetXaxis().SetTitle("#tau_{21}^{DDT}")
			h.GetYaxis().SetTitle("Double CSV")
			h.GetZaxis().SetTitle("Signal Efficiency")
			h.Draw("colz text")
			c.SaveAs("/uscms/home/dryu/DAZSLE/data/LimitSetting/figures/Optimization/{}.pdf".format(c.GetName()))
