import os
import sys
import ROOT
from ROOT import *
import DAZSLE.ZPrimePlusJet.xbb_config as config
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/DAZSLE/CMSSW_7_4_7/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

for jet_type in ["AK8", "CA15"]:
	print jet_type
	f = TFile("/uscms/home/dryu/DAZSLE/data/LimitSetting/Xbb_inputs/histograms_muCR_{}.root".format(jet_type), "READ")
	pass_yields = {}
	fail_yields = {}
	total_yields = {}
	h_pass = {}
	h_fail = {}
	h_total = {}
	for signal in config.simulated_signal_names:
		h_pass[signal] = f.Get("{}_pass".format(signal)).ProjectionX()
		h_fail[signal] = f.Get("{}_fail".format(signal)).ProjectionX()
		h_total[signal] = h_pass[signal].Clone()
		h_total[signal].SetName("h_{}_total")
		h_total[signal].Add(h_fail[signal])
		pass_yields[signal] = h_pass[signal].Integral()
		fail_yields[signal] = h_fail[signal].Integral()
		total_yields[signal] = h_total[signal].Integral()
	print "Signal\tPass\tFail\tTotal"
	for signal in config.simulated_signal_names:
		print "{}\t{}\t{}\t{}".format(signal, pass_yields[signal], fail_yields[signal], total_yields[signal])

	c = TCanvas("c_muonCR_msd_{}".format(jet_type), "c_muonCR_msd_{}".format(jet_type), 700, 500)
	frame = TH1D("frame", "frame", 100, 0., 500.)
	frame.SetMinimum(1.e-3)
	frame.SetMaximum(10.)
	frame.Draw()
	c.SetLogy()
	for style_counter, signal in enumerate(config.simulated_signal_names):
		if "PS" in signal:
			continue
		h_total[signal].SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(config.simulated_signal_names) / 2 + 3))
		h_total[signal].SetMarkerColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(config.simulated_signal_names) / 2 + 3))
		h_total[signal].SetMarkerStyle(20+style_counter)
		h_total[signal].Draw("hist same")
	c.SaveAs("/uscms/home/dryu/DAZSLE/data/EventSelection/figures/{}.pdf".format(c.GetName()))
