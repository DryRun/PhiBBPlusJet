# MC plots for getting a basic idea of the shapes of distributions

import os
import sys
from array import array
import re
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

def Project3DHist(hist, msd_range, pt_range):
	print "[debug] 3D integral = " + str(hist.Integral())
	ybin_min = 1
	ybin_max = hist.GetNbinsY()
	for ybin in xrange(1, hist.GetNbinsY() + 1):
		this_msd = hist.GetYaxis().GetBinCenter(ybin)
		if this_msd < msd_range[0]:
			ybin_min = ybin
		if this_msd > msd_range[1]:
			ybin_max = ybin

	zbin_min = 1
	zbin_max = hist.GetNbinsZ()
	for zbin in xrange(1, hist.GetNbinsZ() + 1):
		this_pt = hist.GetZaxis().GetBinCenter(zbin)
		if this_pt < pt_range[0]:
			zbin_min = zbin
		if this_pt > pt_range[1]:
			zbin_max = zbin
	hist_proj = hist.ProjectionX(hist.GetName() + "_projx", ybin_min, ybin_max, zbin_min, zbin_max	)
	print "[debug] Projection {}, {} and {}, {} integral = {}".format(ybin_min, ybin_max, zbin_min, zbin_max, hist_proj.Integral())
	return hist_proj

def NormalizedPlots():
	xvars = ["dcsv"] # dcsv_vs_msd_vs_pt
	backgrounds = ["qcd","tqq","wqq","zqq","hbb","stqq","vvqq"]
	signals = {"AK8":["Sbb50", "Sbb100", "Sbb125"], "CA15":["Sbb125", "Sbb200", "Sbb300"]}
	for xvar in xvars:
		for jet_type in ["AK8", "CA15"]:
			histogram_file = TFile(config.get_histogram_file("Preselection", jet_type))
			background_histograms = {}
			total_bkgd_histogram = None
			for background in backgrounds:
				if background == "hbb":
					for subbackground in ["hqq125", "tthqq125", "vbfhqq125", "whqq125", "zhqq125"]:
						hname = subbackground + "_" + xvar + "_vs_msd_vs_pt"
						if not background in background_histograms:
							this_3d_histogram = histogram_file.Get(hname)
							this_3d_histogram.SetDirectory(0)
						else:
							this_3d_histogram.Add(histogram_file.Get(hname))
				else:
					hname = background + "_" + xvar  + "_vs_msd_vs_pt"
					this_3d_histogram = histogram_file.Get(hname)
					this_3d_histogram.SetDirectory(0)
				this_histogram = Project3DHist(this_3d_histogram, [40, 500], [450, 1000])
				this_histogram.SetDirectory(0)
				background_histograms[background] = this_histogram
				if not total_bkgd_histogram:
					total_bkgd_histogram = background_histograms[background].Clone()
				else:
					total_bkgd_histogram.Add(background_histograms[background])

			# Normalize to unity
			for background in backgrounds:
				background_histograms[background].Scale(1. / total_bkgd_histogram.Integral())
			total_bkgd_histogram.Scale(1. / total_bkgd_histogram.Integral())

			signal_histograms = {}
			for signal in signals[jet_type]:
				hname = "{}_{}".format(signal, xvar)
				signal_histograms[signal] = histogram_file.Get(hname)
				signal_histograms[signal].SetDirectory(0)
				signal_histograms[signal].Scale(1. / signal_histograms[signal].Integral())

			# Style and stack
			for background in backgrounds:
				background_histograms[background].SetFillColor(style.background_colors[background])
				background_histograms[background].SetLineColor(1)
			bkgds_sorted = backgrounds
			bkgds_sorted.sort(key=lambda x: background_histograms[x].Integral())
			bkgd_stack = THStack("bkgd_stack", "bkgd_stack")
			for bkgd in bkgds_sorted:
				bkgd_stack.Add(background_histograms[bkgd])
			total_bkgd_histogram.SetLineWidth(2)
			total_bkgd_histogram.SetLineColor(1)

			for i, signal in enumerate(signals[jet_type]):
				signal_histograms[signal].SetLineColor(seaborn.GetColorRoot("default", i))
				signal_histograms[signal].SetLineWidth(2)
				signal_histograms[signal].SetLineStyle(2)
				signal_histograms[signal].SetMarkerStyle(20)
				signal_histograms[signal].SetMarkerSize(1)

			# Legend
			l = TLegend(0.25, 0.55, 0.5, 0.85)
			l.SetFillColor(0)
			l.SetBorderSize(0)
			for bkgd in reversed(bkgds_sorted):
				legend_entry = bkgd
				if bkgd in style.background_legends:
					legend_entry = style.background_legends[bkgd]
				l.AddEntry(background_histograms[bkgd], legend_entry, "f")
			for signal in signals[jet_type]:
				l.AddEntry(signal_histograms[signal], signal, "l")

			# Drawing
			xmin = total_bkgd_histogram.GetXaxis().GetXmin()
			xmax = total_bkgd_histogram.GetXaxis().GetXmax()
			ymin = min([x.GetMinimum() for x in signal_histograms.values()] + [total_bkgd_histogram.GetMinimum()])
			ymax = min([x.GetMaximum() for x in signal_histograms.values()] + [total_bkgd_histogram.GetMaximum()])
			ymax *= 2

			cname = "c_normalizedmc_{}_{}".format(xvar, jet_type)
			c = TCanvas(cname, xvar, 800, 600)
			bkgd_stack.Draw("hist")
			bkgd_stack.SetMinimum(ymin)
			bkgd_stack.SetMaximum(ymax)
			bkgd_stack.GetXaxis().SetTitle(style.axis_titles[xvar])
			bkgd_stack.GetYaxis().SetTitleOffset(1.1)
			bkgd_stack.GetYaxis().SetTitle("Normalized to unit area")
			bkgd_stack.Draw("hist")
			total_bkgd_histogram.Draw("hist same")
			for signal in signals[jet_type]:
				signal_histograms[signal].Draw("hist same")
			l.Draw()
			c.SaveAs("/uscms/home/dryu/DAZSLE/data/EventSelection/figures/{}.pdf".format(cname))

			ROOT.SetOwnership(bkgd_stack, False)
			ROOT.SetOwnership(total_bkgd_histogram, False)
			for h in background_histograms.values():
				ROOT.SetOwnership(h, False)
			ROOT.SetOwnership(c, False)

if __name__ == "__main__":
	NormalizedPlots()



