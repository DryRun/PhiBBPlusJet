import os
import sys
from DAZSLE.PhiBBPlusJet.combine_helper.combine_project import CombineProject
from DAZSLE.PhiBBPlusJet.combine_helper.region import Region
from DAZSLE.PhiBBPlusJet.combine_helper.rhalphabet_region import RhalphabetRegion
from ROOT import TFile, TGraph, TColor, gROOT, kOrange, kGreen, TCanvas, TLegend, TH1D, kBlack, kBlue, gStyle
gROOT.SetBatch(True)
import math
from array import array

from MyTools.RootUtils.seaborn_colors import SeabornColors
seaborn_colors = SeabornColors()
seaborn_colors.load_palette("Blues_d", palette_dir=os.path.expandvars("$CMSSW_BASE/python/MyTools/RootUtils/seaborn_palettes"))
seaborn_colors.load_palette("Reds_d", palette_dir=os.path.expandvars("$CMSSW_BASE/python/MyTools/RootUtils/seaborn_palettes"))
seaborn_colors.load_palette("Oranges_d", palette_dir=os.path.expandvars("$CMSSW_BASE/python/MyTools/RootUtils/seaborn_palettes"))
seaborn_colors.load_palette("Greens_d", palette_dir=os.path.expandvars("$CMSSW_BASE/python/MyTools/RootUtils/seaborn_palettes"))
seaborn_colors.load_palette("Purples_d", palette_dir=os.path.expandvars("$CMSSW_BASE/python/MyTools/RootUtils/seaborn_palettes"))
gStyle.SetOptStat(False) 
gStyle.SetOptTitle(0)

signals = ["Sbb50", "Sbb125", "Sbb200", "Sbb300"]
signals.extend(["ZPrime75", "ZPrime125", "ZPrime175", "ZPrime300"])

optimization_wps = []
n2ddt_wps = [0.05, 0.15, 0.26]
dbtag_wps= [0.7, 0.8, 0.9]
for dbtag_wp in dbtag_wps:
	for n2wp_dbpass in n2ddt_wps:
		for n2wp_dbfail in n2ddt_wps:
			opt_wp_string = "dbtag{}_n2wpdbpass{}_n2wpdbfail{}".format(dbtag_wp, n2wp_dbpass, n2wp_dbfail)
			optimization_wps.append(opt_wp_string)
# Also non-b-tagged
for n2wp in n2ddt_wps:
	opt_wp_string = "dbtagnone_n2wp{}".format(n2wp)
	optimization_wps.append(opt_wp_string)

exp_limits = {}
best_exp_limits = {}
for jet_type in ["AK8", "CA15"]:
	exp_limits[jet_type] = {}
	best_exp_limits[jet_type] = {}
	for signal in signals:
		exp_limits[jet_type][signal] = {}

		for opt_wp_name in optimization_wps:

			limit_dir = "/uscms/home/dryu/PhiBB2017/data/datacards/optimization/{}/{}/{}".format(opt_wp_name, jet_type, signal)
			limit_file = TFile("{}/higgsCombineTest.Asymptotic.mH120.root".format(limit_dir), "READ")
			if not limit_file.IsOpen():
				print "{} \t {} \t {} = [No limit file]".format(jet_type, signal, opt_wp_name)
				exp_limits[jet_type][signal][opt_wp_name] = 1.e20
				continue
			limit_tree = limit_file.Get("limit")
			if not limit_tree:
				print "{} \t {} \t {} = No limit tree. See {}".format(jet_type, signal, opt_wp_name, limit_file.GetPath())
				exp_limits[jet_type][signal][opt_wp_name] = 1.e20
				continue
			limit_tree.GetEntry(2)
			exp_limits[jet_type][signal][opt_wp_name] = limit_tree.GetLeaf("limit").GetValue(0)
			print "{} \t {} \t {} = {}".format(jet_type, signal, opt_wp_name, exp_limits[jet_type][signal][opt_wp_name])
		best_exp_limits[jet_type][signal] = min(exp_limits[jet_type][signal].values())

# Limit plots
for jet_type in ["AK8", "CA15"]:
	x = array('d', [75, 125, 175, 300])
	y_opt = array('d', [])
	for mZp in x:
		y_opt.append(best_exp_limits[jet_type]["ZPrime{}".format(int(mZp))])
	y_ref = array('d', [])
	for mZp in x:
		y_ref.append(exp_limits[jet_type]["ZPrime{}".format(int(mZp))]["dbtagnone_n2wp0.05"])
	tg_limits_opt = TGraph(len(x), x, y_opt)
	tg_limits_ref = TGraph(len(x), x, y_ref)
	c = TCanvas("c_limitgain_{}".format(jet_type), "c", 800, 600)
	c.SetLogy()
	frame = TH1D("frame", "frame", 100, 50., 400.)
	frame.GetXaxis().SetTitle("m_{Z'} [GeV]")
	frame.GetYaxis().SetTitle("Exp. limit on #sigma BR #epsilon [pb]")
	frame.SetMinimum(1.0)
	frame.SetMaximum(50.)
	frame.Draw()

	tg_limits_ref.SetLineWidth(2)
	tg_limits_ref.SetLineStyle(2)
	tg_limits_ref.SetLineColor(seaborn_colors.get_root_color("Reds_d", 2))
	tg_limits_ref.Draw("l")

	tg_limits_opt.SetLineWidth(2)
	tg_limits_opt.SetLineStyle(2)
	tg_limits_opt.SetLineColor(seaborn_colors.get_root_color("Greens_d", 2))
	tg_limits_opt.Draw("l")

	l = TLegend(0.7, 0.7, 0.85, 0.85)
	l.SetBorderSize(0)
	l.SetFillColor(0)
	l.AddEntry(tg_limits_ref, "w/o b-tag SR", "l")
	l.AddEntry(tg_limits_opt, "w/ b-tag SR ", "l")
	l.Draw()
	c.SaveAs(os.path.expandvars("$HOME/PhiBB2017/data/limits/optimization/{}.pdf".format(c.GetName())))

