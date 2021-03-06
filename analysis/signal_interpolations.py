import os
import sys
from math import sqrt, floor
from ROOT import *
from DAZSLE.ZPrimePlusJet.histogram_interpolator import HistogramInterpolator 
#import DAZSLE.PhiBBPlusJet.analysis_configuration as config
import DAZSLE.ZPrimePlusJet.xbb_config as config
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsRootUtils.so"))
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()

# 0 = normalization
# 1 = mean
# 2 = gaussian mean
# 3 = lorentzian mean
# r = something related to the approximation... default is 4. Should be fixed.
def fvoigt(x, p, r=4):
	return p[0] * TMath.Voigt(x[0] - p[1], p[2], p[3], r)

def make_validation_plots(interpolated_histogram, simulated_histogram, save_tag, adjacent_histograms=None, x_range=None, do_fits=False):
	# 2D ratio plot
	ratio_histogram = interpolated_histogram.Clone()
	ratio_histogram.Divide(simulated_histogram)
	c_ratio = TCanvas("c_intp_validation_ratio_{}".format(save_tag), save_tag, 700, 500)
	c_ratio.SetRightMargin(0.2)
	ratio_histogram.Draw("colz")
	c_ratio.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.pdf".format(c_ratio.GetName()))
	c_ratio.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.eps".format(c_ratio.GetName()))
	c_ratio.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.C".format(c_ratio.GetName()))

	# 2D pull plot
	pull_histogram = interpolated_histogram.Clone()
	for xbin in xrange(1, pull_histogram.GetNbinsX() + 1):
		for ybin in xrange(1, pull_histogram.GetNbinsY() + 1):
			num = interpolated_histogram.GetBinContent(xbin, ybin)
			num_err = interpolated_histogram.GetBinError(xbin, ybin)
			den = simulated_histogram.GetBinContent(xbin, ybin)
			den_err = simulated_histogram.GetBinError(xbin, ybin)
			total_err = sqrt(num_err**2 + den_err**2)
			if total_err > 0:
				pull = (num - den) / total_err
			else:
				pull = 0.
			pull_histogram.SetBinContent(xbin, ybin, pull)
			pull_histogram.SetBinError(xbin, ybin, 0)

	c_pull = TCanvas("c_intp_validation_pull_{}".format(save_tag), save_tag, 700, 500)
	c_pull.SetRightMargin(0.2)
	pull_histogram.Draw("colz")
	c_pull.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.pdf".format(c_pull.GetName()))
	c_pull.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.eps".format(c_pull.GetName()))
	c_pull.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.C".format(c_pull.GetName()))

	return_data = {}

	# 1D plots
	for ptbin in xrange(0, interpolated_histogram.GetNbinsY() + 1):
		print "[debug] Starting ptbin {}".format(ptbin)
		return_data[ptbin] = {}
		if ptbin == 0:
			interpolated_histogram_1D = interpolated_histogram.ProjectionX("{}_int_ptbinAll".format(interpolated_histogram.GetName()))
			simulated_histogram_1D = simulated_histogram.ProjectionX("{}_sim_ptbinAll".format(simulated_histogram.GetName()))
		else:
			interpolated_histogram_1D = interpolated_histogram.ProjectionX("{}_int_ptbin{}".format(interpolated_histogram.GetName(), ptbin), ptbin, ptbin)
			simulated_histogram_1D = simulated_histogram.ProjectionX("{}_sim_ptbin{}".format(simulated_histogram.GetName(), ptbin), ptbin, ptbin)
		ratio_histogram_1D = interpolated_histogram_1D.Clone()
		ratio_histogram_1D.Divide(simulated_histogram_1D)
		if ptbin == 0:
			c_comparison = TCanvas("c_intp_validation_{}_ptbinAll".format(save_tag), save_tag, 800, 1000)
		else:
			c_comparison = TCanvas("c_intp_validation_{}_ptbin{}".format(save_tag, ptbin), save_tag, 800, 1000)
		l_comparison = TLegend(0.6, 0.6, 0.88, 0.8)
		l_comparison.SetFillColor(0)
		l_comparison.SetBorderSize(0)
		top = TPad("top", "top", 0., 0.5, 1., 1.)
		top.SetBottomMargin(0.02)
		top.Draw()
		top.cd()

		ymax = -1
		if interpolated_histogram_1D.GetMaximum() > ymax:
			ymax = interpolated_histogram_1D.GetMaximum()
		if simulated_histogram_1D.GetMaximum() > ymax:
			ymax = simulated_histogram_1D.GetMaximum()
		if adjacent_histograms:
			for adj_mass, adj_hist in adjacent_histograms.iteritems():
				if adj_hist.GetMaximum() > ymax:
					ymax = adj_hist.GetMaximum()
		ymax = ymax * 1.3

		interpolated_histogram_1D.SetMarkerStyle(24)
		interpolated_histogram_1D.SetMarkerColor(seaborn.GetColorRoot("default", 2))
		interpolated_histogram_1D.SetLineColor(seaborn.GetColorRoot("default", 2))
		interpolated_histogram_1D.SetLineWidth(2)
		interpolated_histogram_1D.GetXaxis().SetTitleSize(0)
		interpolated_histogram_1D.GetXaxis().SetLabelSize(0)
		interpolated_histogram_1D.SetMaximum(ymax)
		if x_range:
			interpolated_histogram_1D.GetXaxis().SetRangeUser(x_range[0], x_range[1])
		interpolated_histogram_1D.Draw()
		simulated_histogram_1D.SetMarkerStyle(20)
		simulated_histogram_1D.SetMarkerColor(seaborn.GetColorRoot("default", 3))
		simulated_histogram_1D.SetLineColor(seaborn.GetColorRoot("default", 3))
		simulated_histogram_1D.SetLineWidth(2)
		simulated_histogram_1D.Draw("same")
		l_comparison.AddEntry(simulated_histogram_1D, "Simulated", "lp")
		l_comparison.AddEntry(interpolated_histogram_1D, "Interpolated", "lp")

		if adjacent_histograms:
			adj_hist_proj = {}
			for adj_mass in sorted(adjacent_histograms.keys()):
				if ptbin == 0:
					adj_hist_proj[adj_mass] = adjacent_histograms[adj_mass].ProjectionX("{}_ptbinAll".format(adjacent_histograms[adj_mass].GetName()))
				else:
					adj_hist_proj[adj_mass] = adjacent_histograms[adj_mass].ProjectionX("{}_ptbin{}".format(adjacent_histograms[adj_mass].GetName(), ptbin), ptbin, ptbin)
				adj_hist_proj[adj_mass].SetMarkerStyle(20)
				adj_hist_proj[adj_mass].SetMarkerSize(0)
				adj_hist_proj[adj_mass].SetLineStyle(3)
				adj_hist_proj[adj_mass].SetLineColor(seaborn.GetColorRoot("pastel", 3))
				adj_hist_proj[adj_mass].SetLineWidth(3)
				adj_hist_proj[adj_mass].Draw("hist same")
				l_comparison.AddEntry(adj_hist_proj[adj_mass], "Sim, m={} GeV".format(adj_mass), "l")

		if do_fits:
			fits = {}
			mean = 0.5 * (interpolated_histogram_1D.GetMean() + simulated_histogram_1D.GetMean())
			rms = 0.5 * (interpolated_histogram_1D.GetRMS() + simulated_histogram_1D.GetRMS())
			for hist_name, hist in {"int":interpolated_histogram_1D, "sim":simulated_histogram_1D}.iteritems():
				return_data[ptbin][hist_name] = {}
				#fits[hist_name] = TF1("voigt_{}".format(hist_name), fvoigt, mean - 1.5*rms, mean + 1.5*rms, 4)
				#fits[hist_name].SetParameter(0, hist.Integral())
				#fits[hist_name].SetParameter(1, mean)
				#fits[hist_name].SetParameter(2, rms/2.)
				#fits[hist_name].SetParameter(3, rms/2.)
				#hist.Fit(fits[hist_name], "QR0")
				#if hist_name == "int":
				#	fits[hist_name].SetLineColor(seaborn.GetColorRoot("pastel", 2))
				#else:
				#	fits[hist_name].SetLineColor(seaborn.GetColorRoot("pastel", 3))
				#fits[hist_name].Draw("same")
				#fl = fits[hist_name].GetParameter(3)
				#dfl = fits[hist_name].GetParError(3)
				#fg = fits[hist_name].GetParameter(2)
				#dfg = fits[hist_name].GetParError(2)

				#fwhm = 0.5346 * fl + sqrt(0.2166 * fl**2 + fg**2)
				#dfwhm = fwhm * sqrt(dfg**2 * (fg / sqrt(0.2166*fl**2+fg**2))**2 + dfl**2 * (0.5346 + 0.2166 * fl / sqrt(0.2166*fl**2+fg**2))**2)

				print "[debug] Fit range {} - {}".format(mean - 1.5*rms, mean + 1.5*rms)
				fits[hist_name] = TF1("gauss_{}".format(hist_name), "gaus(0)", mean - 2*rms, mean + 2*rms)
				fits[hist_name].SetParameter(0, hist.Integral())
				fits[hist_name].SetParameter(1, mean)
				fits[hist_name].SetParameter(2, rms)
				fits[hist_name].SetParLimits(2, rms/2., 2. * rms)
				hist.Fit(fits[hist_name], "R0")
				if hist_name == "int":
					fits[hist_name].SetLineColor(seaborn.GetColorRoot("pastel", 2))
				else:
					fits[hist_name].SetLineColor(seaborn.GetColorRoot("pastel", 3))
				fits[hist_name].Draw("same")

				fwhm = fits[hist_name].GetParameter(2) * sqrt(2.)
				dfwhm = fits[hist_name].GetParError(2) * sqrt(2.)
				l_comparison.AddEntry(fits[hist_name], "Fit {}, FWHM={:0.2f}#pm{:0.2f}, x_{{0}}={:0.2f}#pm{:0.2f})".format(hist_name, fwhm, dfwhm, fits[hist_name].GetParameter(1), fits[hist_name].GetParError(1)), "l")

				return_data[ptbin][hist_name]["fwhm"] = fwhm
				return_data[ptbin][hist_name]["dfwhm"] = dfwhm
				return_data[ptbin][hist_name]["x0"] = fits[hist_name].GetParameter(1)
				return_data[ptbin][hist_name]["dx0"] = fits[hist_name].GetParError(1)
				return_data[ptbin][hist_name]["rms"] = hist.GetRMS()
				return_data[ptbin][hist_name]["norm"] = hist.Integral()
				return_data[ptbin][hist_name]["mean"] = hist.GetMean()

		l_comparison.Draw()

		c_comparison.cd()
		bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
		bottom.SetTopMargin(0.02)
		bottom.SetBottomMargin(0.2)
		bottom.Draw()
		bottom.cd()
		if x_range:
			ratio_histogram_1D.GetXaxis().SetRangeUser(x_range[0], x_range[1])
		ratio_histogram_1D.Draw()

		c_comparison.cd()
		c_comparison.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.pdf".format(c_comparison.GetName()))

		SetOwnership(top, False)
		SetOwnership(bottom, False)
		SetOwnership(c_comparison, False)
		top.IsA().Destructor(top)
		bottom.IsA().Destructor(bottom)
		c_comparison.IsA().Destructor(c_comparison)

	return return_data

def make_summary_plot(sim_hists, int_hists, save_tag):
	c = TCanvas("c_intp_summary_{}".format(save_tag), "c_intp_summary_{}".format(save_tag), 800, 600)
	l = TLegend(0.6, 0.2, 0.88, 0.8)
	l.SetFillColor(0)
	l.SetBorderSize(0)

	xmin = 0.
	xmax = 500.
	ymin = 0.
	ymax = max([h.GetMaximum() for h in sim_hists.values() + int_hists.values()]) * 1.3
	frame = TH1D("frame", "frame", 100, 0., 500.)
	frame.SetMinimum(ymin)
	frame.SetMaximum(ymax)
	frame.GetXaxis().SetTitle("m_{SD}^{PUPPI} [GeV]")

	for sim_mass in sorted(sim_hists.keys()):
		color = seaborn.GetColorRoot("cubehelixlarge", int(floor(30 * sim_mass / 500.)), 30)
		sim_hists[sim_mass].SetLineColor(color)
		sim_hists[sim_mass].SetLineStyle(1)
		sim_hists[sim_mass].SetLineWidth(2)
		sim_hists[sim_mass].SetMarkerStyle(20)
		sim_hists[sim_mass].SetMarkerSize(0)
		sim_hists[sim_mass].SetMarkerColor(color)
		sim_hists[sim_mass].Draw("hist same")
		#sim_hists[sim_mass].Draw("same")
		l.AddEntry(sim_hists[sim_mass], "m_{{X}}={} GeV".format(sim_mass), "l")

	for int_mass in sorted(int_hists.keys()):
		color = seaborn.GetColorRoot("cubehelixlarge", int(floor(30 * int_mass / 500.)), 30)
		int_hists[int_mass].SetLineColor(color)
		int_hists[int_mass].SetLineStyle(2)
		int_hists[int_mass].SetLineWidth(2)
		int_hists[int_mass].SetMarkerStyle(20)
		int_hists[int_mass].SetMarkerSize(0)
		int_hists[int_mass].SetMarkerColor(color)
		int_hists[int_mass].Draw("hist same")
		#int_hists[int_mass].Draw("same")
		l.AddEntry(int_hists[int_mass], "m_{{X}}={} GeV".format(int_mass), "l")

	l.Draw()
	c.SaveAs("/uscms/home/dryu/DAZSLE/data/Validation/figures/{}.pdf".format(c.GetName()))

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Perform signal interpolations on pT vs. mSD histograms.')
	parser.add_argument('--jet_type', type=str, help="AK8 or CA15")
	parser.add_argument('--interpolate', action='store_true', help="Perform signal interpolations")
	parser.add_argument('--interpolate_muCR', action='store_true', help="Perform signal interpolations")
	parser.add_argument('--validate', type=int, help="Validate interpolation at specified mass value (must be a simulated point with adjacent simulated points)")
	parser.add_argument('--plots', action='store_true', help="Validate interpolation at specified mass value (must be a simulated point with adjacent simulated points)")
	parser.add_argument('--input_masses', type=str, default="50,100,125,200,300,350,400,500", help="List of input masses (comma-separated)")
	parser.add_argument('--region', type=str, default="SR", help="SR or N2CR or muCR")
	output_group = parser.add_mutually_exclusive_group() 
	output_group.add_argument('--output_masses', type=str, help="List of output masses (comma-separated)")
	output_group.add_argument('--output_range', type=str, help="Range of output masses (e.g. 100,425,25 for 100-400 GeV in 25 GeV steps)")
	args = parser.parse_args()

	# Create list of input and output masses
	input_masses = [int(x) for x in args.input_masses.split(",")]
	if args.output_masses:
		output_masses = (int(x) for x in args.output_masses.split(","))
	elif args.output_range:
		output_args = [int(x) for x in args.output_range.split(",")]
		output_masses = []
		print output_args
		for x in xrange(output_args[0], output_args[1], output_args[2]):
			if not x in input_masses:
				output_masses.append(x)
	else:
		print "[signal_interpolations] ERROR : Must specify output_masses or output_range"
		sys.exit(1)

	#models = ["Sbb", "PSbb", "ZPrime"]
	models = ["Sbb", "PSbb"]

	if args.interpolate:
		# Input and output files (uses David's configuration. Replace if you are not David.)
		input_file = TFile(config.get_histogram_file(args.region, args.jet_type), "READ")
		output_file = TFile(config.get_interpolation_file(args.region, args.jet_type), "RECREATE")

		# Top-level loop
		for region in ["pass", "fail"]:
			for model in models:
				for syst in ["", "_JESUp", "_JERUp", "_PUUp", "_TriggerUp", "_JESDown", "_JERDown", "_PUDown", "_TriggerDown"]:
					input_histograms = {}
					output_histograms_1D = {}
					model_histogram = None
					for ptbin in xrange(1, 7):
						input_histograms[ptbin] = {}
						output_histograms_1D[ptbin] = {}
						for mass in input_masses:
							if model == "ZPrime" and mass > 300:
								continue
							if not model_histogram:
								model_histogram = input_file.Get("{}{}_{}{}".format(model, mass, region, syst))
								if not model_histogram:
									print "[signal_interpolations] ERROR : Couldn't find histogram {} in file {}".format("{}{}_{}{}".format(model, mass, region, syst), input_file.GetPath())
								model_histogram.SetName("model")
								model_histogram.SetDirectory(0)
								model_histogram.Reset()
								if model_histogram.GetNbinsY() != 6:
									print "ERROR : The model histogram appears to have !=6 y bins, which this code assumes."
									sys.exit(1)
							input_histograms[ptbin][mass] = input_file.Get("{}{}_{}{}".format(model, mass, region, syst)).ProjectionX("{}{}_{}_ptbin{}".format(model, mass, region, ptbin), ptbin, ptbin)
							input_histograms[ptbin][mass].SetDirectory(0)
							# Save a copy to the output file
							output_file.cd()
							input_file.Get("{}{}_{}{}".format(model, mass, region, syst)).Write()
						print input_histograms[ptbin]
						interpolator = HistogramInterpolator(input_histograms[ptbin])
						for mass in output_masses:
							if model == "ZPrime" and mass > 300:
								continue
							output_histograms_1D[ptbin][mass] = interpolator.make_interpolation(mass)
							output_histograms_1D[ptbin][mass].SetName("{}{}_{}{}_ptbin{}".format(model, mass, region, syst, ptbin))
							#output_file.cd()
							#output_histograms_1D[ptbin][mass].Write()

					# Put 1D histograms back together into 2D
					output_histograms = {}
					for mass in output_masses:
						if model == "ZPrime" and mass > 300:
							continue
						output_histograms[mass] = model_histogram.Clone()
						output_histograms[mass].SetName("{}{}_{}{}".format(model, mass, region, syst))
						for msdbin in xrange(1, model_histogram.GetNbinsX() + 1):
							for ptbin in xrange(1, model_histogram.GetNbinsY() + 1):
								output_histograms[mass].SetBinContent(msdbin, ptbin, output_histograms_1D[ptbin][mass].GetBinContent(msdbin))
								output_histograms[mass].SetBinError(msdbin, ptbin, output_histograms_1D[ptbin][mass].GetBinContent(msdbin))
						output_file.cd()
						output_histograms[mass].Write()
		input_file.Close()
		output_file.Close()

	if args.interpolate_muCR:
		# Input and output files (uses David's configuration. Replace if you are not David.)
		input_file = TFile(config.get_histogram_file("muCR", args.jet_type), "READ")
		output_file = TFile(config.get_interpolation_file("muCR", args.jet_type), "RECREATE")
		# Top-level loop
		for region in ["pass", "fail"]:
			for model in models:
				for syst in ['','_JERUp','_JERDown','_JESUp','_JESDown','_MuTriggerUp','_MuTriggerDown','_MuIDUp','_MuIDDown','_MuIsoUp','_MuIsoDown','_PUUp','_PUDown']:
					input_histograms = {}
					output_histograms = {}
					for mass in input_masses:
						if model == "ZPrime" and mass > 300:
							continue
						print "[debug] Trying to get " + "{}{}_{}{}".format(model, mass, region, syst) + " from file " + input_file.GetPath()
						input_histograms[mass] = input_file.Get("{}{}_{}{}".format(model, mass, region, syst)).Clone()
						if not input_histograms[mass]:
							print "ERROR : Couldn't find histogram {}{}_{}{} in file {}".format(model, mass, region, syst, input_file.GetPath())
							sys.exit(1)
						input_histograms[mass].SetDirectory(0)
						# Save a copy of the input to the output file
						output_file.cd()
						input_histograms[mass].Write()
					print input_histograms
					interpolator = HistogramInterpolator(input_histograms)
					for mass in output_masses:
						if model == "ZPrime" and mass > 300:
							continue
						output_histograms[mass] = interpolator.make_interpolation(mass)
						output_histograms[mass].SetName("{}{}_{}{}".format(model, mass, region, syst))
						output_file.cd()
						output_histograms[mass].Write()
		print "Saved interpolations to {}".format(output_file.GetPath())
		input_file.Close()
		output_file.Close()

	if args.validate:
		validation_mass = args.validate
		if not validation_mass in input_masses:
			print "ERROR : The validation mass must be in the list of input masses."
			sys.exit(1)
		input_file = TFile(config.get_histogram_file(args.region, args.jet_type), "READ")

		# Find the neighboring simulated masses
		left_mass = -1
		right_mass = -1
		input_masses_minus_validation = [x for x in input_masses if x != validation_mass]
		for i in xrange(len(input_masses_minus_validation) - 1):
			if input_masses_minus_validation[i] < validation_mass and validation_mass < input_masses_minus_validation[i+1]:
				left_mass = input_masses_minus_validation[i]
				right_mass = input_masses_minus_validation[i+1]
				break
		if left_mass == -1:
			print "ERROR : Couldn't find adjacent simulated masses."
			print "Validation mass = {}".format(validation_mass)
			print "Input masses = ",
			print input_masses_minus_validation
			sys.exit(1)

		# Top-level loop
		total_return_data = {}
		for region in ["pass", "fail"]:
			total_return_data[region] = {}
			for model in models:

				input_histograms = {}
				output_histograms_1D = {}
				model_histogram = None
				for ptbin in xrange(1, 7):
					input_histograms[ptbin] = {}
					for mass in [left_mass, right_mass]:
						if not model_histogram:
							model_histogram = input_file.Get("{}{}_{}".format(model, mass, region))
							model_histogram.SetName("model")
							model_histogram.SetDirectory(0)
							model_histogram.Reset()
							if model_histogram.GetNbinsY() != 6:
								print "ERROR : The model histogram appears to have !=6 y bins, which this code assumes."
								sys.exit(1)
						input_histograms[ptbin][mass] = input_file.Get("{}{}_{}".format(model, mass, region)).ProjectionX("{}{}_{}_ptbin{}".format(model, mass, region, ptbin), ptbin, ptbin)
						print "[debug] Name = {}".format(input_histograms[ptbin][mass].GetName())
					interpolator = HistogramInterpolator(input_histograms[ptbin])
					output_histograms_1D[ptbin] = interpolator.make_interpolation(validation_mass)
					output_histograms_1D[ptbin].SetName("{}{}_{}_ptbin{}".format(model, validation_mass, region, ptbin))

				# Put 1D histograms back together into 2D
				interpolated_histogram = model_histogram.Clone()
				interpolated_histogram.SetName("{}{}_{}_interpolated".format(model, validation_mass, region))
				for msdbin in xrange(1, model_histogram.GetNbinsX() + 1):
					for ptbin in xrange(1, model_histogram.GetNbinsY() + 1):
						interpolated_histogram.SetBinContent(msdbin, ptbin, output_histograms_1D[ptbin].GetBinContent(msdbin))
						interpolated_histogram.SetBinError(msdbin, ptbin, output_histograms_1D[ptbin].GetBinError(msdbin))

				# Get simulated histograms
				simulated_histogram = input_file.Get("{}{}_{}".format(model, validation_mass, region))
				adjacent_histograms = {}
				for adj_mass in [left_mass, right_mass]:
					adjacent_histograms[adj_mass] = input_file.Get("{}{}_{}".format(model, adj_mass, region))

				save_tag = "{}_{}_{}_{}".format(args.jet_type, region, model, validation_mass)
				total_return_data[region][model] = make_validation_plots(interpolated_histogram, simulated_histogram, save_tag, adjacent_histograms=adjacent_histograms, x_range=[left_mass-75.,right_mass+75.], do_fits=True)
				with open("peakdata_{}_{}_{}{}.txt".format(args.jet_type, region, model, validation_mass), 'w') as f:
					f.write("rms[\"{}\"][\"{}\"][\"int\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["int"]["rms"]))
					f.write("rms[\"{}\"][\"{}\"][\"sim\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["sim"]["rms"]))
					f.write("fwhm[\"{}\"][\"{}\"][\"int\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["int"]["fwhm"]))
					f.write("fwhm[\"{}\"][\"{}\"][\"sim\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["sim"]["fwhm"]))
					f.write("dfwhm[\"{}\"][\"{}\"][\"int\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["int"]["dfwhm"]))
					f.write("dfwhm[\"{}\"][\"{}\"][\"sim\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["sim"]["dfwhm"]))
					f.write("norm[\"{}\"][\"{}\"][\"int\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["int"]["norm"]))
					f.write("norm[\"{}\"][\"{}\"][\"sim\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["sim"]["norm"]))
					f.write("mean[\"{}\"][\"{}\"][\"int\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["int"]["mean"]))
					f.write("mean[\"{}\"][\"{}\"][\"sim\"][{}]={}\n".format(args.jet_type, region, validation_mass, total_return_data[region][model][0]["sim"]["mean"]))
		#print total_return_data
	if args.plots:
		sim_file = TFile(config.get_histogram_file(args.region, args.jet_type), "READ")
		int_file = TFile(config.get_interpolation_file(args.region, args.jet_type), "READ")
		for region in ["pass", "fail"]:
			for model in models:
				for ptbin in xrange(0, 7):
					sim_hists = {}
					int_hists = {}
					for sim_mass in input_masses:
						if ptbin == 0:
							sim_hists[sim_mass] = sim_file.Get("{}{}_{}".format(model, sim_mass, region)).ProjectionX("{}{}_{}_ptbinAll".format(model, sim_mass, region))
						else:
							sim_hists[sim_mass] = sim_file.Get("{}{}_{}".format(model, sim_mass, region)).ProjectionX("{}{}_{}_ptbin{}".format(model, sim_mass, region, ptbin), ptbin, ptbin)
						sim_hists[sim_mass].SetDirectory(0)
					for int_mass in output_masses:
						if ptbin == 0:
							int_hists[int_mass] = int_file.Get("{}{}_{}".format(model, int_mass, region)).ProjectionX("{}{}_{}_ptbinAll".format(model, sim_mass, region))
						else:
							int_hists[int_mass] = int_file.Get("{}{}_{}".format(model, int_mass, region)).ProjectionX("{}{}_{}_ptbin{}".format(model, sim_mass, region, ptbin), ptbin, ptbin)
						int_hists[int_mass].SetDirectory(0)
					if ptbin == 0:
						save_tag = "{}_{}_{}_ptbinAll".format(region, model, args.jet_type)
					else:
						save_tag = "{}_{}_{}_ptbin{}".format(region, model, args.jet_type, ptbin)
					make_summary_plot(sim_hists, int_hists, save_tag)


