import os
import sys
from DAZSLE.PhiBBPlusJet.combine_helper.combine_project import CombineProject
from DAZSLE.PhiBBPlusJet.combine_helper.region import Region
from ROOT import TFile
import math

# Load histograms
backgrounds = ["qcd", "tqq", "zqq", "wqq"]
signals = []
for mass in [50,100,125,200,300,350,400,500]:
	signals.append("Sbb{}".format(mass))
for mass in [75, 125, 150, 175, 225, 250, 300, 400]:
	signals.append("ZPrime{}".format(mass))

datacard_topdir = "/uscms/home/dryu/PhiBB2017/data/datacards/simple"
histogram_dir = "/uscms/home/dryu/PhiBB2017/data/histograms"

default_pt_categories = [[450., 500.], 	[500., 550.], [550., 600.], [600., 675.], [675., 800.], [800., 1000.]]

# Slice up the pt vs msd histogram
def SliceSignalHistogram(hist2d, pt_categories=default_pt_categories):
	# Make sure the requested pt bins correspond to boundaries
	pt_axis = hist2d.GetYaxis()
	histogram_pt_boundaries = []
	for bin in xrange(1, pt_axis.GetNbins() + 1):
		histogram_pt_boundaries.append(pt_axis.GetBinLowEdge(bin))
	histogram_pt_boundaries.append(pt_axis.GetBinUpEdge(pt_axis.GetNbins()))

	for pt_category in pt_categories:
		if not pt_category[0] in histogram_pt_boundaries:
			print "[run_histograms::SliceSignalHistogram] ERROR : Bin boundary {} does not correspond to a histogram bin boundary.".format(pt_category[0])
			print histogram_pt_boundaries
			sys.exit(1)
		if not pt_category[1] in histogram_pt_boundaries:
			print "[run_histograms::SliceSignalHistogram] ERROR : Bin boundary {} does not correspond to a histogram bin boundary.".format(pt_category[1])
			print histogram_pt_boundaries
			sys.exit(1)

	histogram_slices = {}
	for islice in xrange(len(pt_categories)):
		ptmin = pt_categories[islice][0]
		ptmax = pt_categories[islice][1]
		binmin = 1e10
		binmax = -1e10
		for bin in xrange(1, hist2d.GetNbinsY() + 1):
			low_edge = hist2d.GetYaxis().GetBinLowEdge(bin)
			up_edge = hist2d.GetYaxis().GetBinUpEdge(bin)
			# Is this bin inside this pt slice (+epsilon)?
			if ptmin - 1.e-5 < low_edge and up_edge < ptmax + 1.e-5:
				if bin < binmin:
					binmin = bin
				if bin > binmax:
					binmax = bin
		histogram_slices[islice] = hist2d.ProjectionX(hist2d.GetName() + "_" + str(islice), binmin, binmax)
	return histogram_slices

for jet_type in ["AK8", "CA15"]:
	for signal in signals:
		histogram_file = TFile("{}/histograms_SR_{}.root".format(histogram_dir, jet_type), "READ")
		datacard_dir = datacard_topdir + "/" + signal
		os.system("mkdir -pv {}".format(datacard_dir))
		project = CombineProject(datacard_dir)
		xvar = project.create_xvar("msd", [40., 600.])

		# Construct regions
		region_containers = {}
		for box in ["pass1", "pass2"]:

			# Load 2D histograms
			h2d = {}
			h2d_slices = {}
			for process in backgrounds + ["data_obs", signal]:
				h2d[process] = histogram_file.Get("{}_{}_{}".format(process, box, "pt_vs_msd"))
				if not h2d[process]:
					print "[simple_limits] ERROR : Couldn't load histogram {} from file {}".format("{}_{}_{}".format(process, box, "pt_vs_msd"), histogram_file.GetPath())
					sys.exit(1)
				h2d_slices[process] = SliceSignalHistogram(h2d[process], pt_categories=default_pt_categories)

			for ptcat in xrange(len(default_pt_categories)):
				# Apply rho boundaries, and set data to zero if the prediction is zero
				if jet_type == "AK8":
					rho_range = [-6.0, -2.1]
				elif jet_type == "CA15":
					rho_range = [-4.7, -1.0]

				for xbin in xrange(1, h2d_slices["data_obs"][ptcat].GetNbinsX()+1):
					bin_mass = h2d_slices[process][ptcat].GetXaxis().GetBinCenter(xbin)
					bin_pt = default_pt_categories[ptcat][0] + (default_pt_categories[ptcat][1]-default_pt_categories[ptcat][0]) * 0.3
					bin_rho = 2. * math.log(bin_mass / bin_pt)
					if bin_rho < rho_range[0] or bin_rho > rho_range[1]:
						for process in backgrounds + ["data_obs", signal]:
							h2d_slices[process][ptcat].SetBinContent(xbin, 0.)
							h2d_slices[process][ptcat].SetBinError(xbin, 0.)

					# Background prediction in this bin
					bin_total_bkgd = 0.
					for process in backgrounds:
						bin_total_bkgd += h2d_slices[process][ptcat].GetBinContent(xbin)
					if bin_total_bkgd == 0:
						h2d_slices["data_obs"][ptcat].SetBinContent(xbin, 0)
						h2d_slices["data_obs"][ptcat].SetBinError(xbin, 0)

				region_name = "{}_ptcat{}".format(box, ptcat)
				region_containers[region_name] = Region(region_name)
				region_containers[region_name].add_xvar(xvar)

				for background in backgrounds:
					# Get pt vs msd histogram
					region_containers[region_name].add_background(background, h2d_slices[background][ptcat])
				region_containers[region_name].add_data("data_obs", h2d_slices["data_obs"][ptcat])
				region_containers[region_name].add_signal(signal, h2d_slices[signal][ptcat])
				project.add_region(region_containers[region_name])
		project.write()

# Sum of pass1 and pass2, i.e. no b-tag categories
datacard_topdir = "/uscms/home/dryu/PhiBB2017/data/datacards/simple_nobtag"
for jet_type in ["AK8", "CA15"]:
	for signal in signals:
		histogram_file = TFile("{}/histograms_SR_{}.root".format(histogram_dir, jet_type), "READ")
		datacard_dir = datacard_topdir + "/" + signal
		os.system("mkdir -pv {}".format(datacard_dir))
		project = CombineProject(datacard_dir)
		xvar = project.create_xvar("msd", [40., 600.])

		# Construct regions
		region_containers = {}

		# Load 2D histograms
		h2d = {}
		h2d_slices = {}
		for process in backgrounds + ["data_obs", signal]:
			h2d[process] = histogram_file.Get("{}_{}_{}".format(process, "pass1", "pt_vs_msd"))
			if not h2d[process]:
				print "[simple_limits] ERROR : Couldn't load histogram {} from file {}".format("{}_{}_{}".format(process, "pass1", "pt_vs_msd"), histogram_file.GetPath())
				sys.exit(1)
			h2d[process].Add(histogram_file.Get("{}_{}_{}".format(process, "pass2", "pt_vs_msd")))
			h2d_slices[process] = SliceSignalHistogram(h2d[process], pt_categories=default_pt_categories)

		for ptcat in xrange(len(default_pt_categories)):
			# Apply rho boundaries, and set data to zero if the prediction is zero
			if jet_type == "AK8":
				rho_range = [-6.0, -2.1]
			elif jet_type == "CA15":
				rho_range = [-4.7, -1.0]

			for xbin in xrange(1, h2d_slices["data_obs"][ptcat].GetNbinsX()+1):
				bin_mass = h2d_slices[process][ptcat].GetXaxis().GetBinCenter(xbin)
				bin_pt = default_pt_categories[ptcat][0] + (default_pt_categories[ptcat][1]-default_pt_categories[ptcat][0]) * 0.3
				bin_rho = 2. * math.log(bin_mass / bin_pt)
				if bin_rho < rho_range[0] or bin_rho > rho_range[1]:
					for process in backgrounds + ["data_obs", signal]:
						h2d_slices[process][ptcat].SetBinContent(xbin, 0.)
						h2d_slices[process][ptcat].SetBinError(xbin, 0.)

				# Background prediction in this bin
				bin_total_bkgd = 0.
				for process in backgrounds:
					bin_total_bkgd += h2d_slices[process][ptcat].GetBinContent(xbin)
				if bin_total_bkgd == 0:
					h2d_slices["data_obs"][ptcat].SetBinContent(xbin, 0)
					h2d_slices["data_obs"][ptcat].SetBinError(xbin, 0)

			region_name = "ptcat{}".format(ptcat)
			region_containers[region_name] = Region(region_name)
			region_containers[region_name].add_xvar(xvar)

			for background in backgrounds:
				# Get pt vs msd histogram
				region_containers[region_name].add_background(background, h2d_slices[background][ptcat])
			region_containers[region_name].add_data("data_obs", h2d_slices["data_obs"][ptcat])
			region_containers[region_name].add_signal(signal, h2d_slices[signal][ptcat])
			project.add_region(region_containers[region_name])
	project.write()


