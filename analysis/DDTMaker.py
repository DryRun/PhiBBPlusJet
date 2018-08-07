import ROOT
from ROOT import *
import os
import sys
import math
from array import array
#import functools
#from multiprocessing import Pool
#import scipy
#sys.path.append('/home/marc/code/python/')
#import PlottingFunctions
#import RootHelperFunctions
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections

input_folder = os.path.expandvars("$HOME/PhiBB2017/data/DDT")
output_folder = os.path.expandvars("$HOME/PhiBB2017/data/DDT/tmp")
ddt_folder = os.path.expandvars("$HOME/PhiBB2017/data/DDT/")


rho_bins = [21, -7., -0.]
pt_bins = [8, 200., 1000.]
rho_bins_array = array('d', [])
for i in xrange(rho_bins[0]+1):
	rho_bins_array.append(rho_bins[1] + (rho_bins[2]-rho_bins[1])/rho_bins[0]*i)
print rho_bins_array
pt_bins_array = array('d', [200., 400., 450., 500., 550., 600., 675., 800., 1000.])

z_bins = {
	"N2":[750, 0., 0.75],
	"dcsv":[100, -1., 1.],
	"dsub":[100, -1., 1.],
}
z_bins_array = {}
for zvar, z_bins_simple in z_bins.iteritems():
	z_bins_array[zvar] = array("d", [])
	for i in xrange(z_bins_simple[0]+1):
		z_bins_array[zvar].append(z_bins_simple[1] + (z_bins_simple[2]-z_bins_simple[1]) / z_bins_simple[0] * i)
	print zvar
	print z_bins_array[zvar]
rho_bins_fine = [280, -7.0, 0.0]
pt_bins_fine = [200, 400., 1000.]

def compute_ddt(name, wp, nPtBins, nRhoBins, H):
	#DDT = TH2F(name, "", nRhoBins, -7.0, -2.0, nPtBins, 200, 800)
	DDT = H.Project3D("yx")
	DDT.SetName(name)
	DDT.Reset()	
	DDT.SetStats(0)
	DDT.SetDirectory(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			if proj.Integral() == 0:
				print "[compute_ddt] WARNING : N2 integral=0 for xbin={}, ybin={}. Setting DDT to -1000".format(x+1, y+1)
				DDT.SetBinContent( x+1, y+1, -1000)
				continue
			p = array('d', [wp])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] )
	return DDT


class DDTPoint:
	def __init__(self,r,p,n,w):
		self.rho = r
		self.pt = p
		self.z = n
		self.weight = w
	def GetSortNearestCrit(self, o):
		Drho = math.fabs(self.rho - o.rho)/5.
		Dpt = math.fabs(self.pt - o.pt)/600.
		return (Drho*Drho) + (Dpt*Dpt)
	def GetSortNearestGlobal(self, pt, rho):
		Drho = math.fabs(self.rho - rho)/5.
		Dpt = math.fabs(self.pt - pt)/600.
		return (Drho*Drho) + (Dpt*Dpt)
def SortByNearest(p1): # Compare to global (recently set) parameter
	return p1.GetSortNearestGlobal(gpT, grho)
def SortByZ(p1):
	return p1.z

def do_ddt_smoothing(rho_pt_points, jet_type, zvar="N2", drho_cut=None, wp=0.26):
	print "[do_ddt_smoothing] INFO : Doing smoothing for {} points, jet type {}, zvar {}, drho_cut {}, wp {}".format(len(rho_pt_points), jet_type, zvar, drho_cut, wp)
	points = []
	points_0mod3 = []
	points_1mod3 = []
	points_2mod3 = []
	for sample in config.samples["qcd"]:
		print "[do_ddt_smoothing] INFO : Processing sample {}".format(sample)
		f = TFile("{}/ddt_ntuple_{}.root".format(input_folder, sample))

		nevents = f.Get("NEvents").Integral()
		xs = cross_sections[sample]
		lumi_weight = 35900. * xs / nevents

		T = f.Get("ddttree")
		containers = {}
		global_vars = ["weight_pu"]
		for tree_var in global_vars:
			containers[tree_var] = array('d', [0.])
			T.SetBranchAddress("{}".format(tree_var), containers[tree_var])
		jet_vars = ["msd", "pt", "rho", "dcsv", "dsub", "N2", "weight_trigger", "weight"]
		for tree_var in jet_vars:
			containers[tree_var] = array('d', [0.])
			T.SetBranchAddress("{}_{}".format(tree_var, args.jet_type), containers[tree_var])

		print_every = int(math.ceil(T.GetEntriesFast() / 20))
		for entry in xrange(T.GetEntriesFast()):
			if entry % print_every == 0:
				print "[do_ddt_smoothing] On entry {} / {}".format(entry, T.GetEntriesFast())
			T.GetEntry(j)
			weight = containers["weight_trigger"][0] * lumi_weight
			if zvar == "N2":
				z = containers["N2"][0]
			elif zvar == "dcsv":
				z = containers["dcsv"][0]
			elif zvar == "dsub":
				z = containers["dsub"][0]
			else:
				print "[do_ddt_smoothing] ERROR : Unknown zvar {}".format(zvar)
				sys.exit(1)
			points.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
			#if j % 3 == 0:
			#	points_0mod3.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
			#if j % 3 == 1:
			#	points_1mod3.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
			#if j % 3 == 2:
			#	points_2mod3.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
	print "[do_ddt_smoothing] INFO : Done loading input data. Starting smoothing."
	points.sort(key=SortByZ)
	#points_0mod3.sort(key=SortByZ)
	#points_1mod3.sort(key=SortByZ)
	#points_2mod3.sort(key=SortByZ)

	# Calculate DDT points
	#from joblib import Parallel, delayed
	#smoothed_ddt_results = Parallel(n_jobs=4)(delayed(do_ddt_point)(rho_pt_point[0], rho_pt_point[1], points, drho_cut=drho_cut, wp=wp) for rho_pt_point in rho_pt_points)
	smoothed_ddt_results = []
	for rho_pt_point in rho_pt_points:
		smoothed_ddt_results.append(do_ddt_point(rho_pt_point[0], rho_pt_point[1], points, drho_cut=drho_cut, wp=wp))

	return smoothed_ddt_results	

def do_ddt_point(rho, pt, points, drho_cut=None, wp=0.26):
	print "[do_ddt_point] INFO : Called with rho={}, pt={}".format(rho, pt)
	# First loop: calculate sum of weights. 
	total_weight = 0.
	for p in points:
		drho = abs(p.rho - rho)
		dpt = abs(p.pt - pt)
		kweight = 1. / (0.001 + math.sqrt((drho/5.)**2 + (dpt/600.)**2))
		if drho_cut != None:
			kweight *= (drho < drho_cut)
		point_weight = p.weight * kweight
		#print "[debug] point_weight = {} * {} = {}".format(p.weight, kweight, point_weight)
		total_weight += point_weight

	if total_weight == 0:
		print "[do_ddt_point] ERROR : total_weight == 0 for rho={}, pt={}, drho_cut={}. Returning -1000.".format(rho, pt, drho_cut)
		return [rho, pt, -1000.]

	# Second loop: determine index of points on either side of total_weight*wp
	sum_weight = 0.
	last_weight = -1.
	last_point_index = -1
	for i, p in enumerate(points):
		kweight = 1. / (0.001 + math.sqrt((drho/5.)**2 + (dpt/600.)**2))
		if drho_cut:
			kweight *= (drho < drho_cut)
		point_weight = p.weight * kweight
		sum_weight += point_weight
		if sum_weight / total_weight < wp:
			continue
		z = p.z
		last_point_index = i
		last_weight = point_weight
		break
	if last_point_index < 0:
		print "ERROR : Didn't find z value corresponding to WP {}".format(wp)
		sys.exit(1)

	p_delta = (sum_weight - wp * total_weight) / last_weight
	this_z = points[last_point_index].z * p_delta + points[last_point_index-1].z * (1. - p_delta)
	return [rho, pt, this_z]

def get_ddt3dhist_path(jet_type, zvar):
	path = "{}/DDT_3D_hists_{}_{}".format(output_folder, jet_type, zvar, wp)
	path += ".root"
	return path

def make_3dhists(jet_type, zvar="N2"):
	if jet_type == "AK8":
		dbtag_var = "dcsv"
	elif jet_type == "CA15":
		dbtag_var = "dsub"

	H3 = TH3F("H3", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H31 = TH3F("H31", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H32 = TH3F("H32", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H33 = TH3F("H33", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])

	H3_dbtag_pass = TH3F("H3_dbtag_pass", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H31_dbtag_pass = TH3F("H31_dbtag_pass", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H32_dbtag_pass = TH3F("H32_dbtag_pass", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H33_dbtag_pass = TH3F("H33_dbtag_pass", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])

	H3_dbtag_fail = TH3F("H3_dbtag_fail", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H31_dbtag_fail = TH3F("H31_dbtag_fail", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H32_dbtag_fail = TH3F("H32_dbtag_fail", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
	H33_dbtag_fail = TH3F("H33_dbtag_fail", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])

	H_nevents = TH2F("H_nevents", ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array)

	H3_samples = {}
	H3_samples_dbtag_pass = {}
	H3_samples_dbtag_fail = {}
	for sample in config.samples["qcd"]:
		print "Processing sample {}".format(sample)
		H3_samples[sample] = TH3F("H3_{}".format(sample), ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
		H3_samples[sample].SetDirectory(0)
		H3_samples_dbtag_pass[sample] = TH3F("H3_{}_dbtag_pass".format(sample), ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
		H3_samples_dbtag_pass[sample].SetDirectory(0)
		H3_samples_dbtag_fail[sample] = TH3F("H3_{}_dbtag_fail".format(sample), ";Jet #rho;Jet p_{T} (GeV)", len(rho_bins_array)-1, rho_bins_array, len(pt_bins_array)-1, pt_bins_array, len(z_bins_array[zvar])-1, z_bins_array[zvar])
		H3_samples_dbtag_fail[sample].SetDirectory(0)
		f = TFile("{}/ddt_ntuple_{}.root".format(input_folder, sample))
		if not f.IsOpen():
			print "ERROR : Can't open file {}".format(f.GetPath())
		nevents = f.Get("NEvents").Integral()
		xs = cross_sections[sample]
		T = f.Get("ddttree")
		containers = {}
		global_vars = ["weight_pu"]
		for tree_var in global_vars:
			containers[tree_var] = array('d', [0.])
			T.SetBranchAddress("{}".format(tree_var), containers[tree_var])
		jet_vars = ["msd", "pt", "rho", "dcsv", "dsub", "N2", "weight_trigger", "weight"]
		for tree_var in jet_vars:
			containers[tree_var] = array('d', [0.])
			T.SetBranchAddress("{}_{}".format(tree_var, args.jet_type), containers[tree_var])
		weight_debug = []
		total_weight_debug = []
		for entry in xrange(T.GetEntriesFast()):
			#for entry in xrange(10000):
			if entry % int(math.ceil(T.GetEntriesFast() / 20)) == 0:
				print "On entry {} / {}".format(entry, T.GetEntriesFast())
			T.GetEntry(entry)
			lumi_weight = 35900. * xs / nevents
			weight_debug.append(containers["weight"][0] * containers["weight_pu"][0] * containers["weight_trigger"][0])
			total_weight = containers["weight"][0] * containers["weight_pu"][0] * containers["weight_trigger"][0] * lumi_weight
			H_nevents.Fill(containers["rho"][0], containers["pt"][0])
			H3.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
			H3_samples[sample].Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
			if entry % 3 == 0:
				H31.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
			if entry % 3 == 1:
				H32.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
			if entry % 3 == 2:
				H33.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)

			if containers[dbtag_var][0] > 0.8:
				H3_dbtag_pass.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				H3_samples_dbtag_pass[sample].Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				if entry % 3 == 0:
					H31_dbtag_pass.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				if entry % 3 == 1:
					H32_dbtag_pass.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				if entry % 3 == 2:
					H33_dbtag_pass.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
			elif containers[dbtag_var][0] < 0.8:
				H3_dbtag_fail.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				H3_samples_dbtag_fail[sample].Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				if entry % 3 == 0:
					H31_dbtag_fail.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				if entry % 3 == 1:
					H32_dbtag_fail.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
				if entry % 3 == 2:
					H33_dbtag_fail.Fill(containers["rho"][0], containers["pt"][0], containers[zvar][0], total_weight)
		weight_debug.sort()
		print weight_debug[-20:]
	# Save big histograms
	f_h3 = TFile(get_ddt3dhist_path(jet_type, zvar), "RECREATE")
	H3.Write()
	H31.Write()
	H32.Write()
	H33.Write()
	H3_dbtag_pass.Write()
	H31_dbtag_pass.Write()
	H32_dbtag_pass.Write()
	H33_dbtag_pass.Write()
	H3_dbtag_fail.Write()
	H31_dbtag_fail.Write()
	H32_dbtag_fail.Write()
	H33_dbtag_fail.Write()
	H_nevents.Write()
	for sample, hist in H3_samples.iteritems():
		hist.Write()
	for sample, hist in H3_samples_dbtag_pass.iteritems():
		hist.Write()
	for sample, hist in H3_samples_dbtag_fail.iteritems():
		hist.Write()
	f_h3.Close()

def get_ddttransf_path(jet_type ,zvar, wp, dbtag_pass=False, dbtag_fail=False):
	path = "{}/DDT_{}_{}_wp{}".format(ddt_folder, zvar, jet_type, wp)
	if dbtag_pass:
		path += "_dbtag_pass"
	elif dbtag_fail:
		path += "_dbtag_fail"
	path += ".root"
	return path


def make_ddt_simple(jet_type, wp, zvar="N2", dbtag_pass=False, dbtag_fail=False):
	f_h3 = TFile(get_ddt3dhist_path(jet_type, zvar), "READ")
	if dbtag_pass:
		hist_suffix = "_dbtag_pass"
	elif dbtag_fail:
		hist_suffix = "_dbtag_fail"
	else:
		hist_suffix = ""
	H3 = f_h3.Get("H3" + hist_suffix)
	# Subtract off samples with large errors
	for sample in ["QCD_HT100to200","QCD_HT200to300","QCD_HT300to500"]:
		print f_h3.Get("H3_{}{}".format(sample, hist_suffix))
		H3.Add(f_h3.Get("H3_{}{}".format(sample, hist_suffix)), -1)
	H31 = f_h3.Get("H31" + hist_suffix)
	H32 = f_h3.Get("H32" + hist_suffix)
	H33 = f_h3.Get("H33" + hist_suffix)
	ddt = compute_ddt("DDT", wp, pt_bins[0], rho_bins[0], H3)
	ddt1 = compute_ddt("DDT1", wp, pt_bins[0], rho_bins[0], H31)
	ddt2 = compute_ddt("DDT2", wp, pt_bins[0], rho_bins[0], H32)
	ddt3 = compute_ddt("DDT3", wp, pt_bins[0], rho_bins[0], H33)

	# Save
	print "Saving to " + get_ddttransf_path(jet_type ,zvar, wp, dbtag_pass=dbtag_pass, dbtag_fail=dbtag_fail)
	output_file = TFile(get_ddttransf_path(jet_type ,zvar, wp, dbtag_pass=dbtag_pass, dbtag_fail=dbtag_fail), "RECREATE")
	ddt.Write()
	ddt1.Write()
	ddt2.Write()
	ddt3.Write()

	# Number of events and GetBinContent**2/GetBinError**2 (effective number of events. content=sum(w), err=sqrt(sum(w^2))). Let w=constant, then content=N*w, err^2=N*w^2)
	h2_rho_pt = H3.Project3D("yx")
	h2_nevents_eff = h2_rho_pt.Clone()
	h2_nevents_eff.Reset()
	for xbin in xrange(1, h2_nevents_eff.GetNbinsX()+1):
		for ybin in xrange(1, h2_nevents_eff.GetNbinsY()+1):
			content = h2_rho_pt.GetBinContent(xbin, ybin)
			err = h2_rho_pt.GetBinError(xbin, ybin)
			if content > 0:
				nevents_eff = content**2 / err**2
			else:
				nevents_eff = 0.
			h2_nevents_eff.SetBinContent(xbin, ybin, nevents_eff)
	h2_nevents_eff.Write()
	output_file.Close()
	f_h3.Close()

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("--jet_type", type=str, help="AK8 or CA15")
	parser.add_argument("--run_hists", action="store_true", help="Run simple binned DDT (step 1: long step to make 3D hists)")
	parser.add_argument("--run_ddt_simple", action="store_true", help="Run simple binned DDT (step 2: short step to compute quantiles)")
	parser.add_argument("--run_smoothing", action="store_true", help="Run smoothed DDT")
	parser.add_argument("--wp", type=float, help="Efficiency working point")
	parser.add_argument("--crun", action="store_true", help="Run on HTCondor")
	parser.add_argument("--ncpu", type=int, default=4, help="Multiprocessing ncpus")
	parser.add_argument("--drho_cut", type=float, help="Add drho cut for smoothing")
	parser.add_argument("--smoothing_subjob", type=int, default=-1, help="For condor, run a specific subjob")
	parser.add_argument("--zvar", type=str, help="N2, dcsv, or dsub")
	parser.add_argument("--dbtag_pass", action="store_true", help="Apply double b-tag pass cut to all events")
	parser.add_argument("--dbtag_fail", action="store_true", help="Apply double b-tag fail cut to all events")
	parser.add_argument("--condor", action="store_true", help="Flag for running on condor")
	args = parser.parse_args()

	if args.condor:
		input_folder = "./"
		output_folder = "./"

	# WPs
	if args.wp:
		wp = args.wp
	elif args.zvar == "N2":
		wp = 0.26
	elif args.zvar == "dcsv" and args.jet_type == "AK8":
		wp = 0.015
	elif args.zvar == "dsub" and args.jet_type == "CA15":
		wp = 0.02
	else:
		print "ERROR : No WP established for var {}, jet type {}".format(args.zvar, args.jet_type)
		sys.exit(1)

	# Setup pipelining parameters
	if args.run_smoothing or args.crun:
		rho_pt_points_all = []
		for rho_bin in xrange(rho_bins_fine[0] + 1):
			for pt_bin in xrange(pt_bins_fine[0] + 1):
				rho = rho_bins_fine[1] + rho_bin * (rho_bins_fine[2] - rho_bins_fine[1]) / rho_bins_fine[0]
				pt = pt_bins_fine[1] + pt_bin * (pt_bins_fine[2] - pt_bins_fine[1]) / pt_bins_fine[0]
				rho_pt_points_all.append((rho, pt))
		points_per_job = 1
		n_subjobs = int(math.ceil(1. * len(rho_pt_points_all) / points_per_job))
		print "[debug] n_subjobs = {}".format(n_subjobs)
		rho_pt_points_subjob = {}
		for isubjob in xrange(n_subjobs):
			rho_pt_points_subjob[isubjob] = []
			for ipoint in xrange(points_per_job):
				j = isubjob * points_per_job + ipoint
				if j >= len(rho_pt_points_all):
					break
				rho_pt_points_subjob[isubjob].append(rho_pt_points_all[j])

	if args.run_hists:
		make_3dhists(args.jet_type, args.zvar)
	if args.run_ddt_simple:
		make_ddt_simple(args.jet_type, wp, args.zvar)
		make_ddt_simple(args.jet_type, wp, args.zvar, dbtag_pass=True)
		make_ddt_simple(args.jet_type, wp, args.zvar, dbtag_fail=True)


	if args.run_smoothing:
		if args.smoothing_subjob >= 0:
			print "[debug] rho_pt points for this subjob: "
			print rho_pt_points_subjob[args.smoothing_subjob]
			if args.drho_cut:
				print "Using rho cut {}".format(args.drho_cut)
				drho_cut = args.drho_cut
			else:
				drho_cut = None
			ddt_smooth_results = do_ddt_smoothing(rho_pt_points_subjob[args.smoothing_subjob], args.jet_type, zvar=args.zvar, drho_cut=drho_cut)

			import pickle
			with open('ddt_output_smooth.pkl', 'w') as f:
				pickle.dump(ddt_smooth_results, f)
		else:
			print "ERROR : --run_smoothing must be accompanied by --smoothing_subjob (takes too long in one go)"
			sys.exit(1)

	if args.crun:
		# Setup the smoothing jobs on condor
		#rho_bins_fine = [280, -7.0, 0.0]
		#pt_bins_fine = [200, 400., 1000.]
		import time
		submission_dir = os.path.expandvars("$HOME/PhiBB2017/data/DDT/condor/smoothing_{}/").format(time.time())
		cwd = os.getcwd()
		os.system("mkdir -pv {}".format(submission_dir))
		os.chdir(submission_dir)

		input_files = []
		for sample in config.samples["qcd"]:
			input_files.append(os.path.expandvars("$HOME/PhiBB2017/data/DDT/ddt_ntuple_{}.root").format(sample))

		job_script_path = "{}/run_csubjob.sh".format(submission_dir)
		job_script = open(job_script_path, 'w')
		job_script.write("#!/bin/bash\n")
		command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/DDTMaker.py --run_smoothing --condor --smoothing_subjob $1 --jet_type {} --var {}".format(args.jet_type, args.zvar)
		if args.drho_cut:
			command += " --drho_cut {}".format(args.drho_cut)
		job_script.write(command + "\n")
		job_script.close()

		submission_command = "csub {} --cmssw --no_retar -n {} -F {} -p".format(job_script_path, n_subjobs, ",".join(input_files))

		# Save csub command for resubmission attempts
		submission_script_path = "{}/csub_command.sh".format(submission_dir)
		submission_script = open(submission_script_path, "w")
		submission_script.write("#!/bin/bash\n")
		submission_script.write(submission_command + "\n")
		submission_script.close()

		os.system(submission_command)


