import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy
#sys.path.append('/home/marc/code/python/')
#import PlottingFunctions
#import RootHelperFunctions
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections

input_folder = "/afs/cern.ch/user/d/dryu/DAZSLE/data/DDT"
output_folder = "/afs/cern.ch/user/d/dryu/DAZSLE/data/DDT/tmp"

rho_bins = [35, -7., -0.]
pt_bins = [8, 200., 1000.]
z_bins = {
	"N2":[750, 0., 0.75],
	"dcsv":[100, -1., 1.],
	"dsub":[100, -1., 1.],
}

rho_bins_fine = [280, -7.0, 0.0]
pt_bins_fine = [200, 400., 1000.]

def compute_ddt(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -7.0, -2.0, nPtBins, 200, 800)
	DDT.SetStats(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			p = array('d', [point])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] );
	return DDT


class DDTPoint:
	def __init__(self,r,p,n,w):
		self.rho = r
		self.pT = p
		self.z = n
		self.weight = w
	def GetSortNearestCrit(self, o):
		Drho = math.fabs(self.rho - o.rho)/5.
		DpT = math.fabs(self.pT - o.pT)/600.
		return (Drho*Drho) + (DpT*DpT)
	def GetSortNearestGlobal(self, pt, rho):
		Drho = math.fabs(self.rho - rho)/5.
		DpT = math.fabs(self.pT - pt)/600.
		return (Drho*Drho) + (DpT*DpT)
def SortByNearest(p1): # Compare to global (recently set) parameter
	return p1.GetSortNearestGlobal(gpT, grho)
def SortByZ(p1):
	return p1.z

def do_ddt_smoothing(rho_pt_points, jet_type, var="N2", drho_cut=None, wp=0.26):
	points = []
	points_0mod3 = []
	points_1mod3 = []
	points_2mod3 = []
	for sample in config.samples["qcd"]:
		f = TFile("{}/ddt_ntuple_{}.root".format(input_folder, sample))
		nevents = f.Get("NEvents").Integral()
		xs = cross_sections[sample]
		lumi_weight = 35900. * xs / nevents
		T = f.Get("ddttree")
		containers = {}
		vars = ["msd", "pt", "dcsv", "dsub", "N2"]
		for var in vars:
			containers[var] = array('d', [0.])
			T.SetBranchAddress("{}_{}".format(var, jet_type), containers[var])
		for entry in xrange(T.GetEntriesFast()):
			T.GetEntry(j)
			weight = containers["weight_trigger"][0] * lumi_weight
			if var == "N2":
				z = containers["N2"][0]
			elif var == "dcsv":
				z = containers["dcsv"][0]
			elif var == "dsub":
				z = containers["dsub"][0]
			points.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
			#if j % 3 == 0:
			#	points_0mod3.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
			#if j % 3 == 1:
			#	points_1mod3.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))
			#if j % 3 == 2:
			#	points_2mod3.append(DDTPoint(containers["rho"][0], containers["pt"][0], z, weight))

	points.sort(key=SortByZ)
	#points_0mod3.sort(key=SortByZ)
	#points_1mod3.sort(key=SortByZ)
	#points_2mod3.sort(key=SortByZ)

	# Calculate DDT points
	from joblib import Parallel, delayed
	smoothed_ddt_results = Parallel(n_jobs=4)(delayed(do_ddt_point)(rho_pt_point[0], rho_pt_point[1], points, drho_cut=drho_cut, wp=wp) for rho_pt_point in rho_pt_points)

	return smoothed_ddt_results	

def do_ddt_point(rho, pt, points, drho_cut=None, wp=0.26):
	# First loop: calculate sum of weights. 
	total_weight = 0.
	for p in points:
		drho = abs(p.rho - rho)
		dpt = abs(p.pt - pt)
		kweight = 1. / (0.001 + math.sqrt((drho/5.)**2 + (dpt/600.)**2))
		if drho_cut:
			kweight *= (drho < drho_cut)
		point_weight = p.weight * kweight
		total_weight += point_weight

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

	p_delta = (sum_weight - zp * total_weight) / last_weight
	this_z = points[last_point_index].z * p_delta + points[last_point_index-1].z * (1. - p_delta)
	return [rho, pt, this_z]

def do_ddt_simple(jet_type):
	H3 = TH3F("H3", ";Jet #rho;Jet p_{T} (GeV)", rho_bins[0], rho_bins[1], rho_bins[2], pt_bins[0], pt_bins[1], pt_bins[2], n2_bins[0], n2_bins[1], n2_bins[2])
	H31 = TH3F("H31", ";Jet #rho;Jet p_{T} (GeV)", rho_bins[0], rho_bins[1], rho_bins[2], pt_bins[0], pt_bins[1], pt_bins[2], n2_bins[0], n2_bins[1], n2_bins[2])
	H32 = TH3F("H32", ";Jet #rho;Jet p_{T} (GeV)", rho_bins[0], rho_bins[1], rho_bins[2], pt_bins[0], pt_bins[1], pt_bins[2], n2_bins[0], n2_bins[1], n2_bins[2])
	H33 = TH3F("H33", ";Jet #rho;Jet p_{T} (GeV)", rho_bins[0], rho_bins[1], rho_bins[2], pt_bins[0], pt_bins[1], pt_bins[2], n2_bins[0], n2_bins[1], n2_bins[2])
	H3.SetStats(0)
	H31.SetStats(0)
	H32.SetStats(0)
	H33.SetStats(0)

	for sample in config.samples["qcd"]:
		f = TFile("{}/ddt_ntuple_{}.root".format(input_folder, sample))
		nevents = f.Get("NEvents").Integral()
		xs = cross_sections[sample]
		T = f.Get("ddttree")
		containers = {}
		vars = ["msd", "pt", "dcsv", "dsub", "N2"]
		for var in vars:
			containers[var] = array('d', [0.])
			T.SetBranchAddress("{}_{}".format(var, args.jet_type), containers[var])

		for entry in xrange(T.GetEntriesFast()):
			T.GetEntry(entry)
			lumi_weight = 35900. * xs / nevents
			H3.Fill(containers["rho"][0], containers["pt"][0], containers["N2"][0], containers["weight"][0] * containers["weight_trigger"][0] * lumi_weight)
			if j % 3 == 0:
				H31.Fill(containers["rho"][0], containers["pt"][0], containers["N2"][0], containers["weight"][0] * containers["weight_trigger"][0] * lumi_weight)
			if j % 3 == 1:
				H32.Fill(containers["rho"][0], containers["pt"][0], containers["N2"][0], containers["weight"][0] * containers["weight_trigger"][0] * lumi_weight)
			if j % 3 == 2:
				H33.Fill(containers["rho"][0], containers["pt"][0], containers["N2"][0], containers["weight"][0] * containers["weight_trigger"][0] * lumi_weight)
	ddt = compute_ddt("N2DDT", Vcut, pt_bins[0], rho_bins[0], H3)
	ddt1 = compute_ddt("N2DDT1", Vcut, pt_bins[0], rho_bins[0], H31)
	ddt2 = compute_ddt("N2DDT2", Vcut, pt_bins[0], rho_bins[0], H32)
	ddt3 = compute_ddt("N2DDT3", Vcut, pt_bins[0], rho_bins[0], H33)

	return ddt, ddt1, ddt2, ddt3

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("--jet_type", type=str, help="AK8 or CA15")
	parser.add_argument("--run", action="store_true", help="Run locally")
	parser.add_argument("--crun", action="store_true", help="Run on HTCondor")
	parser.add_argument("--ncpu", type=int, default=4, help="Multiprocessing ncpus")
	parser.add_argument("--drho_cut", type=float, help="Add drho cut for smoothing")
	parser.add_argument("--smoothing_subjob", type=int, help="For condor, run a specific subjob")
	parser.add_argument("--var", type=str, help="N2, dcsv, or dsub")
	parser.add_argument("--condor", action="store_true", help="Flag for running on condor")
	args = parser.parse_args()

	if args.condor:
		input_folder = "./"
		output_folder = "./"

	if args.run:
		ddt, ddt1, ddt2, ddt3 = do_ddt_simple(args.jet_type)

		# WPs
		if args.var == "N2":
			wp = 0.26
		elif args.var == "dcsv" and args.jet_type == "AK8":
			wp = 0.015
		elif args.var == "dsub" and args.jet_type == "CA15":
			wp = 0.02
		else:
			print "ERROR : No WP established for var {}, jet type {}".format(args.var, args.jet_type)

		if args.smoothing_subjob:
			rho_pt_points_all = []
			for rho_bin in xrange(rho_bins[0] + 1):
				for pt_bin in xrange(pt_bins[0] + 1):
					rho = rho_bins[1] + rho_bin * (rho_bins[2] - rho_bins[1]) / rho_bins[0]
					pt = pt_bins[1] + pt_bin * (pt_bins[2] - pt_bins[1]) / pt_bins[0]
					rho_pt_points_all.append((rho, pt))
			points_per_job = 500
			rho_pt_points_subjob = []
			for i in xrange(points_per_job):
				j = args.smoothing_subjob*points_per_job + i
				rho_pt_points_subjob.append(rho_pt_points_all[j])
			if args.drho_cut:
				drho_cut = args.drho_cut
			else:
				drho_cut = None
			ddt_smooth_results = do_ddt_smoothing(rho_pt_points_subjob, args.jet_type, var=args.var, drho_cut=drho_cut)

		# Save
		output_file = TFile("{}/ddt_output.root".format(output_folder), "RECREATE")
		ddt.Write()
		ddt1.Write()
		ddt2.Write()
		ddt3.Write()

		if args.smoothing_subjob:
			import pickle
			with open('ddt_output_smooth.pkl', 'w') as f:
				pickle.dump(ddt_smooth_results, f)

	if args.crun:
		# Setup the smoothing jobs on condor
		#rho_bins_fine = [280, -7.0, 0.0]
		#pt_bins_fine = [200, 400., 1000.]
		import time
		submission_dir = "/afs/cern.ch/user/d/dryu/DAZSLE/data/DDT/condor/smoothing_{}/".format(time.time())
		cwd = os.getcwd()
		os.system("mkdir -pv {}".format(submission_dir))
		os.chdir(submission_dir)

		njobs = int(math.ceil(1. * len(rho_pt_points) / points_per_job))
		rho_pt_strings = []
		for ijob in njobs:
			job_rho_pt_strings = []
			for ipoint in xrange(points_per_job):
				i = ijob * points_per_job + ipoint
				job_rho_pt_strings.append(str(rho_pt_points[i][0]) + ":" + str(rho_pt_points[i][1]))
			job_rho_pt_string = ",".join(job_rho_pt_strings)
			rho_pt_strings.append(job_rho_pt_string)

		job_script_path = "{}/run_csubjob.sh".format(submission_dir)
		job_script = open(job_script_path, 'w')
		job_script.write("#!/bin/bash\n")
		command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/DDTMaker.py --run --condor --smoothing_subjob $1 --jet_type {} --var {}".format(args.jet_type, args.var)
		if args.drho_cut:
			command += " --drho_cut {}".format(args.drho_cut)
		job_script.write(command + "\n")
		job_script.close()

		submission_command = "csub {} --cmssw --no_retar -n {}".format(job_script_path, njobs)

		# Save csub command for resubmission attempts
		submission_script_path = "{}/csub_command.sh".format(submission_dir)
		submission_script = open(submission_script_path, "w")
		submission_script.write("#!/bin/bash\n")
		submission_script.write(submission_command + "\n")
		submission_script.close()

		os.system(submission_command)


