# Make a simple ntuple from QCD MC baconbits, for input to DDTMaker.py
import os
import sys
import ROOT
from DAZSLE.PhiBBPlusJet.analysis_base import AnalysisBase
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.PhiBBPlusJet.bacon_event_selector import *
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
from math import ceil, sqrt,floor
import array

import ROOT
from ROOT import *
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/HistogramManager.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsRootUtils.so"))
gInterpreter.Declare("#include \"MyTools/AnalysisTools/interface/EventSelector.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsAnalysisTools.so"))
ROOT.gInterpreter.Declare("#include \"DAZSLE/PhiBBPlusJet/interface/BaconData.h\"")
ROOT.gInterpreter.Declare("#include \"DAZSLE/PhiBBPlusJet/interface/BaconEventCutFunctions.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
#from ROOT import gInterpreter, gSystem, gROOT, gStyle, Root, TCanvas, TLegend, TH1F, TFile, TGraphErrors
gROOT.SetBatch(ROOT.kTRUE);

class DDTNtupler(AnalysisBase):
	def __init__(self, sample_name, tree_name="otree", output_filename=None):
		super(DDTNtupler, self).__init__(tree_name=tree_name)
		self._data = BaconData(self._chain)
		self._output_file = TFile(output_filename, "RECREATE")

		self._sample_name = sample_name
		self._input_nevents = 0

	def add_file(self, filename):
		super(DDTNtupler, self).add_file(filename)
		f = ROOT.TFile.Open(filename, "READ")
		if f.Get("NEvents").Integral() == 0:
			print "[DDTNtupler::add_file] ERROR : NEvents.Integral() == 0 for file " + filename
			sys.exit(1)
		self._input_nevents += f.Get("NEvents").Integral()
		f.Close()

	def start(self):
		self._processed_events = 0
		self._output_file.cd()
		self._output_tree = TTree("ddttree", "ddttree")
		self._containers = {}
		branches_double = ["rho", "pt", "msd", "N2", "weight_trigger", "weight"] # kfNLO
		for jet_type in ["AK8", "CA15"]:
			self._containers[jet_type] = {}
			for branch in branches_double:
				self._containers[jet_type][branch] = array.array("d", [0.])
				self._output_tree.Branch(branch + "_" + jet_type, self._containers[jet_type][branch], branch + "_" + jet_type + "/D")

		branches_int = []
		for jet_type in ["AK8", "CA15"]:
			self._containers[jet_type] = {}
			for branch in branches_int:
				self._containers[jet_type][branch] = array.array("i", [0])
				self._output_tree.Branch(branch + "_" + jet_type, self._containers[jet_type][branch], branch + "_" + jet_type + "/I")

		branches_global_double = ["pfmet", "weight_pu"]
		for branch in branches_global_double:
			self._containers[branch] = array.array("d", [0.])
			self._output_tree.Branch(branch, self._containers[branch], branch + "/D")

		branches_global_int = ["n_el", "n_mu", "n_tau"]
		for branch in branches_global_double:
			self._containers[branch] = array.array("i", [0])
			self._output_tree.Branch(branch, self._containers[branch], branch + "/I")

		f_pu = TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/puWeights_All.root", "read")
		self._h_pu_weight = f_pu.Get("puw")
		self._h_pu_weight.SetDirectory(0)
		f_pu.Close()

		self._trig_eff = {}
		for jet_type in ["AK8", "CA15"]:
			# Trigger efficiency weight stuff
			if jet_type == "AK8":
				f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_AK8_SingleMuon_Run2016_V2p1_v03.root", "read")
				trig_den = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
				trig_num = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
			elif jet_type == "CA15":
				f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_CA15_SingleMuon_Run2016_V2p4_v08.root", "read")
				trig_den = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
				trig_num = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
			trig_den.SetDirectory(0)
			trig_num.SetDirectory(0)
			trig_den.RebinX(2)
			trig_num.RebinX(2)
			trig_den.RebinY(5)
			trig_num.RebinY(5)
			self._trig_eff[jet_type] = ROOT.TEfficiency()
			if (ROOT.TEfficiency.CheckConsistency(trig_num, trig_den)):
				self._trig_eff[jet_type] = ROOT.TEfficiency(trig_num, trig_den)
				self._trig_eff[jet_type].SetDirectory(0)
			f_trig.Close()


	def run(self, max_nevents=-1, first_event=0):
		if max_nevents > 0:
			limit_nevents = min(max_nevents, self._chain.GetEntries())
		else:
			limit_nevents = self._chain.GetEntries()

		n_checkpoints = 20
		print_every = int(ceil(1. * limit_nevents / n_checkpoints))

		print "[DDTNtupler::run] INFO : Running loop over tree from event {} to {}".format(first_event, limit_nevents - 1)

		self.start_timer()
		for entry in xrange(first_event, limit_nevents):
			self.print_progress(entry, first_event, limit_nevents, print_every)
			self._data.GetEntry(entry)

			npu = min(self._data.npu, 49.5)
			print self._h_pu_weight.FindBin(npu)
			print self._h_pu_weight.GetBinContent(self._h_pu_weight.FindBin(npu))
			self._containers["weight_pu"][0] = self._h_pu_weight.GetBinContent(self._h_pu_weight.FindBin(npu))
			#self._containers["kfNLO"] = 1.
			self._containers["pfmet"][0] = self._data.pfmet

			print "[debug] self._data.AK8Puppijet0_msd_puppi = " + str(self._data.AK8Puppijet0_msd_puppi)
			self._containers["AK8"]["msd"][0] = self._data.AK8Puppijet0_msd_puppi
			self._containers["AK8"]["pt"][0] = self._data.AK8Puppijet0_pt
			self._containers["AK8"]["rho"][0] = self._data.AK8Puppijet0_rho
			self._containers["AK8"]["N2"][0] = self._data.AK8Puppijet0_N2sdb1
			self._containers["AK8"]["dcsv"][0] = self._data.AK8Puppijet0_doublecsv
			self._containers["AK8"]["dsub"][0] = self._data.AK8Puppijet0_doublesub


			trigger_mass_AK8 = min(self._data.AK8Puppijet0_msd, 300.)
			trigger_pt_AK8 = max(200., min(self._data.AK8Puppijet0_pt, 1000.))
			self._containers["AK8"]["weight_trigger"][0] = self._trig_eff["AK8"].GetEfficiency(self._trig_eff["AK8"].FindFixBin(trigger_mass_AK8, trigger_pt_AK8))
			self._containers["AK8"]["weight"][0] = 1.

			self._containers["CA15"]["msd"] = self._data.CA15Puppijet0_msd_puppi
			self._containers["CA15"]["pt"] = self._data.CA15Puppijet0_pt
			self._containers["CA15"]["rho"] = self._data.CA15Puppijet0_rho
			self._containers["CA15"]["N2"] = self._data.CA15Puppijet0_N2sdb1
			self._containers["CA15"]["dcsv"][0] = self._data.CA15Puppijet0_doublecsv
			self._containers["CA15"]["dsub"][0] = self._data.CA15Puppijet0_doublesub

			trigger_mass_CA15 = min(self._data.CA15Puppijet0_msd, 300.)
			trigger_pt_CA15 = max(200., min(self._data.CA15Puppijet0_pt, 1000.))
			self._containers["CA15"]["weight_trigger"] = self._trig_eff["CA15"].GetEfficiency(self._trig_eff["CA15"].FindFixBin(trigger_mass_CA15, trigger_pt_CA15))
			self._containers["CA15"]["weight"][0] = 1.

			self._output_tree.Fill()

	def finish(self):
		self._output_file.cd()
		self._output_tree.Write()
		self._output_file.Close()

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Produce and plot ieta-iphi histograms to look for buggy events')
	input_group = parser.add_mutually_exclusive_group() 
	input_group.add_argument('--all', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--samples', type=str, help="Sample name(s), comma separated. Must be a key in analysis_configuration.skims.")
	input_group.add_argument('--files', type=str, help="Input file name(s), comma separated")
	parser.add_argument('--label', type=str, help="If running with --files, need to specify a label manually, in lieu of the sample names, for the output file naming.")
	action_group = parser.add_mutually_exclusive_group() 
	action_group.add_argument('--run', action="store_true", help="Run on condor")
	action_group.add_argument('--condor_run', action="store_true", help="Run on condor")
	parser.add_argument('--output_folder', type=str, help="Output folder")
	args = parser.parse_args()

	if args.run or args.condor_run:
		samples = []
		sample_files = {} # Dictionary is sample : [list of files in sample]
		if args.all:
			samples.extend(config.samples["qcd"])
			for sample in config.samples["qcd"]:
				sample_files[sample] = config.skims[sample]
		elif args.samples:
			samples = args.samples.split(",")
			for sample in samples:
				sample_files[sample] = config.skims[sample]
		elif args.files:
			files = args.files.split(",")
			for filename in files:
				if args.label:
					this_sample = args.label
				else:
					print "[event_selection_histograms] ERROR : When running with --files option, you must specify a label for the output!"
					sys.exit(1)
				if not this_sample in sample_files:
					sample_files[this_sample] = []
				sample_files[this_sample].append(filename)
			samples = sample_files.keys()
		print "List of input samples: ",
		print samples
		print "List of samples and files: ",
		print sample_files

	if args.run:
		#from joblib import Parallel
		#from joblib import delayed
		for sample in samples:
			print "\n *** Running sample {}".format(sample)
			tree_name = "Events"
			# Sanity check: make sure tree exists in file
			for filename in sample_files[sample]:
				print "[event_selection_histograms] INFO : Checking contents of file {}".format(filename)
				f = ROOT.TFile.Open(filename, "READ")
				t = f.Get(tree_name)
				if not t:
					if tree_name == "otree":
						backup_tree_name = "Events"
					else:
						backup_tree_name = "otree"
					t_backup = f.Get(backup_tree_name)
					if t_backup:
						print "[setup_limits] WARNING : Didn't find tree {} in input file, but did find {}. Changing the tree name, but try to fix this.".format(tree_name, backup_tree_name)
						tree_name = backup_tree_name
					else:
						print "[setup_limits] ERROR : Didn't find tree {} in input file, nor {}. Quitting!".format(tree_name, backup_tree_name)
						sys.exit(1)
				# Check that the "NEvents" histogram is present
				h_NEvents = f.Get("NEvents")
				if not h_NEvents:
					if "data" in sample:
						print "[setup_limits] ERROR : NEvents histogram in not in this file! It is probably corrupt. This is data, so this problem is fatal."
						sys.exit(1)
					else:
						print "[setup_limits] WARNING : NEvents histogram in not in this file! It is probably corrupt. This is MC, so I am skipping the file. But, you probably want to remove from the input list."
						sample_files[sample].remove(filename)
			if args.output_folder:
				output_filename = "{}/ddt_ntuple_{}.root".format(args.output_folder, sample)
			else:
				output_filename = os.path.expandvars("$HOME/DAZSLE/data/DDT/tmp/ddt_ntuple_{}.root".format(sample))
			ddt_ntupler = DDTNtupler(sample, tree_name=tree_name, output_filename=output_filename)
			for filename in sample_files[sample]:
				ddt_ntupler.add_file(filename)
			ddt_ntupler.start()
			ddt_ntupler.run(max_nevents=100)
			ddt_ntupler.finish()

	# Setup pipeline jobs on HTCondor
	if args.condor_run:
		import time
		hadd_scripts = []
		for sample in samples:
			start_directory = os.getcwd()
			job_tag = "job_{}_{}".format(sample, int(floor(time.time())))
			submission_directory = os.path.expandvars("$HOME/DAZSLE/data/DDT/condor/{}".format(job_tag))
			os.system("mkdir -pv {}".format(submission_directory))
			os.chdir(submission_directory)

			files_per_job = 1
			if "QCD_HT500to700" in sample:
				files_per_job = 5
			elif "QCD_HT700to1000" in sample:
				files_per_job = 5
			elif "QCD" in sample:
				files_per_job = 10
			n_jobs = int(math.ceil(1. * len(sample_files[sample]) / files_per_job))

			job_script_path = "{}/run_csubjob.sh".format(submission_directory)
			job_script = open(job_script_path, 'w')
			job_script.write("#!/bin/bash\n")
			job_script.write("which python\n")
			job_script.write("python --version\n")
			job_script.write("ll /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/\n")
			job_script.write("input_files=( " + " ".join(sample_files[sample]) + " )\n")
			job_script.write("files_per_job=" + str(files_per_job) + "\n")
			job_script.write("first_file_index=$(($1*$files_per_job))\n")
			job_script.write("max_file_index=$((${#input_files[@]}-1))\n")
			job_script.write("if [ $(($first_file_index+$files_per_job-1)) -gt $max_file_index ]; then\n")
			job_script.write("	files_per_job=$(($max_file_index-$first_file_index+1))\n")
			job_script.write("fi\n")
			job_script.write("declare -a this_input_files=(${input_files[@]:$first_file_index:$files_per_job})\n")
			job_script.write("function join { local IFS=\"$1\"; shift; echo \"$*\"; }\n")
			job_script.write("this_input_files_string=\"$(join , ${this_input_files[@]})\"\n")
			job_script.write("echo \"Input files:\"\n")
			job_script.write("echo $this_input_files_string\n")

			job_command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/ddt_ntupler.py --files $this_input_files_string --label {}_csubjob$1 --output_folder . --run ".format(sample)
			job_command += " 2>&1\n"
			job_script.write(job_command)

			# Check if the output file exists
			job_script.write("for f in ./ddt_ntuple*csubjob$1*root; do\n")
			job_script.write("\t[ -e \"$f\" ] && echo \"1\" > jobstatus_csubjob$1.txt || echo \"0\" > jobstatus_csubjob$1.txt \n")
			job_script.write("\tbreak\n")
			job_script.write("done\n")
			job_script.close()
			submission_command = "csub {} --cmssw --no_retar -n {}".format(job_script_path, n_jobs)
			print submission_command

			# Save csub command for resubmission attempts
			submission_script_path = "{}/csub_command.sh".format(submission_directory)
			submission_script = open(submission_script_path, "w")
			submission_script.write("#!/bin/bash\n")
			submission_script.write(submission_command + "\n")
			submission_script.close()

			# Submit jobs
			os.system(submission_command)

			hadd_scripts.append("{}/hadd.sh".format(submission_directory))
			hadd_script = open("{}/hadd.sh".format(submission_directory), "w")
			hadd_script.write("#!/bin/bash\n")
			hadd_script.write("for f in ./jobstatus_csubjob*.txt; do\n")
			hadd_script.write("\tif grep -Fxq \"0\" $f; then\n")
			hadd_script.write("\t\techo \"Subjob failure in $f\"\n")
			hadd_script.write("\tfi\n")
			hadd_script.write("done\n")
			hadd_script.write(os.path.expandvars("hadd $HOME/DAZSLE/data/DDT/ddt_ntuple_{}_{}.root {}/InputHistograms*csubjob*root\n".format(sample, submission_directory)))
			hadd_script.close()
			os.chdir(start_directory)
		# One hadd script to rule them all
		master_hadd_script_path = os.path.expandvars("$HOME/DAZSLE/data/DDT/condor/master_hadd")
		if not args.all:
			master_hadd_script_path += "_" + str(int(floor(time.time())))
		master_hadd_script_path += ".sh"
		master_hadd_script = open(master_hadd_script_path, "w")
		master_hadd_script.write("#!/bin/bash\n")
		for hadd_script_path in hadd_scripts:
			master_hadd_script.write("source " + hadd_script_path + "\n")
		master_hadd_script.close()		
