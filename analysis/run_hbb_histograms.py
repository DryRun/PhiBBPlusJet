# Run histogram jobs on condor
# This script handles submission, hadding, and normalization of the histograms.

import os
import sys
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
import math
from math import floor, ceil
import ROOT

def MergeHistograms(var, selection, box, supersample, use_Vmatched_histograms, use_loose_template):
	return_hist = None
	first = True

	for sample in config.samples[supersample]:
		input_histogram_filename = "$HOME/PhiBB2017/data/histograms/InputHistograms_{}_{}.root".format(sample, args.jet_type)
		print "Opening {}".format(input_histogram_filename)
		input_file = ROOT.TFile(input_histogram_filename, "READ")

		# Main signal histogram
		hname = "h_{}_{}_{}".format(selection, box, var)

		if use_Vmatched_histograms:
			hname.replace(selection, selection + "_matched")
		if use_loose_template:
			# Save normalization before replacing name
			normalization = input_file.Get(hname).Integral()
			hname.replace("pass1", "pass1loose") # Only valid for pass1 (=pass N2, pass dcsv)
		this_histogram = input_file.Get(hname)
		if not this_histogram:
			print "[MergeHistograms] WARNING : Histogram {} doesn't exist in file {}".format(hname, input_file.GetPath())
		if use_loose_template:
			this_histogram.Scale(normalization / this_histogram.Integral())

		# Normalize MC histograms
		if supersample in config.background_names or supersample in config.simulated_signal_names:
			n_input_events = input_file.Get("h_input_nevents").Integral()
			print "\tSample input events = {}".format(n_input_events)
			print "\tSample processed events = {}".format(input_file.Get("h_processed_nevents").Integral())
			print "\tScaled nevents ({} pb-1) = {}".format(luminosity, luminosity * cross_sections[sample])
			if input_file.Get("h_processed_nevents").Integral() == 0:
				print "[setup_limits] ERROR : Processed zero events for sample {}. This is fatal, fix it!"
				sys.exit(1)

			# Normalize histograms
			if "Spin0" in sample or "Sbb" in sample or "ZPrime" in sample:
				# Normalize to visible cross section of 1 pb
				#print "\tNormalizing signal sample {} to visible cross section of 1 pb".format(sample)
				#if this_pass_histogram.GetEntries():
				#	lumi_sf = luminosity / this_pass_histogram.GetEntries()
				#	print "\tLuminosity scale factor = {}".format(lumi_sf)
				#else:
				#	print "[setup_limits] WARNING : Found zero input events for sample {}.".format(sample)
				#	lumi_sf = 0.
				
				# Actually, maybe it's easier to normalize to xs*BR*A(filter)=1pb
				print "\tNormalizing signal sample {} to xs*BR*A=1pb"
				if n_input_events > 0:
					print sample
					lumi_sf = luminosity * cross_sections[sample] / n_input_events
					print "\tLuminosity scale factor = {}".format(lumi_sf)
				else:
					print "[setup_limits] WARNING : Found zero input events for sample {}. Something went wrong in an earlier step. I'll continue, but you need to fix this.".format(sample)
					lumi_sf = 0.
			else:
				if n_input_events > 0:
					print sample
					lumi_sf = luminosity * cross_sections[sample] / n_input_events
					print "\tLuminosity scale factor = {}".format(lumi_sf)
				else:
					print "[setup_limits] WARNING : Found zero input events for sample {}. Something went wrong in an earlier step. I'll continue, but you need to fix this.".format(sample)
					lumi_sf = 0.
			this_histogram.Scale(lumi_sf)

		if first:
			return_hist = this_histogram.Clone()
			return_hist.SetDirectory(0)
			return_hist.SetName("{}_{}_{}".format(supersample, box, var))
			first = False
		else:
			return_hist.Add(this_histogram)
		input_file.Close()
	# End loop over samples in supersample
	return return_hist


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Submit histogram jobs on condor")
	input_group = parser.add_mutually_exclusive_group() 
	input_group.add_argument('--all', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--all_lxplus', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--all_cmslpc', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--supersamples', type=str, help="Supersample name(s), comma separated. Must correspond to something in analysis_configuration.(background_names, signal_names, or data_names).")
	input_group.add_argument('--samples', type=str, help="Sample name(s), comma separated. Must be a key in analysis_configuration.skims.")

	action_group = parser.add_mutually_exclusive_group() 
	action_group.add_argument('--crun', action="store_true", help="Run on condor")
	action_group.add_argument('--combine_outputs', action="store_true", help="Compile results into one file for next step (buildRhalphabet). Also applies luminosity weights to MC.")
	parser.add_argument('--output_folder', type=str, help="Output folder")
	parser.add_argument('--label', type=str, help="If running with --files, need to specify a label manually, in lieu of the sample names, for the output file naming.")
	parser.add_argument('--luminosity', type=float, default=35900, help="Luminosity in pb^-1")
	parser.add_argument('--jet_type', type=str, default="AK8", help="AK8 or CA15")
	parser.add_argument('--skim_inputs', action='store_true', help="Run over skim inputs")
	args = parser.parse_args()

	if args.crun:
		# Make a list of input samples and files
		samples = []
		sample_files = {} # Dictionary is sample : [list of files in sample]
		if args.all or args.all_lxplus or args.all_cmslpc:
			if args.all:
				supersamples = config.supersamples
			elif args.all_lxplus:
				# lxplus: JetHT, SingleMuon, QCD, signal
				supersamples = ["data_obs", "data_singlemu", "qcd"] #"data_obs_ps10", "data_singlemu_ps10"
				args.skim_inputs = True
			elif args.all_cmslpc:
				supersamples = ["stqq", "tqq", "wqq", "zqq", "zll", "wlnu", "vvqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]
				supersamples.extend([x for x in config.simulated_signal_names])
				args.skim_inputs = False
			samples = [] 
			for supersample in supersamples:
				samples.extend(config.samples[supersample])
				for sample in config.samples[supersample]:
					if args.skim_inputs:
						sample_files[sample] = config.skims[sample]
					else:
						sample_files[sample] = config.sklims[sample]
		elif args.supersamples:
			supersamples = args.supersamples.split(",")
			samples = [] 
			for supersample in supersamples:
				samples.extend(config.samples[supersample])
				for sample in config.samples[supersample]:
					if args.skim_inputs:
						sample_files[sample] = config.skims[sample]
					else:
						sample_files[sample] = config.sklims[sample]
		elif args.samples:
			samples = args.samples.split(",")
			for sample in samples:
				if args.skim_inputs:
					sample_files[sample] = config.skims[sample]
				else:
					sample_files[sample] = config.sklims[sample]
		elif args.files:
			files = args.files.split(",")
			for filename in files:
				if args.label:
					this_sample = args.label
				else:
					print "[run_histograms] ERROR : When running with --files option, you must specify a label for the output!"
					sys.exit(1)
				if not this_sample in sample_files:
					sample_files[this_sample] = []
				sample_files[this_sample].append(filename)
			samples = sample_files.keys()
		print "List of input samples: ",
		print samples
		print "List of samples and files: ",
		print sample_files

		import time
		hadd_scripts = []
		for sample in samples:
			start_directory = os.getcwd()
			job_tag = "job_{}_{}_{}".format(sample, args.jet_type, int(floor(time.time())))
			submission_directory = os.path.expandvars("$HOME/PhiBB2017/data/histograms/condor/{}".format(job_tag))
			os.system("mkdir -pv {}".format(submission_directory))
			os.chdir(submission_directory)

			files_per_job = 1
			if args.skim_inputs:
				if "JetHTRun2016" in sample:
					if "ps10" in sample:
						files_per_job = 100
					else:
						files_per_job = 20
				elif "SingleMuRun2016" in sample:
					if "ps10" in sample:
						files_per_job = 100
					else:
						files_per_job = 10
				elif "QCD_HT500to700" in sample:
					files_per_job = 3
				elif "QCD_HT700to1000" in sample:
					files_per_job = 3
				elif "QCD_HT1000to1500" in sample:
					files_per_job = 1
				elif "QCD" in sample:
					files_per_job = 5
				elif "Spin0" in sample or "Sbb" in sample or "ZPrime" in sample:
					files_per_job = 3
			n_jobs = int(ceil(1. * len(sample_files[sample]) / files_per_job))

			job_script_path = "{}/run_csubjob.sh".format(submission_directory)
			job_script = open(job_script_path, 'w')
			job_script.write("#!/bin/bash\n")
			job_script.write("which python\n")
			job_script.write("python --version\n")
			job_script.write("ls -lrth /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/\n")
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
			job_command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/histograms.py --jet_type {} --files $this_input_files_string --label {}_csubjob$1 --output_folder .".format(args.jet_type, sample)
			if args.skim_inputs or args.all_lxplus:
				job_command += " --skim_inputs "

			job_command += " 2>&1\n"
			job_script.write(job_command)

			# Check if the output file exists
			job_script.write("for f in ./InputHistograms*_csubjob$1*root; do\n")
			job_script.write("\t[ -e \"$f\" ] && echo \"1\" > jobstatus_csubjob$1.txt || echo \"0\" > jobstatus_csubjob$1.txt \n")
			job_script.write("\tbreak\n")
			job_script.write("done\n")

			#job_script.write("if ls InputHistograms*{}_csubjob$1*root 1> /dev/null 2&>1; then\n".format(sample))
			#job_script.write("\techo \"1\" > jobstatus_csubjob$1.txt\n")
			#job_script.write("else\n")
			#job_script.write("\techo\"0\" > jobstatus_csubjob$1.txt\n")
			#job_script.write("fi\n")

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
			hadd_script.write(os.path.expandvars("hadd $HOME/PhiBB2017/data/histograms/InputHistograms_{}_{}.root {}/InputHistograms*csubjob*root\n".format(sample, args.jet_type, submission_directory)))
			hadd_script.close()
			os.chdir(start_directory)
		# One hadd script to rule them all
		master_hadd_script_path = os.path.expandvars("$HOME/PhiBB2017/data/histograms/condor/master_hadd_{}".format(args.jet_type))
		if not args.all:
			master_hadd_script_path += "_" + str(int(floor(time.time())))
		master_hadd_script_path += ".sh"
		master_hadd_script = open(master_hadd_script_path, "w")
		master_hadd_script.write("#!/bin/bash\n")
		for hadd_script_path in hadd_scripts:
			master_hadd_script.write("source " + hadd_script_path + "\n")
		master_hadd_script.close()


	if args.combine_outputs:
		luminosity = args.luminosity
		from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections

		selections = ["SR", "Preselection"] # N2CR
		boxes = ["all", "pass1", "pass2", "fail1", "fail2"]
		weight_systematics = {
			"SR":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"Preselection":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
		}
		jet_systematics = ["JESUp", "JESDown", "JERUp", "JERDown"]
		#for systematic in weight_systematics["SR"] + jet_systematics:
		#		selections.append("SR_{}".format(systematic))
		vars = ["pt_vs_msd", "pfmet", "dcsv", "n2ddt", "n2", "pt", "eta", "rho"]

		for selection in selections:
			if selection == "muCR":
				supersamples = ["data_obs", "data_singlemu", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq", "zll", "wlnu"]
			else:
				supersamples = ["data_obs", "data_singlemu", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq"]
			supersamples.extend(config.simulated_signal_names)

			output_file = ROOT.TFile("$HOME/PhiBB2017/data/histograms/histograms_{}_{}.root".format(selection, args.jet_type), "RECREATE")

			for box in boxes:
				for supersample in supersamples:
					for var in vars:
						use_Vmatched_histograms = (supersample in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]) or ("Sbb" in supersample) or ("ZPrime" in supersample)
						use_loose_template = (supersample in ["wqq", "zqq"]) # Use looser DCSV cut for pass shape, to improve statistics
						merged_histogram = MergeHistograms(var=var, selection=selection, box=box, supersample=supersample, use_Vmatched_histograms=use_Vmatched_histograms, use_loose_template=use_loose_template)
						output_file.cd()

						# For muCR, project to 1D
						if "muCR" in selection:
							old_name = merged_histogram.GetName()
							merged_histogram.RebinY(merged_histogram.GetNbinsY())
							merged_histogram.SetName(old_name)

						output_file.cd()
						merged_histogram.Write()

						# Systematics
						if selection == "SR" and var == "pt_vs_msd":
							for systematic in weight_systematics["SR"] + jet_systematics:
								merged_histogram_syst = MergeHistograms(var=var, selection="SR_{}".format(systematic), box=box, supersample=supersample, use_Vmatched_histograms=use_Vmatched_histograms, use_loose_template=use_loose_template)
								output_file.cd()
								merged_histogram_syst.Write()
			output_file.Close()
