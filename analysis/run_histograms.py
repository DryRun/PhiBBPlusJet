# Run histogram jobs on condor
# This script handles submission, hadding, and normalization of the histograms.

import os
import sys
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.DAZSLECommon.baconbits import baconbits
import math
from math import floor, ceil
import ROOT

sys.path.append(".")
from histograms import Histograms

# If box_loose is specified, the shape will be taken from box_loose, and the normalization from box. 
def MergeHistograms(var, selection, box, supersample, use_Vmatched_histograms, wp_string, wp_string_loose=None, systematic=None):
	return_hist = None
	first = True
	original_selection = "" + selection

	for sample in config.samples[supersample]:
		input_histogram_filename = "$HOME/DAZSLE/data/histograms/InputHistograms_{}_{}.root".format(sample, args.jet_type)
		#print "Opening {}".format(input_histogram_filename)
		input_file = ROOT.TFile(input_histogram_filename, "READ")

		# Modify the selection string for systematics and matching
		selection = "" + original_selection
		if use_Vmatched_histograms:
			selection += "_matched"
		if systematic:
			selection += "_" + systematic

		# Main signal histogram
		if wp_string_loose:
			hname = "h_{}_{}_{}_{}".format(selection, box, wp_string_loose, var)
		else:
			hname = "h_{}_{}_{}_{}".format(selection, box, wp_string, var)
		this_histogram = input_file.Get(hname)
		
		# HACK: The JES and JER histograms were generated with the wrong name, switched matched and syst name
		if not this_histogram and use_Vmatched_histograms and systematic:
			selection = "{}_{}_matched".format(original_selection, systematic)
			if wp_string_loose:
				hname = "h_{}_{}_{}_{}".format(selection, box, wp_string_loose, var)
			else:
				hname = "h_{}_{}_{}_{}".format(selection, box, wp_string, var)
			this_histogram = input_file.Get(hname)
		if not this_histogram:
			print "[MergeHistograms] WARNING : Histogram {} doesn't exist in file {}".format(hname, input_file.GetPath())
			sys.exit(1)

		if wp_string_loose:
			hname_normalization = "h_{}_{}_{}_{}".format(selection, box, wp_string, var)
			normalization = input_file.Get(hname_normalization).Integral()
			if this_histogram.Integral() > 0:
				this_histogram.Scale(normalization / this_histogram.Integral())

		# Normalize MC histograms
		if supersample in config.background_names or supersample in config.simulated_signal_names:
			n_input_events = input_file.Get("h_input_nevents").Integral()
			print "\tSample input events = {}".format(n_input_events)
			print "\tSample processed events = {}".format(input_file.Get("h_processed_nevents").Integral())
			print "\tScaled nevents ({} pb-1) = {}".format(luminosity, luminosity * cross_sections[sample])
			if input_file.Get("h_processed_nevents").Integral() == 0:
				print "[setup_limits] ERROR : Processed zero events for sample {} {} {} {}. This is fatal, fix it!".format(supersample, var, selection, box)
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
			if systematic:
				hname_out = "{}_{}_{}_{}".format(supersample, box, var, systematic)
			else:
				hname_out = "{}_{}_{}".format(supersample, box, var)
			return_hist.SetName(hname_out)
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
	action_group.add_argument('--run', action="store_true", help="Run")
	action_group.add_argument('--crun', action="store_true", help="Submit condor jobs")
	action_group.add_argument('--combine_outputs', action="store_true", help="Compile results into one file for next step (buildRhalphabet). Also applies luminosity weights to MC.")

	parser.add_argument('--subjob', type=int, nargs=2, help="Run subjob @0 out of @1")
	parser.add_argument('--year', type=int, help="2016, 2017, or 2018")
	parser.add_argument('--output_folder', type=str, help="Output folder")
	parser.add_argument('--label', type=str, help="If running with --files, need to specify a label manually, in lieu of the sample names, for the output file naming.")
	parser.add_argument('--luminosity', type=float, default=35900, help="Luminosity in pb^-1")
	parser.add_argument('--jet_type', type=str, default="AK8", help="AK8 or CA15")
	parser.add_argument('--skim_inputs', action='store_true', help="Run over skim inputs")
	parser.add_argument('--do_optimization', action='store_true', help="For merging: merge optimization histograms")
	parser.add_argument('--do_ps_weights', action='store_true', help="Fill PS weight hists")
	args = parser.parse_args()

	if not args.year in [2016, 2017, 2018]:
		print "[run_histograms] ERROR : --year must be 2016, 2017, or 2018"
		sys.exit(1)

	if args.subjob:
		if not args.samples:
			print "[run_histograms] ERROR : --subjob is only valid for running over a single sample"
			sys.exit(1)
		elif len(args.samples.split(",")) > 1:
			print "[run_histograms] ERROR : --subjob is only valid for running over a single sample"
			sys.exit(1)

	# Determine input files to run over

	if args.run or args.crun:
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
					sample_files[sample] = baconbits[args.year][sample]
		elif args.supersamples:
			supersamples = args.supersamples.split(",")
			samples = [] 
			for supersample in supersamples:
				samples.extend(config.samples[supersample])
				for sample in config.samples[supersample]:
					sample_files[sample] = baconbits[args.year][sample]
		elif args.samples:
			samples = args.samples.split(",")
			for sample in samples:
				sample_files[sample] = baconbits[args.year][sample]
		#elif args.files:
		#	files = args.files.split(",")
		#	for filename in files:
		#		if args.label:
		#			this_sample = args.label
		#		else:
		#			print "[run_histograms] ERROR : When running with --files option, you must specify a label for the output!"
		#			sys.exit(1)
		#		if not this_sample in sample_files:
		#			sample_files[this_sample] = []
		#		sample_files[this_sample].append(filename)
		#	samples = sample_files.keys()
		print "List of input samples: ",
		print samples
		print "List of samples and files: ",
		print sample_files

		if args.subjob:
			subjob_index = args.subjob[0]
			n_subjobs = args.subjob[0]
			if len(sample_files[sample]) < n_subjobs:
				print "[run_histograms] ERROR : Number of subjobs ({}) exceeds number of files ({})!".format(len(sample_files[sample]), n_subjobs)
				sys.exit(1)
			files_per_job = len(sample_files[sample]) / n_subjobs + 1
			i0 = files_per_job * subjob_index
			i1 = min(files_per_job * (subjob_index + 1) - 1, len(sample_files[sample]) - 1)
			files_to_run = sample_file[sample][i0:i1]
			print "Running over files {}-{}:".format(i0, i1)
			print files_to_run
		else:
			files_to_run = sample_files[sample]

	if args.run:
		for sample in samples:
			print "\n *** Running sample {}".format(sample)

			if "Sbb" in sample or args.skim_inputs or "ZPrime" in sample:
				tree_name = "Events"
			else:
				tree_name = "otree"

			# Sanity check: make sure tree exists in file
			for filename in files_to_run:
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
				
			limit_histogrammer = Histograms(sample, tree_name=tree_name, jet_type=args.jet_type)
			output_file_basename ="InputHistograms_{}_{}_{}.root".format(sample, args.jet_type, args.year) 
			if args.output_folder:
				limit_histogrammer.set_output_path("{}/{}".format(args.output_folder, output_file_basename))
			else:
				limit_histogrammer.set_output_path("/uscms/home/dryu/DAZSLE/data/LimitSetting/{}".format(output_file_basename))
			for filename in files_to_run:
				print "Input file {}".format(filename)
				limit_histogrammer.add_file(filename)
			#limit_histogrammer.set_jet_type(args.jet_type)
			if "JetHT" in sample or "SingleMu" in sample:
				limit_histogrammer.set_data_source("data")
			else:
				limit_histogrammer.set_data_source("simulation")
			if "ps10" in sample:
				limit_histogrammer.set_prescale(10)
			if args.do_ps_weights:
				limit_histogrammer.do_ps_weights(True)
			limit_histogrammer.start()
			limit_histogrammer.run()
			limit_histogrammer.finish()

	elif args.crun:
		import time
		hadd_scripts = []
		for sample in samples:
			start_directory = os.getcwd()
			job_tag = "job_{}_{}_{}".format(sample, args.jet_type, int(floor(time.time())))
			submission_directory = os.path.expandvars("$HOME/DAZSLE/data/histograms/condor/{}".format(job_tag))
			os.system("mkdir -pv {}".format(submission_directory))
			os.chdir(submission_directory)

			files_per_job = 1
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
			elif "WJetsToQQ_HT" in sample or "ZJetsToQQ_HT" in sample:
				files_per_job = 10
			n_jobs = int(ceil(1. * len(sample_files[sample]) / files_per_job))

			job_script_path = "{}/run_csubjob.sh".format(submission_directory)
			job_script = open(job_script_path, 'w')
			job_script.write("#!/bin/bash\n")
			job_script.write("which python\n")
			job_script.write("python --version\n")
			job_script.write("ls -lrth /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/\n")
			#job_script.write("input_files=( " + " ".join(sample_files[sample]) + " )\n")
			#job_script.write("files_per_job=" + str(files_per_job) + "\n")
			#job_script.write("first_file_index=$(($1*$files_per_job))\n")
			#job_script.write("max_file_index=$((${#input_files[@]}-1))\n")
			#job_script.write("if [ $(($first_file_index+$files_per_job-1)) -gt $max_file_index ]; then\n")
			#job_script.write("	files_per_job=$(($max_file_index-$first_file_index+1))\n")
			#job_script.write("fi\n")
			#job_script.write("declare -a this_input_files=(${input_files[@]:$first_file_index:$files_per_job})\n")
			#job_script.write("function join { local IFS=\"$1\"; shift; echo \"$*\"; }\n")
			#job_script.write("this_input_files_string=\"$(join , ${this_input_files[@]})\"\n")
			#job_script.write("echo \"Input files:\"\n")
			#job_script.write("echo $this_input_files_string\n")
			job_command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/run_histograms.py --year {} --jet_type {} --files $this_input_files_string --label {}_csubjob$1 --output_folder . --subjob $1 {}".format(args.year, args.jet_type, sample, n_jobs)
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

	elif args.combine_outputs:
		luminosity = args.luminosity
		from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections

		selections = ["SR", "Preselection"] # N2CR

		if args.do_optimization:
			n2ddt_wps = [0.05, 0.15, 0.26]
			dbtag_wps= [0.7, 0.8, 0.9]
		else:
			n2ddt_wps = [0.26]
			dbtag_wps= [0.9]

		weight_systematics = {
			"SR":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"Preselection":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
		}
		jet_systematics = ["JESUp", "JESDown", "JERUp", "JERDown"]
		#for systematic in weight_systematics["SR"] + jet_systematics:
		#		selections.append("SR_{}".format(systematic))
		vars = ["pt_vs_msd", "pfmet", "dcsv", "n2ddt", "n2", "pt", "eta", "rho"]

		boxes = []
		boxes.append("passn2_passdbtag")
		boxes.append("failn2_passdbtag")
		boxes.append("passn2_faildbtag")
		boxes.append("failn2_faildbtag")

		for dbtag_wp in dbtag_wps:
			for n2wp_dbpass in n2ddt_wps:
				for n2wp_dbfail in n2ddt_wps:
					for selection in selections:
						if selection == "muCR":
							supersamples = ["data_obs", "data_singlemu", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq", "zll", "wlnu"]
						else:
							supersamples = ["data_obs", "data_singlemu", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq"]
						supersamples.extend(config.simulated_signal_names)

						if args.do_optimization:
							opt_wp_string = "_dbtag{}_n2wpdbpass{}_n2wpdbfail{}".format(dbtag_wp, n2wp_dbpass, n2wp_dbfail)
							output_file = ROOT.TFile("$HOME/PhiBB2017/data/histograms/optimization/histograms_{}_{}{}.root".format(selection, args.jet_type, opt_wp_string), "RECREATE")
						else:
							output_file = ROOT.TFile("$HOME/PhiBB2017/data/histograms/histograms_{}_{}{}.root".format(selection, args.jet_type), "RECREATE")

						for box in boxes:
							for supersample in supersamples:
								for var in vars:
									use_Vmatched_histograms = ((supersample in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]) or ("Sbb" in supersample) or ("ZPrime" in supersample)) and (selection == "SR")
									#use_loose_template = (supersample in ["wqq", "zqq"]) # Use looser DCSV cut for pass shape, to improve statistics

									if "passdbtag" in box:
										wp_string = "n2wp{}_dbtag{}".format(n2wp_dbpass, dbtag_wp)
									elif "faildbtag" in box:
										wp_string = "n2wp{}_dbtag{}".format(n2wp_dbfail, dbtag_wp)

									# When to use loose templates?
									# - Resonant EW backgrounds (W/Z/H)
									# - Pass-pass box only. pass-fail and fail-pass should have enough stats.
									if supersample in ["wqq", "zqq", "hhqq125"] and box == "passn2_passdbtag":
										wp_string_loose = "n2wp{}_dbtag0.7".format(n2wp_dbpass)
									else:
										wp_string_loose = None

									merged_histogram = MergeHistograms(var=var, selection=selection, box=box, supersample=supersample, use_Vmatched_histograms=use_Vmatched_histograms, wp_string=wp_string, wp_string_loose=wp_string_loose)
									output_file.cd()

									# For muCR, project to 1D
									if "muCR" in selection:
										old_name = merged_histogram.GetName()
										merged_histogram.RebinY(merged_histogram.GetNbinsY())
										merged_histogram.SetName(old_name)

									output_file.cd()
									merged_histogram.Write()

									# Systematics
									if selection == "SR" and var == "pt_vs_msd" and not "data" in supersample:
										for systematic in weight_systematics["SR"] + jet_systematics:
											merged_histogram_syst = MergeHistograms(var=var, selection=selection, box=box, supersample=supersample, use_Vmatched_histograms=use_Vmatched_histograms, wp_string=wp_string, wp_string_loose=wp_string_loose, systematic=systematic)
											output_file.cd()
											merged_histogram_syst.Write()
						output_file.Close()
