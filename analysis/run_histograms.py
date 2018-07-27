# Run histogram jobs on condor
# This script handles submission, hadding, and normalization of the histograms.

import os
import sys
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
import math
from math import floor, ceil

# Slice up the pt vs msd histogram
def SliceSignalHistogram(hist2d, pt_bins):
	# Make sure the requested pt bins correspond to boundaries
	histogram_pt_bins = hist2d.GetYaxis().GetBins()
	for pt_boundary in pt_bins:
		if not pt_boundary in histogram_pt_bins:
			print "[run_histograms::SliceSignalHistogram] ERROR : Bin boundary {} does not correspond to a histogram bin boundary."
			print histogram_pt_bins
			sys.exit(1)
	histogram_slices = {}
	for islice in xrange(len(pt_bins) - 1):
		ptmin = pt_bins[islice]
		ptmax = pt_bins[islice+1]
		binmin = 1e10
		binmax = -1e10
		for bin in xrange(1, hist2d.GetNbinsY() + 1):
			low_edge = hist2d.GetYaxis().GetBinLowEdge(bin)
			up_edge = hist2d.GetYaxis().GetBinUpEdge(bin)
			# Is this bin inside this pt slice (+epsilon)?
			if ptmin - 1.e-5 < low_edge and high_edge < ptmax + 1.e-5:
				if bin < binmin:
					binmin = bin
				if bin > binmax:
					binmax = bin
		histogram_slices[islice] = hist2d.ProjectionX(hist2d.GetName() + "_" + str(islice), binmin, binmax)
	return histogram_slices


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
				supersamples = ["data_obs", "data_singlemu", "data_obs_ps10", "data_singlemu_ps10", "qcd"]
				args.skim_inputs = True
			elif args.all_cmslpc:
				supersamples = ["stqq", "tqq", "wqq", "zqq", "zll", "wlnu", "vvqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]
				supersamples.extend([x for x in config.simulated_signal_names if not "ZPrime" in x])
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
			submission_directory = os.path.expandvars("$HOME/PhiBB2017/data/Histograms/condor/{}".format(job_tag))
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
			hadd_script.write(os.path.expandvars("hadd $HOME/PhiBB2017/data/Histograms/InputHistograms_{}_{}.root {}/InputHistograms*csubjob*root\n".format(sample, args.jet_type, submission_directory)))
			hadd_script.close()
			os.chdir(start_directory)
		# One hadd script to rule them all
		master_hadd_script_path = os.path.expandvars("$HOME/PhiBB2017/data/Histograms/condor/master_hadd_{}".format(args.jet_type))
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
		extra_vars = ["pfmet", "dcsv", "n2ddt", "n2", "pt", "eta", "rho", "n2ddt_vs_msd_vs_pt", "dcsv_vs_msd_vs_pt"]
		output_file = ROOT.TFile("$HOME/PhiBBPlusJet/data/Histograms/histograms_{}_{}.root".format(selection, args.jet_type), "RECREATE")

		box_histograms = {}

		if selection == "muCR":
			supersamples = ["data_obs", "data_singlemu", "data_obs_ps10", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq", "zll", "wlnu"]
		else:
			supersamples = ["data_obs", "data_singlemu", "data_obs_ps10", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq"]
		supersamples.extend(config.simulated_signal_names)
		for supersample in supersamples:
			first = True
			box_histograms[supersample] = {}
			use_Vmatched_histograms = (supersample in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]) or ("Sbb" in supersample) or ("ZPrime" in supersample)
			use_loose_template = False # (supersample in ["wqq", "zqq"]) # Use looser DCSV cut for pass shape, to improve statistics
			if use_loose_template:
				pass_histograms_syst[supersample + "_normalization"] = {}

			for sample in config.samples[supersample]:
				input_histogram_filename = "$HOME/DAZSLE/data/Histograms/InputHistograms_{}_{}.root".format(sample, args.jet_type)
				print "Opening {}".format(input_histogram_filename)
				input_file = ROOT.TFile(input_histogram_filename, "READ")
				if selection in selection_tau21s:
					pass_histogram_name = "h_{}_tau21ddt{}_{}_pass_dcsv{}".format(selection_prefix, selection_tau21s[selection], args.jet_type, selection_dcsvs[selection])
					fail_histogram_name = "h_{}_tau21ddt{}_{}_fail_dcsv{}".format(selection_prefix, selection_tau21s[selection], args.jet_type, selection_dcsvs[selection])
					nevents_histogram_name = "h_{}_tau21ddt{}_{}_pass_nevents".format(selection_prefix, selection_tau21s[selection], args.jet_type)
				elif selection == "N2SR":
					pass_histogram_name = "h_{}_{}_pass".format(selection_prefix, args.jet_type)
					fail_histogram_name = "h_{}_{}_fail".format(selection_prefix, args.jet_type)
					nevents_histogram_name = "h_{}_{}_pass_nevents".format(selection_prefix, args.jet_type)
					if use_Vmatched_histograms:
						pass_histogram_name += "_matched"
						fail_histogram_name += "_matched"
				else:
					pass_histogram_name = "h_{}_{}_pass_dcsv{}".format(selection_prefix, args.jet_type, params[args.jet_type]["DCSV"])
					fail_histogram_name = "h_{}_{}_fail_dcsv{}".format(selection_prefix, args.jet_type, params[args.jet_type]["DCSV"])
					nevents_histogram_name = "h_{}_{}_pass_nevents".format(selection_prefix, args.jet_type)
					if use_Vmatched_histograms:
						pass_histogram_name += "_matched"
						fail_histogram_name += "_matched"
				if use_loose_template:
					pass_histogram_name_normalization = pass_histogram_name
					if selection == "N2SR":
						pass_histogram_name.replace("N2SR", "N2SR_loose")
						fail_histogram_name.replace("N2SR", "N2SR_loose")
					else:
						pass_histogram_name = pass_histogram_name.replace("dcsv{}".format(params[args.jet_type]["DCSV"]), "dcsv{}".format(params[args.jet_type]["DCSV_LOOSE"]))
				this_pass_histogram = input_file.Get(pass_histogram_name)
				if not this_pass_histogram:
					print "[run_histograms] ERROR : Couldn't find histogram {} in file {}".format(pass_histogram_name, input_file.GetPath())
				this_fail_histogram = input_file.Get(fail_histogram_name)
				if use_loose_template:
					this_pass_histogram_normalization = input_file.Get(pass_histogram_name_normalization)
				this_pass_histogram_syst = {}
				this_fail_histogram_syst = {}
				for systematic in systematics[selection]:
					if selection in selection_tau21s:
						pass_histogram_name = "h_{}_tau21ddt{}_{}_pass_{}_dcsv{}".format(selection_prefix, selection_tau21s[selection], args.jet_type, systematic, selection_dcsvs[selection])
						fail_histogram_name = "h_{}_tau21ddt{}_{}_fail_{}_dcsv{}".format(selection_prefix, selection_tau21s[selection], args.jet_type, systematic, selection_dcsvs[selection])
					else:
						pass_histogram_name = "h_{}_{}_pass_{}".format(selection, args.jet_type, systematic)
						fail_histogram_name = "h_{}_{}_fail_{}".format(selection, args.jet_type, systematic)
						if supersample in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"] or ("Sbb" in supersample) or ("ZPrime" in supersample):
							pass_histogram_name += "_matched"
							fail_histogram_name += "_matched"
					if use_loose_template:
						pass_histogram_name_normalization = pass_histogram_name
						if selection == "N2SR":
							pass_histogram_name.replace("N2SR", "N2SR_loose")
							fail_histogram_name.replace("N2SR", "N2SR_loose")
						else:
							pass_histogram_name = pass_histogram_name.replace("pass", "passloose")
					this_pass_histogram_syst[systematic] = input_file.Get(pass_histogram_name)
					this_fail_histogram_syst[systematic] = input_file.Get(fail_histogram_name)
					if use_loose_template:
						this_pass_histogram_syst[systematic + "_normalization"] = input_file.Get(pass_histogram_name_normalization)
				if supersample in config.background_names or supersample in config.simulated_signal_names:
					n_input_events = input_file.Get("h_input_nevents").Integral()
					print "\tSample input events = {}".format(n_input_events)
					print "\tSample processed events = {}".format(input_file.Get("h_processed_nevents").Integral())
					print "\tSample pass events = {}".format(input_file.Get(nevents_histogram_name).Integral())
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
					this_pass_histogram.Scale(lumi_sf)
					this_fail_histogram.Scale(lumi_sf)
					for systematic in systematics[selection]:
						this_pass_histogram_syst[systematic].Scale(lumi_sf)
						this_fail_histogram_syst[systematic].Scale(lumi_sf)
					if use_loose_template:
						this_pass_histogram_normalization.Scale(lumi_sf)
						for systematic in systematics[selection]:
							this_pass_histogram_syst[systematic + "_normalization"].Scale(lumi_sf)

				if first:
					pass_histograms[supersample] = this_pass_histogram.Clone()
					pass_histograms[supersample].SetDirectory(0)
					pass_histograms[supersample].SetName("{}_pass".format(supersample))
					fail_histograms[supersample] = this_fail_histogram.Clone()
					fail_histograms[supersample].SetDirectory(0)
					fail_histograms[supersample].SetName("{}_fail".format(supersample))
					for systematic in systematics[selection]:
						pass_histograms_syst[supersample][systematic] = this_pass_histogram_syst[systematic].Clone()
						pass_histograms_syst[supersample][systematic].SetDirectory(0)
						pass_histograms_syst[supersample][systematic].SetName("{}_pass_{}".format(supersample, systematic))
						fail_histograms_syst[supersample][systematic] = this_fail_histogram_syst[systematic].Clone()
						fail_histograms_syst[supersample][systematic].SetDirectory(0)
						fail_histograms_syst[supersample][systematic].SetName("{}_fail_{}".format(supersample, systematic))
					if use_loose_template:
						pass_histograms[supersample + "_normalization"] = this_pass_histogram_normalization.Clone()
						pass_histograms[supersample + "_normalization"].SetDirectory(0)
						pass_histograms[supersample + "_normalization"].SetName("{}_pass_normalization".format(supersample))
						for systematic in systematics[selection]:
							pass_histograms_syst[supersample + "_normalization"][systematic] = this_pass_histogram_syst[systematic + "_normalization"].Clone()
							pass_histograms_syst[supersample + "_normalization"][systematic].SetDirectory(0)
							pass_histograms_syst[supersample + "_normalization"][systematic].SetName("{}_pass_{}_normalization".format(supersample, systematic))
					first = False
				else:
					pass_histograms[supersample].Add(this_pass_histogram)
					fail_histograms[supersample].Add(this_fail_histogram)
					for systematic in systematics[selection]:
						pass_histograms_syst[supersample][systematic].Add(this_pass_histogram_syst[systematic])
						fail_histograms_syst[supersample][systematic].Add(this_fail_histogram_syst[systematic])
					if use_loose_template:
						pass_histograms[supersample + "_normalization"].Add(this_pass_histogram_normalization)
						for systematic in systematics[selection]:
							pass_histograms_syst[supersample + "_normalization"][systematic].Add(this_pass_histogram_syst[systematic + "_normalization"])
				#if sample in cross_sections:
				#	n_input_events += input_file.Get("h_input_nevents").Integral()
				input_file.Close()
			output_file.cd()
			if use_loose_template:
				if pass_histograms[supersample].Integral():
					pass_histograms[supersample].Scale(pass_histograms[supersample + "_normalization"].Integral() / pass_histograms[supersample].Integral())
				for systematic in systematics[selection]:
					if pass_histograms_syst[supersample][systematic].Integral():
						pass_histograms_syst[supersample][systematic].Scale(pass_histograms_syst[supersample + "_normalization"][systematic].Integral() / pass_histograms_syst[supersample][systematic].Integral())

			# For muCR, project to 1D
			if "muCR" in selection:
				old_name = pass_histograms[supersample].GetName()
				pass_histograms[supersample].RebinY(pass_histograms[supersample].GetNbinsY())
				pass_histograms[supersample].SetName(old_name)
				old_name = fail_histograms[supersample].GetName()
				fail_histograms[supersample].RebinY(fail_histograms[supersample].GetNbinsY())
				fail_histograms[supersample].SetName(old_name)
				for systematic in systematics[selection]:
					old_name = pass_histograms_syst[supersample][systematic].GetName()
					pass_histograms_syst[supersample][systematic].RebinY(pass_histograms_syst[supersample][systematic].GetNbinsY())
					pass_histograms_syst[supersample][systematic].SetName(old_name)
					old_name = fail_histograms_syst[supersample][systematic].GetName()
					fail_histograms_syst[supersample][systematic].RebinY(fail_histograms_syst[supersample][systematic].GetNbinsY())
					fail_histograms_syst[supersample][systematic].SetName(old_name)
				if use_loose_template:
					old_name = pass_histograms[supersample + "_normalization"].GetName()
					pass_histograms[supersample + "_normalization"].RebinY(pass_histograms[supersample + "_normalization"].GetNbinsY())
					pass_histograms[supersample + "_normalization"].SetName(old_name)
					for systematic in systematics[selection]:
						old_name = pass_histograms_syst[supersample + "_normalization"][systematic].GetName()
						pass_histograms_syst[supersample + "_normalization"][systematic].RebinY(pass_histograms_syst[supersample + "_normalization"][systematic].GetNbinsY())
						pass_histograms_syst[supersample + "_normalization"][systematic].SetName(old_name)

			pass_histograms[supersample].Write()
			fail_histograms[supersample].Write()
			for systematic in systematics[selection]:
				pass_histograms_syst[supersample][systematic].Write()
				fail_histograms_syst[supersample][systematic].Write()
			if use_loose_template:
				pass_histograms[supersample + "_normalization"].Write()
				for systematic in systematics[selection]:
					pass_histograms_syst[supersample + "_normalization"][systematic].Write()

			# Now do the extra histograms for plots
			if "SR" in selection or selection in ["Preselection", "muCR", "N2CR"]:
				extra_histograms = {}
				extra_histograms_pass = {}
				extra_histograms_fail = {}
				for var in extra_vars:
					first = True
					for sample in config.samples[supersample]:
						input_histogram_filename = "$HOME/DAZSLE/data/Histograms/InputHistograms_{}_{}.root".format(sample, args.jet_type)
						print "Opening {}".format(input_histogram_filename)
						input_file = TFile(input_histogram_filename, "READ")
						this_histogram = input_file.Get("h_{}_{}_{}".format(selection, args.jet_type, var))
						if not this_histogram:
							print "ERROR : Couldn't find histogram {} in file {}".format("h_{}_{}_{}".format(selection, args.jet_type, var), input_histogram_filename)
						this_histogram_pass = input_file.Get("h_{}_{}_pass_{}".format(selection, args.jet_type, var))
						this_histogram_fail = input_file.Get("h_{}_{}_fail_{}".format(selection, args.jet_type, var))
						# Normalize histograms
						if supersample in config.background_names or supersample in config.simulated_signal_names:
							n_input_events = input_file.Get("h_input_nevents").Integral()
							if "Spin0" in sample or "Sbb" in sample or "ZPrime" in sample:
								# Normalize to visible cross section of 1 pb
								#print "\tNormalizing signal sample {} to visible cross section of 1 pb".format(sample)
								#pass_events = input_file.Get("h_SR_{}_pass".format(args.jet_type)).Integral()
								#if pass_events:
								#	lumi_sf = luminosity / pass_events
								#	print "\tLuminosity scale factor = {}".format(lumi_sf)
								#else:
								#	print "[setup_limits] WARNING : Found zero input events for sample {}.".format(sample)
								#	lumi_sf = 0.
								
								# Actually, maybe it's easier to normalize to xs*BR*A(filter)=1pb
								print "\tNormalizing signal sample {} to xs*BR*A=1pb".format(supersample)
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
									print "\t{} luminosity scale factor = {}*{}/{}={}".format(sample, luminosity, cross_sections[sample], n_input_events, lumi_sf)
								else:
									print "[setup_limits] WARNING : Found zero input events for sample {}. Something went wrong in an earlier step. I'll continue, but you need to fix this.".format(sample)
									lumi_sf = 0.
							this_histogram.Scale(lumi_sf)
							if this_histogram_pass:
								this_histogram_pass.Scale(lumi_sf)
							if this_histogram_fail:
								this_histogram_fail.Scale(lumi_sf)

						# Add up
						if first:
							first = False
							extra_histograms[var] = this_histogram.Clone()
							extra_histograms[var].SetDirectory(0)
							extra_histograms[var].SetName(supersample + "_" + var)
							if this_histogram_pass:
								extra_histograms_pass[var] = this_histogram_pass.Clone()
								extra_histograms_pass[var].SetDirectory(0)
								extra_histograms_pass[var].SetName(supersample + "_" + var + "_pass")
							if this_histogram_fail:
								extra_histograms_fail[var] = this_histogram_fail.Clone()
								extra_histograms_fail[var].SetDirectory(0)
								extra_histograms_fail[var].SetName(supersample + "_" + var + "_fail")
						else:
							extra_histograms[var].Add(this_histogram)
							if this_histogram_pass:
								extra_histograms_pass[var].Add(this_histogram_pass)
							if this_histogram_fail:
								extra_histograms_fail[var].Add(this_histogram_fail)
						input_file.Close()
					output_file.cd()
					extra_histograms[var].Write()
					if var in extra_histograms_pass:
						extra_histograms_pass[var].Write()
					if var in extra_histograms_fail:
						extra_histograms_fail[var].Write()
				# End loop over extra vars
			# End if SR or muCR
			
			# For matched histograms, also save the matched and unmatched histograms
		# End loop over supersamples
		output_file.Close()
