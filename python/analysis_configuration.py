import os

paths = {
	"data":os.path.expandvars("$HOME/DAZSLE/data/"),
	"LimitSetting":os.path.expandvars("$HOME/DAZSLE/data/LimitSetting")
}

# Signal, background, data names
background_names = [
	"qcd",
	"stqq",
	"tqq",
	"wqq",
	"zqq",
	"zll",
	"wlnu",
	"vvqq",
	"hqq125",
	"tthqq125",
	"vbfhqq125",
	"whqq125",
	"zhqq125",
]
# First 12.05 signal processing. Some samples are missing.
signal_names = []
simulated_signal_names = []
signal_models = ["Sbb", "PSbb"] # ZPrime
signal_model_masses = {
	"Sbb":[50,100,125,200,300,350,400,500], 
	"PSbb":[50,100,125,200,300,350,400,500], 
	"ZPrime":[75, 125, 150, 175, 225, 250, 300, 400],
}

signal_masses = {}
#signal_masses = [25,50,75,100,125,150,200,250,300,350,400,500,600,800]
#signal_masses = [50,75,100,125,150,200,250,300,350,400,500,600,800,1000]
#signal_masses = [50,75,100,125,150,200,250,300,400,500]
for model in signal_models:
	for mass in signal_model_masses[model]:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		simulated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass
data_names = ["data_obs", "data_singlemu", "data_obs_ps10", "data_singlemu_ps10"]
supersamples = []
supersamples.extend(background_names)
supersamples.extend(signal_names)
supersamples.extend(data_names)

interpolated_signal_masses = {}
interpolated_signal_names = []
for model in ["Sbb", "PSbb", "ZPrime"]:
	interpolated_signal_masses[model] = [x for x in range(50, 505, 5) if not x in signal_model_masses[model]]
	for mass in interpolated_signal_masses[model]:
		this_signal_name = "{}{}".format(model, mass)
		signal_names.append(this_signal_name)
		interpolated_signal_names.append(this_signal_name)
		signal_masses[this_signal_name] = mass
#for mass in range(50, 325, 25):
#	this_signal_name = "{}{}".format("ZPrime", mass)
#	signal_names.append(this_signal_name)
#	interpolated_signal_names.append(this_signal_name)
#	signal_masses[this_signal_name] = mass


# Sample names. Dictionary is [signal/background/data name]:[list of samples] 
samples = {
	"qcd":["QCD_HT100to200","QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"],
	"stqq":["ST_t_antitop","ST_t_top","ST_tW_antitop","ST_tW_top"],
	"tqq":["tqq"],
	"wqq":["wqq"],
	"zqq":["zqq"],
	"zll":["zll"],
	"wlnu":["wlnu"],
	"wlnu":['wlnu_HT_100To200','wlnu_HT_200To400','wlnu_HT_400To600','wlnu_HT_600To800','wlnu_HT_800To1200','wlnu_HT_1200To2500','wlnu_HT_2500ToInf'],
	"vvqq":["WWTo4Q", "WZ", "ZZ"],
	#"vvqq":["vvqq"],
	"hqq125":["gghbb"],
	"tthqq125":["tthbb"],
	"vbfhqq125":["vbfhbb"],
	"whqq125":["wmhbb", "wphbb"],
	"zhqq125":["zqqhbb", "znunuhbb", "ggzqqhbb", "ggznunuhbb"],
	#"hbb":["gghbb", "vbfhbb", "wmhbb", "wphbb", "tthbb", "zqqhbb", "znunuhbb", "ggzqqhbb", "ggznunuhbb"],#,"vbfhqq125","zhqq125","whqq125","tthqq125"
	"data_obs":["JetHTRun2016B","JetHTRun2016C","JetHTRun2016D","JetHTRun2016E","JetHTRun2016F","JetHTRun2016G","JetHTRun2016H"],
	"data_singlemu":["SingleMuRun2016B","SingleMuRun2016C","SingleMuRun2016D","SingleMuRun2016E","SingleMuRun2016F","SingleMuRun2016G","SingleMuRun2016H"],
	"DMSbb50":["DMSbb50"],
	"DMSbb75":["DMSbb75"],
	"DMSbb100":["DMSbb100"],
	"DMSbb125":["DMSbb125"],
	"DMSbb150":["DMSbb150"],
	"DMSbb200":["DMSbb200"],
	"DMSbb250":["DMSbb250"],
	"DMSbb300":["DMSbb300"],
	#"DMSbb350":["DMSbb350"],
	"DMSbb400":["DMSbb400"],
	"DMSbb500":["DMSbb500"],
	"Sbb10":["Sbb10"],
	"Sbb20":["Sbb20"],
	"Sbb50":["Sbb50"],
	"Sbb75":["Sbb75"],
	"Sbb100":["Sbb100"],
	"Sbb125":["Sbb125"],
	"Sbb150":["Sbb150"],
	"Sbb200":["Sbb200"],
	"Sbb250":["Sbb250"],
	"Sbb300":["Sbb300"],
	"Sbb350":["Sbb350"],
	"Sbb400":["Sbb400"],
	"Sbb500":["Sbb500"],
	"PSbb10":["PSbb10"],
	"PSbb20":["PSbb20"],
	"PSbb50":["PSbb50"],
	"PSbb75":["PSbb75"],
	"PSbb100":["PSbb100"],
	"PSbb125":["PSbb125"],
	"PSbb150":["PSbb150"],
	"PSbb200":["PSbb200"],
	"PSbb250":["PSbb250"],
	"PSbb300":["PSbb300"],
	"PSbb350":["PSbb350"],
	"PSbb400":["PSbb400"],
	"PSbb500":["PSbb500"],
	"ZPrime50":["ZPrime50"],
	"ZPrime75":["ZPrime75"],
	"ZPrime100":["ZPrime100"],
	"ZPrime125":["ZPrime125"],
	"ZPrime150":["ZPrime150"],
	"ZPrime175":["ZPrime175"],
	"ZPrime200":["ZPrime200"],
	"ZPrime225":["ZPrime225"],
	"ZPrime250":["ZPrime250"],
	"ZPrime300":["ZPrime300"],
	"ZPrime400":["ZPrime400"],
	"ZPrime500":["ZPrime500"],
}
# 10% of the data aliases
samples["data_obs_ps10"] = [x + "_ps10" for x in samples["data_obs"]]
samples["data_singlemu_ps10"] = [x + "_ps10" for x in samples["data_singlemu"]]


#for mass in signal_masses:
#	for spin in ["Scalar", "PseudoScalar"]:
#		samples["Pbb_{}_{}".format(mass, spin)] = ["Spin0_ggPhibb1j_{}_{}".format(mass, spin)]

# Skims. Dictionary is [sample name]:[path to skim].
skims = {}
skims["JetHTRun2016B"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016B.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016B_PromptReco_v2_resub.root"]
skims["JetHTRun2016C"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016C.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016C_PromptReco_v2.root"]
skims["JetHTRun2016D"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016D.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016D_PromptReco_v2.root"]
skims["JetHTRun2016E"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016E.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016E_PromptReco_v2.root"]
skims["JetHTRun2016F"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016F.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016F_PromptReco_v1.root"]
skims["JetHTRun2016G"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016G.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016G_PromptReco_v1.root"]
skims["JetHTRun2016H"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016H.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016H_PromptReco_v2.root"]
skims["SingleMuRun2016B"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016B.txt"))] 
skims["SingleMuRun2016C"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016C.txt"))] 
skims["SingleMuRun2016D"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016D.txt"))] 
skims["SingleMuRun2016E"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016E.txt"))] 
skims["SingleMuRun2016F"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016F.txt"))] 
skims["SingleMuRun2016G"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016G.txt"))] 
skims["SingleMuRun2016H"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016H.txt"))] 
skims["JetHTRun2016B_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016B.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016B_PromptReco_v2_resub.root"]
skims["JetHTRun2016C_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016C.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016C_PromptReco_v2.root"]
skims["JetHTRun2016D_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016D.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016D_PromptReco_v2.root"]
skims["JetHTRun2016E_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016E.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016E_PromptReco_v2.root"]
skims["JetHTRun2016F_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016F.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016F_PromptReco_v1.root"]
skims["JetHTRun2016G_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016G.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016G_PromptReco_v1.root"]
skims["JetHTRun2016H_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/JetHTRun2016H.txt"))] #["root://cmsxrootd-site.fnal.gov//store/user/jduarte1/zprimebits-v11.062/JetHTRun2016H_PromptReco_v2.root"]
skims["SingleMuRun2016B_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016B.txt"))] 
skims["SingleMuRun2016C_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016C.txt"))] 
skims["SingleMuRun2016D_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016D.txt"))] 
skims["SingleMuRun2016E_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016E.txt"))] 
skims["SingleMuRun2016F_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016F.txt"))] 
skims["SingleMuRun2016G_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016G.txt"))] 
skims["SingleMuRun2016H_ps10"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/SingleMuRun2016H.txt"))] 
skims["QCD_HT100to200"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT100to200_13TeV.txt"), "r")]
skims["QCD_HT200to300"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT200to300_13TeV.txt"), "r")]
skims["QCD_HT300to500"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT300to500_13TeV.txt"), "r")]
skims["QCD_HT500to700"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT500to700_13TeV.txt"), "r")]
skims["QCD_HT700to1000"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT700to1000_13TeV.txt"), "r")]
skims["QCD_HT1000to1500"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT1000to1500_13TeV.txt"), "r")]
skims["QCD_HT1500to2000"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT1500to2000_13TeV.txt"), "r")]
skims["QCD_HT2000toInf"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/QCD_HT2000toInf_13TeV.txt"), "r")]
#skims["ST_t_antitop"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin.txt"), "r")]
#skims["ST_t_top"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin.txt"), "r")]
#skims["ST_tW_antitop"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4.txt"), "r")]
#skims["ST_tW_top"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4.txt"), "r")]
#skims["wqq"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/WJetsToQQ_HT180_13TeV.txt"), "r")]
#skims["zqq"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/DYJetsToQQ_HT180_13TeV.txt"), "r")]
#skims["tqq"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/TT_powheg.txt"), "r")]
#skims["WWTo4Q"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/WWTo4Q_13TeV_powheg.txt"), "r")]
#skims["WZ"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/WZ_13TeV_pythia8.txt"), "r")]
#skims["ZZ"] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/ZZ_13TeV_pythia8.txt"), "r")]

for mass in signal_model_masses["Sbb"]:
	skims["Sbb{}".format(mass)] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/Spin0_ggPhibb1j_g1_{}_Scalar.txt".format(mass)))]
for mass in signal_model_masses["PSbb"]:
	skims["PSbb{}".format(mass)] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/Spin0_ggPhibb1j_g1_{}_PseudoScalar.txt".format(mass)))]
#for mass in signal_model_masses["ZPrime"]:
#	skims["ZPrime{}".format(mass)] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/VectorDiJet1Jet_{}_13TeV_madgraph.txt".format(mass)))]

# Function to infer the sample from a file path... removing for now, because this doesn't work well on the batch system. E.g. if you have subfiles with generic names (e.g. Output_subjob1.root), there is no way to get the sample.
#def get_sample_from_skim(skim):
#	found_sample = ""
#
#	for sample, filelist in skims.iteritems():
#		for filename in filelist:
#			# Match whole path
#			if skim == filename:
#				found_sample = sample
#				break
#
#			# Match using os.path.samefile
#			if os.path.samefile(filename, skim):
#				found_sample = sample
#				break
#
#			# Match skim subjobs, with generic filenames (Output.root_job...)
#			if "Output.root_job" in skim:
#				file_sample_tag = os.path.basename(os.path.dirname(filename)) # Sample tag is the parent folder
#				if file_sample_tag in skim:
#
#			# Match using basename matching
#			if os.path.basename(filename) == os.path.basename(skim):
#				found_sample = sample
#				break
#			# 
#
#		if found_sample != "":
#			break
#
#	# Search using basename matching
#	if found_sample == "":
#		for sample, filelist
#	return found_sample

# Sklims
sklim_directory = "root://cmsxrootd-site.fnal.gov//store/user/lpchbb/zprimebits-v12.04/sklim2/"
lpchbb_directory = "root://cmsxrootd-site.fnal.gov//store/user/lpchbb"
sklims = { 
	"gghbb":[lpchbb_directory + "/zprimebits-v12.05/GluGluHToBB_M125_13TeV_powheg_pythia8_CKKW_1000pb_weighted.root"],
	"vbfhbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_all_1000pb_weighted.root"],
	"wmhbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root"],
	"wphbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root"],
	"tthbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/ttHTobb_M125_13TeV_powheg_pythia8_1000pb_weighted.root"],
	"zqqhbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root"],
	"znunuhbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_ext_1000pb_weighted.root"],
	"ggzqqhbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root"],
	"ggznunuhbb":[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_1000pb_weighted.root"],
	#'hqq125'            : [lpchbb_directory+'/GluGluHToBB_M125_13TeV_powheg_pythia8_all_1000pb_weighted.root'], 
	#'vbfhqq125'         : [lpchbb_directory+'/VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_all_1000pb_weighted.root'], 
	#'zhqq125'           : [lpchbb_directory+'/ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root', lpchbb_directory+'/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_1000pb_weighted.root', lpchbb_directory+'/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8_ext_1000pb_weighted.root', lpchbb_directory+'/ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root'], 
	#'whqq125'           : [lpchbb_directory+'/WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root', lpchbb_directory+'/WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8_1000pb_weighted.root'], 
	#'tthqq125'          : [lpchbb_directory+'/ttHTobb_M125_13TeV_powheg_pythia8_1000pb_weighted.root'],#ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8_1000pb_weighted.root'],
	'WWTo4Q'            : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/WWTo4Q_13TeV_powheg_1000pb_weighted.root'],
	'ZZ'                : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/ZZ_13TeV_pythia8_1000pb_weighted.root'],
	'WZ'                : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/WZ_13TeV_pythia8_1000pb_weighted.root'], 
	'zqq'               : [lpchbb_directory+'/zprimebits-v12.05/DYJetsToQQ_HT180_13TeV_1000pb_weighted.root'],
	'zll'               : [lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/DYJetsToLL_M_50_13TeV_ext_1000pb_weighted.root"],
	'ST_t_antitop'      : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin_1000pb_weighted.root'],
	'ST_t_top'          : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin_1000pb_weighted.root'],
	'ST_tW_antitop'     : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4_1000pb_weighted.root'],
	'ST_tW_top'         : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4_1000pb_weighted.root'], 
	#'W'                : [lpchbb_directory+'/WJetsToQQ_HT_600ToInf_13TeV_1000pb_weighted.root'], 
	'wqq'               : [lpchbb_directory+'/zprimebits-v12.05/WJetsToQQ_HT180_13TeV_1000pb_weighted.root'], 
	#'wlnu'              : [lpchbb_directory+'WJetsToLNu_HT_70To100_13TeV.root', lpchbb_directory+'WJetsToLNu_HT_100To200_13TeV.root', lpchbb_directory+'/WJetsToLNu_HT_200To400_13TeV.root', lpchbb_directory+'/WJetsToLNu_HT_400To600_13TeV.root', lpchbb_directory+'/WJetsToLNu_HT_600To800_13TeV.root'], 
	'wlnu_HT_100To200':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_100To200_13TeV_1000pb_weighted.root"],
	'wlnu_HT_200To400':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_200To400_13TeV_1000pb_weighted.root"],
	'wlnu_HT_400To600':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_400To600_13TeV_1000pb_weighted.root"],
	'wlnu_HT_600To800':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_600To800_13TeV_1000pb_weighted.root"],
	'wlnu_HT_800To1200':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_800To1200_13TeV_1000pb_weighted.root"],
	'wlnu_HT_1200To2500':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_1200To2500_13TeV_1000pb_weighted.root"],
	'wlnu_HT_2500ToInf':[lpchbb_directory + "/zprimebits-v12.04/norm2/cvernier/WJetsToLNu_HT_2500ToInf_13TeV_1000pb_weighted.root"],

	#'TTbar'            : [lpchbb_directory+'/TTJets_13TeV_1000pb_weighted.root'], #MadGraph is the old default 
	'tqq'             : [lpchbb_directory+'/zprimebits-v12.04/norm2/cvernier/TT_powheg_1000pb_weighted.root'], #Powheg is the new default 
	#"QCD_HT100to200"    : [lpchbb_directory+'/QCD_HT100to200_13TeV_1000pb_weighted.root'],
	#"QCD_HT200to300"    : [lpchbb_directory+'/QCD_HT200to300_13TeV_all_1000pb_weighted.root'],
	#"QCD_HT300to500"    : [lpchbb_directory+'/QCD_HT300to500_13TeV_all_1000pb_weighted.root'],
	#"QCD_HT500to700"    : [lpchbb_directory+'/QCD_HT500to700_13TeV_ext_1000pb_weighted.root'],
	#"QCD_HT700to1000"   : [lpchbb_directory+'/QCD_HT700to1000_13TeV_ext_1000pb_weighted.root'],
	#"QCD_HT1000to1500"  : [lpchbb_directory+'/QCD_HT1000to1500_13TeV_all_1000pb_weighted.root'],
	#"QCD_HT1500to2000"  : [lpchbb_directory+'/QCD_HT1500to2000_13TeV_all_1000pb_weighted.root'],
	#"QCD_HT2000toInf"   : [lpchbb_directory+'/QCD_HT2000toInf_13TeV_1000pb_weighted.root'],
	#"JetHTRun2016B":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016B.txt"))],
	#"JetHTRun2016C":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016C.txt"))],
	#"JetHTRun2016D":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016D.txt"))],
	#"JetHTRun2016E":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016E.txt"))],
	#"JetHTRun2016F":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016F.txt"))],
	#"JetHTRun2016G":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016G.txt"))],
	#"JetHTRun2016H":[x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/lpchbb_directory/cmslpc/JetHTRun2016H.txt"))],
	#'SingleMuRun2016B': [lpchbb_directory+'/SingleMuonRun2016B_03Feb2017_ver1_v1_fixtrig.root', lpchbb_directory+'/SingleMuonRun2016B_03Feb2017_ver2_v2_fixtrig.root'], 
	#'SingleMuRun2016C': [lpchbb_directory+'/SingleMuonRun2016C_03Feb2017_v1_fixtrig.root'], 
	#'SingleMuRun2016D': [lpchbb_directory+'/SingleMuonRun2016D_03Feb2017_v1_fixtrig.root'],
	#'SingleMuRun2016E': [lpchbb_directory+'/SingleMuonRun2016E_03Feb2017_v1_fixtrig.root'],
	#'SingleMuRun2016F': [lpchbb_directory+'/SingleMuonRun2016F_03Feb2017_v1_fixtrig.root'], 
	#'SingleMuRun2016G': [lpchbb_directory+'/SingleMuonRun2016G_03Feb2017_v1_fixtrig.root'], 
	#'SingleMuRun2016H': [lpchbb_directory+'/SingleMuonRun2016H_03Feb2017_ver2_v1_fixtrig.root', lpchbb_directory+'/SingleMuonRun2016H_03Feb2017_ver3_v1_fixtrig.root'],
}
for mass in signal_model_masses["Sbb"]:
	sklims["Sbb{}".format(mass)] = [lpchbb_directory + "/zprimebits-v12.05/Spin0_ggPhibb1j_g1_{}_Scalar_1000pb_weighted.root".format(mass)]
for mass in signal_model_masses["PSbb"]:
	sklims["PSbb{}".format(mass)] = [lpchbb_directory + "/zprimebits-v12.05/Spin0_ggPhibb1j_g1_{}_PseudoScalar_1000pb_weighted.root".format(mass)]
for mass in signal_model_masses["ZPrime"]:
	sklims["ZPrime{}".format(mass)] = ["root://cmseos.fnal.gov//store/user/lpcbacon/dazsle/zprimebits-v12.07/norm/VectorDiJet1Jet_madgraph_Mphi{}Mchi3000_13TeV_noPF_8_1000pb_weighted.root".format(mass)]

for mass in [50,75,100,125,150,200,250,300,400,500]:
	sklims["DMSbb{}".format(mass)] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/cmslpc/DMSpin0_ggPhibb1j_{}.txt".format(mass)))]
	#sklims["Sbb{}".format(mass)] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/Spin0_ggPhi12j_g1_{}_Scalar_13TeV_madgraph.txt".format(mass)))]
	#sklims["PSbb{}".format(mass)] = [x.strip() for x in open(os.path.expandvars("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/data/skim_directory/lxplus/Spin0_ggPhi12j_g1_{}_PseudoScalar_13TeV_madgraph.txt".format(mass)))]

# Need to re-make these sklims for v12! The processing barely did any events.
#for mass in signal_masses:
#	for spin in ["Scalar", "PseudoScalar"]:
#		sklims["Spin0_ggPhibb1j_{}_{}".format(mass, spin)] = ["root://cmsxrootd-site.fnal.gov//store/user/lpchbb/zprimebits-v12.02/norm/Spin0_ggPhi12j_g1_{}_{}_13TeV_madgraph_1000pb_weighted.root".format(mass, spin)]


def get_sample_from_sklim(sklim):
	found_sample = ""
	for sample, filename in sklims.iteritems():
		if os.path.basename(filename) == os.path.basename(sklim):
			found_sample = sample
			break
	return found_sample

def get_histogram_file(selection, jet_type):
	return paths["LimitSetting"] + "/Xbb_inputs/histograms_{}_{}.root".format(selection, jet_type)

def get_interpolation_file(jet_type):
	return paths["LimitSetting"] + "/Xbb_inputs/interpolations_{}.root".format(jet_type)

