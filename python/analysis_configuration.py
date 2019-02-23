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
	"qcd":["QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"], # QCD_HT100to200
	"stqq":["ST_t_antitop","ST_t_top","ST_tW_antitop","ST_tW_top"],
	"tqq":["tqq"],
	"wqq":["wqq"],
	"zqq":["zqq"],
	"wqq2018":["WJetsToQQ_HT400to600", "WJetsToQQ_HT600to800", "WJetsToQQ_HT800toInf"],
	"zqq2018":["ZJetsToQQ_HT400to600", "ZJetsToQQ_HT600to800", "ZJetsToQQ_HT800toInf"],
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

