import os
import sys
import ROOT
from DAZSLE.PhiBBPlusJet.analysis_base import AnalysisBase
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
import DAZSLE.PhiBBPlusJet.event_selections as event_selections
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
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()

class EventSelectionHistograms(AnalysisBase):
	def __init__(self, sample_name, tree_name="otree", jet_type="AK8"):
		super(EventSelectionHistograms, self).__init__(tree_name=tree_name)
		self._data = BaconData(self._chain)
		self._output_path = ""
		self._sample_name = sample_name
		self._input_nevents = 0
		self._n2_ddt_cut = 0.
		self._dcsv_cut = params[jet_type]["DCSV"]
		self._dcsv_cut_loose = 0.8
		self._dcsv_min = -999.
		self._jet_type = jet_type
		self._selections = ["Preselection", "SR", "muCR", "N2CR", "N2SR"]
		self._do_optimization = False
		self._data_source = "data"
		self._prescale = -1

		# Weight systematics: these only affect the weights used to fill histograms, so can easily be filled in normal running
		self._weight_systematics = {
			"SR":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"N2SR":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"N2CR":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"Preselection":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"muCR":["MuTriggerUp", "MuTriggerDown", "MuIDUp", "MuIDDown", "MuIsoUp", "MuIsoDown", "PUUp", "PUDown"]
		}
		# Jet systematics: these affect the jet pT, so modify the event selection
		self._jet_systematics = ["JESUp", "JESDown", "JERUp", "JERDown"]

	def set_prescale(self, prescale):
		self._prescale = prescale

	def do_optimization(self, do_opt=True):
		self._do_optimization = do_opt

	def set_cut(self, cut_name, cut_value):
		if cut_name == "n2_ddt":
			self._n2_ddt_cut = cut_value
		elif cut_name == "dcsv":
			self._dcsv_cut = cut_value
		else:
			print "[EventSelectionHistograms::set_cut] ERROR : Unknown customizable cut {}. Exiting.".format(cut_name)
			sys.exit(1)

	# Commented out on 17/10/24: setting this later can lead to inconsistencies, so instead specify in CONSTRUCTOR, and never change again.
	#def set_jet_type(self, jet_type):
	#	if jet_type != "AK8" and jet_type != "CA15":
	#		print "[EventSelectionHistograms::set_jet_type] ERROR : Unknown jet type {}. Exiting.".format(jet_type)
	#		sys.exit(1)
	#	self._jet_type = jet_type 

	def set_data_source(self, data_source):
		print "[EventSelectionHistograms::set_data_source] INFO : Setting data source to " + data_source
		if not data_source in ["data", "simulation"]:
			print "[EventSelectionHistograms] ERROR : Data source must be data or simulation."
			sys.exit(1)
		self._data_source = data_source

	# Overload add_file to extract the number of input events to the skims, stored in histogram NEvents in the same file as the trees
	def add_file(self, filename):
		super(EventSelectionHistograms, self).add_file(filename)
		f = ROOT.TFile.Open(filename, "READ")
		if f.Get("NEvents").Integral() == 0:
			print "[EventSelectionHistograms::add_file] ERROR : NEvents.Integral() == 0 for file " + filename
			sys.exit(1)
		self._input_nevents += f.Get("NEvents").Integral()
		f.Close()

	def set_output_path(self, output_path):
		self._output_path = output_path
		os.system("mkdir -pv {}".format(os.path.dirname(self._output_path)))

	def start(self):
		self._processed_events = 0

		# Histograms
		self._pt_bins = array.array("d", [450., 500.,550.,600.,675.,800.,1000.])

		self._histograms = ROOT.Root.HistogramManager()
		self._histograms.AddPrefix("h_")
		self._histograms.AddTH1F("input_nevents", "input_nevents", "", 1, -0.5, 0.5)
		self._histograms.GetTH1F("input_nevents").SetBinContent(1, self._input_nevents)
		self._histograms.AddTH1D("processed_nevents", "processed_nevents", "", 1, -0.5, 0.5)

		# Inclusive distributions
		self._histograms.AddTH1D("inclusive_pt", "inclusive_pt", "p_{T} [GeV]", 100, 0., 1000.)
		self._histograms.AddTH1D("inclusive_eta", "inclusive_eta", "#eta", 100, -5., 5.)
		self._histograms.AddTH1D("inclusive_msd", "inclusive_msd", "", 100, 0., 1000.)
		self._histograms.AddTH1D("inclusive_n2ddt", "inclusive_n2ddt", "", 40, -2., 2.)
		self._histograms.AddTH1D("inclusive_dcsv", "inclusive_dcsv", "", 200, -1., 1.)


		# tau21 optimization
		if self._do_optimization:
			for tau21_ddt_cut in [0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7]:
				self._selections.append("SR_tau21ddt{}".format(tau21_ddt_cut))
				self._weight_systematics["SR_tau21ddt{}".format(tau21_ddt_cut)] = ["TriggerUp", "TriggerDown", "PUUp", "PUDown"]

			# dcsv optimization
			self._dcsv_cuts = [0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]

		# Histograms for each event selection
		self._selection_histograms = {}
		for selection in self._selections:
			self._selection_histograms[selection] = ROOT.Root.HistogramManager()
			self._selection_histograms[selection].AddPrefix("h_{}_{}_".format(selection, self._jet_type))

			self._selection_histograms[selection].AddTH1D("pass_nevents", "pass_nevents", "", 1, -0.5, 0.5)
			self._selection_histograms[selection].AddTH1D("pass_nevents_weighted", "pass_nevents_weighted", "", 1, -0.5, 0.5)
			self._selection_histograms[selection].AddTH2D("pt_dcsv", "pt_dcsv", "p_{T} [GeV]", 200, 0., 2000., "Double b-tag", 20, -1., 1.)

			self._selection_histograms[selection].AddTH2D("pass", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
			self._selection_histograms[selection].AddTH2D("pass_unweighted", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
			self._selection_histograms[selection].AddTH2D("passloose", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
			self._selection_histograms[selection].AddTH2D("fail", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
			self._selection_histograms[selection].AddTH2D("fail_unweighted", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)

			self._selection_histograms[selection].AddTH1D("pfmet", "PF MET", "PF MET [GeV]", 200, 0., 1000.)
			self._selection_histograms[selection].AddTH1D("pass_pfmet", "PF MET", "PF MET [GeV]", 200, 0., 1000.)
			self._selection_histograms[selection].AddTH1D("fail_pfmet", "PF MET", "PF MET [GeV]", 200, 0., 1000.)

			self._selection_histograms[selection].AddTH1D("dcsv", "dcsv", "dcsv", 200, -1., 1.)
			self._selection_histograms[selection].AddTH1D("pass_dcsv", "dcsv", "dcsv", 200, -1., 1.)
			self._selection_histograms[selection].AddTH1D("fail_dcsv", "dcsv", "dcsv", 200, -1., 1.)

			self._selection_histograms[selection].AddTH1D("n2ddt", "n2ddt", "n2ddt", 20, -0.5, 0.5)
			self._selection_histograms[selection].AddTH1D("pass_n2ddt", "n2ddt", "n2ddt", 20, -0.5, 0.5)
			self._selection_histograms[selection].AddTH1D("fail_n2ddt", "n2ddt", "n2ddt", 20, -0.5, 0.5)

			self._selection_histograms[selection].AddTH1D("pt", "pt", "pt", 400, 0., 2000.)
			self._selection_histograms[selection].AddTH1D("pass_pt", "pt", "pt", 400, 0., 2000.)
			self._selection_histograms[selection].AddTH1D("fail_pt", "pt", "pt", 400, 0., 2000.)

			self._selection_histograms[selection].AddTH1D("msd", "msd", "msd", 85, 5, 600)
			self._selection_histograms[selection].AddTH1D("pass_msd", "msd", "msd", 85, 5, 600)
			self._selection_histograms[selection].AddTH1D("fail_msd", "msd", "msd", 85, 5, 600)

			self._selection_histograms[selection].AddTH1D("eta", "eta", "eta", 60, -3., 3.)
			self._selection_histograms[selection].AddTH1D("pass_eta", "eta", "eta", 60, -3., 3.)
			self._selection_histograms[selection].AddTH1D("fail_eta", "eta", "eta", 60, -3., 3.)

			self._selection_histograms[selection].AddTH1D("rho", "rho", "rho", 80, -8., 0.)
			self._selection_histograms[selection].AddTH1D("pass_rho", "rho", "rho", 80, -8., 0.)
			self._selection_histograms[selection].AddTH1D("fail_rho", "rho", "rho", 80, -8., 0.)

			self._selection_histograms[selection].AddTH2D("dcsv_vs_rho", "dcsv_vs_rho", 
				"Double-b", 110, -1.1, 1.1,
				"#rho", 40, -8., 0.)
			self._selection_histograms[selection].AddTH2D("dcsv_vs_pt", "dcsv_vs_pt", 
				"Double-b", 110, -1.1, 1.1,
				"p_{T} [GeV]", 50, 0., 1000.)
			self._selection_histograms[selection].AddTH2D("dcsv_vs_msd", "dcsv_vs_msd", 
				"Double-b", 110, -1.1, 1.1,
				"m_{SD} [GeV]", 100, 40, 740)
			self._selection_histograms[selection].AddTH3D("dcsv_vs_msd_vs_pt", "dcsv_vs_msd_vs_pt", 
				"Double-b", 110, -1.1, 1.1,
				"m_{SD} [GeV]", 80, 40, 600,
				"p_{T} [GeV]", 12, 400, 1000)

			if self._jet_type == "CA15":
				self._selection_histograms[selection].AddTH2D("dcsvalt_vs_rho", "dcsvalt_vs_rho", 
					"Double-b (alternate)", 110, -1.1, 1.1,
					"#rho", 40, -8., 0.)
				self._selection_histograms[selection].AddTH2D("dcsvalt_vs_pt", "dcsvalt_vs_pt", 
					"Double-b (alternate)", 110, -1.1, 1.1,
					"p_{T} [GeV]", 50, 0., 1000.)
				self._selection_histograms[selection].AddTH2D("dcsvalt_vs_msd", "dcsvalt_vs_msd", 
					"Double-b (alternate)", 110, -1.1, 1.1,
					"m_{SD} [GeV]", 100, 40, 740)
				self._selection_histograms[selection].AddTH3D("dcsvalt_vs_msd_vs_pt", "dcsvalt_vs_msd_vs_pt", 
					"Double-b (alternate)", 110, -1.1, 1.1,
					"m_{SD} [GeV]", 100, 40, 740,
					"p_{T} [GeV]", 6, 400, 1000)

			self._selection_histograms[selection].AddTH3D("n2ddt_vs_msd_vs_pt", "n2ddt_vs_msd_vs_pt", 
				"N_{2}^{DDT}", 40, -0.5, 0.5,
				"m_{SD} [GeV]", 40, 40, 600,
				"p_{T} [GeV]", 12, 400, 1000)

			if self._do_optimization:
				for dcsv_cut in self._dcsv_cuts:
					self._selection_histograms[selection].AddTH2D("pass_dcsv{}".format(dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
					self._selection_histograms[selection].AddTH2D("pass_unweighted_dcsv{}".format(dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
					self._selection_histograms[selection].AddTH2D("fail_dcsv{}".format(dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
					self._selection_histograms[selection].AddTH2D("fail_unweighted_dcsv{}".format(dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)

			for systematic in self._weight_systematics[selection] + self._jet_systematics:
				self._selection_histograms[selection].AddTH2D("pass_{}".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("passloose_{}".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("fail_{}".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				if self._do_optimization:
					for dcsv_cut in self._dcsv_cuts:
						self._selection_histograms[selection].AddTH2D("pass_{}_dcsv{}".format(systematic, dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
						self._selection_histograms[selection].AddTH2D("pass_{}_unweighted_dcsv{}".format(systematic, dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)					

						self._selection_histograms[selection].AddTH2D("fail_{}_dcsv{}".format(systematic, dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
						self._selection_histograms[selection].AddTH2D("fail_{}_unweighted_dcsv{}".format(systematic, dcsv_cut), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)

			# Pass/fail histograms for matched V in simulation
			if self._data_source == "simulation":
				self._selection_histograms[selection].AddTH2D("pass_matched", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("pass_unmatched", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("passloose_matched", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("passloose_unmatched", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("fail_matched", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				self._selection_histograms[selection].AddTH2D("fail_unmatched", "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
				for systematic in self._weight_systematics[selection] + self._jet_systematics:
					self._selection_histograms[selection].AddTH2D("pass_{}_matched".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
					self._selection_histograms[selection].AddTH2D("pass_{}_unmatched".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)					
					self._selection_histograms[selection].AddTH2D("passloose_{}_matched".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
					self._selection_histograms[selection].AddTH2D("passloose_{}_unmatched".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)					
					self._selection_histograms[selection].AddTH2D("fail_{}_matched".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)
					self._selection_histograms[selection].AddTH2D("fail_{}_unmatched".format(systematic), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", len(self._pt_bins) - 1, self._pt_bins)					

		# Event selections
		self._event_selectors = {}
		self._event_selectors["SR"] = event_selections.MakeSRSelector(self._jet_type)
		self._event_selectors["N2SR"] = event_selections.MakeN2SRSelector(self._jet_type, dbtag_cut=self._dcsv_cut)
		self._event_selectors["N2CR"] = event_selections.MakeN2CRSelector(self._jet_type)
		self._event_selectors["Preselection"] = event_selections.MakePreselectionSelector(self._jet_type)
		self._event_selectors["muCR"] = event_selections.MakeMuCRSelector(self._jet_type)
		self._event_selectors_syst = {"SR":{}, "muCR":{}, "Preselection":{}, "N2CR":{}, "N2SR":{}}
		for systematic in self._jet_systematics:
			self._event_selectors_syst["SR"][systematic] = event_selections.MakeSRSelector(self._jet_type, jet_systematic=systematic)
			self._event_selectors_syst["N2SR"][systematic] = event_selections.MakeN2SRSelector(self._jet_type, dbtag_cut=self._dcsv_cut, jet_systematic=systematic)
			self._event_selectors_syst["N2CR"][systematic] = event_selections.MakeN2CRSelector(self._jet_type, jet_systematic=systematic)
			self._event_selectors_syst["Preselection"][systematic] = event_selections.MakePreselectionSelector(self._jet_type, jet_systematic=systematic)
			self._event_selectors_syst["muCR"][systematic] = event_selections.MakeMuCRSelector(self._jet_type, jet_systematic=systematic)

		if self._do_optimization:
			for tau21_ddt_cut in [0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7]:
				selection_name = "SR_tau21ddt{}".format(tau21_ddt_cut)
				self._event_selectors[selection_name] = event_selections.MakeSRSelector(self._jet_type, n2_ddt_cut=None, tau21_ddt_cut=tau21_ddt_cut, tag="tau21ddt{}".format(tau21_ddt_cut))
				self._event_selectors_syst[selection_name] = {}
				for systematic in self._jet_systematics:
					self._event_selectors_syst[selection_name][systematic] = event_selections.MakeSRSelector(self._jet_type, jet_systematic=systematic, n2_ddt_cut=None, tau21_ddt_cut=tau21_ddt_cut, tag="tau21ddt{}".format(tau21_ddt_cut))

		# Pileup weight stuff
		f_pu = TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/puWeights_All.root", "read")
		self._h_pu_weight = f_pu.Get("puw")
		self._h_pu_weight.SetDirectory(0)
		self._h_pu_weight_up = f_pu.Get("puw_p")
		self._h_pu_weight_up.SetDirectory(0)
		self._h_pu_weight_down = f_pu.Get("puw_m")
		self._h_pu_weight_down.SetDirectory(0)
		f_pu.Close()

		# Trigger efficiency weight stuff
		if self._jet_type == "AK8":
			f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_AK8_SingleMuon_Run2016_V2p1_v03.root", "read")
			self._trig_den = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
			self._trig_num = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
		elif self._jet_type == "CA15":
			f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_CA15_SingleMuon_Run2016_V2p4_v08.root", "read")
			self._trig_den = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
			self._trig_num = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
		self._trig_den.SetDirectory(0)
		self._trig_num.SetDirectory(0)
		self._trig_den.RebinX(2)
		self._trig_num.RebinX(2)
		self._trig_den.RebinY(5)
		self._trig_num.RebinY(5)
		self._trig_eff = ROOT.TEfficiency()
		if (ROOT.TEfficiency.CheckConsistency(self._trig_num, self._trig_den)):
			self._trig_eff = ROOT.TEfficiency(self._trig_num, self._trig_den)
			self._trig_eff.SetDirectory(0)
		f_trig.Close()

		# get muon trigger efficiency object

		lumi_GH = 16.146
		lumi_BCDEF = 19.721
		lumi_total = lumi_GH + lumi_BCDEF

		f_mutrig_GH = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_Period4.root", "read")
		self._mutrig_eff_GH = f_mutrig_GH.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
		self._mutrig_eff_GH.Sumw2()
		self._mutrig_eff_GH.SetDirectory(0)
		f_mutrig_GH.Close()

		f_mutrig_BCDEF = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_RunBtoF.root", "read")
		self._mutrig_eff_BCDEF = f_mutrig_BCDEF.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
		self._mutrig_eff_BCDEF.Sumw2()
		self._mutrig_eff_BCDEF.SetDirectory(0)
		f_mutrig_BCDEF.Close()

		self._mutrig_eff = self._mutrig_eff_GH.Clone('pt_abseta_DATA_mutrig_ave')
		self._mutrig_eff.Scale(lumi_GH / lumi_total)
		self._mutrig_eff.Add(self._mutrig_eff_BCDEF, lumi_BCDEF / lumi_total)

		# get muon ID efficiency object

		f_muid_GH = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_GH.root", "read")
		self._muid_eff_GH = f_muid_GH.Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muid_eff_GH.Sumw2()
		self._muid_eff_GH.SetDirectory(0)
		f_muid_GH.Close()

		f_muid_BCDEF = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_BCDEF.root", "read")
		self._muid_eff_BCDEF = f_muid_BCDEF.Get(
			"MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muid_eff_BCDEF.Sumw2()
		self._muid_eff_BCDEF.SetDirectory(0)
		f_muid_BCDEF.Close()

		self._muid_eff = self._muid_eff_GH.Clone('pt_abseta_DATA_muid_ave')
		self._muid_eff.Scale(lumi_GH / lumi_total)
		self._muid_eff.Add(self._muid_eff_BCDEF, lumi_BCDEF / lumi_total)

		# get muon ISO efficiency object

		f_muiso_GH = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_ISO_GH.root", "read")
		self._muiso_eff_GH = f_muiso_GH.Get("LooseISO_LooseID_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muiso_eff_GH.Sumw2()
		self._muiso_eff_GH.SetDirectory(0)
		f_muiso_GH.Close()

		f_muiso_BCDEF = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_ISO_BCDEF.root", "read")
		self._muiso_eff_BCDEF = f_muiso_BCDEF.Get("LooseISO_LooseID_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muiso_eff_BCDEF.Sumw2()
		self._muiso_eff_BCDEF.SetDirectory(0)
		f_muiso_BCDEF.Close()

		self._muiso_eff = self._muiso_eff_GH.Clone('pt_abseta_DATA_muiso_ave')
		self._muiso_eff.Scale(lumi_GH / lumi_total)
		self._muiso_eff.Add(self._muiso_eff_BCDEF, lumi_BCDEF / lumi_total)


	def run(self, max_nevents=-1, first_event=0):
		if max_nevents > 0:
			limit_nevents = min(max_nevents, self._chain.GetEntries())
		else:
			limit_nevents = self._chain.GetEntries()

		n_checkpoints = 20
		print_every = int(ceil(1. * limit_nevents / n_checkpoints))

		print "[EventSelectionHistograms::run] INFO : Running loop over tree from event {} to {}".format(first_event, limit_nevents - 1)

		self.start_timer()
		for entry in xrange(first_event, limit_nevents):
			self.print_progress(entry, first_event, limit_nevents, print_every)
			self._data.GetEntry(entry)

			# Prescale before anything
			if self._prescale > 0:
				if self._data.evtNum % self._prescale != 0:
					continue

			self._histograms.GetTH1D("processed_nevents").Fill(0)
			self._processed_events += 1

			npu = min(self._data.npu, 49.5)
			pu_weight = self._h_pu_weight.GetBinContent(self._h_pu_weight.FindBin(npu))
			pu_weight_up = self._h_pu_weight_up.GetBinContent(self._h_pu_weight_up.FindBin(npu))
			pu_weight_down = self._h_pu_weight_down.GetBinContent(self._h_pu_weight_down.FindBin(npu))

			k_vjets = 1.
			k_ttbar = 1.
			w_scale = {
				(0, 500):1.0,
				(500, 600):1.0,
				(600, 700):1.0,
				(700, 800):1.2,
				(800, 900):1.25,
				(900, 1000):1.25,
				(1000, 3000):1.0
			}
			if self._sample_name == 'wqq' or self._sample_name == 'W':
				k_vjets = self._data.kfactor * 1.35  # ==1 for not V+jets events
				for pt_range, w_sf in w_scale.iteritems():
					if pt_range[0] < self._data.genVPt < pt_range[1]:
						k_vjets *= w_sf
			elif self._sample_name == 'zqq' or self._sample_name == 'DY':
				k_vjets = self._data.kfactor * 1.45  # ==1 for not V+jets events
			elif self._sample_name == 'tqq':
				k_ttbar = self._data.topPtWeight

			for selection in self._selections:
				# Get weights
				if self._data_source == "data":
					event_weight = 1.
					event_weight_syst = {}
					if "SR" in selection or "Preselection" in selection or "N2CR" in selection:
						event_weight_syst["TriggerUp"] = 1.
						event_weight_syst["TriggerDown"] = 1.
						event_weight_syst["PUUp"] = 1.
						event_weight_syst["PUDown"] = 1.
					elif "muCR" in selection:
						event_weight_syst["MuTriggerUp"] = 1.
						event_weight_syst["MuTriggerDown"] = 1.
						event_weight_syst["MuIDUp"] = 1.
						event_weight_syst["MuIDDown"] = 1.
						event_weight_syst["MuIsoUp"] = 1.
						event_weight_syst["MuIsoDown"] = 1.
						event_weight_syst["PUUp"] = 1.
						event_weight_syst["PUDown"] = 1.
					else:
						print "[event_selection_histograms::run] ERROR : selection {} is not known. Fix!".format(selection)
						sys.exit(1)
				else:
					if "SR" in selection or "Preselection" in selection or "N2CR" in selection:
						if self._jet_type == "AK8":
							trigger_mass = min(self._data.AK8Puppijet0_msd, 300.)
							trigger_pt = max(200., min(self._data.AK8Puppijet0_pt, 1000.))
						elif self._jet_type == "CA15":
							trigger_mass = min(self._data.CA15Puppijet0_msd, 300.)
							trigger_pt = max(200., min(self._data.CA15Puppijet0_pt, 1000.))
						trigger_weight = self._trig_eff.GetEfficiency(self._trig_eff.FindFixBin(trigger_mass, trigger_pt))
						trigger_weight_up = trigger_weight + self._trig_eff.GetEfficiencyErrorUp(self._trig_eff.FindFixBin(trigger_mass, trigger_pt))
						trigger_weight_down = trigger_weight - self._trig_eff.GetEfficiencyErrorLow(
							self._trig_eff.FindFixBin(trigger_mass, trigger_pt))
						if trigger_weight <= 0 or trigger_weight_down <= 0 or trigger_weight_up <= 0:
							#print 'trigger_weights are %f, %f, %f, setting all to 1' % (trigger_weight, trigger_weight_up, trigger_weight_down)
							trigger_weight = 1
							trigger_weight_down = 1
							trigger_weight_up = 1

						event_weight = pu_weight * k_vjets * trigger_weight * k_ttbar
						event_weight_syst = {}
						event_weight_syst["TriggerUp"] = pu_weight * k_vjets * k_ttbar *trigger_weight_up
						event_weight_syst["TriggerDown"] = pu_weight * k_vjets * k_ttbar * trigger_weight_down
						event_weight_syst["PUUp"] = pu_weight_up * k_vjets * k_ttbar * trigger_weight
						event_weight_syst["PUDown"] = pu_weight_down * k_vjets * k_ttbar * trigger_weight

					elif "muCR" in selection:
						mutrigweight = 1
						mutrigweightDown = 1
						mutrigweightUp = 1
						if self._data.nmuLoose > 0:
							muPtForTrig = max(52., min(self._data.vmuoLoose0_pt, 700.))
							muEtaForTrig = min(abs(self._data.vmuoLoose0_eta), 2.3)
							mutrigweight = self._mutrig_eff.GetBinContent(self._mutrig_eff.FindBin(muPtForTrig, muEtaForTrig))
							mutrigweightUp = mutrigweight + self._mutrig_eff.GetBinError(
								self._mutrig_eff.FindBin(muPtForTrig, muEtaForTrig))
							mutrigweightDown = mutrigweight - self._mutrig_eff.GetBinError(
								self._mutrig_eff.FindBin(muPtForTrig, muEtaForTrig))
							if mutrigweight <= 0 or mutrigweightDown <= 0 or mutrigweightUp <= 0:
								print 'mutrigweights are %f, %f, %f, setting all to 1' % (
								mutrigweight, mutrigweightUp, mutrigweightDown)
								mutrigweight = 1
								mutrigweightDown = 1
								mutrigweightUp = 1

						muidweight = 1
						muidweightDown = 1
						muidweightUp = 1
						if self._data.nmuLoose > 0:
							muPtForId = max(20., min(self._data.vmuoLoose0_pt, 100.))
							muEtaForId = min(abs(self._data.vmuoLoose0_eta), 2.3)
							muidweight = self._muid_eff.GetBinContent(self._muid_eff.FindBin(muPtForId, muEtaForId))
							muidweightUp = muidweight + self._muid_eff.GetBinError(self._muid_eff.FindBin(muPtForId, muEtaForId))
							muidweightDown = muidweight - self._muid_eff.GetBinError(self._muid_eff.FindBin(muPtForId, muEtaForId))
							if muidweight <= 0 or muidweightDown <= 0 or muidweightUp <= 0:
								print 'muidweights are %f, %f, %f, setting all to 1' % (muidweight, muidweightUp, muidweightDown)
								muidweight = 1
								muidweightDown = 1
								muidweightUp = 1

						muisoweight = 1
						muisoweightDown = 1
						muisoweightUp = 1
						if self._data.nmuLoose > 0:
							muPtForIso = max(20., min(self._data.vmuoLoose0_pt, 100.))
							muEtaForIso = min(abs(self._data.vmuoLoose0_eta), 2.3)
							muisoweight = self._muiso_eff.GetBinContent(self._muiso_eff.FindBin(muPtForIso, muEtaForIso))
							muisoweightUp = muisoweight + self._muiso_eff.GetBinError(
								self._muiso_eff.FindBin(muPtForIso, muEtaForIso))
							muisoweightDown = muisoweight - self._muiso_eff.GetBinError(
								self._muiso_eff.FindBin(muPtForIso, muEtaForIso))
							if muisoweight <= 0 or muisoweightDown <= 0 or muisoweightUp <= 0:
								print 'muisoweights are %f, %f, %f, setting all to 1' % (
								muisoweight, muisoweightUp, muisoweightDown)
								muisoweight = 1
								muisoweightDown = 1
								muisoweightUp = 1

						event_weight = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweight
						event_weight_syst = {}
						event_weight_syst["MuTriggerUp"] = pu_weight * k_vjets * k_ttbar * mutrigweightUp * muidweight * muisoweight
						event_weight_syst["MuTriggerDown"] = pu_weight * k_vjets * k_ttbar * mutrigweightDown * muidweight * muisoweight
						event_weight_syst["MuIDUp"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweightUp * muisoweight
						event_weight_syst["MuIDDown"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweightDown * muisoweight
						event_weight_syst["MuIsoUp"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweightUp
						event_weight_syst["MuIsoDown"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweightDown
						event_weight_syst["PUUp"] = pu_weight_up * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweight
						event_weight_syst["PUDown"] = pu_weight_down * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweight

				# Pick up AK8 or CA15 event variables here, to avoid mistakes later
				if self._jet_type == "AK8":
					fatjet_pt = self._data.AK8Puppijet0_pt
					fatjet_eta = self._data.AK8Puppijet0_eta
					fatjet_msd = self._data.AK8Puppijet0_msd_puppi
					fatjet_dcsv = self._data.AK8Puppijet0_doublecsv
					fatjet_dcsv_alt = 0.
					fatjet_n2ddt = self._data.AK8Puppijet0_N2DDT
					fatjet_rho = self._data.AK8Puppijet0_rho
					fatjet_phi = self._data.AK8Puppijet0_phi
				elif self._jet_type == "CA15":
					fatjet_pt = self._data.CA15Puppijet0_pt
					fatjet_eta = self._data.CA15Puppijet0_eta
					fatjet_msd = self._data.CA15Puppijet0_msd_puppi
					fatjet_dcsv = self._data.CA15Puppijet0_doublesub
					fatjet_dcsv_alt = self._data.CA15Puppijet0_doublecsv
					fatjet_n2ddt = self._data.CA15Puppijet0_N2DDT
					fatjet_rho = self._data.CA15Puppijet0_rho
					fatjet_phi = self._data.CA15Puppijet0_phi

				# Inclusive histograms
				self._histograms.GetTH1D("inclusive_pt").Fill(fatjet_pt, event_weight)
				self._histograms.GetTH1D("inclusive_eta").Fill(fatjet_eta, event_weight)
				self._histograms.GetTH1D("inclusive_msd").Fill(fatjet_msd, event_weight)
				self._histograms.GetTH1D("inclusive_n2ddt").Fill(fatjet_n2ddt, event_weight)
				self._histograms.GetTH1D("inclusive_dcsv").Fill(fatjet_dcsv, event_weight)

				# Run selection and fill histograms
				self._event_selectors[selection].process_event(self._data, event_weight)
				if selection == "N2SR":
					pass_box = (fatjet_n2ddt < 0.)
					fail_box = (fatjet_n2ddt >= 0.)
					pass_loose_box = (fatjet_n2ddt < 0.)
				else:
					pass_box = (fatjet_dcsv > self._dcsv_cut)
					fail_box = (fatjet_dcsv > self._dcsv_min)
					pass_loose_box = (fatjet_dcsv > self._dcsv_cut_loose)
				if self._event_selectors[selection].event_pass():
					self._selection_histograms[selection].GetTH1D("pass_nevents").Fill(0)
					self._selection_histograms[selection].GetTH1D("pass_nevents_weighted").Fill(0, event_weight)

					# For simulated V(bb), match fat jet to parent truth particle
					if self._data_source == "simulation":
						vmatched = False
						if self._data.genVPt > 0 and self._data.genVMass > 0:
							matching_dphi = abs(math.acos(math.cos(self._data.genVPhi - fatjet_phi)))
							matching_dpt = abs(self._data.genVPt - fatjet_pt) / self._data.genVPt
							matching_dmass = abs(self._data.genVMass - fatjet_msd) / self._data.genVMass
							vmatched = matching_dphi < 0.8 and matching_dpt < 0.5 and matching_dmass < 0.3

					self._selection_histograms[selection].GetTH2D("pt_dcsv").Fill(fatjet_pt, fatjet_dcsv, event_weight)
					self._selection_histograms[selection].GetTH1D("pfmet").Fill(self._data.pfmet, event_weight)
					self._selection_histograms[selection].GetTH1D("dcsv").Fill(fatjet_dcsv, event_weight)
					self._selection_histograms[selection].GetTH1D("n2ddt").Fill(fatjet_n2ddt, event_weight)
					self._selection_histograms[selection].GetTH1D("pt").Fill(fatjet_pt, event_weight)
					self._selection_histograms[selection].GetTH1D("msd").Fill(fatjet_msd, event_weight)
					self._selection_histograms[selection].GetTH1D("eta").Fill(fatjet_eta, event_weight)
					self._selection_histograms[selection].GetTH1D("rho").Fill(fatjet_rho, event_weight)
					self._selection_histograms[selection].GetTH2D("dcsv_vs_rho").Fill(fatjet_dcsv, fatjet_rho, event_weight)
					self._selection_histograms[selection].GetTH2D("dcsv_vs_pt").Fill(fatjet_dcsv, fatjet_pt, event_weight)
					self._selection_histograms[selection].GetTH2D("dcsv_vs_msd").Fill(fatjet_dcsv, fatjet_msd, event_weight)
					self._selection_histograms[selection].GetTH3D("dcsv_vs_msd_vs_pt").Fill(fatjet_dcsv, fatjet_msd, fatjet_pt, event_weight)

					if self._jet_type == "CA15":
						self._selection_histograms[selection].GetTH2D("dcsvalt_vs_rho").Fill(self._data.CA15Puppijet0_doublecsv, fatjet_rho, event_weight)
						self._selection_histograms[selection].GetTH2D("dcsvalt_vs_pt").Fill(self._data.CA15Puppijet0_doublecsv, fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH2D("dcsvalt_vs_msd").Fill(self._data.CA15Puppijet0_doublecsv, fatjet_msd, event_weight)
						self._selection_histograms[selection].GetTH3D("dcsvalt_vs_msd_vs_pt").Fill(self._data.CA15Puppijet0_doublecsv, fatjet_msd, fatjet_pt, event_weight)

					# Pass and fail histograms
					if pass_box:
						self._selection_histograms[selection].GetTH2D("pass").Fill(fatjet_msd, fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH2D("pass_unweighted").Fill(fatjet_msd, fatjet_pt)
						if self._data_source == "simulation":
							if vmatched:
								self._selection_histograms[selection].GetTH2D("pass_matched").Fill(fatjet_msd, fatjet_pt, event_weight)
							else:
								self._selection_histograms[selection].GetTH2D("pass_unmatched").Fill(fatjet_msd, fatjet_pt, event_weight)
						for systematic in self._weight_systematics[selection]:
							self._selection_histograms[selection].GetTH2D("pass_{}".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
							if self._data_source == "simulation":
								if vmatched:
									self._selection_histograms[selection].GetTH2D("pass_{}_matched".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
								else:
									self._selection_histograms[selection].GetTH2D("pass_{}_unmatched".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
						self._selection_histograms[selection].GetTH1D("pass_pfmet").Fill(self._data.pfmet, event_weight)
						self._selection_histograms[selection].GetTH1D("pass_dcsv").Fill(fatjet_dcsv, event_weight)
						self._selection_histograms[selection].GetTH1D("pass_n2ddt").Fill(fatjet_n2ddt, event_weight)
						self._selection_histograms[selection].GetTH1D("pass_pt").Fill(fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH1D("pass_msd").Fill(fatjet_msd, event_weight)
						self._selection_histograms[selection].GetTH1D("pass_eta").Fill(fatjet_eta, event_weight)
						self._selection_histograms[selection].GetTH1D("pass_rho").Fill(fatjet_rho, event_weight)

					elif fail_box:
						self._selection_histograms[selection].GetTH2D("fail").Fill(fatjet_msd, fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH2D("fail_unweighted").Fill(fatjet_msd, fatjet_pt)
						if self._data_source == "simulation":
							if vmatched:
								self._selection_histograms[selection].GetTH2D("fail_matched").Fill(fatjet_msd, fatjet_pt, event_weight)
							else:
								self._selection_histograms[selection].GetTH2D("fail_unmatched").Fill(fatjet_msd, fatjet_pt, event_weight)
						for systematic in self._weight_systematics[selection]:
							if not systematic in event_weight_syst:
								print "[debug] ERROR : Systematic {} is not in event_weight_syst.".format(systematic)
								print "[debug] Data source = " + self._data_source
								print event_weight_syst
								sys.exit(1)
							self._selection_histograms[selection].GetTH2D("fail_{}".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
							if self._data_source == "simulation":
								if vmatched:
									self._selection_histograms[selection].GetTH2D("fail_{}_matched".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
								else:
									self._selection_histograms[selection].GetTH2D("fail_{}_unmatched".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
						self._selection_histograms[selection].GetTH1D("fail_pfmet").Fill(self._data.pfmet, event_weight)
						self._selection_histograms[selection].GetTH1D("fail_dcsv").Fill(fatjet_dcsv, event_weight)
						self._selection_histograms[selection].GetTH1D("fail_n2ddt").Fill(fatjet_n2ddt, event_weight)
						self._selection_histograms[selection].GetTH1D("fail_pt").Fill(fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH1D("fail_msd").Fill(fatjet_msd, event_weight)
						self._selection_histograms[selection].GetTH1D("fail_eta").Fill(fatjet_eta, event_weight)
						self._selection_histograms[selection].GetTH1D("fail_rho").Fill(fatjet_rho, event_weight)
					if fatjet_dcsv > self._dcsv_cut_loose:
						self._selection_histograms[selection].GetTH2D("passloose").Fill(fatjet_msd, fatjet_pt, event_weight)
						if self._data_source == "simulation":
							if vmatched:
								self._selection_histograms[selection].GetTH2D("passloose_matched").Fill(fatjet_msd, fatjet_pt, event_weight)
							else:
								self._selection_histograms[selection].GetTH2D("passloose_unmatched").Fill(fatjet_msd, fatjet_pt, event_weight)

						for systematic in self._weight_systematics[selection]:

							self._selection_histograms[selection].GetTH2D("passloose_{}".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
							if self._data_source == "simulation":
								if vmatched:
									self._selection_histograms[selection].GetTH2D("passloose_{}_matched".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])
								else:
									self._selection_histograms[selection].GetTH2D("passloose_{}_unmatched".format(systematic)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])

					if self._do_optimization:
						for dcsv_cut in self._dcsv_cuts:
							if fatjet_dcsv > dcsv_cut:
								self._selection_histograms[selection].GetTH2D("pass_dcsv{}".format(dcsv_cut)).Fill(fatjet_msd, fatjet_pt, event_weight)
								self._selection_histograms[selection].GetTH2D("pass_unweighted_dcsv{}".format(dcsv_cut)).Fill(fatjet_msd, fatjet_pt)
							elif fatjet_dcsv > self._dcsv_min:
								self._selection_histograms[selection].GetTH2D("fail_dcsv{}".format(dcsv_cut)).Fill(fatjet_msd, fatjet_pt, event_weight)
								self._selection_histograms[selection].GetTH2D("fail_unweighted_dcsv{}".format(dcsv_cut)).Fill(fatjet_msd, fatjet_pt)



				# Run systematics that affect event selection
				for systematic in self._jet_systematics:
					self._event_selectors_syst[selection][systematic].process_event(self._data, event_weight)
					if self._event_selectors_syst[selection][systematic].event_pass():
						if self._jet_type == "AK8":
							if systematic == "JESUp":
								fatjet_pt_syst = self._data.AK8Puppijet0_pt_JESUp
							elif systematic == "JESDown":
								fatjet_pt_syst = self._data.AK8Puppijet0_pt_JESDown
							elif systematic == "JERUp":
								fatjet_pt_syst = self._data.AK8Puppijet0_pt_JERUp
							elif systematic == "JERDown":
								fatjet_pt_syst = self._data.AK8Puppijet0_pt_JERDown
							else:
								print "ERROR : Systematic not recognized: " + systematic
								sys.exit(1)
						elif self._jet_type == "CA15":
							if systematic == "JESUp":
								fatjet_pt_syst = self._data.CA15Puppijet0_pt_JESUp
							elif systematic == "JESDown":
								fatjet_pt_syst = self._data.CA15Puppijet0_pt_JESDown
							elif systematic == "JERUp":
								fatjet_pt_syst = self._data.CA15Puppijet0_pt_JERUp
							elif systematic == "JERDown":
								fatjet_pt_syst = self._data.CA15Puppijet0_pt_JERDown
							else:
								print "ERROR : Systematic not recognized: " + systematic
								sys.exit(1)

						# For simulated V(bb), match fat jet to parent truth particle
						if self._data_source == "simulation":
							vmatched = False
							if self._data.genVPt > 0 and self._data.genVMass > 0:
								matching_dphi = abs(math.acos(math.cos(self._data.genVPhi - fatjet_phi)))
								matching_dpt = abs(self._data.genVPt - fatjet_pt_syst) / self._data.genVPt
								matching_dmass = abs(self._data.genVMass - fatjet_msd) / self._data.genVMass
								vmatched = matching_dphi < 0.8 and matching_dpt < 0.5 and matching_dmass < 0.3

						if pass_box:
							self._selection_histograms[selection].GetTH2D("pass_{}".format(systematic)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
							if self._data_source == "simulation":
								if vmatched:
									self._selection_histograms[selection].GetTH2D("pass_{}_matched".format(systematic)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
								else:
									self._selection_histograms[selection].GetTH2D("pass_{}_unmatched".format(systematic)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
						elif fail_box:
							self._selection_histograms[selection].GetTH2D("fail_{}".format(systematic)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
							if self._data_source == "simulation":
								if vmatched:
									self._selection_histograms[selection].GetTH2D("fail_{}_matched".format(systematic)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
								else:
									self._selection_histograms[selection].GetTH2D("fail_{}_unmatched".format(systematic)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)

						if self._do_optimization:
							for dcsv_cut in self._dcsv_cuts:
								if fatjet_dcsv > dcsv_cut:
									self._selection_histograms[selection].GetTH2D("pass_{}_dcsv{}".format(systematic, dcsv_cut)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
									self._selection_histograms[selection].GetTH2D("pass_{}_unweighted_dcsv{}".format(systematic, dcsv_cut)).Fill(fatjet_msd, fatjet_pt_syst)
								elif fatjet_dcsv > self._dcsv_min:
									self._selection_histograms[selection].GetTH2D("fail_{}_dcsv{}".format(systematic, dcsv_cut)).Fill(fatjet_msd, fatjet_pt_syst, event_weight)
									self._selection_histograms[selection].GetTH2D("fail_{}_unweighted_dcsv{}".format(systematic, dcsv_cut)).Fill(fatjet_msd, fatjet_pt_syst)


	def finish(self):
		if self._output_path == "":
			self._output_path = "/uscms/home/dryu/DAZSLE/data/LimitSetting/InputHistograms_{}.root".format(time.time)
			print "[SignalCutflow::finish] WARNING : Output path was not provided! Saving to {}".format(self._output_path)
		print "[SignalCutflow::finish] INFO : Saving histograms to {}".format(self._output_path)
		f_out = ROOT.TFile(self._output_path, "RECREATE")
		self._histograms.SaveAll(f_out)
		for selection, histogrammer in self._selection_histograms.iteritems():
			histogrammer.SaveAll(f_out)
		for selection, selector in self._event_selectors.iteritems():
			selector.print_cutflow()
			selector.make_cutflow_histograms(f_out)
			selector.save_nminusone_histograms(f_out)
		f_out.Close()

# Not using this right now! It was intended for joblib parallel processing, but you weren't able to figure out joblib on condor. 
#def RunSingleFile(args, samples, filename, label="", tree_name="Events"):
#	print "Input file {}".format(filename)
#	limit_histogrammer = EventSelectionHistograms(label, tree_name=tree_name)
#	if args.output_folder:
#		limit_histogrammer.set_output_path("{}/InputHistograms_{}.root".format(args.output_folder, label))
#	else:
#		limit_histogrammer.set_output_path("/uscms/home/dryu/DAZSLE/data/LimitSetting/InputHistograms_{}.root".format(label))
#	limit_histogrammer.add_file(filename)
#	limit_histogrammer.start()
#	limit_histogrammer.run()
#	limit_histogrammer.finish()

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Produce and plot ieta-iphi histograms to look for buggy events')
	input_group = parser.add_mutually_exclusive_group() 
	input_group.add_argument('--all', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--all_lxplus', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--all_cmslpc', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--supersamples', type=str, help="Supersample name(s), comma separated. Must correspond to something in analysis_configuration.(background_names, signal_names, or data_names).")
	input_group.add_argument('--samples', type=str, help="Sample name(s), comma separated. Must be a key in analysis_configuration.skims.")
	input_group.add_argument('--files', type=str, help="Input file name(s), comma separated")
	parser.add_argument('--n_jobs', type=int, default=4, help="For --run, specify the number of parallel jobs.")
	action_group = parser.add_mutually_exclusive_group() 
	action_group.add_argument('--combine_outputs', action="store_true", help="Compile results into one file for next step (buildRhalphabet). Also applies luminosity weights to MC.")
	action_group.add_argument('--run', action="store_true", help="Run")
	action_group.add_argument('--condor_run', action="store_true", help="Run on condor")
	#action_group.add_argument('--rhalphabet', action="store_true", help="Run rhalpabet and create workspaces for combine")
	#action_group.add_argument('--datacards', action="store_true", help="Create datacards for combine")
	parser.add_argument('--output_folder', type=str, help="Output folder")
	parser.add_argument('--label', type=str, help="If running with --files, need to specify a label manually, in lieu of the sample names, for the output file naming.")
	parser.add_argument('--luminosity', type=float, default=35900, help="Luminosity in pb^-1")
	parser.add_argument('--jet_type', type=str, default="AK8", help="AK8 or CA15")
	parser.add_argument('--data_source', type=str, default="data", help="data or simulation")
	parser.add_argument('--skim_inputs', action='store_true', help="Run over skim inputs")
	parser.add_argument('--do_optimization', action='store_true', help="Make tau21DDT opt plots")
	#parser.add_argument('--prescale', type=int, help="'Prescale' the input events (using evtNum % ps == 0)") # Currently this is set from the sample name, i.e. JetHTRun2016B_ps10. 
	args = parser.parse_args()

	if args.run or args.condor_run:
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
		#print "List of sample => input files:",
		#print sample_files


	if args.run:
		#from joblib import Parallel
		#from joblib import delayed
		for sample in samples:
			print "\n *** Running sample {}".format(sample)
			if "Sbb" in sample or args.skim_inputs or "ZPrime" in sample:
				tree_name = "Events"
			else:
				tree_name = "otree"

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
				
			limit_histogrammer = EventSelectionHistograms(sample, tree_name=tree_name, jet_type=args.jet_type)
			if args.do_optimization:
				limit_histogrammer.do_optimization()
			output_file_basename ="InputHistograms_{}_{}.root".format(sample, args.jet_type) 
			if args.output_folder:
				limit_histogrammer.set_output_path("{}/{}".format(args.output_folder, output_file_basename))
			else:
				limit_histogrammer.set_output_path("/uscms/home/dryu/DAZSLE/data/LimitSetting/{}".format(output_file_basename))
			for filename in sample_files[sample]:
				print "Input file {}".format(filename)
				limit_histogrammer.add_file(filename)
			#limit_histogrammer.set_jet_type(args.jet_type)
			if "JetHTRun2016" in sample or "SingleMuRun2016" in sample:
				limit_histogrammer.set_data_source("data")
			else:
				limit_histogrammer.set_data_source("simulation")
			if "ps10" in sample:
				limit_histogrammer.set_prescale(10)
			limit_histogrammer.start()
			limit_histogrammer.run()
			limit_histogrammer.finish()
			#subjob_index = 0
			#for subjob_index, filename in enumerate(sample_files[sample]):
			#	RunSingleFile(args, samples, filename, label=sample)
			# joblib doesn't work on condor! Yet.
			#if len(sample_files[sample]) >= 2:
			#	Parallel(n_jobs=4)(
			#		delayed(RunSingleFile)(
			#			args, 
			#			samples, 
			#			sample_files[sample][i], 
			#			label="{}_subjob{}".format(sample, i)
			#		) for i in xrange(len(sample_files[sample]))
			#	)
			#else:
			#	RunSingleFile(args, samples, sample_files[sample][0], label=sample)

	if args.condor_run:
		import time
		hadd_scripts = []
		for sample in samples:
			start_directory = os.getcwd()
			job_tag = "job_{}_{}_{}".format(sample, args.jet_type, int(floor(time.time())))
			submission_directory = os.path.expandvars("$HOME/DAZSLE/data/LimitSetting/condor/{}".format(job_tag))
			os.system("mkdir -pv {}".format(submission_directory))
			os.chdir(submission_directory)

			files_per_job = 1
			if args.skim_inputs:
				if "JetHTRun2016" in sample:
					if "ps10" in sample:
						files_per_job = 100
					else:
						files_per_job = 30
				elif "SingleMuRun2016" in sample:
					if "ps10" in sample:
						files_per_job = 100
					else:
						files_per_job = 15
				elif "QCD_HT500to700" in sample:
					files_per_job = 30
				elif "QCD_HT700to1000" in sample:
					files_per_job = 5
				elif "QCD" in sample:
					files_per_job = 10
				elif "Spin0" in sample or "Sbb" in sample or "ZPrime" in sample:
					files_per_job = 3
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

			job_command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/event_selection_histograms.py --jet_type {} --files $this_input_files_string --label {}_csubjob$1 --output_folder . --run ".format(args.jet_type, sample)
			if args.skim_inputs or args.all_lxplus:
				job_command += " --skim_inputs "
			if args.do_optimization:
				job_command += " --do_optimization "

			job_command += " 2>&1\n"
			job_script.write(job_command)

			# Check if the output file exists
			job_script.write("if ls InputHistograms*{}_csubjob$1*root 1> /dev/null 2&>1; then\n".format(sample))
			job_script.write("\techo \"1\" > jobstatus_csubjob$1.txt\n")
			job_script.write("else\n")
			job_script.write("\techo\"0\" > jobstatus_csubjob$1.txt\n")
			job_script.write("fi\n")

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
			hadd_script.write("for f in jobstatus_csubjob*.txt; do\n")
			hadd_script.write("\tif grep -Fxq \"0\" $f; then\n")
			hadd_script.write("\t\techo \"Subjob failure in $f\"\n")
			hadd_script.write("\tfi\n")
			hadd_script.write("done\n")
			hadd_script.write(os.path.expandvars("hadd $HOME/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root {}/InputHistograms*csubjob*root\n".format(sample, args.jet_type, submission_directory)))
			hadd_script.close()
			os.chdir(start_directory)
		# One hadd script to rule them all
		master_hadd_script_path = os.path.expandvars("$HOME/DAZSLE/data/LimitSetting/condor/master_hadd_{}".format(args.jet_type))
		if not args.all:
			master_hadd_script_path += "_" + str(int(floor(time.time())))
		master_hadd_script_path += ".sh"
		master_hadd_script = open(master_hadd_script_path, "w")
		master_hadd_script.write("#!/bin/bash\n")
		for hadd_script_path in hadd_scripts:
			master_hadd_script.write("source " + hadd_script_path + "\n")
		master_hadd_script.close()		

	if args.combine_outputs:
		luminosity = args.luminosity # in pb^-1
		from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections
		systematics = {
			"SR":["JESUp", "JESDown", "JERUp", "JERDown", "TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"N2SR":["JESUp", "JESDown", "JERUp", "JERDown", "TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"N2CR":["JESUp", "JESDown", "JERUp", "JERDown", "TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"Preselection":["JESUp", "JESDown", "JERUp", "JERDown", "TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"muCR":["JESUp", "JESDown", "JERUp", "JERDown", "MuTriggerUp", "MuTriggerDown", "MuIDUp", "MuIDDown", "MuIsoUp", "MuIsoDown", "PUUp", "PUDown"]
		}
		selections = ["SR", "muCR", "Preselection", "N2CR", "N2SR"] # N2CR
		extra_vars = ["pfmet", "dcsv", "n2ddt", "pt", "eta", "rho"]
		selection_tau21s = {}
		selection_dcsvs = {}
		if args.do_optimization:
			for tau21_ddt_cut in [0.4, 0.45, 0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7]:
				for dcsv_cut in [0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]:
					selections.append("SR_tau21ddt{}_dcsv{}".format(tau21_ddt_cut, dcsv_cut))
					systematics["SR_tau21ddt{}_dcsv{}".format(tau21_ddt_cut, dcsv_cut)] = ["JESUp", "JESDown", "JERUp", "JERDown", "TriggerUp", "TriggerDown", "PUUp", "PUDown"]
					selection_tau21s["SR_tau21ddt{}_dcsv{}".format(tau21_ddt_cut, dcsv_cut)] = tau21_ddt_cut
					selection_dcsvs["SR_tau21ddt{}_dcsv{}".format(tau21_ddt_cut, dcsv_cut)] = dcsv_cut
		for selection in selections:
			if "SR" in selection:
				selection_prefix = "SR"
			if "N2CR" in selection:
				selection_prefix = "N2CR"
			if "N2SR" in selection:
				selection_prefix = "N2SR"
			elif "muCR" in selection:
				selection_prefix = "muCR"
			elif "Preselection" in selection:
				selection_prefix = "Preselection"
			output_file = ROOT.TFile("/uscms/home/dryu/DAZSLE/data/LimitSetting/histograms_{}_{}.root".format(selection, args.jet_type), "RECREATE")
			pass_histograms = {}
			pass_histograms_syst = {}
			fail_histograms = {}
			fail_histograms_syst = {}
			# data_obs, data_singlemu - not ready yet
			# "zll", "wlnu", "vvqq" - you need to find the cross sections, and split into appropriate samples
			if selection == "muCR":
				supersamples = ["data_obs", "data_singlemu", "data_obs_ps10", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq", "zll", "wlnu"]
			else:
				supersamples = ["data_obs", "data_singlemu", "data_obs_ps10", "qcd", "tqq", "wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125", "stqq", "vvqq"]
			supersamples.extend(config.simulated_signal_names)
			for supersample in supersamples:
				first = True
				pass_histograms_syst[supersample] = {}
				fail_histograms_syst[supersample] = {}
				use_Vmatched_histograms = (supersample in ["wqq", "zqq", "hqq125","tthqq125","vbfhqq125","whqq125","zhqq125"]) or ("Sbb" in supersample) or ("ZPrime" in supersample)
				use_loose_template = (supersample in ["wqq", "zqq"]) # Use looser DCSV cut for pass shape, to improve statistics
				if use_loose_template:
					pass_histograms_syst[supersample + "_normalization"] = {}

				for sample in config.samples[supersample]:
					input_histogram_filename = "/uscms/home/dryu/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root".format(sample, args.jet_type)
					print "Opening {}".format(input_histogram_filename)
					input_file = ROOT.TFile(input_histogram_filename, "READ")
					if selection in selection_tau21s:
						pass_histogram_name = "h_{}_tau21ddt{}_{}_pass_dcsv{}".format(selection_prefix, selection_tau21s[selection], args.jet_type, selection_dcsvs[selection])
						fail_histogram_name = "h_{}_tau21ddt{}_{}_fail_dcsv{}".format(selection_prefix, selection_tau21s[selection], args.jet_type, selection_dcsvs[selection])
						nevents_histogram_name = "h_{}_tau21ddt{}_{}_pass_nevents".format(selection_prefix, selection_tau21s[selection], args.jet_type)
					else:
						pass_histogram_name = "h_{}_{}_pass".format(selection_prefix, args.jet_type)
						fail_histogram_name = "h_{}_{}_fail".format(selection_prefix, args.jet_type)
						nevents_histogram_name = "h_{}_{}_pass_nevents".format(selection_prefix, args.jet_type)
						if use_Vmatched_histograms:
							pass_histogram_name += "_matched"
							fail_histogram_name += "_matched"
					if use_loose_template:
						pass_histogram_name_normalization = pass_histogram_name
						pass_histogram_name = pass_histogram_name.replace("pass", "passloose")
					this_pass_histogram = input_file.Get(pass_histogram_name)
					if not this_pass_histogram:
						print "[event_selection_histograms] ERROR : Couldn't find histogram {} in file {}".format(pass_histogram_name, input_file.GetPath())
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
							input_histogram_filename = "/uscms/home/dryu/DAZSLE/data/LimitSetting/InputHistograms_{}_{}.root".format(sample, args.jet_type)
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
								this_histogram_pass.Scale(lumi_sf)
								this_histogram_fail.Scale(lumi_sf)

							# Add up
							if first:
								first = False
								extra_histograms[var] = this_histogram.Clone()
								extra_histograms[var].SetDirectory(0)
								extra_histograms[var].SetName(supersample + "_" + var)
								extra_histograms_pass[var] = this_histogram_pass.Clone()
								extra_histograms_pass[var].SetDirectory(0)
								extra_histograms_pass[var].SetName(supersample + "_" + var + "_pass")
								extra_histograms_fail[var] = this_histogram_fail.Clone()
								extra_histograms_fail[var].SetDirectory(0)
								extra_histograms_fail[var].SetName(supersample + "_" + var + "_fail")
							else:
								extra_histograms[var].Add(this_histogram)
								extra_histograms_pass[var].Add(this_histogram_pass)
								extra_histograms_fail[var].Add(this_histogram_fail)
							input_file.Close()
						output_file.cd()
						extra_histograms[var].Write()
						extra_histograms_pass[var].Write()
						extra_histograms_fail[var].Write()
					# End loop over extra vars
				# End if SR or muCR
				
				# For matched histograms, also save the matched and unmatched histograms
			# End loop over supersamples
			output_file.Close()

	#if args.rhalphabet:
	#	top_directory = "/uscms/home/dryu/DAZSLE/data/LimitSetting/"
	#	os.system("mkdir -pv {}/combine".format(top_directory))
	#	rhalphabet_command = "python /uscms_data/d3/dryu/DAZSLE/CMSSW_7_4_7/src/DAZSLE/PhiBBPlusJet/fitting/PbbJet/buildRhalphabetPbb.py -i {}/hists_1D.root -b -o {}/combine --pseudo".format(top_directory, top_directory)
	#	print rhalphabet_command
	#	os.system(rhalphabet_command)

	#if args.datacards:
	#	os.system("python PbbJet/makeCardsPbb.py -i /uscms/home/dryu/DAZSLE/data/LimitSetting/combine -o /uscms/home/dryu/DAZSLE/data/LimitSetting/combine -b --pseudo")
