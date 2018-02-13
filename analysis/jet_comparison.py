import os
import sys
import ROOT
from DAZSLE.PhiBBPlusJet.analysis_base import AnalysisBase
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
import DAZSLE.PhiBBPlusJet.event_selections as event_selections
from DAZSLE.PhiBBPlusJet.bacon_event_selector import *
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
from math import ceil, sqrt, floor, cos, acos
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

class JetComparison(AnalysisBase):
	def __init__(self, sample_name, tree_name="otree"):
		super(JetComparison, self).__init__(tree_name=tree_name)
		self._data = BaconData(self._chain)
		self._output_path = ""
		self._sample_name = sample_name
		self._input_nevents = 0
		self._n2_ddt_cut = 0.
		self._dcsv_cut = {"AK8":0.9, "CA15":0.9}
		self._n2_ddt_cut = 0.
		self._data_source = "data"
		jet_types = ["AK8", "CA15"]
		self._prescale = -1

	def set_data_source(self, data_source):
		print "[JetComparison::set_data_source] INFO : Setting data source to " + data_source
		if not data_source in ["data", "simulation"]:
			print "[JetComparison] ERROR : Data source must be data or simulation."
			sys.exit(1)
		self._data_source = data_source

	# Overload add_file to extract the number of input events to the skims, stored in histogram NEvents in the same file as the trees
	def add_file(self, filename):
		super(JetComparison, self).add_file(filename)
		f = ROOT.TFile.Open(filename, "READ")
		if f.Get("NEvents").Integral() == 0:
			print "[JetComparison::add_file] ERROR : NEvents.Integral() == 0 for file " + filename
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
		self._histograms.AddTH1D("AK8_pass_nevents", "pass_nevents", "", 1, -0.5, 0.5)
		self._histograms.AddTH1D("AK8_pass_nevents_weighted", "pass_nevents_weighted", "", 1, -0.5, 0.5)
		self._histograms.AddTH1D("AK8_pt", "AK8_pt", "AK8 p_{T} [GeV]", 100, 0., 1000.)
		self._histograms.AddTH1D("AK8_eta", "AK8_eta", "AK8 #eta", 100, -5., 5.)
		self._histograms.AddTH1D("AK8_msd", "AK8_msd", "AK8 m_{SD} [GeV]", 85, 5., 600.)
		self._histograms.AddTH1D("AK8_n2ddt", "AK8_n2ddt", "AK8 N_{2}^{DDT}", 40, -2., 2.)
		self._histograms.AddTH1D("AK8_dcsv", "AK8_dcsv", "AK8 dcsv", 200, -1., 1.)
		self._histograms.AddTH1D("AK8_rho", "AK8_rho", "AK8 #rho", 80, -8., 0.)

		self._histograms.AddTH1D("CA15_pass_nevents", "pass_nevents", "", 1, -0.5, 0.5)
		self._histograms.AddTH1D("CA15_pass_nevents_weighted", "pass_nevents_weighted", "", 1, -0.5, 0.5)
		self._histograms.AddTH1D("CA15_pt", "inclusive_pt", "CA15 p_{T} [GeV]", 100, 0., 1000.)
		self._histograms.AddTH1D("CA15_eta", "inclusive_eta", "CA15 #eta", 100, -5., 5.)
		self._histograms.AddTH1D("CA15_msd", "inclusive_msd", "CA15 m_{SD} [GeV]", 85, 5., 600.)
		self._histograms.AddTH1D("CA15_n2ddt", "inclusive_n2ddt", "N_{2}^{DDT}", 40, -2., 2.)
		self._histograms.AddTH1D("CA15_dcsv", "inclusive_dcsv", "CA15 dsub", 200, -1., 1.)
		self._histograms.AddTH1D("CA15_rho", "rho", "CA15 #rho", 80, -8., 0.)

		self._histograms.AddTH2D("mAK8_vs_mCA15", "mAK8_vs_mCA15", "AK8 m_{SD} [GeV]", 85, 5, 600, "CA15 m_{SD} [GeV]", 85, 5, 600)
		self._histograms.AddTH2D("mCA15_over_mAK8_vs_mAK8", "mCA15_over_mAK8_vs_mAK8", "m_{CA15}/m_{AK8}", 40, 0., 4., "AK8 m_{SD} [GeV]", 80, 40, 600)
		self._histograms.AddTH2D("ptAK8_vs_ptCA15", "ptAK8_vs_ptCA15", "AK8 p_{T} [GeV]", 100, 0., 1000., "CA15 p_{T} [GeV]", 100, 0., 1000.)

		self._histograms.AddTH2D("mAK8_vs_mCA15_dR1p5", "mAK8_vs_mCA15", "AK8 m_{SD} [GeV]", 85, 5, 600, "CA15 m_{SD} [GeV]", 85, 5, 600)
		self._histograms.AddTH2D("mCA15_over_mAK8_vs_mAK8_dR1p5", "mCA15_over_mAK8_vs_mAK8", "m_{CA15}/m_{AK8}", 40, 0., 4., "AK8 m_{SD} [GeV]", 80, 40, 600)
		self._histograms.AddTH2D("ptAK8_vs_ptCA15_dR1p5", "ptAK8_vs_ptCA15", "AK8 p_{T} [GeV]", 100, 0., 1000., "CA15 p_{T} [GeV]", 100, 0., 1000.)


		# Event selections
		self._event_selectors = {}
		self._event_selectors["AK8"] = 	BaconEventSelector("EventSelector_SR_AK8")
		self._event_selectors["AK8"].add_cut("Min_AK8Puppijet0_pt", {"Min_AK8Puppijet0_pt":450., "systematic":"nominal"})
		#self._event_selectors["AK8"].add_cut("Min_AK8Puppijet0_msd_puppi", 40.)
		self._event_selectors["AK8"].add_cut("AK8Puppijet0_isTightVJet")
		self._event_selectors["AK8"].add_cut("Max_neleLoose", 0)
		self._event_selectors["AK8"].add_cut("Max_nmuLoose", 0)
		self._event_selectors["AK8"].add_cut("Max_ntau", 0)
		self._event_selectors["AK8"].add_cut("Max_pfmet", {"Max_pfmet":140., "systematic":"nominal"})
		self._event_selectors["AK8"].add_cut("Max_AK8Puppijet0_N2DDT", 0.)
		self._event_selectors["AK8"].add_cut("Min_AK8Puppijet0_doublecsv", 0.9)

		self._event_selectors["CA15"] = BaconEventSelector("EventSelector_SR_CA15")
		self._event_selectors["CA15"].add_cut("Min_CA15Puppijet0_pt", {"Min_CA15Puppijet0_pt":450., "systematic":"nominal"})
		#self._event_selectors["CA15"].add_cut("Min_CA15Puppijet0_msd_puppi", 40.)
		self._event_selectors["CA15"].add_cut("CA15Puppijet0_isTightVJet")
		self._event_selectors["CA15"].add_cut("Max_neleLoose", 0)
		self._event_selectors["CA15"].add_cut("Max_nmuLoose", 0)
		self._event_selectors["CA15"].add_cut("Max_ntau", 0)
		self._event_selectors["CA15"].add_cut("Max_pfmet", {"Max_pfmet":140., "systematic":"nominal"})
		self._event_selectors["CA15"].add_cut("Max_CA15Puppijet0_N2DDT", 0.)
		self._event_selectors["CA15"].add_cut("Min_CA15Puppijet0_doublesub", 0.9)

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
		self._trig_eff = {}
		for jet_type in ["AK8", "CA15"]:
			if jet_type == "AK8":
				f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_AK8_SingleMuon_Run2016_V2p1_v03.root", "read")
				trig_den = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
				trig_num = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
			elif jet_type == "CA15":
				f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_CA15_SingleMuon_Run2016_V2p4_v08.root", "read")
				trig_den = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
				trig_num = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
			trig_den.RebinX(2)
			trig_num.RebinX(2)
			trig_den.RebinY(5)
			trig_num.RebinY(5)
			self._trig_eff[jet_type] = ROOT.TEfficiency()
			if (ROOT.TEfficiency.CheckConsistency(trig_num, trig_den)):
				self._trig_eff[jet_type] = ROOT.TEfficiency(trig_num, trig_den)
				self._trig_eff[jet_type].SetDirectory(0)
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

		print "[JetComparison::run] INFO : Running loop over tree from event {} to {}".format(first_event, limit_nevents - 1)

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

			# Weights: should be product of AK8*CA15, I think (event must pass both to be considered)
			event_weights = {}
			for jet_type  in ["AK8", "CA15"]:
				# Get weights
				if self._data_source == "data":
					event_weights[jet_type] = 1.
				else:
					if jet_type == "AK8":
						trigger_mass = min(self._data.AK8Puppijet0_msd, 300.)
						trigger_pt = max(200., min(self._data.AK8Puppijet0_pt, 1000.))
					elif jet_type == "CA15":
						trigger_mass = min(self._data.CA15Puppijet0_msd, 300.)
						trigger_pt = max(200., min(self._data.CA15Puppijet0_pt, 1000.))
					trigger_weight = self._trig_eff[jet_type].GetEfficiency(self._trig_eff[jet_type].FindFixBin(trigger_mass, trigger_pt))
					trigger_weight_up = trigger_weight + self._trig_eff[jet_type].GetEfficiencyErrorUp(self._trig_eff[jet_type].FindFixBin(trigger_mass, trigger_pt))
					trigger_weight_down = trigger_weight - self._trig_eff[jet_type].GetEfficiencyErrorLow(
						self._trig_eff[jet_type].FindFixBin(trigger_mass, trigger_pt))
					if trigger_weight <= 0 or trigger_weight_down <= 0 or trigger_weight_up <= 0:
						#print 'trigger_weights are %f, %f, %f, setting all to 1' % (trigger_weight, trigger_weight_up, trigger_weight_down)
						trigger_weight = 1
						trigger_weight_down = 1
						trigger_weight_up = 1
					event_weights[jet_type] = pu_weight * k_vjets * trigger_weight * k_ttbar
			event_weight = event_weights["AK8"] * event_weights["CA15"]

			# Pick up AK8 and CA15 event variables here, to avoid mistakes later
			AK8_pt = self._data.AK8Puppijet0_pt
			AK8_eta = self._data.AK8Puppijet0_eta
			AK8_msd = self._data.AK8Puppijet0_msd_puppi
			AK8_dcsv = self._data.AK8Puppijet0_doublecsv
			AK8_n2ddt = self._data.AK8Puppijet0_N2DDT
			AK8_rho = self._data.AK8Puppijet0_rho
			AK8_phi = self._data.AK8Puppijet0_phi

			CA15_pt = self._data.CA15Puppijet0_pt
			CA15_eta = self._data.CA15Puppijet0_eta
			CA15_msd = self._data.CA15Puppijet0_msd_puppi
			CA15_dcsv = self._data.CA15Puppijet0_doublesub
			CA15_n2ddt = self._data.CA15Puppijet0_N2DDT
			CA15_rho = self._data.CA15Puppijet0_rho
			CA15_phi = self._data.CA15Puppijet0_phi

			dR = sqrt((AK8_eta - CA15_eta)**2 + acos(cos(AK8_phi - CA15_phi))**2);

			event_pass = {}
			self._event_selectors["AK8"].process_event(self._data, event_weight)
			self._event_selectors["CA15"].process_event(self._data, event_weight)			
			event_pass["AK8"] = self._event_selectors["AK8"].event_pass()
			event_pass["CA15"] = self._event_selectors["CA15"].event_pass()

			if event_pass["AK8"]:
				self._histograms.GetTH1D("AK8_pass_nevents").Fill(0)
				self._histograms.GetTH1D("AK8_pass_nevents_weighted").Fill(0, event_weights["AK8"])
				self._histograms.GetTH1D("AK8_pt").Fill(AK8_pt, event_weights["AK8"])
				self._histograms.GetTH1D("AK8_eta").Fill(AK8_eta, event_weights["AK8"])
				self._histograms.GetTH1D("AK8_msd").Fill(AK8_msd, event_weights["AK8"])
				self._histograms.GetTH1D("AK8_n2ddt").Fill(AK8_n2ddt, event_weights["AK8"])
				self._histograms.GetTH1D("AK8_dcsv").Fill(AK8_dcsv, event_weights["AK8"])
				self._histograms.GetTH1D("AK8_rho").Fill(AK8_rho, event_weights["AK8"])
			if event_pass["CA15"]:
				self._histograms.GetTH1D("CA15_pass_nevents").Fill(0)
				self._histograms.GetTH1D("CA15_pass_nevents_weighted").Fill(0, event_weights["CA15"])
				self._histograms.GetTH1D("CA15_pt").Fill(CA15_pt, event_weights["CA15"])
				self._histograms.GetTH1D("CA15_eta").Fill(CA15_eta, event_weights["CA15"])
				self._histograms.GetTH1D("CA15_msd").Fill(CA15_msd, event_weights["CA15"])
				self._histograms.GetTH1D("CA15_n2ddt").Fill(CA15_n2ddt, event_weights["CA15"])
				self._histograms.GetTH1D("CA15_dcsv").Fill(CA15_dcsv, event_weights["CA15"])
				self._histograms.GetTH1D("CA15_rho").Fill(CA15_rho, event_weights["CA15"])

			if event_pass["AK8"] and event_pass["CA15"]:
				self._histograms.GetTH2D("mAK8_vs_mCA15").Fill(AK8_msd, CA15_msd, event_weights["AK8"] * event_weights["CA15"])
				self._histograms.GetTH2D("ptAK8_vs_ptCA15").Fill(AK8_pt, CA15_pt, event_weights["AK8"] * event_weights["CA15"])
				if AK8_msd != 0. and CA15_msd != 0.:
					self._histograms.GetTH2D("mCA15_over_mAK8_vs_mAK8").Fill(CA15_msd / AK8_msd, AK8_msd, event_weights["AK8"] * event_weights["CA15"])

				if dR < 1.5:
					self._histograms.GetTH2D("mAK8_vs_mCA15_dR1p5").Fill(AK8_msd, CA15_msd, event_weights["AK8"] * event_weights["CA15"])
					self._histograms.GetTH2D("ptAK8_vs_ptCA15_dR1p5").Fill(AK8_pt, CA15_pt, event_weights["AK8"] * event_weights["CA15"])
					if AK8_msd != 0. and CA15_msd != 0.:
						self._histograms.GetTH2D("mCA15_over_mAK8_vs_mAK8_dR1p5").Fill(CA15_msd / AK8_msd, AK8_msd, event_weights["AK8"] * event_weights["CA15"])

			# If events fail one or the other selection, fill that jet type with -1
			if event_pass["AK8"] and not event_pass["CA15"]:
				self._histograms.GetTH2D("mAK8_vs_mCA15").Fill(AK8_msd, -1, event_weights["AK8"] * event_weights["CA15"])
				self._histograms.GetTH2D("ptAK8_vs_ptCA15").Fill(AK8_pt, -1, event_weights["AK8"] * event_weights["CA15"])
			elif not event_pass["AK8"] and event_pass["CA15"]:
				self._histograms.GetTH2D("mAK8_vs_mCA15").Fill(-1, CA15_msd, event_weights["AK8"] * event_weights["CA15"])
				self._histograms.GetTH2D("ptAK8_vs_ptCA15").Fill(-1, CA15_pt, event_weights["AK8"] * event_weights["CA15"])
			if dR < 1.5:
				if event_pass["AK8"] and not event_pass["CA15"]:
					self._histograms.GetTH2D("mAK8_vs_mCA15_dR1p5").Fill(AK8_msd, -1, event_weights["AK8"] * event_weights["CA15"])
					self._histograms.GetTH2D("ptAK8_vs_ptCA15_dR1p5").Fill(AK8_pt, -1, event_weights["AK8"] * event_weights["CA15"])
				elif not event_pass["AK8"] and event_pass["CA15"]:
					self._histograms.GetTH2D("mAK8_vs_mCA15_dR1p5").Fill(-1, CA15_msd, event_weights["AK8"] * event_weights["CA15"])
					self._histograms.GetTH2D("ptAK8_vs_ptCA15_dR1p5").Fill(-1, CA15_pt, event_weights["AK8"] * event_weights["CA15"])

	def finish(self):
		if self._output_path == "":
			self._output_path = "/uscms/home/dryu/DAZSLE/data/LimitSetting/JetComparisons_{}.root".format(time.time)
			print "[SignalCutflow::finish] WARNING : Output path was not provided! Saving to {}".format(self._output_path)
		print "[SignalCutflow::finish] INFO : Saving histograms to {}".format(self._output_path)
		f_out = ROOT.TFile(self._output_path, "RECREATE")
		self._histograms.SaveAll(f_out)
		for selection, selector in self._event_selectors.iteritems():
			selector.print_cutflow()
			selector.make_cutflow_histograms(f_out)
			selector.save_nminusone_histograms(f_out)
		f_out.Close()

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Compare AK8 vs CA15')
	input_group = parser.add_mutually_exclusive_group() 
	input_group.add_argument('--all', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--all_lxplus', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--all_cmslpc', action="store_true", help="Run over all supersamples")
	input_group.add_argument('--supersamples', type=str, help="Supersample name(s), comma separated. Must correspond to something in analysis_configuration.(background_names, signal_names, or data_names).")
	input_group.add_argument('--samples', type=str, help="Sample name(s), comma separated. Must be a key in analysis_configuration.skims.")
	input_group.add_argument('--files', type=str, help="Input file name(s), comma separated")
	parser.add_argument('--max_nevents', type=int, help="Max nevents (for development)")
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
				
			jet_comparisoner = JetComparison(sample, tree_name=tree_name)
			if args.do_optimization:
				jet_comparisoner.do_optimization()
			output_file_basename ="JetComparison_{}.root".format(sample) 
			if args.output_folder:
				jet_comparisoner.set_output_path("{}/{}".format(args.output_folder, output_file_basename))
			else:
				jet_comparisoner.set_output_path("/afs/cern.ch/user/d/dryu/DAZSLE/data/JetComparison/{}".format(output_file_basename))
			for filename in sample_files[sample]:
				print "Input file {}".format(filename)
				jet_comparisoner.add_file(filename)
			if "JetHTRun2016" in sample or "SingleMuRun2016" in sample:
				jet_comparisoner.set_data_source("data")
			else:
				jet_comparisoner.set_data_source("simulation")
			if "ps10" in sample:
				jet_comparisoner.set_prescale(10)
			jet_comparisoner.start()
			if args.max_nevents:
				jet_comparisoner.run(max_nevents=args.max_nevents)
			else:
				jet_comparisoner.run()
			jet_comparisoner.finish()

	if args.condor_run:
		import time
		hadd_scripts = []
		for sample in samples:
			start_directory = os.getcwd()
			job_tag = "job_{}_{}".format(sample, int(floor(time.time())))
			submission_directory = os.path.expandvars("$HOME/DAZSLE/data/JetComparison/condor/{}".format(job_tag))
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
					files_per_job = 5
				elif "QCD_HT700to1000" in sample:
					files_per_job = 5
				elif "QCD_HT1000to1500" in sample:
					files_per_job = 1
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

			job_command = "python $CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/analysis/jet_comparison.py --files $this_input_files_string --label {}_csubjob$1 --output_folder . --run ".format(sample)
			if args.skim_inputs or args.all_lxplus:
				job_command += " --skim_inputs "
			if args.do_optimization:
				job_command += " --do_optimization "

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
			hadd_script.write(os.path.expandvars("hadd $HOME/DAZSLE/data/LimitSetting/InputHistograms_{}.root {}/InputHistograms*csubjob*root\n".format(sample, submission_directory)))
			hadd_script.close()
			os.chdir(start_directory)
		# One hadd script to rule them all
		master_hadd_script_path = os.path.expandvars("$HOME/DAZSLE/data/LimitSetting/condor/master_hadd")
		if not args.all:
			master_hadd_script_path += "_" + str(int(floor(time.time())))
		master_hadd_script_path += ".sh"
		master_hadd_script = open(master_hadd_script_path, "w")
		master_hadd_script.write("#!/bin/bash\n")
		for hadd_script_path in hadd_scripts:
			master_hadd_script.write("source " + hadd_script_path + "\n")
		master_hadd_script.close()		
