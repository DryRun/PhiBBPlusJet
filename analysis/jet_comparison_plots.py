import os
import sys
import ROOT
from DAZSLE.PhiBBPlusJet.analysis_base import AnalysisBase
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
import DAZSLE.PhiBBPlusJet.event_selections as event_selections
from DAZSLE.PhiBBPlusJet.bacon_event_selector import *
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
from DAZSLE.ZPrimePlusJet.cross_sections import cross_sections
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


h_num = None
h_den = None
for sample in ["QCD_HT2000toInf", "QCD_HT1500to2000", "QCD_HT100to200", "QCD_HT200to300", "QCD_HT300to500", "QCD_HT1000to1500", "QCD_HT700to1000", "QCD_HT500to700"]:
	f = TFile("/afs/cern.ch/user/d/dryu/DAZSLE/data/JetComparison/JetComparison_{}.root".format(sample), "READ")
	nevents = f.Get("h_input_nevents").Integral()
	this_h_num = f.Get("h_CA15_pt_msd_trigAK8")
	this_h_den = f.Get("h_CA15_pt_msd")
	this_h_num.Scale(cross_sections[sample] * 35700. / nevents)
	this_h_den.Scale(cross_sections[sample] * 35700. / nevents)
	if not h_num:
		h_num = this_h_num.Clone()
		h_num.SetDirectory(0)
		h_den = this_h_den.Clone()
		h_den.SetDirectory(0)
	else:
		h_num.Add(this_h_num)
		h_den.Add(this_h_den)

h_ratio = h_num.Clone()
h_ratio.Divide(h_num, h_den, 1., 1., "B")
c = TCanvas("c_jetcomparison_CA15trigAK8", "c", 800, 600)
h_ratio.Draw("colz")
c.SaveAs("/afs/cern.ch/user/d/dryu/DAZSLE/data/JetComparison/{}.pdf".format(c.GetName()))

