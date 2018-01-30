import ROOT
from ROOT import *
import os
import sys
import math
from array import array
#import functools
#from multiprocessing import Pool
#import scipy
#sys.path.append('/home/marc/code/python/')
#import PlottingFunctions
#import RootHelperFunctions
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
from DAZSLE.PhiBBPlusJet.cross_sections import cross_sections
os.path.append(".")
from DDTMaker import get_ddttransf_path

input_folder = "/afs/cern.ch/user/d/dryu/DAZSLE/data/DDT"
output_folder = "/afs/cern.ch/user/d/dryu/DAZSLE/data/DDT/tmp"

rho_bins = [35, -7., -0.]
pt_bins = [16, 200., 1000.]
z_bins = {
	"N2":[750, 0., 0.75],
	"dcsv":[100, -1., 1.],
	"dsub":[100, -1., 1.],
}

rho_bins_fine = [280, -7.0, 0.0]
pt_bins_fine = [200, 400., 1000.]

