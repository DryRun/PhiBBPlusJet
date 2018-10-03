import os
import sys
import ROOT
from DAZSLE.PhiBBPlusJet.combine_helper.roofit_hist_container import RoofitHistContainer

# Stores roofit objects 
class RhalphabetRegion():
	def __init__(self, name):
		self._region_name = name
		self._background_names = []
		self._backgrounds_pass = {}
		self._backgrounds_fail = {}
		self._signal_name = None
		self._signal_pass = None
		self._signal_fail = None
		self._data_pass = None
		self._data_fail = None
		self._xvar = None
		self._systematics_norm = {}

		self._save_vars = [] # List of RooFit variables, to ensure they don't get garbage collected early 

	def name(self):
		return self._region_name

	def add_xvar(self, xvar):
		print "[Region::add_xvar] DEBUG : Adding xvar: ",
		print xvar
		self._xvar = xvar

	def add_simple_background(self, bkgd_name, pass_hist, fail_hist, normalization_var=None):
		if not self._xvar:
			print "[Region::add_background] ERROR : Need to add an independent variable with add_xvar(xvar) first."
			sys.exit(1)
		if hist.Integral() == 0:
			print "[Region::add_background] WARNING : Histogram from background {}, region {} has zero norm. Ignoring.".format(bkgd_name, self._region_name)
			return
		self._background_names.append(bkgd_name)
		if pass_hist:
			self._backgrounds_pass[bkgd_name] = RoofitHistContainer(pass_hist, self._xvar, name="{}_pass_{}".format(bkgd_name, self._region_name), normalization_var=normalization_var)
		if fail_hist:
			self._backgrounds_fail[bkgd_name] = RoofitHistContainer(fail_hist, self._xvar, name="{}_fail_{}".format(bkgd_name, self._region_name), normalization_var=normalization_var)

	def add_signal(self, signal_name, pass_hist, fail_hist, normalization_var=None):
		if not self._xvar:
			print "[Region::add_signal] ERROR : Need to add an independent variable with add_xvar(xvar) first."
			sys.exit(1)
		self._signal_name = signal_name
		if pass_hist:
			self._signal_pass = RoofitHistContainer(pass_hist, self._xvar, name="{}_{}".format(signal_name, self._region_name), normalization_var=normalization_var)
		else:
			print "[rhalphabet_region::add_signal] ERROR : pass_hist=None, but this is required."
			sys.exit(1)
		if fail_hist:
			self._signal_fail = RoofitHistContainer(fail_hist, self._xvar, name="{}_{}".format(signal_name, self._region_name), normalization_var=normalization_var)

	def add_data(self, data_name, pass_hist, fail_hist):
		if not self._xvar:
			print "[Region::add_data] ERROR : Need to add an independent variable with add_xvar(xvar) first."
			sys.exit(1)
		self._data_pass = RoofitHistContainer(pass_hist, self._xvar, name="{}_{}".format(data_name, self._region_name))
		self._data_fail = RoofitHistContainer(fail_hist, self._xvar, name="{}_{}".format(data_name, self._region_name))

	def add_syst_norm(self, name, syst_dict):
		self._systematics_norm[name] = syst_dict

	def write(self, datacard_path, ws_path):
		# Create datacard
		with open(datacard_path, 'w') as datacard:
			datacard.write("imax 1 number of channels\n")
			datacard.write("jmax * number of processes minus 1\n")
			datacard.write("kmax * number of nuisance parameters\n")
			datacard.write("---------------\n")
			# shapes * bin_name file wsname:histogramname($PROCESS, $SYSTEMATIC)
			datacard.write("shapes * bin0 {} {}:$PROCESS_{}\n".format(os.path.basename(ws_path), ws_name, self._region_name))
			datacard.write("---------------\n")
			datacard.write("bin bin0\n")
			datacard.write("observation -1.0\n")
			datacard.write("---------------\n")
			datacard.write("bin {}\n".format(" bin0"*(len(self._background_names)+1)))
			datacard.write("process {} {}\n".format(self._signal_name, " ".join(self._background_names)))

			all_processes = [self._signal_name]
			all_processes.extend(self._background_names)

			process_index_string = "process "
			for i in xrange(len(all_processes)):
				process_index_string += " " + str(i)
			datacard.write(process_index_string + "\n")
			rate_string = "rate "
			for i in all_processes:
				rate_string += " -1"
			datacard.write(rate_string + "\n")
			datacard.write("---------------\n")

			for systematic in self._systematics_norm.iteritems():
				this_syst_string = systematic + " "
				for process in all_processes:
					if process in self._systematics_norm[systematic]:
						this_syst_string += " {}".format(self._systematics_norm[systematic][process])
					else:
						this_syst_string += " -"
				datacard.write(this_syst_string + "\n")

		# Create workspace
		ws_name = "ws_pass_{}".format(self._region_name)
		ws_pass = ROOT.RooWorkspace(ws_name)
		for background_name in self._background_names:
			self._backgrounds_pass[background_name].RooDataHist().SetName("{}_pass_{}".format(background_name, self._region_name))
			getattr(ws, "import")(self._backgrounds_pass[background_name].RooDataHist())

		self._signal_pass.RooDataHist().SetName("{}_pass_{}".format(self._signal_name, self._region_name))
		getattr(ws, "import")(self._signal_pass.RooDataHist())

		self._data_pass.RooDataHist().SetName("{}_pass_{}".format("data_obs", self._region_name))
		getattr(ws, "import")(self._data_pass.RooDataHist())

		self._rhbackground_pass.SetName("{}_pass_{}".format(self._rhbackground_name, self._region_name))
		ws_pass.writeToFile(ws_path)		

		ws_name = "ws_fail_{}".format(self._region_name)
		ws_fail = ROOT.RooWorkspace(ws_name)
		for background_name in self._background_names:
			self._backgrounds_fail[background_name].RooDataHist().SetName("{}_fail_{}".format(background_name, self._region_name))
			getattr(ws, "import")(self._backgrounds_fail[background_name].RooDataHist())

		self._signal_fail.RooDataHist().SetName("{}_fail_{}".format(self._signal_name, self._region_name))
		getattr(ws, "import")(self._signal_fail.RooDataHist())

		self._data_fail.RooDataHist().SetName("{}_fail_{}".format("data_obs", self._region_name))
		getattr(ws, "import")(self._data_fail.RooDataHist())

		self._rhbackground_fail.SetName("{}_fail_{}".format(self._rhbackground_name, self._region_name))
		ws_fail.writeToFile(ws_path)		


	# y_value = y value to use in polynomial
	def add_rhbackground(self, rhbackground_name, y_value, n_x=6, n_y=5, x_range=[-7., 0.], y_range=[400., 1000.]):
		self._rhbackground_name = rhbackground_name

		# Make polynomial
		# - Coefficients
		self._rh_coefficients = []
		inital_tf = self._data_pass.Integral() / self._data_fail.Integral()
		initial_coeff_value = inital_tf / (n_x * n_y) # Initial values of Bernstein polynomial coefficients corresponding to f(x)=average_tf
		for i in xrange(n_x+1):
			self._rh_coefficients.append([])
			for j in xrange(n_y+1):
				coeff_name = "poly_coeff_x{}_y{}".format(i, j)
				self._rh_coefficients[i].append(RooRealVar(coeff_name, coeff_name, initial_coeff_value, -30., 30.))

		# Make background templates
		# Fail: each bin gets a (unconstrained) RooRealVar for the background prediction
		# Pass: RooFormulaVar, pass=fail*TF(rho, pt)
		self._rh_fail_vars = RooArgList()
		self._rh_pass_vars = RooArgList()

		for xbin in xrange(1, self._data_fail.GetNbinsX() + 1):
			# Independent x variable: rho scaled to the interval [0, 1]
			poly_xvar = RooConstVar(
				"poly_x_xbin{}_{}".format(self._xbin, self._region_name),
				"poly_x_xbin{}_{}".format(self._xbin, self._region_name),
				(self._data_fail.GetXaxis().GetBinCenter(xbin) - x_range[0]) / (x_range[1] - x_range[0])
			)

			# Independent y variable: pT scaled to the interval [0, 1]
			poly_yvar = RooConstVar(
				"poly_y_xbin{}_{}".format(self._xbin, self._region_name), 
				"poly_y_xbin{}_{}".format(self._xbin, self._region_name), 
				(y_value - y_range[0]) / (y_range[1] - y_range[0])
			)

			# Polynomial strings
			polystr_y = self.generate_bernstein_string(n_y)
			polystr_x = self.generate_bernstein_string(n_x)
			
			# Make polynomials
			# - First, make a list of polynomials as a function of y (the x polynomial will be built from these)
			xpoly_variable_list = RooArgList() # List corresponds to the [@0, @1, ...] in the polynomial string
			for i in xrange(n_x+1):
				ypoly_variable_list = RooArgList() # List corresponds to the [@0, @1, ...] in the polynomial string
				for j in xrange(j_y+1):
					ypoly_variable_list.add(self._rh_coefficients[i][j])
				ypoly_variable_list.add(poly_yvar)
				y_poly_name = "subpolynomial_y_xbin{}_{}".format(xbin, self._region_name)
				y_polynomial = RooFormulaVar(y_poly_name, y_poly_name, polystr_y, ypoly_variable_list)
				xpoly_variable_list.add(y_polynomial)
				self._save_vars.append(y_polynomial)

			# - Make x polynomial out of y polynomials
			xpoly_variable_list.add(poly_xvar)
			x_poly_name = "polynomial_xbin{}_{}".format(xbin, self._region_name)
			x_polynomial = RooFormulaVar(x_poly_name, x_poly_name, polystr_x, xpoly_variable_list)
			self._save_vars.append(x_polynomial)


			# Make variable corresponding to background counts in pass and fail regions
			# - Fail is just a RooRealVar (unconstrained). Initialize to the data counts minus the MC backgrounds.
			fail_bin_content = self._data_fail.TH1().GetBinContent(xbin) - sum([h.TH1().GetBinContent(xbin) for h in self._backgrounds_fail])
			fail_bin_content = max(0., fail_bin_content) # Ensure non-negative
			fail_bin_unc = math.sqrt(fail_bin_content) * 50. + 10.
			fail_var_name = "{}_fail_{}_xbin{}".format(rhbackground_name, self._region_name, xbin)
			fail_var = RooRealVar(fail_var_name, fail_var_name, fail_bin_content, 0., fail_bin_content + fail_bin_unc)
			self._save_vars.append(fail_var)

			# - Pass is a RooFormulaVar representing pass=fail*TF(rho, pt)
			pass_var_name = "{}_pass_{}_xbin{}".format(rhbackground_name, self._region_name, xbin)
			# Safety check: if the fail bin has too few events, set fail=pass=0 to constant vars
			if fail_bin_content == 0:
				fail_bin_var.setConstant(True)
				pass_var = RooRealVar(pass_var_name, pass_var_name, 0., 0., 0.)
				pass_var.setConstant(True)
			else:
				pass_var = RooFormulaVar(pass_var_name, pass_var_name, "@0*max(@1, 0.000001)", RooArgList(fail_var, x_polynomial))
			self._save_vars.append(pass_var)

			self._rh_pass_vars.add(pass_var)
			self._rh_fail_vars.add(fail_var)

		# End loop over x bins
		 
		# Make histogram objects out of the individual bin vars
		# - Note the class, RooParametricHist, is a combine class, not RooFit.
		# - Constructor:  RooParametricHist (const char *name, const char *title, RooAbsReal& _x, RooArgList& _pars, const TH1 &_shape);
		# - I think arg _shape is just used for the x-axis bins
		self._rhbackground_pass = RooParametricHist(
			"{}_pass_{}".format(rhbackground_name, self._region_name), 
			"{}_pass_{}".format(rhbackground_name, self._region_name), 
			self._xvar, 
			self._rh_pass_vars, 
			self._data_fail.TH1()
		)
		self._rhbackground_fail = RooParametricHist(
			"{}_fail_{}".format(rhbackground_name, self._region_name), 
			"{}_fail_{}".format(rhbackground_name, self._region_name), 
			self._xvar, 
			self._rh_fail_vars, 
			self._data_fail.TH1()
		)

		self._rh_pass_norm = RooAddition(rhalph_bkgd_name + "_pass_" + category + "_norm",
								  rhalph_bkgd_name + "_pass_" + category + "_norm", pass_bins)
		self._rh_fail_norm = RooAddition(rhalph_bkgd_name + "_fail_" + category + "_norm",
								  rhalph_bkgd_name + "_fail_" + category + "_norm", fail_bins)
