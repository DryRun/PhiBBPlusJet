import os
import sys
import math
import cmath
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import seaborn as sns
darks = sns.color_palette("dark")
pastels = sns.color_palette("pastel")

# Constants
# - Masses
masses = {
	"u":2.3e-3,
	"d":4.8e-3,
	"s":95e-3,
	"c":1.275,
	"b":4.18,
	"t":173.2,
	"e":0.511e-3,
	"mu":0.1066e-3,
	"tau":1.777
}

#fermions = ["u", "d", "s", "c", "b", "t", "e", "mu", "tau"]
quarks = ["u", "d", "s", "c", "b", "t"]

# - Width formulas
# -- Width to SM fermions
def partial_width(qflavor, mzp, gq, v_type):
	if masses[qflavor] > mzp / 2.:
		return 0.
	else:
		beta = math.sqrt(1. - (4. * masses[qflavor]**2 / mzp**2))
		if v_type == "vector":
			return 3. * gq**2 * mzp / (12. * math.pi) * (1 + 2. * masses[qflavor]**2 / mzp**2) * beta
		elif v_type == "axialvector":
			return 3. * gq**2 * mzp / (12. * math.pi) * beta**3

def chi_width(mzp, mchi, gchi, v_type):
	if mchi > mzp / 2.:
		return 0.
	else:
		beta = math.sqrt(1. - (4. * mchi**2 / mzp**2))
		if v_type == "vector":
			return gchi**2 * mzp / (12. * pi) * (1. + 2. * mchi**2 / mzp**2) * beta
		elif v_type == "axialvector":
			return gchi**2 * mzp / (12. * pi) * beta**3

def total_width(mzp, gq, mchi, gchi, v_type):
	total_width = 0.
	for qflavor in quarks:
		total_width += partial_width(qflavor, mzp, gq, v_type)
	total_width += chi_width(mzp, mchi, gchi, v_type)
	return total_width

def branching_fraction(qflavor, mzp, gq, mchi, gchi, v_type):
	return partial_width(qflavor, mzp, gq, v_type) / total_width(mzp, gq, mchi, gchi, v_type)

def branching_fraction_gg(mzp, gq, mchi, gchi, v_type):
	return gg_width(mzp, gq, v_type) / total_width(mzp, gq, mchi, gchi, v_type)

def plot_all_types(gq, mass_range=[1., 800.], mchi=1.e10, gchi=0.):
	print "plot()"
	mzps = np.arange(mass_range[0], mass_range[1], 1.)

	branching_fractions_vector = {}
	for qflavor in quarks:
		branching_fractions_vector[qflavor] = np.array([branching_fraction(qflavor, mzp, gq, mchi, gchi, "vector") for mzp in mzps])
	branching_fractions_vector["udsc"] = branching_fractions_vector["u"] + branching_fractions_vector["d"] + branching_fractions_vector["s"] + branching_fractions_vector["c"]

	branching_fractions_axialvector = {}
	for qflavor in quarks:
		branching_fractions_axialvector[qflavor] = np.array([branching_fraction(qflavor, mzp, gq, mchi, gchi, "axialvector") for mzp in mzps])
	branching_fractions_axialvector["udsc"] = branching_fractions_axialvector["u"] + branching_fractions_axialvector["d"] + branching_fractions_axialvector["s"] + branching_fractions_axialvector["c"]


	fig, ax1 = plt.subplots()
	ax1.set_xlabel(r"$m_{Z'}$ [GeV]")

	ax1.set_yscale("linear")
	ax1.set_xlim(0., mass_range[1])
	ax1.set_ylim(0., 1.)
	ax1.set_ylabel("Branching fraction")

	plotobj_vector_t, = ax1.plot(mzps, branching_fractions_vector["t"], color=darks[2], linewidth=1.0, linestyle="-", label=r"$t$")
	plotobj_vector_b, = ax1.plot(mzps, branching_fractions_vector["b"], color=darks[0], linewidth=1.0, linestyle="-", label=r"$b$")
	plotobj_vector_udsc, = ax1.plot(mzps, branching_fractions_vector["udsc"], color=darks[5], linewidth=1.0, linestyle="-", label=r"$u/d/s/c$")
	vector_legend = plt.legend(handles=[plotobj_vector_t, plotobj_vector_b, plotobj_vector_udsc], loc=(0.1, 0.33), title="Vector")
	ax1.add_artist(vector_legend)

	plotobj_axialvector_t, = ax1.plot(mzps, branching_fractions_axialvector["t"], color=darks[2], linewidth=1.0, linestyle="--", label=r"$t$")
	plotobj_axialvector_b, = ax1.plot(mzps, branching_fractions_axialvector["b"], color=darks[0], linewidth=1.0, linestyle="--", label=r"$b$")
	plotobj_axialvector_udsc, = ax1.plot(mzps, branching_fractions_axialvector["udsc"], color=darks[5], linewidth=1.0, linestyle="--", label=r"$u/d/s/c$")
	axialvector_legend = plt.legend(handles=[plotobj_axialvector_t, plotobj_axialvector_b, plotobj_axialvector_udsc], loc=(0.3, 0.33), title="Axial vector")
	ax1.add_artist(axialvector_legend)

	plt.savefig(os.path.expandvars("$HOME/DAZSLE/data/Signal/figures/brs_zp_both.png"))

def plot(gq, v_type, mass_range=[1., 800.], mchi=1.e10, gchi=0.):
	print "plot()"
	mzps = np.arange(mass_range[0], mass_range[1], 1.)
	quark_widths = {}
	for qflavor in quarks:
		quark_widths[qflavor] = np.array([partial_width(qflavor, mzp, gq, v_type) for mzp in mzps])
	quark_widths["udsc"] = quark_widths["u"] + quark_widths["d"] + quark_widths["s"] + quark_widths["c"]
	width_total = np.array([total_width(mzp, gq, mchi, gchi, v_type) for mzp in mzps])

	branching_fractions = {}
	for qflavor in quarks:
		branching_fractions[qflavor] = np.array([branching_fraction(qflavor, mzp, gq, mchi, gchi, v_type) for mzp in mzps])
	branching_fractions["udsc"] = branching_fractions["u"] + branching_fractions["d"] + branching_fractions["s"] + branching_fractions["c"]


	fig, ax1 = plt.subplots()
	if v_type == "vector":
		ax1.set_xlabel(r"$m_{Z'_{v}}$ [GeV]")
	elif v_type == "axialvector":
		ax1.set_xlabel(r"$m_{Z'_{av}}$ [GeV]")


	ax1.set_yscale("linear")
	ax1.set_xlim(0., mass_range[1])
	ax1.set_ylim(0., 1.)
	ax1.set_ylabel("Branching fraction")
	ax1.plot(mzps, branching_fractions["t"], color=darks[2], linewidth=1.0, linestyle="-", label=r"$t$")
	ax1.plot(mzps, branching_fractions["b"], color=darks[0], linewidth=1.0, linestyle="-", label=r"$b$")
	ax1.plot(mzps, branching_fractions["udsc"], color=darks[5], linewidth=1.0, linestyle="-", label=r"$u/d/s/c$")
	ax1.legend(loc="center left")
	plt.savefig(os.path.expandvars("$HOME/DAZSLE/data/Signal/figures/brs_zp_{}.png".format(v_type)))

	ax2 = ax1.twinx()
	ax2.set_ylabel(r"$\Gamma$ [GeV]")
	ax2.set_yscale("log")
	ax2.set_xlim(0., mass_range[1])
	ax2.set_ylim(0.001, 100.)
	ax2.plot(mzps, width_total, color="black", linewidth=2.0, linestyle="--", label="Total")
	ax2.plot(mzps, quark_widths["t"], color=pastels[2], linewidth=1.0, linestyle="--", label=r"$t$")
	ax2.plot(mzps, quark_widths["b"], color=pastels[0], linewidth=1.0, linestyle="--", label=r"$b$")
	ax2.plot(mzps, quark_widths["udsc"], color=pastels[5], linewidth=1.0, linestyle="--", label=r"$u/d/s/c$")
	ax2.legend(loc="center right")

	#plt.show()
	plt.savefig(os.path.expandvars("$HOME/DAZSLE/data/Signal/figures/brsandwidths_zp_{}.png".format(v_type)))
	print "Done with plot()"

if __name__ == "__main__":
	#for i in xrange(1, 100):
	#	tau = 2. / 100. * i
	#	print "form_factor({}) = {}".format(tau, form_factor(tau, "vector"))
	#print branching_fraction("b", 300., 1., 1.e10, 0., "vector")
	print "Welcome to branching_fractions_zprime"
	plot_all_types(1.)
	#plot(1., "vector")
	#plot(1., "axialvector")