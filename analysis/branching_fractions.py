import os
import sys
import math
import cmath
import numpy as np
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
hvev = 246.
alpha_s = 0.118

# - Yukawa couplings
y = {}
for f in masses.keys():
	y[f] = masses[f] * hvev / math.sqrt(2)

fermions = ["u", "d", "s", "c", "b", "t", "e", "mu", "tau"]
quarks = ["u", "d", "s", "c", "b", "t"]

# - Width formulas
# -- Width to SM fermions
def partial_width(f, mphi, gqphi, phi_type):
	if masses[f] > mphi / 2.:
		return 0.
	else:
		if f in quarks:
			nc = 3.
		else:
			nc = 1.
		if phi_type == "scalar":
			exponent = 1.5
		elif phi_type == "pseudoscalar":
			exponent = 0.5
		return nc * y[f]**2 * gqphi**2 * mphi / (16 * math.pi) * (1. - (2. * masses[f] / mphi)**2)**exponent

def chi_width(mphi, mchi, gchi, phi_type):
	if mchi > mphi / 2.:
		return 0.
	else:
		return gchi**2 * mphi / (8. * math.pi) * (1. - (2. * mchi / mphi)**2)

def form_factor(tau, phi_type):
	if phi_type == "scalar":
		if tau == 1.:
			# Return lim(tau --> 1.)
			return 1.
		else:
			return tau * (1. + (1. - tau) * cmath.atan(1. / cmath.sqrt(tau - 1.))**2)
	elif phi_type == "pseudoscalar":
		if tau == 1.:
			return (math.pi/2.)**2
		else:
			return tau * cmath.atan(1. / cmath.sqrt(tau - 1.))**2

def gg_width(mphi, gqphi, phi_type):
	return alpha_s**2 * y["t"]**2 * gqphi**2 * mphi**3 / (32 * math.pi**3 * hvev**2) * abs(form_factor(4*masses["t"]**2/mphi**2, phi_type))**2

def total_width(mphi, gqphi, mchi, gchi, phi_type):
	total_width = 0.
	for f in fermions:
		total_width += partial_width(f, mphi, gqphi, phi_type)
	total_width += chi_width(mphi, mchi, gchi, phi_type)
	total_width += gg_width(mphi, gqphi, phi_type)
	return total_width

def branching_fraction(f, mphi, gqphi, mchi, gchi, phi_type):
	return partial_width(f, mphi, gqphi, phi_type) / total_width(mphi, gqphi, mchi, gchi, phi_type)

def branching_fraction_gg(mphi, gqphi, mchi, gchi, phi_type):
	return gg_width(mphi, gqphi, phi_type) / total_width(mphi, gqphi, mchi, gchi, phi_type)

def plot_all_types(gqphi, mass_range=[1., 500.], mchi=1.e10, gchi=0.):
	mphis = np.arange(mass_range[0], mass_range[1], 1.)
	width_scalar_gg = np.array([gg_width(mphi, gqphi, "scalar") for mphi in mphis])
	width_scalar_tt = np.array([partial_width("t", mphi, gqphi, "scalar") for mphi in mphis])
	width_scalar_bb = np.array([partial_width("b", mphi, gqphi, "scalar") for mphi in mphis])
	width_scalar_cc = np.array([partial_width("c", mphi, gqphi, "scalar") for mphi in mphis])
	width_scalar_total = np.array([total_width(mphi, gqphi, mchi, gchi, "scalar") for mphi in mphis])

	width_pseudoscalar_gg = np.array([gg_width(mphi, gqphi, "pseudoscalar") for mphi in mphis])
	width_pseudoscalar_tt = np.array([partial_width("t", mphi, gqphi, "pseudoscalar") for mphi in mphis])
	width_pseudoscalar_bb = np.array([partial_width("b", mphi, gqphi, "pseudoscalar") for mphi in mphis])
	width_pseudoscalar_cc = np.array([partial_width("c", mphi, gqphi, "pseudoscalar") for mphi in mphis])
	width_pseudoscalar_total = np.array([total_width(mphi, gqphi, mchi, gchi, "pseudoscalar") for mphi in mphis])

	br_scalar_gg = np.array([branching_fraction_gg(mphi, gqphi, mchi, gchi, "scalar") for mphi in mphis])
	br_scalar_tt = np.array([branching_fraction("t", mphi, gqphi, mchi, gchi, "scalar") for mphi in mphis])
	br_scalar_bb = np.array([branching_fraction("b", mphi, gqphi, mchi, gchi, "scalar") for mphi in mphis])
	br_scalar_cc = np.array([branching_fraction("c", mphi, gqphi, mchi, gchi, "scalar") for mphi in mphis])

	br_pseudoscalar_gg = np.array([branching_fraction_gg(mphi, gqphi, mchi, gchi, "pseudoscalar") for mphi in mphis])
	br_pseudoscalar_tt = np.array([branching_fraction("t", mphi, gqphi, mchi, gchi, "pseudoscalar") for mphi in mphis])
	br_pseudoscalar_bb = np.array([branching_fraction("b", mphi, gqphi, mchi, gchi, "pseudoscalar") for mphi in mphis])
	br_pseudoscalar_cc = np.array([branching_fraction("c", mphi, gqphi, mchi, gchi, "pseudoscalar") for mphi in mphis])

	fig, ax1 = plt.subplots()
	ax1.set_xlabel(r"$m_{\Phi/A}$ [GeV]")

	ax1.set_yscale("linear")
	ax1.set_ylim(0., 1.)
	ax1.set_ylabel("Branching fraction")

	# Plot scalar
	plotobj_br_scalar_tt, = ax1.plot(mphis, br_scalar_tt, color=darks[2], linewidth=1.0, linestyle="-", label=r"$Scalar t\bar{t}$")
	plotobj_br_scalar_bb, = ax1.plot(mphis, br_scalar_bb, color=darks[0], linewidth=1.0, linestyle="-", label=r"$Scalar b\bar{b}$")
	plotobj_br_scalar_cc, = ax1.plot(mphis, br_scalar_cc, color=darks[5], linewidth=1.0, linestyle="-", label=r"$Scalar c\bar{c}$")
	plotobj_br_scalar_gg, = ax1.plot(mphis, br_scalar_gg, color=darks[1], linewidth=1.0, linestyle="-", label=r"$Scalar gg$")
	scalar_legend = plt.legend(handles=[plotobj_br_scalar_tt, plotobj_br_scalar_bb, plotobj_br_scalar_cc, plotobj_br_scalar_gg, ], loc="center left")
	ax1.add_artist(scalar_legend)

	# Plot pseudoscalar
	plotobj_br_pseudoscalar_tt, = ax1.plot(mphis, br_pseudoscalar_tt, color=darks[2], linewidth=1.0, linestyle="--", label=r"$Pseudoscalar t\bar{t}$")
	plotobj_br_pseudoscalar_bb, = ax1.plot(mphis, br_pseudoscalar_bb, color=darks[0], linewidth=1.0, linestyle="--", label=r"$Pseudoscalar b\bar{b}$")
	plotobj_br_pseudoscalar_cc, = ax1.plot(mphis, br_pseudoscalar_cc, color=darks[5], linewidth=1.0, linestyle="--", label=r"$Pseudoscalar c\bar{c}$")
	plotobj_br_pseudoscalar_gg, = ax1.plot(mphis, br_pseudoscalar_gg, color=darks[1], linewidth=1.0, linestyle="--", label=r"$Pseudoscalar gg$")
	pseudoscalar_legend = plt.legend(handles=[plotobj_br_pseudoscalar_tt, plotobj_br_pseudoscalar_bb, plotobj_br_pseudoscalar_cc, plotobj_br_pseudoscalar_gg, ], loc="center right")
	ax1.add_artist(pseudoscalar_legend)

	plt.savefig(os.path.expandvars("$HOME/DAZSLE/data/Signal/figures/brs_spin0_both.png"))


def plot(gqphi, phi_type, mass_range=[1., 500.], mchi=1.e10, gchi=0.):
	mphis = np.arange(mass_range[0], mass_range[1], 1.)
	width_gg = np.array([gg_width(mphi, gqphi, phi_type) for mphi in mphis])
	width_tt = np.array([partial_width("t", mphi, gqphi, phi_type) for mphi in mphis])
	width_bb = np.array([partial_width("b", mphi, gqphi, phi_type) for mphi in mphis])
	width_cc = np.array([partial_width("c", mphi, gqphi, phi_type) for mphi in mphis])
	width_total = np.array([total_width(mphi, gqphi, mchi, gchi, phi_type) for mphi in mphis])

	br_gg = np.array([branching_fraction_gg(mphi, gqphi, mchi, gchi, phi_type) for mphi in mphis])
	br_tt = np.array([branching_fraction("t", mphi, gqphi, mchi, gchi, phi_type) for mphi in mphis])
	br_bb = np.array([branching_fraction("b", mphi, gqphi, mchi, gchi, phi_type) for mphi in mphis])
	br_cc = np.array([branching_fraction("c", mphi, gqphi, mchi, gchi, phi_type) for mphi in mphis])

	fig, ax1 = plt.subplots()
	if phi_type == "scalar":
		ax1.set_xlabel(r"$m_{\Phi}$ [GeV]")
	elif phi_type == "pseudoscalar":
		ax1.set_xlabel(r"$m_{A}$ [GeV]")


	ax1.set_yscale("linear")
	ax1.set_ylim(0., 1.)
	ax1.set_ylabel("Branching fraction")
	ax1.plot(mphis, br_tt, color=darks[2], linewidth=1.0, linestyle="-", label=r"$t\bar{t}$")
	ax1.plot(mphis, br_bb, color=darks[0], linewidth=1.0, linestyle="-", label=r"$b\bar{b}$")
	ax1.plot(mphis, br_cc, color=darks[5], linewidth=1.0, linestyle="-", label=r"$c\bar{c}$")
	ax1.plot(mphis, br_gg, color=darks[1], linewidth=1.0, linestyle="-", label=r"$gg$")
	ax1.legend(loc="center right")
	plt.savefig(os.path.expandvars("$HOME/DAZSLE/data/Signal/figures/brs_spin0_{}.png".format(phi_type)))

	ax2 = ax1.twinx()
	ax2.set_ylabel(r"$\Gamma$ [GeV]")
	ax2.set_yscale("log")
	ax2.set_xlim(0., 500.)
	ax2.set_ylim(0.001, 100.)
	ax2.plot(mphis, width_total, color="black", linewidth=2.0, linestyle="--", label="Total")
	ax2.plot(mphis, width_tt, color=pastels[2], linewidth=1.0, linestyle="--", label=r"$t\bar{t}$")
	ax2.plot(mphis, width_bb, color=pastels[0], linewidth=1.0, linestyle="--", label=r"$b\bar{b}$")
	ax2.plot(mphis, width_cc, color=pastels[5], linewidth=1.0, linestyle="--", label=r"$c\bar{c}$")
	ax2.plot(mphis, width_gg, color=pastels[1], linewidth=1.0, linestyle="--", label=r"$gg$")
	ax2.legend(loc="center right")

	#plt.show()
	plt.savefig(os.path.expandvars("$HOME/DAZSLE/data/Signal/figures/brsandwidths_{}.png".format(phi_type)))


if __name__ == "__main__":
	#for i in xrange(1, 100):
	#	tau = 2. / 100. * i
	#	print "form_factor({}) = {}".format(tau, form_factor(tau, "scalar"))
	#print branching_fraction("b", 300., 1., 1.e10, 0., "scalar")

	plot_all_types(1.)
	#plot(1., "scalar")
	#plot(1., "pseudoscalar")