import os
import sys
from ROOT import *

class LimitPlot():
	_limit_graph = None
	_limit_graphs_exp = {}
	_canvas = None

	def load_limit_graph(self, limit_graph):
		self._limit_graph = limit_graph

	def load_limit_graph_exp(self, sigma, limit_graph):
		self._limit_graphs_exp[sigma] = limit_graph

	def draw(self, name,
		canvas_x=800,
		canvas_y=600,
		legend_coords=[0.6, 0.6, 0.88, 0.8],
		legend_header=None,
		logx=False,
		logy=True,
		x_range=None,
		y_range=None,
		x_title=None,
		y_title=None,
		):
		self._canvas = TCanvas(name, name, canvas_x, canvas_y)
		self._legend = TLegend(legend_coords[0], legend_coords[1], legend_coords[2], legend_coords[3])
		self._legend.SetBorderSize(0)
		self._legend.SetFillColor(0)
		if legend_header:
			self._legend.SetHeader(legend_header)

		if x_range:
			x_min = x_range[0]
			x_max = x_range[1]
		else:
			x_min = min(self._limit_graph.GetX())
			x_max = max(self._limit_graph.GetX())
			if logx:
				x_min = x_min / 5.
				x_max = x_max * 5.
			else:
				delta_x = x_max - x_min
				x_min = x_min - 0.1 * delta_x
				x_max = x_max + 0.1 * delta_x
		if y_range:
			y_min = y_range[0]
			y_max = y_range[1]
		else:
			y_min = min(self._limit_graph.GetY())
			y_max = max(self._limit_graph.GetY())
			if logy:
				y_min = y_min / 10.
				y_max = y_max * 10.
			else:
				delta_y = y_max - y_min
				y_min = y_min - 0.1 * delta_y
				y_max = y_max + 0.1 * delta_y

		self._frame = TH1F("frame", "frame", 100, x_min, x_max)
		self._frame.SetMinimum(y_min)
		self._frame.SetMaximum(y_max)
		if x_title:
			self._frame.GetXaxis().SetTitle(x_title)
		if y_title:
			self._frame.GetYaxis().SetTitle(y_title)
		self._frame.Draw()

		self._limit_graph.SetLineColor(1)
		self._limit_graph.SetLineWidth(1)
		self._limit_graph.SetLineStyle(1)

		self._limit_graphs_exp[0].SetLineColor(0)
		self._limit_graphs_exp[0].SetLineWidth(1)
		self._limit_graphs_exp[0].SetLineStyle(2)

		# Make fills
		self._limit_fills_exp = {}
		for pmsigma in [1, 2]:
			self._limit_fills_exp[pmsigma] = TGraph(self._limit_graphs_exp[pmsigma].GetN() * 2)
			for i in xrange(self._limit_graphs_exp[pmsigma].GetN()):
				self._limit_fills_exp[pmsigma].SetPoint(i, self._limit_graphs_exp[-1*pmsigma].GetX()[i], self._limit_graphs_exp[-1*pmsigma].GetY()[i])

			for i in xrange(self._limit_graphs_exp[pmsigma].GetN()-1, -1, -1):
				j = self._limit_graphs_exp[pmsigma].GetN() * 2 - i - 1
				self._limit_fills_exp[pmsigma].SetPoint(j, self._limit_graphs_exp[pmsigma].GetX()[i], self._limit_graphs_exp[pmsigma].GetY()[i])

			if pmsigma == 1:
				self._limit_fills_exp[pmsigma].SetFillColor(kGreen)
				self._limit_fills_exp[pmsigma].SetFillStyle(1001)
			elif pmsigma == 2:
				self._limit_fills_exp[pmsigma].SetFillColor(kOrange)
				self._limit_fills_exp[pmsigma].SetFillStyle(1001)

		# Draw stuff
		self._limit_fills_exp[2].Draw("f")
		self._limit_fills_exp[1].Draw("f")
		self._limit_graphs_exp[0].Draw("l")
		self._limit_graph.Draw("l")

		# Legend
		self._legend.AddEntry(self._limit_graph, "Observed", "l")
		self._legend.AddEntry(self._limit_graphs_exp[0], "Expected", "l")
		self._legend.AddEntry(self._limit_fills_exp[1], "Expected 68%", "l")
		self._legend.AddEntry(self._limit_fills_exp[2], "Expected 95%", "l")
		self._legend.Draw()

	def save(self, directory, formats=["pdf"]):
		for file_format in formats:
			self._canvas.SaveAs("{}/{}.{}".format(directory, self._canvas.name, file_format))