import os
import sys
from ROOT import *
from cutflow import Cutflow

# Specialization of Cutflow to event selection. 
class EventSelector(Cutflow):
	def __init__(self, name="UnnamedEventSelector", event=None):
		super(EventSelector, self).__init__(name)
		self._event = event
		self._object_selectors = {}
		self._event_pass = False
		self._event_pass_nminusone = {}
		self._cut_results = {}
		self._object_name = "Event"

	def add_object_selector(self, object_name, object_selector):
		self._object_selectors[object_name] = object_selector

	def get_object_pass(self, object_name, i):
		return self._object_selectors[object_name].get_object_pass(i)

	def event(self):
		return self._event

	# Set pointer to the event data. If the event instance is constant, better to set this in the constructor; this method enables changing the pointer, if it happens to be reallocated.
	def set_event(self, event):
		self._event = event

	def process_event(self, weight=1):
		self.reset()
		self._pass_calls += 1
		self._pass_calls_weighted += weight

		for object_selector in self._object_selectors:
			object_selector.process_event(weight)

		# Run cuts
		this_pass = True
		for cut in self._cut_list:
			self._cut_results[cut] = getattr(self, cut)() # Or self._cut_functions[cut](self._event)?
			if not self._cut_results[cut]:
				if this_pass:
					this_pass = False
					self._cutflow_counter[cut] += 1
					self._cutflow_counter_weighted[cut] += weight
				self._cut_counter[cut] += 1
				self._cut_counter_weighted[cut] += weight
			if this_pass:
				self._pass_counter[cut] += 1
				self._pass_counter_weighted[cut] += weight
		self._event_pass = this_pass

		# N-1 histograms
		for cut in self._cut_list:
			# Calculate N-1
			self._event_pass_nminusone[cut] = True
			for cut2 in self._cut_list:
				if cut == cut2:
					continue
				if not self._cut_results[cut2]:
					self._event_pass_nminusone[cut] = False
					break

			# Fill histograms
			if self._event_pass_nminusone[cut]:
				if cut in self._nminusone_histograms:
					for variable_name, histogram in self._nminusone_histograms[cut].iteritems():
						histogram.Fill(self._return_data[variable_name])

	def event_pass(self):
		return self._event_pass

	def event_pass_nminusone(self, cut_name):
		return self._event_pass_nminusone[cut_name]

	def event_passes_cut(self, cut_name):
		return self._cut_results[cut_name]

	def reset(self):
		self._event = None
		self._cut_results.clear()
		self._event_pass = False
		self._event_pass_nminusone.clear()

	def add_cut(self, cut_name, cut_logic, return_data=None):
		self._cut_list.append(cut_name)

		# Create the cut function
		exec("""
def {}_{}(obj, event):
	{}
""".format(self._name, cut_name, cut_logic))

		# Attach cut function to this instance
		exec("self.{} = MethodType({}, self, EventSelector".format(cut_name, ))
