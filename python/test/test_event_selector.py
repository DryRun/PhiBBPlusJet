from DAZSLE.PhiBBPlusJet.event_selector import EventSelector, add_cut, add_nm1_hist

es1 = EventSelector("e")

@add_cut(es1)
@add_nm1_hist(es1, "pt", "p_{T} [GeV]", 100, 0., 1000.)
def min_pt(self, event):
  self._return_data["pt"] = event.pt
  return event.pt > 450.

es2 = EventSelector("e")

@add_cut(es2)
@add_nm1_hist(es2, "pt", "p_{T} [GeV]", 100, 0., 1000.)
def min_pt(self, event):
  self._return_data["pt"] = event.pt
  return event.pt > 500.

class Event:
  def __init__(self):
    self.pt = 0.

e1 = Event()
e1.pt = 100.
es1.min_pt(e1)
es2.min_pt(e1)

e2 = Event()
e2.pt = 1000.
es1.min_pt(e2)
es2.min_pt(e2)

e3 = Event()
e3.pt = 475.
es1.min_pt(e3)
es2.min_pt(e3)