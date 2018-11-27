#ifndef BaconEventCutFunctions_h
#define BaconEventCutFunctions_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TMath.h"
#include "TPython.h"

#include "MyTools/RootUtils/interface/Constants.h"
#include "MyTools/AnalysisTools/interface/EventSelector.h"
#include "DAZSLE/DAZSLECommon/interface/BaconData.h"

template class EventSelector<BaconData>;

namespace BaconEventCutFunctions {
	bool Min_dphi_mu_jet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_neleLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_neleLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_nmuLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_nmuLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_vmuoLoose0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_vmuoLoose0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_vmuoLoose0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_vmuoLoose0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_ntau(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_ntau(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_nphoLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_nphoLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool AK8Puppijet0_isTightVJet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_pfmet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_pfmet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_puppet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_puppet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_nAK4PuppijetsPt30dR08_0(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_nAK4PuppijetsMPt50dR08_0(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);

	bool Min_AK8Puppijet0_msd(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_msd(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8Puppijet0_msd_puppi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_msd_puppi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_AK8Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_AK8Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);

	bool CA15Puppijet0_isTightVJet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15CHSjet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15CHSjet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15CHSjet1_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15CHSjet1_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet1_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet1_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet1_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet1_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet1_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet1_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet1_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet1_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15CHSjet2_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15CHSjet2_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet2_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet2_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet2_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet2_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet2_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet2_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet2_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet2_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Min_CA15Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);
	bool Max_CA15Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector);

	void Configure(EventSelector<BaconData>* p_event_selector, TString name = "");

}

#endif