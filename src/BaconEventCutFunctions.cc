#ifndef BaconEventCutFunctions_cxx
#define BaconEventCutFunctions_cxx

#include "DAZSLE/PhiBBPlusJet/interface/BaconEventCutFunctions.h"

namespace BaconEventCutFunctions {
	bool Min_dphi_mu_jet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double phi_mu = p_data.vmuoLoose0_phi;
		double phi_jet = 0.;
		if (p_event_selector->GetCutDescriptors("Min_dphi_mu_jet")[0].EqualTo("AK8")) {
			phi_jet = p_data.AK8Puppijet0_phi;
		} else if (p_event_selector->GetCutDescriptors("Min_dphi_mu_jet")[0].EqualTo("CA15")) {
			phi_jet = p_data.CA15Puppijet0_phi;
		} else {
			std::cerr << "[BaconEventCutFunctions::Min_dphi_mu_jet] ERROR : Cut descriptor [0] must be AK8 or CA15." << std::endl;
			exit(1);
		}
		double delta_phi = TMath::ACos(TMath::Cos(phi_mu - phi_jet));
		p_event_selector->SetReturnData("", delta_phi);
		return (delta_phi > p_event_selector->GetCutParameters("Min_dphi_mu_jet")[0]);
	}

	bool Min_AK8Puppijet0_msd(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8Puppijet0_msd", p_data.AK8Puppijet0_msd);
		return (p_data.AK8Puppijet0_msd >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_msd")[0]);
	}

	bool Max_AK8Puppijet0_msd(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8Puppijet0_msd", p_data.AK8Puppijet0_msd);
		return (p_data.AK8Puppijet0_msd <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_msd")[0]);
	}

	bool Min_AK8Puppijet0_msd_puppi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8Puppijet0_msd_puppi", p_data.AK8Puppijet0_msd_puppi);
		return (p_data.AK8Puppijet0_msd_puppi >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_msd_puppi")[0]);
	}

	bool Max_AK8Puppijet0_msd_puppi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8Puppijet0_msd_puppi", p_data.AK8Puppijet0_msd_puppi);
		return (p_data.AK8Puppijet0_msd_puppi <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_msd_puppi")[0]);
	}

	bool Min_AK8Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double jet_pt = 0.;
		if (p_event_selector->GetCutDescriptors("Min_AK8Puppijet0_pt").size() > 0) {
			TString systematic = (p_event_selector->GetCutDescriptors("Min_AK8Puppijet0_pt")[0]);
			systematic.ToLower();
			if (systematic.EqualTo("nominal")) {
				jet_pt = p_data.AK8Puppijet0_pt;
			} else if (systematic.EqualTo("jesup")) {
				jet_pt = p_data.AK8Puppijet0_pt_JESUp;
			} else if (systematic.EqualTo("jesdown")) {
				jet_pt = p_data.AK8Puppijet0_pt_JESDown;
			} else if (systematic.EqualTo("jerup")) {
				jet_pt = p_data.AK8Puppijet0_pt_JERUp;
			} else if (systematic.EqualTo("jerdown")) {
				jet_pt = p_data.AK8Puppijet0_pt_JERDown;
			} else {
				std::cerr << "[Min_AK8Puppijet0_pt] ERROR : Systematic " << systematic << " not known." << std::endl;
				exit(1);
			}
		} else {
			jet_pt = p_data.AK8Puppijet0_pt;
		}
		p_event_selector->SetReturnData("Min_AK8Puppijet0_pt", jet_pt);
		return (jet_pt >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_pt")[0]);
	}

	bool Max_AK8Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double jet_pt = 0.;
		if (p_event_selector->GetCutDescriptors("Max_AK8Puppijet0_pt").size() > 0) {
			TString systematic = (p_event_selector->GetCutDescriptors("Max_AK8Puppijet0_pt")[0]);
			systematic.ToLower();
			if (systematic.EqualTo("nominal")) {
				jet_pt = p_data.AK8Puppijet0_pt;
			} else if (systematic.EqualTo("jesup")) {
				jet_pt = p_data.AK8Puppijet0_pt_JESUp;
			} else if (systematic.EqualTo("jesdown")) {
				jet_pt = p_data.AK8Puppijet0_pt_JESDown;
			} else if (systematic.EqualTo("jerup")) {
				jet_pt = p_data.AK8Puppijet0_pt_JERUp;
			} else if (systematic.EqualTo("jerdown")) {
				jet_pt = p_data.AK8Puppijet0_pt_JERDown;
			} else {
				std::cerr << "[Max_AK8Puppijet0_pt] ERROR : Systematic " << systematic << " not known." << std::endl;
				exit(1);
			}
		} else {
			jet_pt = p_data.AK8Puppijet0_pt;
		}
		p_event_selector->SetReturnData("Max_AK8Puppijet0_pt", jet_pt);
		return (jet_pt <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_pt")[0]);
	}

	bool Min_AK8CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8CHSjet0_doublecsv", p_data.AK8CHSjet0_doublecsv);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		p_event_selector->SetReturnData("AK8Puppijet0_tau21", p_data.AK8Puppijet0_tau21);
		p_event_selector->SetReturnData("AK8Puppijet0_tau32", p_data.AK8Puppijet0_tau32);
		return (p_data.AK8CHSjet0_doublecsv >= p_event_selector->GetCutParameters("Min_AK8CHSjet0_doublecsv")[0]);
	}

	bool Max_AK8CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8CHSjet0_doublecsv", p_data.AK8CHSjet0_doublecsv);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		p_event_selector->SetReturnData("AK8Puppijet0_tau21", p_data.AK8Puppijet0_tau21);
		p_event_selector->SetReturnData("AK8Puppijet0_tau32", p_data.AK8Puppijet0_tau32);
		return (p_data.AK8CHSjet0_doublecsv <= p_event_selector->GetCutParameters("Max_AK8CHSjet0_doublecsv")[0]);
	}

	bool Min_nmuLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_nmuLoose", p_data.nmuLoose);
		return (p_data.nmuLoose >= (int)p_event_selector->GetCutParameters("Min_nmuLoose")[0]);
	}

	bool Max_nmuLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_nmuLoose", p_data.nmuLoose);
		return (p_data.nmuLoose <= (int)p_event_selector->GetCutParameters("Max_nmuLoose")[0]);
	}

	bool Min_vmuoLoose0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_vmuoLoose0_pt", p_data.vmuoLoose0_pt);
		return (p_data.vmuoLoose0_pt >= p_event_selector->GetCutParameters("Min_vmuoLoose0_pt")[0]);
	}

	bool Max_vmuoLoose0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_vmuoLoose0_pt", p_data.vmuoLoose0_pt);
		return (p_data.vmuoLoose0_pt <= p_event_selector->GetCutParameters("Max_vmuoLoose0_pt")[0]);
	}

	bool Min_vmuoLoose0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_vmuoLoose0_abseta", p_data.vmuoLoose0_eta);
		return (fabs(p_data.vmuoLoose0_eta) >= p_event_selector->GetCutParameters("Min_vmuoLoose0_abseta")[0]);
	}

	bool Max_vmuoLoose0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_vmuoLoose0_abseta", p_data.vmuoLoose0_eta);
		return (fabs(p_data.vmuoLoose0_eta) <= p_event_selector->GetCutParameters("Max_vmuoLoose0_abseta")[0]);
	}

	bool Min_neleLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_neleLoose", p_data.neleLoose);
		return (p_data.neleLoose >= (int)p_event_selector->GetCutParameters("Min_neleLoose")[0]);
	}

	bool Max_neleLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_neleLoose", p_data.neleLoose);
		return (p_data.neleLoose <= (int)p_event_selector->GetCutParameters("Max_neleLoose")[0]);
	}

	bool Min_ntau(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_ntau", p_data.neleLoose);
		return (p_data.ntau >= (int)p_event_selector->GetCutParameters("Min_ntau")[0]);
	}

	bool Max_ntau(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_ntau", p_data.neleLoose);
		return (p_data.ntau <= (int)p_event_selector->GetCutParameters("Max_ntau")[0]);
	}

	bool Min_nphoLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_nphoLoose", p_data.neleLoose);
		return (p_data.nphoLoose >= (int)p_event_selector->GetCutParameters("Min_nphoLoose")[0]);
	}

	bool Max_nphoLoose(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_nphoLoose", p_data.neleLoose);
		return (p_data.nphoLoose <= (int)p_event_selector->GetCutParameters("Max_nphoLoose")[0]);
	}

	bool AK8Puppijet0_isTightVJet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("AK8Puppijet0_isTightVJet", p_data.AK8Puppijet0_isTightVJet);
		return (p_data.AK8Puppijet0_isTightVJet == 1);
	}

	bool Min_pfmet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_pfmet", p_data.pfmet);
		return (p_data.pfmet >= p_event_selector->GetCutParameters("Min_pfmet")[0]);
	}

	bool Max_pfmet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_pfmet", p_data.pfmet);
		return (p_data.pfmet <= p_event_selector->GetCutParameters("Max_pfmet")[0]);
	}

	bool Min_puppet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double this_puppet = 0.;
		if (p_event_selector->GetCutDescriptors("Min_puppet").size() > 0) {
			TString systematic = (p_event_selector->GetCutDescriptors("Min_puppet")[0]);
			systematic.ToLower();
			if (systematic.EqualTo("nominal")) {
				this_puppet = p_data.puppet;
			} else if (systematic.EqualTo("jesup")) {
				this_puppet = p_data.puppet_JESUp;
			} else if (systematic.EqualTo("jesdown")) {
				this_puppet = p_data.puppet_JESDown;
			} else if (systematic.EqualTo("jerup")) {
				this_puppet = p_data.puppet_JERUp;
			} else if (systematic.EqualTo("jerdown")) {
				this_puppet = p_data.puppet_JERDown;
			} else {
				std::cerr << "[Min_puppet] ERROR : Systematic " << systematic << " not known." << std::endl;
				exit(1);
			}
		} else {
			this_puppet = p_data.puppet;
		}
		p_event_selector->SetReturnData("Min_puppet", this_puppet);
		return (this_puppet >= p_event_selector->GetCutParameters("Min_puppet")[0]);
	}

	bool Max_puppet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double this_puppet = 0.;
		if (p_event_selector->GetCutDescriptors("Max_puppet").size() > 0) {
			TString systematic = (p_event_selector->GetCutDescriptors("Max_puppet")[0]);
			systematic.ToLower();
			if (systematic.EqualTo("nominal")) {
				this_puppet = p_data.puppet;
			} else if (systematic.EqualTo("jesup")) {
				this_puppet = p_data.puppet_JESUp;
			} else if (systematic.EqualTo("jesdown")) {
				this_puppet = p_data.puppet_JESDown;
			} else if (systematic.EqualTo("jerup")) {
				this_puppet = p_data.puppet_JERUp;
			} else if (systematic.EqualTo("jerdown")) {
				this_puppet = p_data.puppet_JERDown;
			} else {
				std::cerr << "[Max_puppet] ERROR : Systematic " << systematic << " not known." << std::endl;
				exit(1);
			}
		} else {
			this_puppet = p_data.puppet;
		}
		p_event_selector->SetReturnData("Max_puppet", this_puppet);
		return (this_puppet <= p_event_selector->GetCutParameters("Max_puppet")[0]);
	}

	bool Max_nAK4PuppijetsPt30dR08_0(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_nAK4PuppijetsPt30dR08_0", p_data.nAK4PuppijetsPt30dR08_0);
		return (p_data.nAK4PuppijetsPt30dR08_0 <= (int)p_event_selector->GetCutParameters("Max_nAK4PuppijetsPt30dR08_0")[0]);
	}

	bool Min_nAK4PuppijetsMPt50dR08_0(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_nAK4PuppijetsMPt50dR08_0", p_data.nAK4PuppijetsMPt50dR08_0);
		return (p_data.nAK4PuppijetsMPt50dR08_0 >= (int)p_event_selector->GetCutParameters("Min_nAK4PuppijetsMPt50dR08_0")[0]);
	}

	bool Min_AK8Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8Puppijet0_tau21", p_data.AK8Puppijet0_tau21);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		return (p_data.AK8Puppijet0_tau21 >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_tau21")[0]);
	}

	bool Max_AK8Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8Puppijet0_tau21", p_data.AK8Puppijet0_tau21);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		p_event_selector->SetReturnData("AK8CHSjet0_doublecsv", p_data.AK8CHSjet0_doublecsv);
		return (p_data.AK8Puppijet0_tau21 <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_tau21")[0]);
	}

	bool Min_AK8Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8Puppijet0_tau21DDT", p_data.AK8Puppijet0_tau21DDT);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		return (p_data.AK8Puppijet0_tau21DDT >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_tau21DDT")[0]);
	}

	bool Max_AK8Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8Puppijet0_tau21DDT", p_data.AK8Puppijet0_tau21DDT);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		p_event_selector->SetReturnData("AK8CHSjet0_doublecsv", p_data.AK8CHSjet0_doublecsv);
		return (p_data.AK8Puppijet0_tau21DDT <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_tau21DDT")[0]);
	}

	bool Min_AK8Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8Puppijet0_N2DDT", p_data.AK8Puppijet0_N2DDT);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		return (p_data.AK8Puppijet0_N2DDT >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_N2DDT")[0]);
	}

	bool Max_AK8Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8Puppijet0_N2DDT", p_data.AK8Puppijet0_N2DDT);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		p_event_selector->SetReturnData("AK8CHSjet0_doublecsv", p_data.AK8CHSjet0_doublecsv);
		return (p_data.AK8Puppijet0_N2DDT <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_N2DDT")[0]);
	}

	bool Min_AK8Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_AK8Puppijet0_tau32", p_data.AK8Puppijet0_tau32);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		p_event_selector->SetReturnData("AK8CHSjet0_doublecsv", p_data.AK8CHSjet0_doublecsv);
		return (p_data.AK8Puppijet0_tau32 >= p_event_selector->GetCutParameters("Min_AK8Puppijet0_tau32")[0]);
	}

	bool Max_AK8Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_AK8Puppijet0_tau32", p_data.AK8Puppijet0_tau32);
		p_event_selector->SetReturnData("AK8Puppijet0_pt", p_data.AK8Puppijet0_pt);
		return (p_data.AK8Puppijet0_tau32 <= p_event_selector->GetCutParameters("Max_AK8Puppijet0_tau32")[0]);
	}

	bool Min_CA15Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_N2DDT", p_data.CA15Puppijet0_N2DDT);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		return (p_data.CA15Puppijet0_N2DDT >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_N2DDT")[0]);
	}

	bool Max_CA15Puppijet0_N2DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_N2DDT", p_data.CA15Puppijet0_N2DDT);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		p_event_selector->SetReturnData("CA15CHSjet0_doublecsv", p_data.CA15CHSjet0_doublecsv);
		return (p_data.CA15Puppijet0_N2DDT <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_N2DDT")[0]);
	}

	bool Min_CA15Puppijet0_msd(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_msd", p_data.CA15Puppijet0_msd);
		return (p_data.CA15Puppijet0_msd >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_msd")[0]);
	}

	bool Max_CA15Puppijet0_msd(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_msd", p_data.CA15Puppijet0_msd);
		return (p_data.CA15Puppijet0_msd <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_msd")[0]);
	}

	bool CA15Puppijet0_isTightVJet(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("CA15Puppijet0_isTightVJet", p_data.CA15Puppijet0_isTightVJet);
		return (p_data.CA15Puppijet0_isTightVJet == 1);
	}

	bool Min_CA15CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15CHSjet0_doublecsv", p_data.CA15CHSjet0_doublecsv);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		p_event_selector->SetReturnData("CA15Puppijet0_tau21", p_data.CA15Puppijet0_tau21);
		p_event_selector->SetReturnData("CA15Puppijet0_tau32", p_data.CA15Puppijet0_tau32);
		return (p_data.CA15CHSjet0_doublecsv >= p_event_selector->GetCutParameters("Min_CA15CHSjet0_doublecsv")[0]);
	}
	bool Max_CA15CHSjet0_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15CHSjet0_doublecsv", p_data.CA15CHSjet0_doublecsv);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		p_event_selector->SetReturnData("CA15Puppijet0_tau21", p_data.CA15Puppijet0_tau21);
		p_event_selector->SetReturnData("CA15Puppijet0_tau32", p_data.CA15Puppijet0_tau32);
		return (p_data.CA15CHSjet0_doublecsv <= p_event_selector->GetCutParameters("Max_CA15CHSjet0_doublecsv")[0]);
	}
	bool Min_CA15Puppijet0_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_doublesub", p_data.CA15Puppijet0_doublesub);
		return (p_data.CA15Puppijet0_doublesub >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_doublesub")[0]);
	}
	bool Max_CA15Puppijet0_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_doublesub", p_data.CA15Puppijet0_doublesub);
		return (p_data.CA15Puppijet0_doublesub <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_doublesub")[0]);
	}
	bool Min_CA15Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		return (p_data.CA15Puppijet0_pt >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_pt")[0]);
	}
	bool Max_CA15Puppijet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		return (p_data.CA15Puppijet0_pt <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_pt")[0]);
	}
	bool Min_CA15CHSjet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double jet_pt = 0.;
		if (p_event_selector->GetCutDescriptors("Min_CA15CHSjet0_pt").size() > 0) {
			TString systematic = (p_event_selector->GetCutDescriptors("Min_CA15CHSjet0_pt")[0]);
			systematic.ToLower();
			if (systematic.EqualTo("nominal")) {
				jet_pt = p_data.CA15CHSjet0_pt;
			//} else if (systematic.EqualTo("jesup")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JESUp;
			//} else if (systematic.EqualTo("jesdown")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JESDown;
			//} else if (systematic.EqualTo("jerup")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JERUp;
			//} else if (systematic.EqualTo("jerdown")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JERDown;
			} else {
				std::cerr << "[Min_CA15CHSjet0_pt] ERROR : Systematic " << systematic << " not known, or not implemented for CA15CHS jets." << std::endl;
				exit(1);
			}
		} else {
			jet_pt = p_data.CA15CHSjet0_pt;
		}
		p_event_selector->SetReturnData("Min_CA15CHSjet0_pt", jet_pt);
		return (jet_pt >= p_event_selector->GetCutParameters("Min_CA15CHSjet0_pt")[0]);
	}
	bool Max_CA15CHSjet0_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		double jet_pt = 0.;
		if (p_event_selector->GetCutDescriptors("Max_CA15CHSjet0_pt").size() > 0) {
			TString systematic = (p_event_selector->GetCutDescriptors("Max_CA15CHSjet0_pt")[0]);
			systematic.ToLower();
			if (systematic.EqualTo("nominal")) {
				jet_pt = p_data.CA15CHSjet0_pt;
			//} else if (systematic.EqualTo("jesup")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JESUp;
			//} else if (systematic.EqualTo("jesdown")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JESDown;
			//} else if (systematic.EqualTo("jerup")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JERUp;
			//} else if (systematic.EqualTo("jerdown")) {
			//	jet_pt = p_data.CA15CHSjet0_pt_JERDown;
			} else {
				std::cerr << "[Max_CA15CHSjet0_pt] ERROR : Systematic " << systematic << " not known, or not implemented for CA15CHS jets." << std::endl;
				exit(1);
			}
		} else {
			jet_pt = p_data.CA15CHSjet0_pt;
		}
		p_event_selector->SetReturnData("Max_CA15CHSjet0_pt", jet_pt);
		return (jet_pt <= p_event_selector->GetCutParameters("Max_CA15CHSjet0_pt")[0]);
	}
	bool Min_CA15Puppijet0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_abseta", p_data.CA15Puppijet0_eta);
		return (TMath::Abs(p_data.CA15Puppijet0_eta) >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_abseta")[0]);
	}
	bool Max_CA15Puppijet0_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_abseta", p_data.CA15Puppijet0_eta);
		return (TMath::Abs(p_data.CA15Puppijet0_eta) <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_abseta")[0]);
	}
	bool Min_CA15Puppijet0_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_phi", p_data.CA15Puppijet0_phi);
		return (p_data.CA15Puppijet0_phi >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_phi")[0]);
	}
	bool Max_CA15Puppijet0_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_phi", p_data.CA15Puppijet0_phi);
		return (p_data.CA15Puppijet0_phi <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_phi")[0]);
	}
	bool Min_CA15CHSjet1_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15CHSjet1_doublecsv", p_data.CA15CHSjet1_doublecsv);
		return (p_data.CA15CHSjet1_doublecsv >= p_event_selector->GetCutParameters("Min_CA15CHSjet1_doublecsv")[0]);
	}
	bool Max_CA15CHSjet1_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15CHSjet1_doublecsv", p_data.CA15CHSjet1_doublecsv);
		return (p_data.CA15CHSjet1_doublecsv <= p_event_selector->GetCutParameters("Max_CA15CHSjet1_doublecsv")[0]);
	}
	bool Min_CA15Puppijet1_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet1_doublesub", p_data.CA15Puppijet1_doublesub);
		return (p_data.CA15Puppijet1_doublesub >= p_event_selector->GetCutParameters("Min_CA15Puppijet1_doublesub")[0]);
	}
	bool Max_CA15Puppijet1_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet1_doublesub", p_data.CA15Puppijet1_doublesub);
		return (p_data.CA15Puppijet1_doublesub <= p_event_selector->GetCutParameters("Max_CA15Puppijet1_doublesub")[0]);
	}
	bool Min_CA15Puppijet1_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet1_pt", p_data.CA15Puppijet1_pt);
		return (p_data.CA15Puppijet1_pt >= p_event_selector->GetCutParameters("Min_CA15Puppijet1_pt")[0]);
	}
	bool Max_CA15Puppijet1_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet1_pt", p_data.CA15Puppijet1_pt);
		return (p_data.CA15Puppijet1_pt <= p_event_selector->GetCutParameters("Max_CA15Puppijet1_pt")[0]);
	}
	bool Min_CA15Puppijet1_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet1_abseta", p_data.CA15Puppijet1_eta);
		return (TMath::Abs(p_data.CA15Puppijet1_eta) >= p_event_selector->GetCutParameters("Min_CA15Puppijet1_abseta")[0]);
	}
	bool Max_CA15Puppijet1_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet1_abseta", p_data.CA15Puppijet1_eta);
		return (TMath::Abs(p_data.CA15Puppijet1_eta) <= p_event_selector->GetCutParameters("Max_CA15Puppijet1_abseta")[0]);
	}
	bool Min_CA15Puppijet1_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet1_phi", p_data.CA15Puppijet1_phi);
		return (p_data.CA15Puppijet1_phi >= p_event_selector->GetCutParameters("Min_CA15Puppijet1_phi")[0]);
	}
	bool Max_CA15Puppijet1_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet1_phi", p_data.CA15Puppijet1_phi);
		return (p_data.CA15Puppijet1_phi <= p_event_selector->GetCutParameters("Max_CA15Puppijet1_phi")[0]);
	}
	bool Min_CA15CHSjet2_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15CHSjet2_doublecsv", p_data.CA15CHSjet2_doublecsv);
		return (p_data.CA15CHSjet2_doublecsv >= p_event_selector->GetCutParameters("Min_CA15CHSjet2_doublecsv")[0]);
	}
	bool Max_CA15CHSjet2_doublecsv(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15CHSjet2_doublecsv", p_data.CA15CHSjet2_doublecsv);
		return (p_data.CA15CHSjet2_doublecsv <= p_event_selector->GetCutParameters("Max_CA15CHSjet2_doublecsv")[0]);
	}
	bool Min_CA15Puppijet2_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet2_doublesub", p_data.CA15Puppijet2_doublesub);
		return (p_data.CA15Puppijet2_doublesub >= p_event_selector->GetCutParameters("Min_CA15Puppijet2_doublesub")[0]);
	}
	bool Max_CA15Puppijet2_doublesub(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet2_doublesub", p_data.CA15Puppijet2_doublesub);
		return (p_data.CA15Puppijet2_doublesub <= p_event_selector->GetCutParameters("Max_CA15Puppijet2_doublesub")[0]);
	}
	bool Min_CA15Puppijet2_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet2_pt", p_data.CA15Puppijet2_pt);
		return (p_data.CA15Puppijet2_pt >= p_event_selector->GetCutParameters("Min_CA15Puppijet2_pt")[0]);
	}
	bool Max_CA15Puppijet2_pt(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet2_pt", p_data.CA15Puppijet2_pt);
		return (p_data.CA15Puppijet2_pt <= p_event_selector->GetCutParameters("Max_CA15Puppijet2_pt")[0]);
	}
	bool Min_CA15Puppijet2_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet2_abseta", p_data.CA15Puppijet2_eta);
		return (TMath::Abs(p_data.CA15Puppijet2_eta) >= p_event_selector->GetCutParameters("Min_CA15Puppijet2_abseta")[0]);
	}
	bool Max_CA15Puppijet2_abseta(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet2_abseta", p_data.CA15Puppijet2_eta);
		return (TMath::Abs(p_data.CA15Puppijet2_eta) <= p_event_selector->GetCutParameters("Max_CA15Puppijet2_abseta")[0]);
	}
	bool Min_CA15Puppijet2_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet2_phi", p_data.CA15Puppijet2_phi);
		return (p_data.CA15Puppijet2_phi >= p_event_selector->GetCutParameters("Min_CA15Puppijet2_phi")[0]);
	}
	bool Max_CA15Puppijet2_phi(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet2_phi", p_data.CA15Puppijet2_phi);
		return (p_data.CA15Puppijet2_phi <= p_event_selector->GetCutParameters("Max_CA15Puppijet2_phi")[0]);
	}

	bool Min_CA15Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_tau21", p_data.CA15Puppijet0_tau21);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		return (p_data.CA15Puppijet0_tau21 >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_tau21")[0]);
	}

	bool Max_CA15Puppijet0_tau21(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_tau21", p_data.CA15Puppijet0_tau21);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		p_event_selector->SetReturnData("CA15CHSjet0_doublecsv", p_data.CA15CHSjet0_doublecsv);
		return (p_data.CA15Puppijet0_tau21 <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_tau21")[0]);
	}

	bool Min_CA15Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_tau21DDT", p_data.CA15Puppijet0_tau21DDT);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		return (p_data.CA15Puppijet0_tau21DDT >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_tau21DDT")[0]);
	}

	bool Max_CA15Puppijet0_tau21DDT(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_tau21DDT", p_data.CA15Puppijet0_tau21DDT);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		p_event_selector->SetReturnData("CA15CHSjet0_doublecsv", p_data.CA15CHSjet0_doublecsv);
		return (p_data.CA15Puppijet0_tau21DDT <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_tau21DDT")[0]);
	}

	bool Min_CA15Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Min_CA15Puppijet0_tau32", p_data.CA15Puppijet0_tau32);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		p_event_selector->SetReturnData("CA15CHSjet0_doublecsv", p_data.CA15CHSjet0_doublecsv);
		return (p_data.CA15Puppijet0_tau32 >= p_event_selector->GetCutParameters("Min_CA15Puppijet0_tau32")[0]);
	}

	bool Max_CA15Puppijet0_tau32(const BaconData& p_data, EventSelector<BaconData>* p_event_selector) {
		p_event_selector->SetReturnData("Max_CA15Puppijet0_tau32", p_data.CA15Puppijet0_tau32);
		p_event_selector->SetReturnData("CA15Puppijet0_pt", p_data.CA15Puppijet0_pt);
		return (p_data.CA15Puppijet0_tau32 <= p_event_selector->GetCutParameters("Max_CA15Puppijet0_tau32")[0]);
	}


	void Configure(EventSelector<BaconData>* p_event_selector, TString name) {
		p_event_selector->AddCutFunction("Min_dphi_mu_jet", &Min_dphi_mu_jet);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_msd", &Min_AK8Puppijet0_msd);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_msd", &Max_AK8Puppijet0_msd);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_msd_puppi", &Min_AK8Puppijet0_msd_puppi);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_msd_puppi", &Max_AK8Puppijet0_msd_puppi);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_pt", &Min_AK8Puppijet0_pt);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_pt", &Max_AK8Puppijet0_pt);
		p_event_selector->AddCutFunction("Min_AK8CHSjet0_doublecsv", &Min_AK8CHSjet0_doublecsv);
		p_event_selector->AddCutFunction("Max_AK8CHSjet0_doublecsv", &Max_AK8CHSjet0_doublecsv);
		p_event_selector->AddCutFunction("Min_nmuLoose", &Min_nmuLoose);
		p_event_selector->AddCutFunction("Max_nmuLoose", &Max_nmuLoose);
		p_event_selector->AddCutFunction("Min_vmuoLoose0_pt", &Min_vmuoLoose0_pt);
		p_event_selector->AddCutFunction("Max_vmuoLoose0_pt", &Max_vmuoLoose0_pt);
		p_event_selector->AddCutFunction("Min_vmuoLoose0_abseta", &Min_vmuoLoose0_abseta);
		p_event_selector->AddCutFunction("Max_vmuoLoose0_abseta", &Max_vmuoLoose0_abseta);
		p_event_selector->AddCutFunction("Min_neleLoose", &Min_neleLoose);
		p_event_selector->AddCutFunction("Max_neleLoose", &Max_neleLoose);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_tau21", &Min_AK8Puppijet0_tau21);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_tau21", &Max_AK8Puppijet0_tau21);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_tau21DDT", &Min_AK8Puppijet0_tau21DDT);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_tau21DDT", &Max_AK8Puppijet0_tau21DDT);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_N2DDT", &Min_AK8Puppijet0_N2DDT);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_N2DDT", &Max_AK8Puppijet0_N2DDT);
		p_event_selector->AddCutFunction("Min_AK8Puppijet0_tau32", &Min_AK8Puppijet0_tau32);
		p_event_selector->AddCutFunction("Max_AK8Puppijet0_tau32", &Max_AK8Puppijet0_tau32);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_msd", &Min_CA15Puppijet0_msd);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_msd", &Max_CA15Puppijet0_msd);
		p_event_selector->AddCutFunction("Min_CA15CHSjet0_pt", &Min_CA15CHSjet0_pt);
		p_event_selector->AddCutFunction("Max_CA15CHSjet0_pt", &Max_CA15CHSjet0_pt);
		p_event_selector->AddCutFunction("Min_CA15CHSjet0_doublecsv", &Min_CA15CHSjet0_doublecsv);
		p_event_selector->AddCutFunction("Max_CA15CHSjet0_doublecsv", &Max_CA15CHSjet0_doublecsv);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_doublesub", &Min_CA15Puppijet0_doublesub);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_doublesub", &Max_CA15Puppijet0_doublesub);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_pt", &Min_CA15Puppijet0_pt);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_pt", &Max_CA15Puppijet0_pt);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_abseta", &Min_CA15Puppijet0_abseta);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_abseta", &Max_CA15Puppijet0_abseta);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_phi", &Min_CA15Puppijet0_phi);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_phi", &Max_CA15Puppijet0_phi);
		p_event_selector->AddCutFunction("Min_CA15CHSjet1_doublecsv", &Min_CA15CHSjet1_doublecsv);
		p_event_selector->AddCutFunction("Max_CA15CHSjet1_doublecsv", &Max_CA15CHSjet1_doublecsv);
		p_event_selector->AddCutFunction("Min_CA15Puppijet1_doublesub", &Min_CA15Puppijet1_doublesub);
		p_event_selector->AddCutFunction("Max_CA15Puppijet1_doublesub", &Max_CA15Puppijet1_doublesub);
		p_event_selector->AddCutFunction("Min_CA15Puppijet1_pt", &Min_CA15Puppijet1_pt);
		p_event_selector->AddCutFunction("Max_CA15Puppijet1_pt", &Max_CA15Puppijet1_pt);
		p_event_selector->AddCutFunction("Min_CA15Puppijet1_abseta", &Min_CA15Puppijet1_abseta);
		p_event_selector->AddCutFunction("Max_CA15Puppijet1_abseta", &Max_CA15Puppijet1_abseta);
		p_event_selector->AddCutFunction("Min_CA15Puppijet1_phi", &Min_CA15Puppijet1_phi);
		p_event_selector->AddCutFunction("Max_CA15Puppijet1_phi", &Max_CA15Puppijet1_phi);
		p_event_selector->AddCutFunction("Min_CA15CHSjet2_doublecsv", &Min_CA15CHSjet2_doublecsv);
		p_event_selector->AddCutFunction("Max_CA15CHSjet2_doublecsv", &Max_CA15CHSjet2_doublecsv);
		p_event_selector->AddCutFunction("Min_CA15Puppijet2_doublesub", &Min_CA15Puppijet2_doublesub);
		p_event_selector->AddCutFunction("Max_CA15Puppijet2_doublesub", &Max_CA15Puppijet2_doublesub);
		p_event_selector->AddCutFunction("Min_CA15Puppijet2_pt", &Min_CA15Puppijet2_pt);
		p_event_selector->AddCutFunction("Max_CA15Puppijet2_pt", &Max_CA15Puppijet2_pt);
		p_event_selector->AddCutFunction("Min_CA15Puppijet2_abseta", &Min_CA15Puppijet2_abseta);
		p_event_selector->AddCutFunction("Max_CA15Puppijet2_abseta", &Max_CA15Puppijet2_abseta);
		p_event_selector->AddCutFunction("Min_CA15Puppijet2_phi", &Min_CA15Puppijet2_phi);
		p_event_selector->AddCutFunction("Max_CA15Puppijet2_phi", &Max_CA15Puppijet2_phi);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_tau21", &Min_CA15Puppijet0_tau21);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_tau21", &Max_CA15Puppijet0_tau21);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_tau21DDT", &Min_CA15Puppijet0_tau21DDT);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_tau21DDT", &Max_CA15Puppijet0_tau21DDT);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_tau32", &Min_CA15Puppijet0_tau32);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_tau32", &Max_CA15Puppijet0_tau32);
		p_event_selector->AddCutFunction("Min_CA15Puppijet0_N2DDT", &Min_CA15Puppijet0_N2DDT);
		p_event_selector->AddCutFunction("Max_CA15Puppijet0_N2DDT", &Max_CA15Puppijet0_N2DDT);
		p_event_selector->AddCutFunction("Min_ntau", &Min_ntau);
		p_event_selector->AddCutFunction("Max_ntau", &Max_ntau);
		p_event_selector->AddCutFunction("Min_nphoLoose", &Min_nphoLoose);
		p_event_selector->AddCutFunction("Max_nphoLoose", &Max_nphoLoose);
		p_event_selector->AddCutFunction("AK8Puppijet0_isTightVJet", &AK8Puppijet0_isTightVJet);
		p_event_selector->AddCutFunction("CA15Puppijet0_isTightVJet", &CA15Puppijet0_isTightVJet);
		p_event_selector->AddCutFunction("Min_pfmet", &Min_pfmet);
		p_event_selector->AddCutFunction("Max_pfmet", &Max_pfmet);
		p_event_selector->AddCutFunction("Min_puppet", &Min_puppet);
		p_event_selector->AddCutFunction("Max_puppet", &Max_puppet);
		p_event_selector->AddCutFunction("Max_nAK4PuppijetsPt30dR08_0", &Max_nAK4PuppijetsPt30dR08_0);
		p_event_selector->AddCutFunction("Min_nAK4PuppijetsMPt50dR08_0", &Min_nAK4PuppijetsMPt50dR08_0);

		p_event_selector->AddNMinusOneHistogram("Min_dphi_mu_jet", "#Delta#phi(#mu, jet)", 100, -2.*TMath::Pi(), 2.*TMath::Pi());
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_msd", "m_{SD} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_msd", "m_{SD} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_msd_puppi", "m_{SD}^{(PUPPI)} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_msd_puppi", "m_{SD} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_AK8CHSjet0_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8CHSjet0_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_nmuLoose", "n_{#mu}", 11, -0.5, 10.5);
		p_event_selector->AddNMinusOneHistogram("Max_nmuLoose", "n_{#mu}", 11, -0.5, 10.5);
		p_event_selector->AddNMinusOneHistogram("Min_vmuoLoose0_pt", "#mu p_{T} [GeV]", 1000, 0., 1000.);
		p_event_selector->AddNMinusOneHistogram("Max_vmuoLoose0_pt", "#mu p_{T} [GeV]", 1000, 0., 1000.);
		p_event_selector->AddNMinusOneHistogram("Min_vmuoLoose0_abseta", "#mu #eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Max_vmuoLoose0_abseta", "#mu #eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Min_neleLoose", "n_{e}", 11, -0.5, 10.5);
		p_event_selector->AddNMinusOneHistogram("Max_neleLoose", "n_{e}", 11, -0.5, 10.5);
		p_event_selector->AddNMinusOneHistogram("Min_ntau", "n_{#tau}", 11, -0.5, 10.5); 
		p_event_selector->AddNMinusOneHistogram("Max_ntau", "n_{#tau}", 11, -0.5, 10.5); 
		p_event_selector->AddNMinusOneHistogram("Min_nphoLoose", "n_{#gamma}", 11, -0.5, 10.5); 
		p_event_selector->AddNMinusOneHistogram("Max_nphoLoose", "n_{#gamma}", 11, -0.5, 10.5); 
		p_event_selector->AddNMinusOneHistogram("AK8Puppijet0_isTightVJet", "Is tight VJet", 2, -0.5, 1.5); 
		p_event_selector->AddNMinusOneHistogram("CA15Puppijet0_isTightVJet", "Is tight VJet", 2, -0.5, 1.5);
		p_event_selector->AddNMinusOneHistogram("Min_pfmet", "MET [GeV]", 100, 0., 1000.); 
		p_event_selector->AddNMinusOneHistogram("Max_pfmet", "MET [GeV]", 100, 0., 1000.); 
		p_event_selector->AddNMinusOneHistogram("Min_puppet", "PUPPET [GeV]", 100, 0., 1000.); 
		p_event_selector->AddNMinusOneHistogram("Max_puppet", "PUPPET [GeV]", 100, 0., 1000.); 
		p_event_selector->AddNMinusOneHistogram("Max_nAK4PuppijetsPt30dR08_0", "Max_nAK4PuppijetsPt30dR08_0", 21, -0.5, 20.5);
		p_event_selector->AddNMinusOneHistogram("Min_nAK4PuppijetsMPt50dR08_0", "Min_nAK4PuppijetsMPt50dR08_0", 21, -0.5, 20.5);
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_tau21", "#tau_{21}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_tau21", "#tau_{21}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_tau21DDT", "#tau_{21}^{DDT}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_tau21DDT", "#tau_{21}^{DDT}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_N2DDT", "#N_{2}^{DDT}", 100, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_N2DDT", "#N_{2}^{DDT}", 100, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Min_AK8Puppijet0_tau32", "#tau_{32}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Max_AK8Puppijet0_tau32", "#tau_{32}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_msd", "m_{SD} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_msd", "m_{SD} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15CHSjet0_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15CHSjet0_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_doublesub", "Double sub", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_doublesub", "Double sub", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15CHSjet0_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15CHSjet0_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_abseta", "#eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_abseta", "#eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_phi", "#phi", 50, -5, 5);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_phi", "#phi", 50, -5, 5);
		p_event_selector->AddNMinusOneHistogram("Min_CA15CHSjet1_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15CHSjet1_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet1_doublesub", "Double sub", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet1_doublesub", "Double sub", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet1_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet1_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet1_abseta", "#eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet1_abseta", "#eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet1_phi", "#phi", 50, -5, 5);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet1_phi", "#phi", 50, -5, 5);
		p_event_selector->AddNMinusOneHistogram("Min_CA15CHSjet2_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15CHSjet2_doublecsv", "Double b-tag", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet2_doublesub", "Double sub", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet2_doublesub", "Double sub", 40, -2., 2.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet2_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet2_pt", "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet2_abseta", "#eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet2_abseta", "#eta", 50, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet2_phi", "#phi", 50, -5, 5);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet2_phi", "#phi", 50, -5, 5);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_tau21", "#tau_{21}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_tau21", "#tau_{21}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_tau21DDT", "#tau_{21}^{DDT}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_tau21DDT", "#tau_{21}^{DDT}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_tau32", "#tau_{32}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_tau32", "#tau_{32}", 50, 0., 1.);
		p_event_selector->AddNMinusOneHistogram("Min_CA15Puppijet0_N2DDT", "#N_{2}^{DDT}", 100, -5., 5.);
		p_event_selector->AddNMinusOneHistogram("Max_CA15Puppijet0_N2DDT", "#N_{2}^{DDT}", 100, -5., 5.);

		p_event_selector->AddNMinusOne2DHistogram("Min_AK8CHSjet0_doublecsv", "AK8Puppijet0_pt", "Double b-tag", 20, -1., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Max_AK8CHSjet0_doublecsv", "AK8Puppijet0_pt", "Double b-tag", 20, -1., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Min_AK8CHSjet0_doublecsv", "AK8Puppijet0_tau21", "Double b-tag", 20, -1., 1., "#tau_{21}", 10, 0., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Max_AK8CHSjet0_doublecsv", "AK8Puppijet0_tau21", "Double b-tag", 20, -1., 1., "#tau_{21}", 10, 0., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Min_AK8CHSjet0_doublecsv", "AK8Puppijet0_tau32", "Double b-tag", 20, -1., 1., "#tau_{32}", 10, 0., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Max_AK8CHSjet0_doublecsv", "AK8Puppijet0_tau32", "Double b-tag", 20, -1., 1., "#tau_{32}", 10, 0., 1.);

		p_event_selector->AddNMinusOne2DHistogram("Min_AK8Puppijet0_tau21", "AK8Puppijet0_pt", "#tau_{21}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Max_AK8Puppijet0_tau21", "AK8Puppijet0_pt", "#tau_{21}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Min_AK8Puppijet0_tau21", "AK8CHSjet0_doublecsv", "#tau_{21}", 10, 0., 1., "CSV", 20, -1., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Max_AK8Puppijet0_tau21", "AK8CHSjet0_doublecsv", "#tau_{21}", 10, 0., 1., "CSV", 20, -1., 1.);

		p_event_selector->AddNMinusOne2DHistogram("Min_AK8Puppijet0_tau32", "AK8Puppijet0_pt", "#tau_{32}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Max_AK8Puppijet0_tau32", "AK8Puppijet0_pt", "#tau_{32}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);

		p_event_selector->AddNMinusOne2DHistogram("Min_CA15CHSjet0_doublecsv", "CA15Puppijet0_pt", "Double b-tag", 20, -1., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Max_CA15CHSjet0_doublecsv", "CA15Puppijet0_pt", "Double b-tag", 20, -1., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Min_CA15CHSjet0_doublecsv", "CA15Puppijet0_tau21", "Double b-tag", 20, -1., 1., "#tau_{21}", 10, 0., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Max_CA15CHSjet0_doublecsv", "CA15Puppijet0_tau21", "Double b-tag", 20, -1., 1., "#tau_{21}", 10, 0., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Min_CA15CHSjet0_doublecsv", "CA15Puppijet0_tau32", "Double b-tag", 20, -1., 1., "#tau_{32}", 10, 0., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Max_CA15CHSjet0_doublecsv", "CA15Puppijet0_tau32", "Double b-tag", 20, -1., 1., "#tau_{32}", 10, 0., 1.);

		p_event_selector->AddNMinusOne2DHistogram("Min_CA15Puppijet0_tau21", "CA15Puppijet0_pt", "#tau_{21}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Max_CA15Puppijet0_tau21", "CA15Puppijet0_pt", "#tau_{21}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Min_CA15Puppijet0_tau21", "CA15CHSjet0_doublecsv", "#tau_{21}", 10, 0., 1., "CSV", 20, -1., 1.);
		p_event_selector->AddNMinusOne2DHistogram("Max_CA15Puppijet0_tau21", "CA15CHSjet0_doublecsv", "#tau_{21}", 10, 0., 1., "CSV", 20, -1., 1.);

		p_event_selector->AddNMinusOne2DHistogram("Min_CA15Puppijet0_tau32", "CA15Puppijet0_pt", "#tau_{32}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);
		p_event_selector->AddNMinusOne2DHistogram("Max_CA15Puppijet0_tau32", "CA15Puppijet0_pt", "#tau_{32}", 10, 0., 1., "p_{T} [GeV]", 200, 0., 2000.);

		if (!name.EqualTo("")) {
			p_event_selector->SetName(name);
		} else {
			p_event_selector->SetName("BaconEventSelector");
		}
		p_event_selector->SetObjectName("Event");
	}
}
#endif