#ifndef BaconData_cc
#define BaconData_cc

#include "DAZSLE/PhiBBPlusJet/interface/BaconData.h"

BaconData::BaconData(TTree *tree) : BaconTree(tree) {
	// Default configs
	_jettype = kAK8;
	_jetordering = kPt;

	// Histogram for N2 DDT
	TFile *f_n2ddt_AK8 = new TFile("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ZqqJet/h3_n2ddt_26eff_36binrho11pt_Spring16.root","read");
	n2_ddt_transformation_AK8_ = (TH1D*)f_n2ddt_AK8->Get("h2ddt");
	n2_ddt_transformation_AK8_->SetName("h2ddt_AK8");
	n2_ddt_transformation_AK8_->SetDirectory(0);
	f_n2ddt_AK8->Close();
	delete f_n2ddt_AK8;

	//TFile *f_n2ddt_CA15 = new TFile("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/PbbJet/h3_n2ddt_CA15.root","read");
	TFile *f_n2ddt_CA15 = new TFile("$CMSSW_BASE/src/DAZSLE/PhiBBPlusJet/data/DDT_N2_CA15_wp0.26.root","read");
	n2_ddt_transformation_CA15_ = (TH1D*)f_n2ddt_CA15->Get("DDT");
	n2_ddt_transformation_CA15_->SetName("h2ddt_CA15");
	n2_ddt_transformation_CA15_->SetDirectory(0);
	f_n2ddt_CA15->Close();
	delete f_n2ddt_CA15;

	// PUPPI weight functions
	// Based on https://github.com/thaarres/PuppiSoftdropMassCorr Summer16
	puppi_corr_gen_ = new TF1("corrGEN", "[0]+[1]*pow(x*[2],-[3])", 200, 3500);
	puppi_corr_gen_->SetParameter(0, 1.00626);
	puppi_corr_gen_->SetParameter(1, -1.06161);
	puppi_corr_gen_->SetParameter(2, 0.0799900);
	puppi_corr_gen_->SetParameter(3, 1.20454);

	puppi_corr_reco_cen_ = new TF1("corrRECO_cen", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)", 200, 3500);
	puppi_corr_reco_cen_->SetParameter(0, 1.09302);
	puppi_corr_reco_cen_->SetParameter(1, -0.000150068);
	puppi_corr_reco_cen_->SetParameter(2, 3.44866e-07);
	puppi_corr_reco_cen_->SetParameter(3, -2.68100e-10);
	puppi_corr_reco_cen_->SetParameter(4, 8.67440e-14);
	puppi_corr_reco_cen_->SetParameter(5, -1.00114e-17);

	puppi_corr_reco_for_ = new TF1("corrRECO_for", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)", 200, 3500);
	puppi_corr_reco_for_->SetParameter(0, 1.27212);
	puppi_corr_reco_for_->SetParameter(1, -0.000571640);
	puppi_corr_reco_for_->SetParameter(2, 8.37289e-07);
	puppi_corr_reco_for_->SetParameter(3, -5.20433e-10);
	puppi_corr_reco_for_->SetParameter(4, 1.45375e-13);
	puppi_corr_reco_for_->SetParameter(5, -1.50389e-17);

	AK8Puppijet0_tau21 = 0.;
	CA15Puppijet0_tau21 = 0.;
	AK8Puppijet0_N2DDT = 0.;
	CA15Puppijet0_N2DDT = 0.;
	AK8Puppijet0_msd_puppi = 0.;
	CA15Puppijet0_msd_puppi = 0.;

	puppet_JESUp = 0.;
	puppet_JESDown = 0.;
	puppet_JERUp = 0.;
	puppet_JERDown = 0.;
	pfmet_JESUp = 0.;
	pfmet_JESDown = 0.;
	pfmet_JERUp = 0.;
	pfmet_JERDown = 0.;
}

BaconData::~BaconData() {
}

Int_t BaconData::GetEntry(Long64_t entry) {
	// Clear previous event

	// Load new event
	Int_t ret = fChain->GetEntry(entry);

	/*** Computed variables ***/
	AK8Puppijet0_msd_puppi = AK8Puppijet0_msd * PUPPIweight(AK8Puppijet0_pt, AK8Puppijet0_eta);
	AK8Puppijet1_msd_puppi = AK8Puppijet1_msd * PUPPIweight(AK8Puppijet1_pt, AK8Puppijet1_eta);
	AK8Puppijet2_msd_puppi = AK8Puppijet2_msd * PUPPIweight(AK8Puppijet2_pt, AK8Puppijet2_eta);

	CA15Puppijet0_msd_puppi = CA15Puppijet0_msd * PUPPIweight(CA15Puppijet0_pt, CA15Puppijet0_eta);
	CA15Puppijet1_msd_puppi = CA15Puppijet1_msd * PUPPIweight(CA15Puppijet1_pt, CA15Puppijet1_eta);
	CA15Puppijet2_msd_puppi = CA15Puppijet2_msd * PUPPIweight(CA15Puppijet2_pt, CA15Puppijet2_eta);

	AK8Puppijet0_tau21DDT = AK8Puppijet0_tau21 + 0.063*TMath::Log(AK8Puppijet0_msd_puppi*AK8Puppijet0_msd_puppi/AK8Puppijet0_pt);
	AK8Puppijet1_tau21DDT = AK8Puppijet1_tau21 + 0.063*TMath::Log(AK8Puppijet1_msd_puppi*AK8Puppijet1_msd_puppi/AK8Puppijet1_pt);
	AK8Puppijet2_tau21DDT = AK8Puppijet2_tau21 + 0.063*TMath::Log(AK8Puppijet2_msd_puppi*AK8Puppijet2_msd_puppi/AK8Puppijet2_pt);

	CA15Puppijet0_tau21DDT = CA15Puppijet0_tau21 + 0.063*TMath::Log(CA15Puppijet0_msd*CA15Puppijet0_msd/CA15Puppijet0_pt);
	CA15Puppijet1_tau21DDT = CA15Puppijet1_tau21 + 0.063*TMath::Log(CA15Puppijet1_msd*CA15Puppijet1_msd/CA15Puppijet1_pt);
	CA15Puppijet2_tau21DDT = CA15Puppijet2_tau21 + 0.063*TMath::Log(CA15Puppijet2_msd*CA15Puppijet2_msd/CA15Puppijet2_pt);

	AK8Puppijet0_rho = 2 * TMath::Log(AK8Puppijet0_msd_puppi/ AK8Puppijet0_pt);
	AK8Puppijet1_rho = 2 * TMath::Log(AK8Puppijet1_msd_puppi/ AK8Puppijet1_pt);
	AK8Puppijet2_rho = 2 * TMath::Log(AK8Puppijet2_msd_puppi/ AK8Puppijet2_pt);

	CA15Puppijet0_rho = 2 * TMath::Log(CA15Puppijet0_msd_puppi/ CA15Puppijet0_pt);
	CA15Puppijet1_rho = 2 * TMath::Log(CA15Puppijet1_msd_puppi/ CA15Puppijet1_pt);
	CA15Puppijet2_rho = 2 * TMath::Log(CA15Puppijet2_msd_puppi/ CA15Puppijet2_pt);

	// AK8Puppijet0_N2DDT
	int rho_index = n2_ddt_transformation_AK8_->GetXaxis()->FindBin(AK8Puppijet0_rho);
	if (rho_index > n2_ddt_transformation_AK8_->GetXaxis()->GetNbins()) {
		rho_index = n2_ddt_transformation_AK8_->GetXaxis()->GetNbins();
	} else if (rho_index <= 0) {
		rho_index = 1;
	}
	int pt_index = n2_ddt_transformation_AK8_->GetYaxis()->FindBin(AK8Puppijet0_pt);
	if (pt_index > n2_ddt_transformation_AK8_->GetYaxis()->GetNbins()) {
		pt_index = n2_ddt_transformation_AK8_->GetYaxis()->GetNbins();
	} else if (pt_index <= 0) {
		pt_index = 1;
	}
	AK8Puppijet0_N2DDT = AK8Puppijet0_N2sdb1 - n2_ddt_transformation_AK8_->GetBinContent(rho_index, pt_index);

	rho_index = n2_ddt_transformation_AK8_->GetXaxis()->FindBin(AK8Puppijet1_rho);
	if (rho_index > n2_ddt_transformation_AK8_->GetXaxis()->GetNbins()) {
		rho_index = n2_ddt_transformation_AK8_->GetXaxis()->GetNbins();
	} else if (rho_index <= 0) {
		rho_index = 1;
	}
	pt_index = n2_ddt_transformation_AK8_->GetYaxis()->FindBin(AK8Puppijet1_pt);
	if (pt_index > n2_ddt_transformation_AK8_->GetYaxis()->GetNbins()) {
		pt_index = n2_ddt_transformation_AK8_->GetYaxis()->GetNbins();
	} else if (pt_index <= 0) {
		pt_index = 1;
	}
	AK8Puppijet1_N2DDT = AK8Puppijet1_N2sdb1 - n2_ddt_transformation_AK8_->GetBinContent(rho_index, pt_index);

	rho_index = n2_ddt_transformation_AK8_->GetXaxis()->FindBin(AK8Puppijet2_rho);
	if (rho_index > n2_ddt_transformation_AK8_->GetXaxis()->GetNbins()) {
		rho_index = n2_ddt_transformation_AK8_->GetXaxis()->GetNbins();
	} else if (rho_index <= 0) {
		rho_index = 1;
	}
	pt_index = n2_ddt_transformation_AK8_->GetYaxis()->FindBin(AK8Puppijet2_pt);
	if (pt_index > n2_ddt_transformation_AK8_->GetYaxis()->GetNbins()) {
		pt_index = n2_ddt_transformation_AK8_->GetYaxis()->GetNbins();
	} else if (pt_index <= 0) {
		pt_index = 1;
	}
	AK8Puppijet2_N2DDT = AK8Puppijet2_N2sdb1 - n2_ddt_transformation_AK8_->GetBinContent(rho_index, pt_index);

	// CA15Puppijet0_N2DDT
	rho_index = n2_ddt_transformation_CA15_->GetXaxis()->FindBin(CA15Puppijet0_rho);
	if (rho_index > n2_ddt_transformation_CA15_->GetXaxis()->GetNbins()) {
		rho_index = n2_ddt_transformation_CA15_->GetXaxis()->GetNbins();
	} else if (rho_index <= 0) {
		rho_index = 1;
	}
	pt_index = n2_ddt_transformation_CA15_->GetYaxis()->FindBin(CA15Puppijet0_pt);
	if (pt_index > n2_ddt_transformation_CA15_->GetYaxis()->GetNbins()) {
		pt_index = n2_ddt_transformation_CA15_->GetYaxis()->GetNbins();
	} else if (pt_index <= 0) {
		pt_index = 1;
	}
	CA15Puppijet0_N2DDT = CA15Puppijet0_N2sdb1 - n2_ddt_transformation_CA15_->GetBinContent(rho_index, pt_index);

	rho_index = n2_ddt_transformation_CA15_->GetXaxis()->FindBin(CA15Puppijet1_rho);
	if (rho_index > n2_ddt_transformation_CA15_->GetXaxis()->GetNbins()) {
		rho_index = n2_ddt_transformation_CA15_->GetXaxis()->GetNbins();
	} else if (rho_index <= 0) {
		rho_index = 1;
	}
	pt_index = n2_ddt_transformation_CA15_->GetYaxis()->FindBin(CA15Puppijet1_pt);
	if (pt_index > n2_ddt_transformation_CA15_->GetYaxis()->GetNbins()) {
		pt_index = n2_ddt_transformation_CA15_->GetYaxis()->GetNbins();
	} else if (pt_index <= 0) {
		pt_index = 1;
	}
	CA15Puppijet1_N2DDT = CA15Puppijet1_N2sdb1 - n2_ddt_transformation_CA15_->GetBinContent(rho_index, pt_index);

	rho_index = n2_ddt_transformation_CA15_->GetXaxis()->FindBin(CA15Puppijet2_rho);
	if (rho_index > n2_ddt_transformation_CA15_->GetXaxis()->GetNbins()) {
		rho_index = n2_ddt_transformation_CA15_->GetXaxis()->GetNbins();
	} else if (rho_index <= 0) {
		rho_index = 1;
	}
	pt_index = n2_ddt_transformation_CA15_->GetYaxis()->FindBin(CA15Puppijet2_pt);
	if (pt_index > n2_ddt_transformation_CA15_->GetYaxis()->GetNbins()) {
		pt_index = n2_ddt_transformation_CA15_->GetYaxis()->GetNbins();
	} else if (pt_index <= 0) {
		pt_index = 1;
	}
	CA15Puppijet2_N2DDT = CA15Puppijet2_N2sdb1 - n2_ddt_transformation_CA15_->GetBinContent(rho_index, pt_index);


	// MET JES/JER
    double puppet_x = puppet * TMath::Cos(puppetphi);
    double puppet_y = puppet * TMath::Sin(puppetphi);
    puppet_JESUp = TMath::Sqrt((puppet_x + MetXCorrjesUp) * (puppet_x + MetXCorrjesUp) + (puppet_y + MetYCorrjesUp) * (puppet_y + MetYCorrjesUp));
    puppet_JESDown = TMath::Sqrt((puppet_x + MetXCorrjesDown) * (puppet_x + MetXCorrjesDown) + (puppet_y + MetYCorrjesDown) * (puppet_y + MetYCorrjesDown));
    puppet_JERUp = TMath::Sqrt((puppet_x + MetXCorrjerUp) * (puppet_x + MetXCorrjerUp) + (puppet_y + MetYCorrjerUp) * (puppet_y + MetYCorrjerUp));
    puppet_JERDown = TMath::Sqrt((puppet_x + MetXCorrjerDown) * (puppet_x + MetXCorrjerDown) + (puppet_y + MetYCorrjerDown) * (puppet_y + MetYCorrjerDown));

    double pfmet_x = pfmet * TMath::Cos(pfmetphi);
    double pfmet_y = pfmet * TMath::Sin(pfmetphi);
    pfmet_JESUp = TMath::Sqrt(TMath::Power(pfmet_x + MetXCorrjesUp, 2) + TMath::Power(pfmet_y + MetYCorrjesUp, 2));
    pfmet_JESDown = TMath::Sqrt(TMath::Power(pfmet_x + MetXCorrjesDown, 2) + TMath::Power(pfmet_y + MetYCorrjesDown, 2));
    pfmet_JERUp = TMath::Sqrt(TMath::Power(pfmet_x + MetXCorrjerUp, 2) + TMath::Power(pfmet_y + MetYCorrjerUp, 2));
    pfmet_JERDown = TMath::Sqrt(TMath::Power(pfmet_x + MetXCorrjerDown, 2) + TMath::Power(pfmet_y + MetYCorrjerDown, 2));


    // Selected jet stuff
    WhichJet_t which_jet = kAK8_0;
    if (_jettype == kAK8) {
    	if (_jetordering == kPt) {
    		double max_pt = 0.;
    		if (AK8Puppijet0_pt > max_pt) {
    			max_pt = AK8Puppijet0_pt;
    			which_jet = kAK8_0;
    		}
    		if (AK8Puppijet1_pt > max_pt) {
    			max_pt = AK8Puppijet1_pt;
    			which_jet = kAK8_1;
    		}
    		if (AK8Puppijet2_pt > max_pt) {
    			max_pt = AK8Puppijet2_pt;
    			which_jet = kAK8_2;
    		}
    	} else if (_jetordering == kDbtag) {
    		double max_dbtag = 0.;
    		if (AK8Puppijet0_doublecsv > max_dbtag) {
    			max_dbtag = AK8Puppijet0_doublecsv;
    			which_jet = kAK8_0;
    		}
    		if (AK8Puppijet1_doublecsv > max_dbtag) {
    			max_dbtag = AK8Puppijet1_doublecsv;
    			which_jet = kAK8_1;
    		}
    		if (AK8Puppijet2_doublecsv > max_dbtag) {
    			max_dbtag = AK8Puppijet2_doublecsv;
    			which_jet = kAK8_2;
    		}
    	} else if (_jetordering == kN2DDT) {
    		double max_n2ddt = 0.;
    		if (AK8Puppijet0_N2DDT > max_n2ddt) {
    			max_n2ddt = AK8Puppijet0_N2DDT;
    			which_jet = kAK8_0;
    		}
    		if (AK8Puppijet1_N2DDT > max_n2ddt) {
    			max_n2ddt = AK8Puppijet1_N2DDT;
    			which_jet = kAK8_1;
    		}
    		if (AK8Puppijet2_N2DDT > max_n2ddt) {
    			max_n2ddt = AK8Puppijet2_N2DDT;
    			which_jet = kAK8_2;
    		}
    	}
    } else if (_jettype == kCA15) {
    	if (_jetordering == kPt) {
    		double max_pt = 0.;
    		if (CA15Puppijet0_pt > max_pt) {
    			max_pt = CA15Puppijet0_pt;
    			which_jet = kCA15_0;
    		}
    		if (CA15Puppijet1_pt > max_pt) {
    			max_pt = CA15Puppijet1_pt;
    			which_jet = kCA15_1;
    		}
    		if (CA15Puppijet2_pt > max_pt) {
    			max_pt = CA15Puppijet2_pt;
    			which_jet = kCA15_2;
    		}
    	} else if (_jetordering == kDbtag) {
    		double max_dbtag = 0.;
    		if (CA15Puppijet0_doublesub > max_dbtag) {
    			max_dbtag = CA15Puppijet0_doublesub;
    			which_jet = kCA15_0;
    		}
    		if (CA15Puppijet1_doublesub > max_dbtag) {
    			max_dbtag = CA15Puppijet1_doublesub;
    			which_jet = kCA15_1;
    		}
    		if (CA15Puppijet2_doublesub > max_dbtag) {
    			max_dbtag = CA15Puppijet2_doublesub;
    			which_jet = kCA15_2;
    		}
    	} else if (_jetordering == kN2DDT) {
    		double max_n2ddt = 0.;
    		if (CA15Puppijet0_N2DDT > max_n2ddt) {
    			max_n2ddt = CA15Puppijet0_N2DDT;
    			which_jet = kCA15_0;
    		}
    		if (CA15Puppijet1_N2DDT > max_n2ddt) {
    			max_n2ddt = CA15Puppijet1_N2DDT;
    			which_jet = kCA15_1;
    		}
    		if (CA15Puppijet2_N2DDT > max_n2ddt) {
    			max_n2ddt = CA15Puppijet2_N2DDT;
    			which_jet = kCA15_2;
    		}
    	}
    }

    // Set selected jet vars
    if (which_jet == kAK8_0) {
		SelectedJet_pt           = AK8Puppijet0_pt;
		SelectedJet_eta          = AK8Puppijet0_eta;
		SelectedJet_phi          = AK8Puppijet0_phi;
		SelectedJet_csv          = AK8Puppijet0_csv;
		SelectedJet_CHF          = AK8Puppijet0_CHF;
		SelectedJet_NHF          = AK8Puppijet0_NHF;
		SelectedJet_NEMF         = AK8Puppijet0_NEMF;
		SelectedJet_tau21        = AK8Puppijet0_tau21;
		SelectedJet_tau32        = AK8Puppijet0_tau32;
		SelectedJet_msd          = AK8Puppijet0_msd;
		SelectedJet_rho          = AK8Puppijet0_rho;
		SelectedJet_minsubcsv    = AK8Puppijet0_minsubcsv;
		SelectedJet_maxsubcsv    = AK8Puppijet0_maxsubcsv;
		SelectedJet_doublecsv    = AK8Puppijet0_doublecsv;
		SelectedJet_doublesub    = AK8Puppijet0_doublesub;
		SelectedJet_ptraw        = AK8Puppijet0_ptraw;
		SelectedJet_genpt        = AK8Puppijet0_genpt;
		SelectedJet_e2_b1        = AK8Puppijet0_e2_b1;
		SelectedJet_e3_b1        = AK8Puppijet0_e3_b1;
		SelectedJet_e3_v1_b1     = AK8Puppijet0_e3_v1_b1;
		SelectedJet_e3_v2_b1     = AK8Puppijet0_e3_v2_b1;
		SelectedJet_e4_v1_b1     = AK8Puppijet0_e4_v1_b1;
		SelectedJet_e4_v2_b1     = AK8Puppijet0_e4_v2_b1;
		SelectedJet_e2_b2        = AK8Puppijet0_e2_b2;
		SelectedJet_e3_b2        = AK8Puppijet0_e3_b2;
		SelectedJet_e3_v1_b2     = AK8Puppijet0_e3_v1_b2;
		SelectedJet_e3_v2_b2     = AK8Puppijet0_e3_v2_b2;
		SelectedJet_e4_v1_b2     = AK8Puppijet0_e4_v1_b2;
		SelectedJet_e4_v2_b2     = AK8Puppijet0_e4_v2_b2;
		SelectedJet_e2_sdb1      = AK8Puppijet0_e2_sdb1;
		SelectedJet_e3_sdb1      = AK8Puppijet0_e3_sdb1;
		SelectedJet_e3_v1_sdb1   = AK8Puppijet0_e3_v1_sdb1;
		SelectedJet_e3_v2_sdb1   = AK8Puppijet0_e3_v2_sdb1;
		SelectedJet_e4_v1_sdb1   = AK8Puppijet0_e4_v1_sdb1;
		SelectedJet_e4_v2_sdb1   = AK8Puppijet0_e4_v2_sdb1;
		SelectedJet_e2_sdb2      = AK8Puppijet0_e2_sdb2;
		SelectedJet_e3_sdb2      = AK8Puppijet0_e3_sdb2;
		SelectedJet_e3_v1_sdb2   = AK8Puppijet0_e3_v1_sdb2;
		SelectedJet_e3_v2_sdb2   = AK8Puppijet0_e3_v2_sdb2;
		SelectedJet_e4_v1_sdb2   = AK8Puppijet0_e4_v1_sdb2;
		SelectedJet_e4_v2_sdb2   = AK8Puppijet0_e4_v2_sdb2;
		SelectedJet_N2sdb1       = AK8Puppijet0_N2sdb1;
		SelectedJet_N2sdb2       = AK8Puppijet0_N2sdb2;
		SelectedJet_M2sdb1       = AK8Puppijet0_M2sdb1;
		SelectedJet_M2sdb2       = AK8Puppijet0_M2sdb2;
		SelectedJet_D2sdb1       = AK8Puppijet0_D2sdb1;
		SelectedJet_D2sdb2       = AK8Puppijet0_D2sdb2;
		SelectedJet_N2b1         = AK8Puppijet0_N2b1;
		SelectedJet_N2b2         = AK8Puppijet0_N2b2;
		SelectedJet_M2b1         = AK8Puppijet0_M2b1;
		SelectedJet_M2b2         = AK8Puppijet0_M2b2;
		SelectedJet_D2b1         = AK8Puppijet0_D2b1;
		SelectedJet_D2b2         = AK8Puppijet0_D2b2;
		SelectedJet_pt_old       = AK8Puppijet0_pt_old;
		SelectedJet_pt_JESUp     = AK8Puppijet0_pt_JESUp;
		SelectedJet_pt_JESDown   = AK8Puppijet0_pt_JESDown;
		SelectedJet_pt_JERUp     = AK8Puppijet0_pt_JERUp;
		SelectedJet_pt_JERDown   = AK8Puppijet0_pt_JERDown;
		SelectedJet_isTightVJet  = AK8Puppijet0_isTightVJet;
		SelectedJet_tau21DDT     = AK8Puppijet0_tau21DDT;
		SelectedJet_rho          = AK8Puppijet0_rho;
		SelectedJet_N2DDT        = AK8Puppijet0_N2DDT;
		SelectedJet_msd_puppi    = AK8Puppijet0_msd_puppi;
		SelectedJet_nParticles = AK8Puppijet0_nParticles;
	} else if (which_jet == kAK8_1) {
		SelectedJet_pt           = AK8Puppijet1_pt;
		SelectedJet_eta          = AK8Puppijet1_eta;
		SelectedJet_phi          = AK8Puppijet1_phi;
		SelectedJet_csv          = AK8Puppijet1_csv;
		SelectedJet_CHF          = AK8Puppijet1_CHF;
		SelectedJet_NHF          = AK8Puppijet1_NHF;
		SelectedJet_NEMF         = AK8Puppijet1_NEMF;
		SelectedJet_tau21        = AK8Puppijet1_tau21;
		SelectedJet_tau32        = AK8Puppijet1_tau32;
		SelectedJet_msd          = AK8Puppijet1_msd;
		SelectedJet_rho          = AK8Puppijet1_rho;
		SelectedJet_minsubcsv    = AK8Puppijet1_minsubcsv;
		SelectedJet_maxsubcsv    = AK8Puppijet1_maxsubcsv;
		SelectedJet_doublecsv    = AK8Puppijet1_doublecsv;
		SelectedJet_doublesub    = AK8Puppijet1_doublesub;
		SelectedJet_ptraw        = AK8Puppijet1_ptraw;
		SelectedJet_genpt        = AK8Puppijet1_genpt;
		SelectedJet_e2_b1        = AK8Puppijet1_e2_b1;
		SelectedJet_e3_b1        = AK8Puppijet1_e3_b1;
		SelectedJet_e3_v1_b1     = AK8Puppijet1_e3_v1_b1;
		SelectedJet_e3_v2_b1     = AK8Puppijet1_e3_v2_b1;
		SelectedJet_e4_v1_b1     = AK8Puppijet1_e4_v1_b1;
		SelectedJet_e4_v2_b1     = AK8Puppijet1_e4_v2_b1;
		SelectedJet_e2_b2        = AK8Puppijet1_e2_b2;
		SelectedJet_e3_b2        = AK8Puppijet1_e3_b2;
		SelectedJet_e3_v1_b2     = AK8Puppijet1_e3_v1_b2;
		SelectedJet_e3_v2_b2     = AK8Puppijet1_e3_v2_b2;
		SelectedJet_e4_v1_b2     = AK8Puppijet1_e4_v1_b2;
		SelectedJet_e4_v2_b2     = AK8Puppijet1_e4_v2_b2;
		SelectedJet_e2_sdb1      = AK8Puppijet1_e2_sdb1;
		SelectedJet_e3_sdb1      = AK8Puppijet1_e3_sdb1;
		SelectedJet_e3_v1_sdb1   = AK8Puppijet1_e3_v1_sdb1;
		SelectedJet_e3_v2_sdb1   = AK8Puppijet1_e3_v2_sdb1;
		SelectedJet_e4_v1_sdb1   = AK8Puppijet1_e4_v1_sdb1;
		SelectedJet_e4_v2_sdb1   = AK8Puppijet1_e4_v2_sdb1;
		SelectedJet_e2_sdb2      = AK8Puppijet1_e2_sdb2;
		SelectedJet_e3_sdb2      = AK8Puppijet1_e3_sdb2;
		SelectedJet_e3_v1_sdb2   = AK8Puppijet1_e3_v1_sdb2;
		SelectedJet_e3_v2_sdb2   = AK8Puppijet1_e3_v2_sdb2;
		SelectedJet_e4_v1_sdb2   = AK8Puppijet1_e4_v1_sdb2;
		SelectedJet_e4_v2_sdb2   = AK8Puppijet1_e4_v2_sdb2;
		SelectedJet_N2sdb1       = AK8Puppijet1_N2sdb1;
		SelectedJet_N2sdb2       = AK8Puppijet1_N2sdb2;
		SelectedJet_M2sdb1       = AK8Puppijet1_M2sdb1;
		SelectedJet_M2sdb2       = AK8Puppijet1_M2sdb2;
		SelectedJet_D2sdb1       = AK8Puppijet1_D2sdb1;
		SelectedJet_D2sdb2       = AK8Puppijet1_D2sdb2;
		SelectedJet_N2b1         = AK8Puppijet1_N2b1;
		SelectedJet_N2b2         = AK8Puppijet1_N2b2;
		SelectedJet_M2b1         = AK8Puppijet1_M2b1;
		SelectedJet_M2b2         = AK8Puppijet1_M2b2;
		SelectedJet_D2b1         = AK8Puppijet1_D2b1;
		SelectedJet_D2b2         = AK8Puppijet1_D2b2;
		SelectedJet_pt_old       = AK8Puppijet1_pt_old;
		SelectedJet_pt_JESUp     = AK8Puppijet1_pt_JESUp;
		SelectedJet_pt_JESDown   = AK8Puppijet1_pt_JESDown;
		SelectedJet_pt_JERUp     = AK8Puppijet1_pt_JERUp;
		SelectedJet_pt_JERDown   = AK8Puppijet1_pt_JERDown;
		SelectedJet_isTightVJet  = AK8Puppijet1_isTightVJet;
		SelectedJet_tau21DDT     = AK8Puppijet1_tau21DDT;
		SelectedJet_rho          = AK8Puppijet1_rho;
		SelectedJet_N2DDT        = AK8Puppijet1_N2DDT;
		SelectedJet_msd_puppi    = AK8Puppijet1_msd_puppi;
		SelectedJet_nParticles = AK8Puppijet0_nParticles;
	} else if (which_jet == kAK8_2) {
		SelectedJet_pt           = AK8Puppijet2_pt;
		SelectedJet_eta          = AK8Puppijet2_eta;
		SelectedJet_phi          = AK8Puppijet2_phi;
		SelectedJet_csv          = AK8Puppijet2_csv;
		SelectedJet_CHF          = AK8Puppijet2_CHF;
		SelectedJet_NHF          = AK8Puppijet2_NHF;
		SelectedJet_NEMF         = AK8Puppijet2_NEMF;
		SelectedJet_tau21        = AK8Puppijet2_tau21;
		SelectedJet_tau32        = AK8Puppijet2_tau32;
		SelectedJet_msd          = AK8Puppijet2_msd;
		SelectedJet_rho          = AK8Puppijet2_rho;
		SelectedJet_minsubcsv    = AK8Puppijet2_minsubcsv;
		SelectedJet_maxsubcsv    = AK8Puppijet2_maxsubcsv;
		SelectedJet_doublecsv    = AK8Puppijet2_doublecsv;
		SelectedJet_doublesub    = AK8Puppijet2_doublesub;
		SelectedJet_ptraw        = AK8Puppijet2_ptraw;
		SelectedJet_genpt        = AK8Puppijet2_genpt;
		SelectedJet_e2_b1        = AK8Puppijet2_e2_b1;
		SelectedJet_e3_b1        = AK8Puppijet2_e3_b1;
		SelectedJet_e3_v1_b1     = AK8Puppijet2_e3_v1_b1;
		SelectedJet_e3_v2_b1     = AK8Puppijet2_e3_v2_b1;
		SelectedJet_e4_v1_b1     = AK8Puppijet2_e4_v1_b1;
		SelectedJet_e4_v2_b1     = AK8Puppijet2_e4_v2_b1;
		SelectedJet_e2_b2        = AK8Puppijet2_e2_b2;
		SelectedJet_e3_b2        = AK8Puppijet2_e3_b2;
		SelectedJet_e3_v1_b2     = AK8Puppijet2_e3_v1_b2;
		SelectedJet_e3_v2_b2     = AK8Puppijet2_e3_v2_b2;
		SelectedJet_e4_v1_b2     = AK8Puppijet2_e4_v1_b2;
		SelectedJet_e4_v2_b2     = AK8Puppijet2_e4_v2_b2;
		SelectedJet_e2_sdb1      = AK8Puppijet2_e2_sdb1;
		SelectedJet_e3_sdb1      = AK8Puppijet2_e3_sdb1;
		SelectedJet_e3_v1_sdb1   = AK8Puppijet2_e3_v1_sdb1;
		SelectedJet_e3_v2_sdb1   = AK8Puppijet2_e3_v2_sdb1;
		SelectedJet_e4_v1_sdb1   = AK8Puppijet2_e4_v1_sdb1;
		SelectedJet_e4_v2_sdb1   = AK8Puppijet2_e4_v2_sdb1;
		SelectedJet_e2_sdb2      = AK8Puppijet2_e2_sdb2;
		SelectedJet_e3_sdb2      = AK8Puppijet2_e3_sdb2;
		SelectedJet_e3_v1_sdb2   = AK8Puppijet2_e3_v1_sdb2;
		SelectedJet_e3_v2_sdb2   = AK8Puppijet2_e3_v2_sdb2;
		SelectedJet_e4_v1_sdb2   = AK8Puppijet2_e4_v1_sdb2;
		SelectedJet_e4_v2_sdb2   = AK8Puppijet2_e4_v2_sdb2;
		SelectedJet_N2sdb1       = AK8Puppijet2_N2sdb1;
		SelectedJet_N2sdb2       = AK8Puppijet2_N2sdb2;
		SelectedJet_M2sdb1       = AK8Puppijet2_M2sdb1;
		SelectedJet_M2sdb2       = AK8Puppijet2_M2sdb2;
		SelectedJet_D2sdb1       = AK8Puppijet2_D2sdb1;
		SelectedJet_D2sdb2       = AK8Puppijet2_D2sdb2;
		SelectedJet_N2b1         = AK8Puppijet2_N2b1;
		SelectedJet_N2b2         = AK8Puppijet2_N2b2;
		SelectedJet_M2b1         = AK8Puppijet2_M2b1;
		SelectedJet_M2b2         = AK8Puppijet2_M2b2;
		SelectedJet_D2b1         = AK8Puppijet2_D2b1;
		SelectedJet_D2b2         = AK8Puppijet2_D2b2;
		SelectedJet_pt_old       = AK8Puppijet2_pt_old;
		SelectedJet_pt_JESUp     = AK8Puppijet2_pt_JESUp;
		SelectedJet_pt_JESDown   = AK8Puppijet2_pt_JESDown;
		SelectedJet_pt_JERUp     = AK8Puppijet2_pt_JERUp;
		SelectedJet_pt_JERDown   = AK8Puppijet2_pt_JERDown;
		SelectedJet_isTightVJet  = AK8Puppijet2_isTightVJet;
		SelectedJet_tau21DDT     = AK8Puppijet2_tau21DDT;
		SelectedJet_rho          = AK8Puppijet2_rho;
		SelectedJet_N2DDT        = AK8Puppijet2_N2DDT;
		SelectedJet_msd_puppi    = AK8Puppijet2_msd_puppi;
		SelectedJet_nParticles = AK8Puppijet0_nParticles;
	} else if (which_jet == kCA15_0) {
		SelectedJet_pt           = CA15Puppijet0_pt;
		SelectedJet_eta          = CA15Puppijet0_eta;
		SelectedJet_phi          = CA15Puppijet0_phi;
		SelectedJet_csv          = CA15Puppijet0_csv;
		SelectedJet_CHF          = CA15Puppijet0_CHF;
		SelectedJet_NHF          = CA15Puppijet0_NHF;
		SelectedJet_NEMF         = CA15Puppijet0_NEMF;
		SelectedJet_tau21        = CA15Puppijet0_tau21;
		SelectedJet_tau32        = CA15Puppijet0_tau32;
		SelectedJet_msd          = CA15Puppijet0_msd;
		SelectedJet_rho          = CA15Puppijet0_rho;
		SelectedJet_minsubcsv    = CA15Puppijet0_minsubcsv;
		SelectedJet_maxsubcsv    = CA15Puppijet0_maxsubcsv;
		SelectedJet_doublecsv    = CA15Puppijet0_doublecsv;
		SelectedJet_doublesub    = CA15Puppijet0_doublesub;
		SelectedJet_ptraw        = CA15Puppijet0_ptraw;
		SelectedJet_genpt        = CA15Puppijet0_genpt;
		SelectedJet_e2_b1        = CA15Puppijet0_e2_b1;
		SelectedJet_e3_b1        = CA15Puppijet0_e3_b1;
		SelectedJet_e3_v1_b1     = CA15Puppijet0_e3_v1_b1;
		SelectedJet_e3_v2_b1     = CA15Puppijet0_e3_v2_b1;
		SelectedJet_e4_v1_b1     = CA15Puppijet0_e4_v1_b1;
		SelectedJet_e4_v2_b1     = CA15Puppijet0_e4_v2_b1;
		SelectedJet_e2_b2        = CA15Puppijet0_e2_b2;
		SelectedJet_e3_b2        = CA15Puppijet0_e3_b2;
		SelectedJet_e3_v1_b2     = CA15Puppijet0_e3_v1_b2;
		SelectedJet_e3_v2_b2     = CA15Puppijet0_e3_v2_b2;
		SelectedJet_e4_v1_b2     = CA15Puppijet0_e4_v1_b2;
		SelectedJet_e4_v2_b2     = CA15Puppijet0_e4_v2_b2;
		SelectedJet_e2_sdb1      = CA15Puppijet0_e2_sdb1;
		SelectedJet_e3_sdb1      = CA15Puppijet0_e3_sdb1;
		SelectedJet_e3_v1_sdb1   = CA15Puppijet0_e3_v1_sdb1;
		SelectedJet_e3_v2_sdb1   = CA15Puppijet0_e3_v2_sdb1;
		SelectedJet_e4_v1_sdb1   = CA15Puppijet0_e4_v1_sdb1;
		SelectedJet_e4_v2_sdb1   = CA15Puppijet0_e4_v2_sdb1;
		SelectedJet_e2_sdb2      = CA15Puppijet0_e2_sdb2;
		SelectedJet_e3_sdb2      = CA15Puppijet0_e3_sdb2;
		SelectedJet_e3_v1_sdb2   = CA15Puppijet0_e3_v1_sdb2;
		SelectedJet_e3_v2_sdb2   = CA15Puppijet0_e3_v2_sdb2;
		SelectedJet_e4_v1_sdb2   = CA15Puppijet0_e4_v1_sdb2;
		SelectedJet_e4_v2_sdb2   = CA15Puppijet0_e4_v2_sdb2;
		SelectedJet_N2sdb1       = CA15Puppijet0_N2sdb1;
		SelectedJet_N2sdb2       = CA15Puppijet0_N2sdb2;
		SelectedJet_M2sdb1       = CA15Puppijet0_M2sdb1;
		SelectedJet_M2sdb2       = CA15Puppijet0_M2sdb2;
		SelectedJet_D2sdb1       = CA15Puppijet0_D2sdb1;
		SelectedJet_D2sdb2       = CA15Puppijet0_D2sdb2;
		SelectedJet_N2b1         = CA15Puppijet0_N2b1;
		SelectedJet_N2b2         = CA15Puppijet0_N2b2;
		SelectedJet_M2b1         = CA15Puppijet0_M2b1;
		SelectedJet_M2b2         = CA15Puppijet0_M2b2;
		SelectedJet_D2b1         = CA15Puppijet0_D2b1;
		SelectedJet_D2b2         = CA15Puppijet0_D2b2;
		SelectedJet_pt_old       = CA15Puppijet0_pt_old;
		SelectedJet_pt_JESUp     = CA15Puppijet0_pt_JESUp;
		SelectedJet_pt_JESDown   = CA15Puppijet0_pt_JESDown;
		SelectedJet_pt_JERUp     = CA15Puppijet0_pt_JERUp;
		SelectedJet_pt_JERDown   = CA15Puppijet0_pt_JERDown;
		SelectedJet_isTightVJet  = CA15Puppijet0_isTightVJet;
		SelectedJet_tau21DDT     = CA15Puppijet0_tau21DDT;
		SelectedJet_rho          = CA15Puppijet0_rho;
		SelectedJet_N2DDT        = CA15Puppijet0_N2DDT;
		SelectedJet_msd_puppi    = CA15Puppijet0_msd_puppi;
		SelectedJet_nParticles = CA15Puppijet0_nParticles;
	} else if (which_jet == kCA15_1) {
		SelectedJet_pt           = CA15Puppijet1_pt;
		SelectedJet_eta          = CA15Puppijet1_eta;
		SelectedJet_phi          = CA15Puppijet1_phi;
		SelectedJet_csv          = CA15Puppijet1_csv;
		SelectedJet_CHF          = CA15Puppijet1_CHF;
		SelectedJet_NHF          = CA15Puppijet1_NHF;
		SelectedJet_NEMF         = CA15Puppijet1_NEMF;
		SelectedJet_tau21        = CA15Puppijet1_tau21;
		SelectedJet_tau32        = CA15Puppijet1_tau32;
		SelectedJet_msd          = CA15Puppijet1_msd;
		SelectedJet_rho          = CA15Puppijet1_rho;
		SelectedJet_minsubcsv    = CA15Puppijet1_minsubcsv;
		SelectedJet_maxsubcsv    = CA15Puppijet1_maxsubcsv;
		SelectedJet_doublecsv    = CA15Puppijet1_doublecsv;
		SelectedJet_doublesub    = CA15Puppijet1_doublesub;
		SelectedJet_ptraw        = CA15Puppijet1_ptraw;
		SelectedJet_genpt        = CA15Puppijet1_genpt;
		SelectedJet_e2_b1        = CA15Puppijet1_e2_b1;
		SelectedJet_e3_b1        = CA15Puppijet1_e3_b1;
		SelectedJet_e3_v1_b1     = CA15Puppijet1_e3_v1_b1;
		SelectedJet_e3_v2_b1     = CA15Puppijet1_e3_v2_b1;
		SelectedJet_e4_v1_b1     = CA15Puppijet1_e4_v1_b1;
		SelectedJet_e4_v2_b1     = CA15Puppijet1_e4_v2_b1;
		SelectedJet_e2_b2        = CA15Puppijet1_e2_b2;
		SelectedJet_e3_b2        = CA15Puppijet1_e3_b2;
		SelectedJet_e3_v1_b2     = CA15Puppijet1_e3_v1_b2;
		SelectedJet_e3_v2_b2     = CA15Puppijet1_e3_v2_b2;
		SelectedJet_e4_v1_b2     = CA15Puppijet1_e4_v1_b2;
		SelectedJet_e4_v2_b2     = CA15Puppijet1_e4_v2_b2;
		SelectedJet_e2_sdb1      = CA15Puppijet1_e2_sdb1;
		SelectedJet_e3_sdb1      = CA15Puppijet1_e3_sdb1;
		SelectedJet_e3_v1_sdb1   = CA15Puppijet1_e3_v1_sdb1;
		SelectedJet_e3_v2_sdb1   = CA15Puppijet1_e3_v2_sdb1;
		SelectedJet_e4_v1_sdb1   = CA15Puppijet1_e4_v1_sdb1;
		SelectedJet_e4_v2_sdb1   = CA15Puppijet1_e4_v2_sdb1;
		SelectedJet_e2_sdb2      = CA15Puppijet1_e2_sdb2;
		SelectedJet_e3_sdb2      = CA15Puppijet1_e3_sdb2;
		SelectedJet_e3_v1_sdb2   = CA15Puppijet1_e3_v1_sdb2;
		SelectedJet_e3_v2_sdb2   = CA15Puppijet1_e3_v2_sdb2;
		SelectedJet_e4_v1_sdb2   = CA15Puppijet1_e4_v1_sdb2;
		SelectedJet_e4_v2_sdb2   = CA15Puppijet1_e4_v2_sdb2;
		SelectedJet_N2sdb1       = CA15Puppijet1_N2sdb1;
		SelectedJet_N2sdb2       = CA15Puppijet1_N2sdb2;
		SelectedJet_M2sdb1       = CA15Puppijet1_M2sdb1;
		SelectedJet_M2sdb2       = CA15Puppijet1_M2sdb2;
		SelectedJet_D2sdb1       = CA15Puppijet1_D2sdb1;
		SelectedJet_D2sdb2       = CA15Puppijet1_D2sdb2;
		SelectedJet_N2b1         = CA15Puppijet1_N2b1;
		SelectedJet_N2b2         = CA15Puppijet1_N2b2;
		SelectedJet_M2b1         = CA15Puppijet1_M2b1;
		SelectedJet_M2b2         = CA15Puppijet1_M2b2;
		SelectedJet_D2b1         = CA15Puppijet1_D2b1;
		SelectedJet_D2b2         = CA15Puppijet1_D2b2;
		SelectedJet_pt_old       = CA15Puppijet1_pt_old;
		SelectedJet_pt_JESUp     = CA15Puppijet1_pt_JESUp;
		SelectedJet_pt_JESDown   = CA15Puppijet1_pt_JESDown;
		SelectedJet_pt_JERUp     = CA15Puppijet1_pt_JERUp;
		SelectedJet_pt_JERDown   = CA15Puppijet1_pt_JERDown;
		SelectedJet_isTightVJet  = CA15Puppijet1_isTightVJet;
		SelectedJet_tau21DDT     = CA15Puppijet1_tau21DDT;
		SelectedJet_rho          = CA15Puppijet1_rho;
		SelectedJet_N2DDT        = CA15Puppijet1_N2DDT;
		SelectedJet_msd_puppi    = CA15Puppijet1_msd_puppi;
		SelectedJet_nParticles = CA15Puppijet0_nParticles;
	} else if (which_jet == kCA15_2) {
		SelectedJet_pt           = CA15Puppijet2_pt;
		SelectedJet_eta          = CA15Puppijet2_eta;
		SelectedJet_phi          = CA15Puppijet2_phi;
		SelectedJet_csv          = CA15Puppijet2_csv;
		SelectedJet_CHF          = CA15Puppijet2_CHF;
		SelectedJet_NHF          = CA15Puppijet2_NHF;
		SelectedJet_NEMF         = CA15Puppijet2_NEMF;
		SelectedJet_tau21        = CA15Puppijet2_tau21;
		SelectedJet_tau32        = CA15Puppijet2_tau32;
		SelectedJet_msd          = CA15Puppijet2_msd;
		SelectedJet_rho          = CA15Puppijet2_rho;
		SelectedJet_minsubcsv    = CA15Puppijet2_minsubcsv;
		SelectedJet_maxsubcsv    = CA15Puppijet2_maxsubcsv;
		SelectedJet_doublecsv    = CA15Puppijet2_doublecsv;
		SelectedJet_doublesub    = CA15Puppijet2_doublesub;
		SelectedJet_ptraw        = CA15Puppijet2_ptraw;
		SelectedJet_genpt        = CA15Puppijet2_genpt;
		SelectedJet_e2_b1        = CA15Puppijet2_e2_b1;
		SelectedJet_e3_b1        = CA15Puppijet2_e3_b1;
		SelectedJet_e3_v1_b1     = CA15Puppijet2_e3_v1_b1;
		SelectedJet_e3_v2_b1     = CA15Puppijet2_e3_v2_b1;
		SelectedJet_e4_v1_b1     = CA15Puppijet2_e4_v1_b1;
		SelectedJet_e4_v2_b1     = CA15Puppijet2_e4_v2_b1;
		SelectedJet_e2_b2        = CA15Puppijet2_e2_b2;
		SelectedJet_e3_b2        = CA15Puppijet2_e3_b2;
		SelectedJet_e3_v1_b2     = CA15Puppijet2_e3_v1_b2;
		SelectedJet_e3_v2_b2     = CA15Puppijet2_e3_v2_b2;
		SelectedJet_e4_v1_b2     = CA15Puppijet2_e4_v1_b2;
		SelectedJet_e4_v2_b2     = CA15Puppijet2_e4_v2_b2;
		SelectedJet_e2_sdb1      = CA15Puppijet2_e2_sdb1;
		SelectedJet_e3_sdb1      = CA15Puppijet2_e3_sdb1;
		SelectedJet_e3_v1_sdb1   = CA15Puppijet2_e3_v1_sdb1;
		SelectedJet_e3_v2_sdb1   = CA15Puppijet2_e3_v2_sdb1;
		SelectedJet_e4_v1_sdb1   = CA15Puppijet2_e4_v1_sdb1;
		SelectedJet_e4_v2_sdb1   = CA15Puppijet2_e4_v2_sdb1;
		SelectedJet_e2_sdb2      = CA15Puppijet2_e2_sdb2;
		SelectedJet_e3_sdb2      = CA15Puppijet2_e3_sdb2;
		SelectedJet_e3_v1_sdb2   = CA15Puppijet2_e3_v1_sdb2;
		SelectedJet_e3_v2_sdb2   = CA15Puppijet2_e3_v2_sdb2;
		SelectedJet_e4_v1_sdb2   = CA15Puppijet2_e4_v1_sdb2;
		SelectedJet_e4_v2_sdb2   = CA15Puppijet2_e4_v2_sdb2;
		SelectedJet_N2sdb1       = CA15Puppijet2_N2sdb1;
		SelectedJet_N2sdb2       = CA15Puppijet2_N2sdb2;
		SelectedJet_M2sdb1       = CA15Puppijet2_M2sdb1;
		SelectedJet_M2sdb2       = CA15Puppijet2_M2sdb2;
		SelectedJet_D2sdb1       = CA15Puppijet2_D2sdb1;
		SelectedJet_D2sdb2       = CA15Puppijet2_D2sdb2;
		SelectedJet_N2b1         = CA15Puppijet2_N2b1;
		SelectedJet_N2b2         = CA15Puppijet2_N2b2;
		SelectedJet_M2b1         = CA15Puppijet2_M2b1;
		SelectedJet_M2b2         = CA15Puppijet2_M2b2;
		SelectedJet_D2b1         = CA15Puppijet2_D2b1;
		SelectedJet_D2b2         = CA15Puppijet2_D2b2;
		SelectedJet_pt_old       = CA15Puppijet2_pt_old;
		SelectedJet_pt_JESUp     = CA15Puppijet2_pt_JESUp;
		SelectedJet_pt_JESDown   = CA15Puppijet2_pt_JESDown;
		SelectedJet_pt_JERUp     = CA15Puppijet2_pt_JERUp;
		SelectedJet_pt_JERDown   = CA15Puppijet2_pt_JERDown;
		SelectedJet_isTightVJet  = CA15Puppijet2_isTightVJet;
		SelectedJet_tau21DDT     = CA15Puppijet2_tau21DDT;
		SelectedJet_rho          = CA15Puppijet2_rho;
		SelectedJet_N2DDT        = CA15Puppijet2_N2DDT;
		SelectedJet_msd_puppi    = CA15Puppijet2_msd_puppi;
		SelectedJet_nParticles = CA15Puppijet0_nParticles;
	}


	return ret;
}

Double_t BaconData::PUPPIweight(double pt, double eta) const {
	return (TMath::Abs(eta) < 1.3 ? 
		puppi_corr_gen_->Eval(pt) * puppi_corr_reco_cen_->Eval(pt) : 
		puppi_corr_gen_->Eval(pt) * puppi_corr_reco_for_->Eval(pt));
}

Bool_t BaconData::IsVMatched(double matching_dR) const {
	double dR = TMath::Sqrt(TMath::Power(SelectedJet_eta - genVEta, 2) + TMath::Power(SelectedJet_phi - genVPhi, 2));
	return dR <= matching_dR;
}


#endif