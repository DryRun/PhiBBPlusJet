#ifndef BaconData_cc
#define BaconData_cc

#include "DAZSLE/PhiBBPlusJet/interface/BaconData.h"

BaconData::BaconData(TTree *tree) : BaconTree(tree) {
	// Histogram for N2 DDT
	TFile *f_n2ddt_AK8 = new TFile("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ZqqJet/h3_n2ddt_26eff_36binrho11pt_Spring16.root","read");
	n2_ddt_transformation_AK8_ = (TH1D*)f_n2ddt_AK8->Get("h2ddt");
	n2_ddt_transformation_AK8_->SetName("h2ddt_AK8");
	n2_ddt_transformation_AK8_->SetDirectory(0);
	f_n2ddt_AK8->Close();
	delete f_n2ddt_AK8;

	TFile *f_n2ddt_CA15 = new TFile("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/PbbJet/h3_n2ddt_CA15.root","read");
	n2_ddt_transformation_CA15_ = (TH1D*)f_n2ddt_CA15->Get("h2ddt");
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

	CA15Puppijet0_msd_puppi = CA15Puppijet0_msd * PUPPIweight(CA15Puppijet0_pt, CA15Puppijet0_eta);

	AK8Puppijet0_tau21DDT = AK8Puppijet0_tau21 + 0.063*TMath::Log(AK8Puppijet0_msd_puppi*AK8Puppijet0_msd_puppi/AK8Puppijet0_pt);

	CA15Puppijet0_tau21DDT = CA15Puppijet0_tau21 + 0.063*TMath::Log(CA15Puppijet0_msd*CA15Puppijet0_msd/CA15Puppijet0_pt);

	AK8Puppijet0_rho = 2 * TMath::Log(AK8Puppijet0_msd_puppi/ AK8Puppijet0_pt);

	CA15Puppijet0_rho = 2 * TMath::Log(CA15Puppijet0_msd_puppi/ CA15Puppijet0_pt);

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

	return ret;
}

Double_t BaconData::PUPPIweight(double pt, double eta) const {
	return (TMath::Abs(eta) < 1.3 ? 
		puppi_corr_gen_->Eval(pt) * puppi_corr_reco_cen_->Eval(pt) : 
		puppi_corr_gen_->Eval(pt) * puppi_corr_reco_for_->Eval(pt));
}


#endif