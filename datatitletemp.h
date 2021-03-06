
DATA HISTO NAMES

EtaCan80[i] = new TCanvas(Form("EtaCan80_%d",i),Form("Jet80 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan80[i]->Divide(3,1);
	EtaCan80[i]->cd(1);
	h_data_Eta80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_etacent_R%d_cent%d",radius,i));
	h_data_Eta80[i]->GetXaxis()->SetTitle("Eta");
	h_data_Eta80[i]->GetXaxis()->CenterTitle();
	h_data_Eta80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_Eta80[i]->GetYaxis()->CenterTitle();
	h_data_Eta80[i]->Draw();
	EtaCan80[i]->cd(2);
	h_data_Phi80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phicent_R%d_cent%d",radius,i));
	h_data_Phi80[i]->GetXaxis()->SetTitle("Phi");
	h_data_Phi80[i]->GetXaxis()->CenterTitle();
	h_data_Phi80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_Phi80[i]->GetYaxis()->CenterTitle();
	h_data_Phi80[i]->Draw();
	EtaCan80[i]->cd(3);
	h_data_Phieta80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phieta_R%d_cent%d",radius,i));
	h_data_Phieta80[i]->GetXaxis()->SetTitle("Eta");
	h_data_Phieta80[i]->GetXaxis()->CenterTitle();
	h_data_Phieta80[i]->GetYaxis()->SetTitle("Phi");
	h_data_Phieta80[i]->GetYaxis()->CenterTitle();
	h_data_Phieta80[i]->Draw("colz,logz");

	EtaCan65[i] = new TCanvas(Form("EtaCan65_%d",i),Form("Jet65 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan65[i]->Divide(3,1);
	EtaCan65[i]->cd(1);
	h_data_Eta65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_etacent_R%d_cent%d",radius,i));
	h_data_Eta65[i]->GetXaxis()->SetTitle("Eta");
	h_data_Eta65[i]->GetXaxis()->CenterTitle();
	h_data_Eta65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_Eta65[i]->GetYaxis()->CenterTitle();
	h_data_Eta65[i]->Draw();
	EtaCan65[i]->cd(2);
	h_data_Phi65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phicent_R%d_cent%d",radius,i));
	h_data_Phi65[i]->GetXaxis()->SetTitle("Phi");
	h_data_Phi65[i]->GetXaxis()->CenterTitle();
	h_data_Phi65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_Phi65[i]->GetYaxis()->CenterTitle();
	h_data_Phi65[i]->Draw();
	EtaCan65[i]->cd(3);
	h_data_Phieta65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phieta_R%d_cent%d",radius,i));
	h_data_Phieta65[i]->GetXaxis()->SetTitle("Eta");
	h_data_Phieta65[i]->GetXaxis()->CenterTitle();
	h_data_Phieta65[i]->GetYaxis()->SetTitle("Phi");
	h_data_Phieta65[i]->GetYaxis()->CenterTitle();
	h_data_Phieta65[i]->Draw("colz,logz");

	EtaCan55[i] = new TCanvas(Form("EtaCan55_%d",i),Form("Jet55 R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaCan55[i]->Divide(3,1);
	EtaCan55[i]->cd(1);
	h_data_Eta55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_etacent_R%d_cent%d",radius,i));
	h_data_Eta55[i]->GetXaxis()->SetTitle("Eta");
	h_data_Eta55[i]->GetXaxis()->CenterTitle();
	h_data_Eta55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_Eta55[i]->GetYaxis()->SetRange(-0.5,2.5);
	h_data_Eta55[i]->GetYaxis()->CenterTitle();
	h_data_Eta55[i]->Draw();
	EtaCan55[i]->cd(2);
	h_data_Phi55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phicent_R%d_cent%d",radius,i));
	h_data_Phi55[i]->GetXaxis()->SetTitle("Phi");
	h_data_Phi55[i]->GetXaxis()->CenterTitle();
	h_data_Phi55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_Phi55[i]->GetYaxis()->CenterTitle();
	h_data_Phi55[i]->Draw();
	EtaCan55[i]->cd(3);
	h_data_Phieta55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phieta_R%d_cent%d",radius,i));
	h_data_Phieta55[i]->GetXaxis()->SetTitle("Eta");
	h_data_Phieta55[i]->GetXaxis()->CenterTitle();
	h_data_Phieta55[i]->GetYaxis()->SetTitle("Phi");
	h_data_Phieta55[i]->GetYaxis()->CenterTitle();
	h_data_Phieta55[i]->Draw("colz,logz");
		
	//PRECUT LOOP

	EtaPrecutCan80[i] = new TCanvas(Form("EtaPrecutCan80_%d",i),Form("Jet80 Precut R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaPrecutCan80[i]->Divide(3,1);
	EtaPrecutCan80[i]->cd(1);
	h_data_EtaPrecut80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_etaprecut_cent_R%d_cent%d",radius,i));
	h_data_EtaPrecut80[i]->GetXaxis()->SetTitle("Eta");
	h_data_EtaPrecut80[i]->GetXaxis()->CenterTitle();
	h_data_EtaPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_EtaPrecut80[i]->GetYaxis()->CenterTitle();
	h_data_EtaPrecut80[i]->Draw();
	EtaPrecutCan80[i]->cd(2);
	h_data_PhiPrecut80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phiprecut_cent_R%d_cent%d",radius,i));
	h_data_PhiPrecut80[i]->GetXaxis()->SetTitle("Phi");
	h_data_PhiPrecut80[i]->GetXaxis()->CenterTitle();
	h_data_PhiPrecut80[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_PhiPrecut80[i]->GetYaxis()->CenterTitle();
	h_data_PhiPrecut80[i]->Draw();
	EtaPrecutCan80[i]->cd(3);
	h_data_PhiEtaPrecut80[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet80_phieta_precut_R%d_cent%d",radius,i));
	h_data_PhiEtaPrecut80[i]->GetXaxis()->SetTitle("Eta");
	h_data_PhiEtaPrecut80[i]->GetXaxis()->CenterTitle();
	h_data_PhiEtaPrecut80[i]->GetYaxis()->SetTitle("Phi");
	h_data_PhiEtaPrecut80[i]->GetYaxis()->CenterTitle();
	h_data_PhiEtaPrecut80[i]->Draw("colz,logz");

	EtaPrecutCan65[i] = new TCanvas(Form("EtaPrecutCan65_%d",i),Form("Jet65 Precut R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaPrecutCan65[i]->Divide(3,1);
	EtaPrecutCan65[i]->cd(1);
	h_data_EtaPrecut65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_etaprecut_cent_R%d_cent%d",radius,i));
	h_data_EtaPrecut65[i]->GetXaxis()->SetTitle("Eta");
	h_data_EtaPrecut65[i]->GetXaxis()->CenterTitle();
	h_data_EtaPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_EtaPrecut65[i]->GetYaxis()->CenterTitle();
	h_data_EtaPrecut65[i]->Draw();
	EtaPrecutCan65[i]->cd(2);
	h_data_PhiPrecut65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phiprecut_cent_R%d_cent%d",radius,i));
	h_data_PhiPrecut65[i]->GetXaxis()->SetTitle("Phi");
	h_data_PhiPrecut65[i]->GetXaxis()->CenterTitle();
	h_data_PhiPrecut65[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_PhiPrecut65[i]->GetYaxis()->CenterTitle();
	h_data_PhiPrecut65[i]->Draw();
	EtaPrecutCan65[i]->cd(3);
	h_data_PhiEtaPrecut65[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet65_phieta_precut_R%d_cent%d",radius,i));
	h_data_PhiEtaPrecut65[i]->GetXaxis()->SetTitle("Eta");
	h_data_PhiEtaPrecut65[i]->GetXaxis()->CenterTitle();
	h_data_PhiEtaPrecut65[i]->GetYaxis()->SetTitle("Phi");
	h_data_PhiEtaPrecut65[i]->GetYaxis()->CenterTitle();
	h_data_PhiEtaPrecut65[i]->Draw("colz,logz");
	
	EtaPrecutCan55[i] = new TCanvas(Form("EtaPrecutCan55_%d",i),Form("Jet55 Precut R%d Eta and Phi plots %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaPrecutCan55[i]->Divide(3,1);
	EtaPrecutCan55[i]->cd(1);
	h_data_EtaPrecut55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_etaprecut_cent_R%d_cent%d",radius,i));
	h_data_EtaPrecut55[i]->GetXaxis()->SetTitle("Eta");
	h_data_EtaPrecut55[i]->GetXaxis()->CenterTitle();
	h_data_EtaPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_EtaPrecut55[i]->GetYaxis()->CenterTitle();
	h_data_EtaPrecut55[i]->Draw();
	EtaPrecutCan55[i]->cd(2);
	h_data_PhiPrecut55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phiprecut_cent_R%d_cent%d",radius,i));
	h_data_PhiPrecut55[i]->GetXaxis()->SetTitle("Phi");
	h_data_PhiPrecut55[i]->GetXaxis()->CenterTitle();
	h_data_PhiPrecut55[i]->GetYaxis()->SetTitle("Multiplicity");
	h_data_PhiPrecut55[i]->GetYaxis()->CenterTitle();
	h_data_PhiPrecut55[i]->Draw();
	EtaPrecutCan55[i]->cd(3);
	h_data_PhiEtaPrecut55[i] = (TH1F*)fdata->Get(Form("hpbpb_Jet55_phieta_precut_R%d_cent%d",radius,i));
	h_data_PhiEtaPrecut55[i]->GetXaxis()->SetTitle("Eta");
	h_data_PhiEtaPrecut55[i]->GetXaxis()->CenterTitle();
	h_data_PhiEtaPrecut55[i]->GetYaxis()->SetTitle("Phi");
	h_data_PhiEtaPrecut55[i]->GetYaxis()->CenterTitle();
	h_data_PhiEtaPrecut55[i]->Draw("colz,logz");
	
	EtaRatCan[i] = new TCanvas(Form("EtaRatios_%d",i),Form("Eta Cut/Precut Ratios R%d %2.0f - %2.0f cent",radius,5*boundaries_cent[i],5*boundaries_cent[i+1]),1750,500);
	EtaRatCan[i]->Divide(3,1);
	EtaRatCan[i]->cd(1);
	EtaRatio80[i] = new TH1F(Form("E80%d",i),Form("Eta 80 Cut/Precut"),30,-2.5,2.5);
	EtaRatio80[i]->Divide(h_data_Eta80[i],h_data_EtaPrecut80[i]);
	EtaRatio80[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio80[i]->GetXaxis()->CenterTitle();
	EtaRatio80[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio80[i]->Draw();
	EtaRatCan[i]->cd(2);
	EtaRatio65[i] = new TH1F(Form("E65%d",i),Form("Eta 65 Cut/Precut"),30,-2.5,2.5);
	EtaRatio65[i]->Divide(h_data_Eta65[i],h_data_EtaPrecut65[i]);
	EtaRatio65[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio65[i]->GetXaxis()->CenterTitle();
	EtaRatio65[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio65[i]->Draw();
	EtaRatCan[i]->cd(3);
	EtaRatio55[i] = new TH1F(Form("E55%d",i),Form("Eta 55 Cut/Precut"),30,-2.5,2.5);
	EtaRatio55[i]->Divide(h_data_Eta55[i],h_data_EtaPrecut55[i]);
	EtaRatio55[i]->GetXaxis()->SetTitle("Eta");
	EtaRatio55[i]->GetXaxis()->CenterTitle();
	EtaRatio55[i]->GetYaxis()->SetRangeUser(0.6,1.4);
	EtaRatio55[i]->Draw();