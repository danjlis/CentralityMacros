/* QA_centrality */
// Author: Daniel Lis - August 15 - 2023
// Brief : This macro class gives a quick QA analysis
//         of a run in sPHENIX
//      To be run on the output of the CentralityReco Module


//void QA_FindCentralities()
void QA_MBDCalibrations(const int runnumber, const int loadRunCalibration, const int doPerRunCalibration);
void QA_MakeChargeSum(const int runnumber, const int loadRunCalibration);

void QA_centrality(
		   const int runnumber,
		   const int loadRunCalibration = 1,
		   const int doPerRunCalibration = 1,
		   const int doFindCentralities = 1,
		   const int doVertexSelection = 1
		   )
{

  QA_MBDCalibrations(runnumber, loadRunCalibration, doPerRunCalibration);  
  //QA_MakeChargeSum(runnumber, loadRunCalibration);
  return;
}

void QA_MakeChargeSum(const int runnumber, const int loadRunCalibration)
{

  int before = 2;
  int start_looking = 10;
  int min_range = 8;
  TFile *calibfile;
  float gain_corr[128];
  float g;
  double tthresh = 16;
  double cthresh = 1;
  if (loadRunCalibration)
    {
      calibfile = new TFile(Form("../calib/calib_mbd_%d.root", runnumber), "r");
        
      if (!calibfile)
	{
	  std::cout << " No calibration File found " <<std::endl;
	  return;
	}
  
      // 64 gain corrections on each side
      //  0 -  63 -- North
      // 64 - 127 -- South
      
    
      TNtuple *s = (TNtuple*) calibfile->Get("mbd_calib");
      s->SetBranchAddress("peak", &g);
      
      for (int i = 0 ; i < 128; i ++)
	{
	  s->GetEntry(i);
	  gain_corr[i] = g;
	}
    }

  int nbin = 2500;
  int maxrange = 2500;
  bool minbias = false;


  TH1D *h_charge_sum = new TH1D("h_charge_sum","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias = new TH1D("h_charge_sum_min_bias","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_30 = new TH1D("h_charge_sum_min_bias_w_vertex_30","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_recalib_sum = new TH1D("h_charge_recalib_sum","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_recalib_sum_min_bias = new TH1D("h_charge_recalib_sum_min_bias","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_recalib_sum_min_bias_w_vertex_30 = new TH1D("h_charge_recalib_sum_min_bias_w_vertex_30","",  nbin, -0.5, (float)maxrange - 0.5);

  TFile *file = new TFile(Form("../output/run%d/trees_%d.root", runnumber, runnumber), "r");
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms

  float zdc_energy_sum_n, zdc_energy_sum_s;
  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];

  float z_vertex;
  
  TTree *t = (TTree*)file->Get("T");
  t->SetBranchAddress("zdc_energy_sum_n",&zdc_energy_sum_n);
  t->SetBranchAddress("zdc_energy_sum_s",&zdc_energy_sum_s);
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge_raw);
  t->SetBranchAddress("mbd_time_raw",mbd_time_raw);
  for (int i = 0; i < t->GetEntries();i++)
    {
      minbias = false;
      z_vertex = 100;      
      if (i%10000 == 0) cout << i << endl;
      t->GetEntry(i);

      double charge_sum = 0;
      double charge_recalib_sum = 0;

      //// vertex
      int hits_n = 0;
      int hits_s = 0;      
      std::vector<float> time_sum_n;
      std::vector<float> time_sum_s;
      float sum_n = 0.;
      float sum_s = 0.;
      for (int ich = 0 ; ich < 64; ich++)
	{
	  
	  if (mbd_charge[ich] > cthresh && mbd_time[ich] < tthresh)
	    {
	      if (ich == 56 || ich == 54) continue;
	      hits_n++;
	      time_sum_n.push_back(mbd_time[ich]);
	      sum_n += mbd_time[ich];
	    }      
	  if (mbd_charge[ich + 64] > cthresh && mbd_time[ich+64] < tthresh)
	    {
	      hits_s++;
	      time_sum_s.push_back(mbd_time[ich+64]);
	      sum_s += mbd_time[ich+64];
	    }

      
	} 
      
      sort(time_sum_n.begin(), time_sum_n.end());
      sort(time_sum_s.begin(), time_sum_s.end());

      if (hits_s >=2 && hits_n >=2){
	minbias = true;
	float median_n;
	if (time_sum_n.size()%2)
	  {
	    median_n = time_sum_n.at(time_sum_n.size()/2);
	  }
	else
	  {
	    median_n = (time_sum_n.at(time_sum_n.size()/2) + time_sum_n.at(time_sum_n.size()/2 - 1))/2.;
	  }
	float median_s;
	if (time_sum_s.size()%2)
	  {
	    median_s = time_sum_s.at(time_sum_s.size()/2);
	  }
	else
	  {
	    median_s = (time_sum_s.at(time_sum_s.size()/2) + time_sum_s.at(time_sum_s.size()/2 - 1))/2.;
	  }

	z_vertex = 15*(median_n - median_s);

      }


      ///
      for (int ich = 0 ; ich < 128; ich++)
	{
	  if (mbd_time[ich] < tthresh && mbd_charge[ich] > cthresh)
	    {
	      charge_sum += mbd_charge[ich];
	    }
	  if (mbd_time[ich] < tthresh && mbd_charge_raw[ich]*gain_corr[ich] > cthresh)
	    {
	      charge_recalib_sum += mbd_charge_raw[ich]*gain_corr[ich];
	    }
	}

      h_charge_sum->Fill(charge_sum);
      h_charge_recalib_sum->Fill(charge_recalib_sum);

      if (zdc_energy_sum_n < 40 || zdc_energy_sum_s < 40) continue;

      h_charge_sum_min_bias->Fill(charge_sum);
      h_charge_recalib_sum_min_bias->Fill(charge_recalib_sum);

      if (TMath::Abs(z_vertex) > 30) continue;

      h_charge_sum_min_bias_w_vertex_30->Fill(charge_sum);
      h_charge_recalib_sum_min_bias_w_vertex_30->Fill(charge_recalib_sum);

    }

  TFile *fout = new TFile(Form("../plots/mbd_charge_sum_%d.root", runnumber), "RECREATE");
  h_charge_sum->Write();
  h_charge_recalib_sum->Write();
  h_charge_sum_min_bias->Write();
  h_charge_recalib_sum_min_bias->Write();
  h_charge_sum_min_bias_w_vertex_30->Write();
  h_charge_recalib_sum_min_bias_w_vertex_30->Write();
  fout->Close();
  
  
}

void QA_MBDCalibrations(const int runnumber, const int loadRunCalibration, const int doPerRunCalibration)
{
  
  // Making the calibration file for the MBD

  int before = 2;
  int start_looking = 10;
  int min_range = 8;
  TFile *calibfile;
  float gain_corr[128];
  float g;

  if (loadRunCalibration)
    {
      calibfile = new TFile("../calib/calib_mbd.root","r");
  
      if (!calibfile)
	{
	  std::cout << " No calibration File found " <<std::endl;
	  return;
	}
  
      // 64 gain corrections on each side
      //  0 -  63 -- North
      // 64 - 127 -- South
      
    
      TNtuple *s = (TNtuple*) calibfile->Get("mbd_calib");
      s->SetBranchAddress("peak", &g);
      
      for (int i = 0 ; i < 128; i ++)
	{
	  s->GetEntry(i);
	  gain_corr[i] = g;
	}
    }

  // Output from the centrality module.

  TFile *file = new TFile(Form("../output/run%d/trees_%d.root", runnumber, runnumber), "r");
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms

  TH1D *h_charge[128];
  TH1D *h_charge_raw[128];
  TH1D *h_charge_recalib[128];
  TH1D *h_time[128];
  TH1D *h_time_raw[128];

  for (int j = 0; j < 128; j++)
    {
      h_charge[j] = new TH1D(Form("h_mbd_charge_ch%d", j), "",500, 0, 10);
      h_charge_raw[j] = new TH1D(Form("h_mbd_charge_raw_ch%d", j),"", 600, -0.5, 599.5);
      h_charge_recalib[j] = new TH1D(Form("h_mbd_charge_recalib_ch%d", j),"", 500, 0, 10);
      h_time[j] = new TH1D(Form("h_mbd_time_ch%d", j), "", 1000, -25, 25);
      h_time_raw[j] = new TH1D(Form("h_mbd_time_raw_ch%d", j), "",1000, 0, 15000);
    }


  float zdc_energy_sum_n, zdc_energy_sum_s;
  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];
  int isMinBias;
  float z_vertex;
  
  TTree *t = (TTree*)file->Get("T");
  t->SetBranchAddress("zdc_energy_sum_n",&zdc_energy_sum_n);
  t->SetBranchAddress("zdc_energy_sum_s",&zdc_energy_sum_s);
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge_raw);
  t->SetBranchAddress("mbd_time_raw",mbd_time_raw);
  t->SetBranchAddress("z_vertex",&z_vertex);
  t->SetBranchAddress("isMinBias",&isMinBias);
  
  for (int i = 0; i < t->GetEntries();i++)
    {
      if (i%10000 == 0) cout << i << endl;
      t->GetEntry(i);
      if (zdc_energy_sum_n < 40 || zdc_energy_sum_s < 40) continue;
      for (int ich = 0 ; ich < 128; ich++)
	{
	  h_charge[ich]->Fill(mbd_charge[ich]);
	  h_time[ich]->Fill(mbd_time[ich]);
	  h_charge_raw[ich]->Fill(mbd_charge_raw[ich]);
	  h_time_raw[ich]->Fill(mbd_time_raw[ich]);
	}
    }

  TH1D *h_peaks = new TH1D("h_peaks","", 64, -0.5, 1.5);
  TH1D *h_peaks_raw = new TH1D("h_peaks_raw","", 64, 0, 500);
  TH1D *h_peaks_recalib = new TH1D("h_peaks_recalib","", 64, 0.5, 1.5);

  TF1 *f_lan_w_exp = new TF1("lan_w_exp","[0]*TMath::Landau(x,[1],[2],3) + expo(3)");
  f_lan_w_exp->SetParLimits(0, 100, 10000000);
  f_lan_w_exp->SetParLimits(1, 30, 400);
  f_lan_w_exp->SetParLimits(2, 0.1, 4000);
  TFile *fcalibout;
  TNtuple *tn;
  fcalibout = new TFile(Form("../calib/calib_mbd_%d.root", runnumber),"RECREATE");
  tn = new TNtuple("mbd_calib","mbd_calib","channel:peak:width");

  for (int ich = 0 ; ich < 128; ich++){

    std::cout << "Raw Fitting :: Channel "<<ich;
    
    double local_min = 0;
    double local_max = 0;
    double max_value = 0;
    int good = 0;
    int low_bin = 0;
    double upper;
    int high_bin = 0;
    TH1D *hsmooth = (TH1D*) h_charge_raw[ich]->Clone();
    hsmooth->Smooth();
    hsmooth->Smooth();
    // Finding the dip
    for (int ic = start_looking; ic < hsmooth->GetNbinsX(); ic++)
      {
	good = 0;
	
	for (int j = -1*min_range; j <= min_range; j++)
	  {
	    if (j == 0) continue;
	    
	    if (hsmooth->GetBinContent(ic+1) <= hsmooth->GetBinContent(ic+1 +j))
	      good++;
	    
	  }
	if (good == 2*min_range)
	  {
	    local_min = hsmooth->GetBinLowEdge(ic);
	    low_bin = ic;
	    break;
	  }
      }
    if (low_bin == 0)
      {
	low_bin = start_looking;
      }

    // Find local Max after the dip
    
    for (int ic = low_bin; ic < hsmooth->GetNbinsX(); ic++)
      {
	if (hsmooth->GetBinContent(ic+1) > max_value)
	  {
	    high_bin = ic+1; 
	    max_value = hsmooth->GetBinContent(high_bin);
	    local_max = hsmooth->GetBinLowEdge(high_bin);
	  }
      }
    
    int tries = 0;
    while (tries < 5)
      {

	

	if (tries == 0)
	  {
	    if (local_max > 600){
	      local_max = 150;
	    }
	    else if (local_max < 40)
	      {
		local_max = 41;
	      }
	    if (local_min > 400) 
	      {
		local_min = 30;
		local_max = 150;
	      }
	    else if (local_min < 20)
	      {
		local_min = 20;
		
	      }

	    f_lan_w_exp->SetParLimits(1, (local_max - 20 < 15?15:local_max - 20), 500);
	    f_lan_w_exp->SetParameter(0, 10000);

	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 10);
	    f_lan_w_exp->SetParameter(3, 1);
	    f_lan_w_exp->SetParameter(4, 0.001);
	    
	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));
	    hsmooth->Fit("lan_w_exp","Q","",local_min, upper);
	    h_charge_raw[ich]->Fit("lan_w_exp","Q","",local_min, upper);
	  }
	else
	  {

	    local_max = hsmooth->GetBinLowEdge(high_bin);
	    f_lan_w_exp->SetParameter(0, 10000);

	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 10);
	    f_lan_w_exp->SetParameter(3, 1);
	    f_lan_w_exp->SetParameter(4, 0.001);

	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));

	    hsmooth->Fit("lan_w_exp","Q","",local_min, upper);

	    h_charge_raw[ich]->Fit("lan_w_exp","Q","",local_min, upper);

	  }
	std::cout <<" try " <<tries <<endl;
	if (f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF() < 5)
	  {
	    break;
	  }
	tries++;
      }
    
    std::cout <<" " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_exp->GetParameter(1)<<endl;
    gain_corr[ich] = 1./((double) f_lan_w_exp->GetParameter(1));
  }    
    
  for (int i = 0; i < t->GetEntries();i++)
    {
      if (i%10000 == 0) cout << i << endl;
      t->GetEntry(i);
      if (zdc_energy_sum_n < 40 || zdc_energy_sum_s < 40) continue;
      for (int ich = 0 ; ich < 128; ich++)
	{
	  h_charge_recalib[ich]->Fill(mbd_charge_raw[ich]*gain_corr[ich]);
	}
    }


  for (int ich = 0 ; ich < 128; ich++){

    std::cout << "Recalib Fitting :: Channel "<<ich;

    f_lan_w_exp->SetParLimits(0, 100, 50000);
    f_lan_w_exp->SetParLimits(1, 0.2, 1.5);
    f_lan_w_exp->SetParLimits(2, 0.01, 0.8);

    
    double local_min = 0;
    double local_max = 0;
    double max_value = 0;
    int good = 0;
    int low_bin = 0;
    int high_bin = 0;
    float upper;
    TH1D *hsmooth = (TH1D*)h_charge_recalib[ich]->Clone();
    hsmooth->Smooth();
    hsmooth->Smooth();
    // Finding the dip
    int tries = 0;
    while (tries < 5)
      {
	for (int ic = start_looking; ic < hsmooth->GetNbinsX(); ic++)
	  {
	    good = 0;
	    for (int j = -1*min_range; j <= min_range; j++)
	      {
		if (j == 0) continue;

		if (hsmooth->GetBinContent(ic+1) <= hsmooth->GetBinContent(ic+1 +j))
		  good++;

	      }
	    if (good == 2*min_range)
	      {
		local_min = hsmooth->GetBinLowEdge(ic - before);
		low_bin = ic - before;
		break;
	      }
	  }

	// Find local Max after the dip

	for (int ic = low_bin; ic < hsmooth->GetNbinsX(); ic++)
	  {
	    if (hsmooth->GetBinContent(ic+1) > max_value)
	      {
		high_bin = ic+1; 
		max_value = hsmooth->GetBinContent(high_bin);
		local_max = hsmooth->GetBinLowEdge(high_bin);
	      }
	  }

	if (local_max > 2){
	  local_max = 1;
	  local_min = 0.2;
	}
	else if (local_max < 0.3){
	  local_max = 1;
	  local_min = 0.2;
	}

	if (tries == 0)
	  {
	    f_lan_w_exp->SetParLimits(1, (local_max - 0.5 < 0.2?0.2:local_max - 0.5), local_max + 0.5);
	    f_lan_w_exp->SetParameter(0, 10000);
	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 0.12);
	    f_lan_w_exp->SetParameter(3, 1);
	    f_lan_w_exp->SetParameter(4, 0.001);
	    
	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));
	  }
	else 
	  {
	    f_lan_w_exp->SetParLimits(1, (local_max - 0.2 < 0.2?0.2:local_max - 0.2), local_max + 0.2);
	    f_lan_w_exp->SetParameter(0, 10000);
	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 0.12);
	    f_lan_w_exp->SetParameter(3, 2);
	    f_lan_w_exp->SetParameter(4, 0.001);
	    
	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));

	  }
	hsmooth->Fit("lan_w_exp","Q","",local_min, upper);
	h_charge_recalib[ich]->Fit("lan_w_exp","Q","",local_min, upper);

	std::cout <<" try " <<tries <<endl;
	if (f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF() < 5)
	  {
	    break;
	  }
	tries++;
      }
    h_peaks_recalib->Fill(f_lan_w_exp->GetParameter(1));

    gain_corr[ich] = gain_corr[ich]*f_lan_w_exp->GetParameter(1);

    std::cout <<" -- low/guess/high = " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_exp->GetParameter(1)<<endl;
    tn->Fill(ich, gain_corr[ich], ((double) f_lan_w_exp->GetParameter(2)));
  }
  fcalibout->Write();
  fcalibout->Close();

  for (int ich = 0 ; ich < 128; ich++){

    std::cout << "Raw Fitting :: Channel "<<ich;

    f_lan_w_exp->SetParLimits(0, 100, 50000);
    f_lan_w_exp->SetParLimits(1, 0.2, 2.0);
    f_lan_w_exp->SetParLimits(2, 0.1, 0.8);

    
    double local_min = 0;
    double local_max = 0;
    double max_value = 0;
    int good = 0;
    int low_bin = 0;
    int high_bin = 0;
    TH1D *hsmooth = (TH1D*)h_charge[ich]->Clone();
    hsmooth->Smooth();
    hsmooth->Smooth();
    // Finding the dip

    for (int ic = start_looking; ic < hsmooth->GetNbinsX(); ic++)
      {
	good = 0;
	for (int j = -1*min_range; j <= min_range; j++)
	  {
	    if (j == 0) continue;

	    if (hsmooth->GetBinContent(ic+1) <= hsmooth->GetBinContent(ic+1 +j))
	      good++;

	  }
	if (good == 2*min_range)
	  {
	    local_min = hsmooth->GetBinLowEdge(ic - before);
	    low_bin = ic - before;
	    break;
	  }
      }

    // Find local Max after the dip

    for (int ic = low_bin; ic < hsmooth->GetNbinsX(); ic++)
      {
	if (hsmooth->GetBinContent(ic+1) > max_value)
	  {
	    high_bin = ic+1; 
	    max_value = hsmooth->GetBinContent(high_bin);
	    local_max = hsmooth->GetBinLowEdge(high_bin);
	  }
      }
    if (local_max > 2){
      local_max = 1;
      local_min = 0.4;
    }
    else if (local_max < 0.3){
      local_max = 0.5;
      local_min = 0.2;
    }

    f_lan_w_exp->SetParameter(0, 10000);
    f_lan_w_exp->SetParameter(1, local_max);
    f_lan_w_exp->SetParameter(2, 0.12);
    f_lan_w_exp->SetParameter(3, 7.4);
    f_lan_w_exp->SetParameter(4, -0.2);

    
    float upper = local_min + 6*(local_max - local_min);
    
    h_charge[ich]->Fit("lan_w_exp","Q","",local_min, upper);

    h_peaks->Fill(f_lan_w_exp->GetParameter(1));

    std::cout <<" -- low/guess/high = " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_exp->GetParameter(1)<<endl;
  }

  
  // As function of number of tubes hit, RMS of the time dists;
  TFile *fout = new TFile(Form("../plots/mbd_calib_plots_%d.root", runnumber), "RECREATE");
  for (int ich = 0 ; ich < 128; ich++)
    {
      h_charge[ich]->Write();
      h_time[ich]->Write();
      h_charge_raw[ich]->Write();
      h_charge_recalib[ich]->Write();
      h_time_raw[ich]->Write();
    }

  h_peaks->Write();
  h_peaks_raw->Write();
  h_peaks_recalib->Write();
  fout->Close();
}
