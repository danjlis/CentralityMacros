

void DrawMBDTime(int runnumber = 21813)
{
  gStyle->SetOptStat(0);
  float thresh = 0.5;
  float timethresh = 15;
  float timethreshs = 12;

  TFile *file = new TFile(Form("../output/run%d/trees_%d.root", runnumber, runnumber), "r");

  // 3 vertex cuts with centralities between them compared.
  
  TH1D *h_time_n = new TH1D("h_time_n","", 100, -25, 25);
  TH1D *h_time_s = new TH1D("h_time_s","", 100, -25, 25);
  TH1D *h_tubes_n = new TH1D("h_tubes_n","", 65, -0.5, 64.5);
  TH1D *h_tubes_s = new TH1D("h_tubes_s","", 65, -0.5, 64.5);
  TProfile *hp_rms_by_hits_n = new TProfile("hp_rms_by_hits_n", "", 63, .5, 64.5);
  TProfile *hp_rms_by_hits_s = new TProfile("hp_rms_by_hits_s", "", 63, .5, 64.5);

  TH1D *distance_from_mean[128];
  for (int i = 0; i < 128; i++)
    {
      distance_from_mean[i] = new TH1D(Form("distance_from_mean_%d", i), "", 100, 0, 25);
    }

  // 0 - mean
  // 1 - median
  
  TProfile *vertex_by_mean = new TProfile("vertex_by_mean","", 65, -0.5, 64.5,"s");
  TProfile *vertex_by_median = new TProfile("vertex_by_median","", 65, -0.5, 64.5,"s");
  TProfile *vertex_by_mean_cover = new TProfile("vertex_by_mean_cover","", 65, -0.5, 64.5,"s");
  TProfile *vertex_by_median_cover = new TProfile("vertex_by_median_cover","", 65, -0.5, 64.5,"s");

  TProfile2D *vertex_by_mean2 = new TProfile2D("vertex_by_mean2","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TProfile2D *vertex_by_median2 = new TProfile2D("vertex_by_median2","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TProfile2D *vertex_by_mean_cover2 = new TProfile2D("vertex_by_mean_cover2","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TProfile2D *vertex_by_median_cover2 = new TProfile2D("vertex_by_median_cover2","", 65, -0.5, 64.5, 65, -0.5, 64.5);

  TProfile *diff_vertex_methods = new TProfile("diff_vertex_methods","", 65, -0.5, 64.5);
  TProfile *diff_vertex_methods_cover = new TProfile("diff_vertex_methods_cover","", 65, -0.5, 64.5);
  TH1D *h_vertex_mean = new TH1D("h_vertex_mean", "", 501, -250, 250);
  TH1D *h_vertex_median = new TH1D("h_vertex_median", "", 501, -250, 250);
  TH1D *h_vertex_mean_cover = new TH1D("h_vertex_mean_cover", "", 501, -250, 250);
  TH1D *h_vertex_median_cover = new TH1D("h_vertex_median_cover", "", 501, -250, 250);

  TH2D *nhit_dist = new TH2D("nhit_dist","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TH2D *nhit_dist_cover = new TH2D("nhit_dist_cover","", 65, -0.5, 64.5, 65, -0.5, 64.5);


  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];
  int isMinBias;
  float z_vertex;
  
  TTree *t = (TTree*)file->Get("T");
  
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge);
  t->SetBranchAddress("mbd_time_raw",mbd_time);
  t->SetBranchAddress("z_vertex",&z_vertex);

  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      h_time_n->Reset();
      h_time_s->Reset();
      int hits_n = 0;
      int hits_n_cover = 0;
      int hits_s = 0;      
      std::vector<float> time_sum_n;
      std::vector<float> time_sum_n_cover;
      std::vector<float> time_sum_s;
      float sum_n = 0.;
      float sum_n_cover = 0.;
      float sum_s = 0.;
      for (int ich = 0 ; ich < 64; ich++)
	{
	  
	  if (mbd_charge[ich] > thresh && mbd_time[ich] < timethresh)
	    {
	      hits_n++;
	      h_time_n->Fill(mbd_time[ich]);
	      time_sum_n.push_back(mbd_time[ich]);
	      sum_n += mbd_time[ich];
	      if (!(ich == 56 || ich == 54))
		{
		  time_sum_n_cover.push_back(mbd_time[ich]);
		  sum_n_cover += mbd_time[ich];
		  hits_n_cover++;
		}
	    }      
	  if (mbd_charge[ich + 64] > thresh && mbd_time[ich+64] < timethreshs)
	    {
	      hits_s++;
	      h_time_s->Fill(mbd_time[ich + 64]);
	      time_sum_s.push_back(mbd_time[ich+64]);
	      sum_s += mbd_time[ich+64];
	    }

       
	  
	} 
      nhit_dist->Fill(hits_s, hits_n);
      h_tubes_s->Fill(hits_s);
      h_tubes_n->Fill(hits_n);

      sort(time_sum_n.begin(), time_sum_n.end());
      sort(time_sum_n_cover.begin(), time_sum_n_cover.end());
      sort(time_sum_s.begin(), time_sum_s.end());

      if (hits_s >=2 && hits_n >=2){
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

	float vertex_median = 15*(median_n - median_s);
	float vertex_mean = 15*(sum_n/static_cast<float>(hits_n) - (sum_s/static_cast<float>(hits_s)));
	vertex_by_median->Fill(min(hits_s, hits_n), vertex_median);
	vertex_by_mean->Fill(min(hits_s, hits_n), vertex_mean);
 	vertex_by_median2->Fill(hits_s, hits_n, vertex_median);
	vertex_by_mean2->Fill(hits_s, hits_n, vertex_mean);

	h_vertex_mean->Fill(vertex_mean);
	h_vertex_median->Fill(vertex_median);
	diff_vertex_methods->Fill(min(hits_s, hits_n), vertex_median - vertex_mean);
      }

      if (hits_s >=2 && hits_n_cover >=2){
	float median_n_cover;
	if (time_sum_n_cover.size()%2)
	  {
	    median_n_cover = time_sum_n_cover.at(time_sum_n_cover.size()/2);
	  }
	else
	  {
	    median_n_cover = (time_sum_n_cover.at(time_sum_n_cover.size()/2) + time_sum_n_cover.at(time_sum_n_cover.size()/2 - 1))/2.;
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

	float vertex_median_cover = 15*(median_n_cover - median_s);
	float vertex_mean_cover = 15*(sum_n_cover/static_cast<float>(hits_n_cover) - (sum_s/static_cast<float>(hits_s)));

	vertex_by_median_cover->Fill(min(hits_s, hits_n_cover), vertex_median_cover);
	vertex_by_mean_cover->Fill(min(hits_s, hits_n_cover), vertex_mean_cover);
	vertex_by_median_cover2->Fill(hits_s, hits_n_cover, vertex_median_cover);
	vertex_by_mean_cover2->Fill(hits_s, hits_n_cover, vertex_mean_cover);
	nhit_dist_cover->Fill(hits_s, hits_n_cover);
	diff_vertex_methods_cover->Fill(min(hits_s, hits_n_cover), vertex_median_cover - vertex_mean_cover);
	h_vertex_mean_cover->Fill(vertex_mean_cover);
	h_vertex_median_cover->Fill(vertex_median_cover);

      }

      for (int ich = 0 ; ich < 64; ich++)
	{
	  
	  if (mbd_charge[ich] > thresh && mbd_time[ich] < timethresh)
	    {
	      distance_from_mean[ich]->Fill(TMath::Abs(h_time_n->GetMean() - mbd_time[ich]));
	    }      
	  if (mbd_charge[ich + 64] > thresh && mbd_time[ich+64] < timethreshs)
	    {
	      distance_from_mean[ich+64]->Fill(TMath::Abs(h_time_n->GetMean() - mbd_time[ich+64]));
	    }
	  
	  
	} 
      
      hp_rms_by_hits_s->Fill(hits_s, h_time_s->GetRMS());
      hp_rms_by_hits_n->Fill(hits_n, h_time_n->GetRMS());
    }

  TCanvas *c1 = new TCanvas("c1", "c1");
 
  hp_rms_by_hits_n->Draw("hist");
  hp_rms_by_hits_s->SetLineColor(kRed);
  hp_rms_by_hits_s->Draw("hist same");


  TCanvas *c2 = new TCanvas("c2", "c2");
  
  h_vertex_mean->Draw("hist");
  gPad->SetLogy();

  h_vertex_median->SetLineColor(kRed);
  h_vertex_median->Draw("hist same");
  h_vertex_mean_cover->SetLineStyle(4);
  h_vertex_mean_cover->Draw("hist same");
  h_vertex_median_cover->SetLineStyle(4);
  h_vertex_median_cover->SetLineColor(kRed);
  h_vertex_median_cover->Draw("hist same");
  
  TH1D *h_half_mast = new TH1D("h_half_mast","", 128, -0.5, 127.5);

  for (int i = 0; i < 128;i++)
    {
      for (int j = 0; j < distance_from_mean[i]->GetNbinsX();j++)
	{
	  if (distance_from_mean[i]->GetBinContent(j+1) < distance_from_mean[i]->GetBinContent( distance_from_mean[i]->GetMaximumBin())/2.)
	    {
	      h_half_mast->Fill(i, distance_from_mean[i]->GetBinLowEdge(j));
	      break;
	    }
	}

    }


  TFile *out = new TFile(Form("mbd_vertex_plots_%d.root", runnumber), "RECREATE");
  hp_rms_by_hits_n->Write();
  hp_rms_by_hits_s->Write();
  vertex_by_mean->Write();
  vertex_by_median->Write();
  vertex_by_mean_cover->Write();
  vertex_by_median_cover->Write();

  vertex_by_mean2->Write();
  vertex_by_median2->Write();
  vertex_by_mean_cover2->Write();
  vertex_by_median_cover2->Write();

  h_vertex_mean->Write();
  h_vertex_median->Write();
  h_vertex_mean_cover->Write();
  h_vertex_median_cover->Write();

  h_tubes_s->Write();
  h_tubes_n->Write();
  diff_vertex_methods->Write();
  nhit_dist->Write();
  nhit_dist_cover->Write();
  h_half_mast->Write();
  for (int i = 0; i < 128; i++)
    {
      
      distance_from_mean[i]->Write();
    }
  out->Close();

}
