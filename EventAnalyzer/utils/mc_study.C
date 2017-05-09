#include "Canvas.h"
#include "Plotter.h"

void mc_study()
{
  const double pt_bins[] = { 50., 150., 500., 1000. };
  const unsigned short num_pt_bins = sizeof( pt_bins )/sizeof( pt_bins[0] )-1;
  double mean_pt_bins[num_pt_bins], err_pt_bins[num_pt_bins];
  TH1D* h_ptdist[num_pt_bins];
  for ( unsigned short i=0; i<num_pt_bins; i++ ) {
    h_ptdist[i] = new TH1D( Form( "ptdist_%d", i ), "", 100, -5., 5. );
    mean_pt_bins[i] = ( pt_bins[i+1]+pt_bins[i] )/2.;
    err_pt_bins[i] = ( pt_bins[i+1]-pt_bins[i] )/2.;
  }

  const double diph_pt_bins[] = { 0., 2., 10., 20., 50. };
  const unsigned short num_diph_pt_bins = sizeof( diph_pt_bins )/sizeof( diph_pt_bins[0] )-1;
  double mean_diph_pt_bins[num_diph_pt_bins], err_diph_pt_bins[num_diph_pt_bins];
  for ( unsigned short i=0; i<num_diph_pt_bins; i++ ) {
    mean_diph_pt_bins[i] = ( diph_pt_bins[i+1]+diph_pt_bins[i] )/2.;
    err_diph_pt_bins[i] = ( diph_pt_bins[i+1]-diph_pt_bins[i] )/2.;
  }
  TH1D* h_diph_ptdist_numer = new TH1D( "diph_ptdist_numer", "", num_diph_pt_bins, diph_pt_bins ),
       *h_diph_ptdist_denom = (TH1D*)h_diph_ptdist_numer->Clone( "diph_ptdist_denom" );

  const float dr_cut = 0.1;

  //const char* file = "output_0000_GammaGammaToGammaGamma_13TeV_fpmc.root";
  const char* file = "output_GammaGammaToGammaGamma_fpmc_v2.root";
  //const char* file = "output_DiPhotonJetsBox_MGG-80toInf_Sherpa_v3.root";
  TFile f( file );
  TTree* t = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  unsigned int num_dipho;
  const unsigned short max_dipho = 50;
  float photon1_pt[max_dipho], photon1_eta[max_dipho], photon1_phi[max_dipho], photon1_energy[max_dipho],
        photon2_pt[max_dipho], photon2_eta[max_dipho], photon2_phi[max_dipho], photon2_energy[max_dipho];
  float gen_photon1_pt[max_dipho], gen_photon1_eta[max_dipho], gen_photon1_phi[max_dipho], gen_photon1_energy[max_dipho],
        gen_photon2_pt[max_dipho], gen_photon2_eta[max_dipho], gen_photon2_phi[max_dipho], gen_photon2_energy[max_dipho];
  float dipho_vx[max_dipho], dipho_vy[max_dipho], dipho_vz[max_dipho];
  t->SetBranchAddress( "num_diphoton", &num_dipho );
  t->SetBranchAddress( "diphoton_pt1", photon1_pt );
  t->SetBranchAddress( "diphoton_eta1", photon1_eta );
  t->SetBranchAddress( "diphoton_phi1", photon1_phi );
  t->SetBranchAddress( "diphoton_energy1", photon1_energy );
  t->SetBranchAddress( "diphoton_pt2", photon2_pt );
  t->SetBranchAddress( "diphoton_eta2", photon2_eta );
  t->SetBranchAddress( "diphoton_phi2", photon2_phi );
  t->SetBranchAddress( "diphoton_energy2", photon2_energy );
  t->SetBranchAddress( "diphoton_genpt1", gen_photon1_pt );
  t->SetBranchAddress( "diphoton_geneta1", gen_photon1_eta );
  t->SetBranchAddress( "diphoton_genphi1", gen_photon1_phi );
  t->SetBranchAddress( "diphoton_genenergy1", gen_photon1_energy );
  t->SetBranchAddress( "diphoton_genpt2", gen_photon2_pt );
  t->SetBranchAddress( "diphoton_geneta2", gen_photon2_eta );
  t->SetBranchAddress( "diphoton_genphi2", gen_photon2_phi );
  t->SetBranchAddress( "diphoton_genenergy2", gen_photon2_energy );
  t->SetBranchAddress( "diphoton_vertex_x", dipho_vx );
  t->SetBranchAddress( "diphoton_vertex_y", dipho_vy );
  t->SetBranchAddress( "diphoton_vertex_z", dipho_vz );

  unsigned int num_genpho;
  const unsigned short max_genpho = 10;
  float genpho_pt[max_genpho], genpho_eta[max_genpho], genpho_phi[max_genpho], genpho_energy[max_genpho];
  t->SetBranchAddress( "num_gen_photon", &num_genpho );
  t->SetBranchAddress( "gen_photon_pt", genpho_pt );
  t->SetBranchAddress( "gen_photon_eta", genpho_eta );
  t->SetBranchAddress( "gen_photon_phi", genpho_phi );
  t->SetBranchAddress( "gen_photon_energy", genpho_energy );

  const unsigned short max_genpart = 10;
  unsigned int num_genpart;
  int genpart_pdgid[max_genpart];
  t->SetBranchAddress( "num_gen_part", &num_genpart );
  t->SetBranchAddress( "gen_part_pdgid", genpart_pdgid );

  float dipho_vtx_smear_x, dipho_vtx_smear_y, dipho_vtx_smear_z;
  t->SetBranchAddress( "diphoton_genvertex_smeared_x", &dipho_vtx_smear_x );
  t->SetBranchAddress( "diphoton_genvertex_smeared_y", &dipho_vtx_smear_y );
  t->SetBranchAddress( "diphoton_genvertex_smeared_z", &dipho_vtx_smear_z );
  float pu_weight;
  t->SetBranchAddress( "pileup_weight", &pu_weight );

  TH1D* h_ptsingle_gen = new TH1D( "pt_gen", "Single photon p_{T}\\Events\\GeV", 100, 0., 1000. ),
       *h_ptsingle_reco = (TH1D*) h_ptsingle_gen->Clone( "pt_reco" );
  TH1D* h_ptpair_gen = new TH1D( "ptpair_gen", "Diphoton p_{T}\\Events\\GeV", 100, 0., 50. ),
       *h_ptpair_reco = (TH1D*) h_ptpair_gen->Clone( "ptpair_reco" );
  TH1D* h_etasingle_gen = new TH1D( "etasingle_gen", "Leading photon #eta\\Events\\?.2f", 50, -2.5, 2.5 ),
       *h_etasingle_reco = (TH1D*) h_etasingle_gen->Clone( "etasingle_reco" );
  TH1D* h_phisingle_gen = new TH1D( "phisingle_gen", "Leading photon #phi\\Events\\?.2f", 20, -TMath::Pi(), TMath::Pi() ),
       *h_phisingle_reco = (TH1D*) h_phisingle_gen->Clone( "phisingle_reco" );
  TH1D* h_energysingle_gen = new TH1D( "energysingle_gen", "Leading photon energy\\Events\\GeV?.2f", 50, 0., 500. ),
       *h_energysingle_reco = (TH1D*) h_energysingle_gen->Clone( "energysingle_reco" );
  TH1D* h_ptdiff = new TH1D( "ptdiff", "p_{T}^{reco} - p_{T}^{gen}\\Events\\GeV?.2f", 100, -25., 25. );
  TH2D* h_ptevol = new TH2D( "ptevol", "Generated leading photon p_{T}\\Reconstructed - generated leading photon p_{T}", 250, 0., 1000., 100, -25., 25. );
  TH1D* h_vtxdist_gen = new TH1D( "vtxdist_gen", "d(generated/reconstructed diphoton vertex)\\Events fraction\\cm?.2f", 200, 0., 20. ),
       *h_vtxdist_gen_purw = (TH1D*) h_vtxdist_gen->Clone( "vtxdist_gen_purw" );
  TH1D* h_match_mass = new TH1D( "match_mass", "(m_{reco}(#gamma#gamma)-m_{gen}(#gamma#gamma))/m_{gen}(#gamma#gamma)\\Events fraction\\?.4f", 200, -5, 5 ),
       *h_match_ptpair = new TH1D( "match_pt", "(p_{T}^{reco}(#gamma#gamma)-p_{T}^{gen}(#gamma#gamma))/p_{T}^{gen}(#gamma#gamma)\\Events fraction", 100, -15., 15. );
  TH1D* h_diphvtxz_gen = new TH1D( "diphvtx_gen", "Diphoton vertex z\\Events\\cm", 60, -30., 30. ),
       *h_diphvtxz_reco = (TH1D*)h_diphvtxz_gen->Clone( "diphvtx_reco" );

  const unsigned int num_entries = t->GetEntries();
  for ( unsigned long long i=0; i<num_entries; ++i ) {
    t->GetEntry( i );
    bool has_pho1_match = false,
         has_pho2_match = false;

    for ( unsigned int j=0; j<num_dipho; ++j ) {
      TLorentzVector pho1, pho2,
                     pho1_gen, pho2_gen;
      TVector3 dipho_vtx( dipho_vx[j], dipho_vy[j], dipho_vz[j] ),
               dipho_gen_vtx( dipho_vtx_smear_x, dipho_vtx_smear_y, dipho_vtx_smear_z );

      /*pho1.SetPtEtaPhiM( photon1_pt[j], photon1_eta[j], photon1_phi[j], 0. );
      pho2.SetPtEtaPhiM( photon2_pt[j], photon2_eta[j], photon2_phi[j], 0. );*///FIXME
      pho1.SetPtEtaPhiE( photon1_pt[j], photon1_eta[j], photon1_phi[j], photon1_energy[j] );
      pho2.SetPtEtaPhiE( photon2_pt[j], photon2_eta[j], photon2_phi[j], photon2_energy[j] );
      pho1_gen.SetPtEtaPhiE( gen_photon1_pt[j], gen_photon1_eta[j], gen_photon1_phi[j], gen_photon1_energy[j] );
      pho2_gen.SetPtEtaPhiE( gen_photon2_pt[j], gen_photon2_eta[j], gen_photon2_phi[j], gen_photon2_energy[j] );

      //cout << pho1.M() << " | " << pho1_gen.M() << "\t" << pho2.M() << " | " << pho2_gen.M() << endl;

      h_ptsingle_reco->Fill( photon1_pt[j] );
      h_ptsingle_gen->Fill( gen_photon1_pt[j] );
      h_etasingle_reco->Fill( photon1_eta[j] );
      h_etasingle_gen->Fill( gen_photon1_eta[j] );
      h_phisingle_reco->Fill( photon1_phi[j] );
      h_phisingle_gen->Fill( gen_photon1_phi[j] );
      h_energysingle_reco->Fill( photon1_energy[j] );
      h_energysingle_gen->Fill( gen_photon1_energy[j] );

      h_ptdiff->Fill( photon1_pt[j]-gen_photon1_pt[j] );
      h_ptevol->Fill( gen_photon1_pt[j], photon1_pt[j]-gen_photon1_pt[j] );

      int bin_id = -1;
      for ( unsigned int l=0; l<num_pt_bins; l++ ) {
        if ( gen_photon1_pt[j]<pt_bins[l] || gen_photon1_pt[j]>pt_bins[l+1] ) continue;
        bin_id = l;
      }
      if ( bin_id>=0 ) {
        h_ptdist[bin_id]->Fill( ( photon1_pt[j]-gen_photon1_pt[j] )/gen_photon1_pt[j] );
      }

      h_diphvtxz_gen->Fill( dipho_gen_vtx.z() );
      h_diphvtxz_reco->Fill( dipho_vtx.z() );

      const double reco_gen_dist = ( dipho_vtx-dipho_gen_vtx ).Mag();
      //const double reco_gen_dist = fabs( dipho_vz[j]-dipho_vtx_smear_z ); //FIXME
//cout << reco_gen_dist << "\t" << dipho_vz[j] << "\t" << dipho_vtx_smear_z << endl;
      h_vtxdist_gen->Fill( reco_gen_dist, 1./num_entries );
      h_vtxdist_gen_purw->Fill( reco_gen_dist, pu_weight/num_entries );

      const TLorentzVector pair_gen( pho1_gen+pho2_gen ),
                           pair_reco( pho1+pho2 );

      h_diph_ptdist_numer->Fill( pair_reco.Pt() );
      h_diph_ptdist_denom->Fill( pair_gen.Pt() );
      
      h_ptpair_gen->Fill( pair_gen.Pt() );
      h_ptpair_reco->Fill( pair_reco.Pt() );
      h_match_mass->Fill( ( pair_reco.M()-pair_gen.M() )/pair_gen.M(), pu_weight/num_entries );
      h_match_ptpair->Fill( ( pair_reco.Pt()-pair_gen.Pt() )/pair_gen.Pt(), pu_weight/num_entries );
    }
  }
  gStyle->SetOptStat( 0 );
  const char* output_dir = "/afs/cern.ch/user/l/lforthom/www/private/twophoton/mc_study",
             *top_title = "CMS Simulation - Excl. #gamma#gamma #rightarrow #gamma#gamma, #sqrt{s} = 13 TeV, 2016 PU cond.";
  {
    Canvas c( "mc_study_ptevol", top_title );
    h_ptevol->Draw( "colz" );
    c.Prettify( h_ptevol );
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_ptresol", top_title );
    double pt_error_mean[num_pt_bins], pt_error_rms[num_pt_bins];
    for ( unsigned short i=0; i<num_pt_bins; i++ ) {
      pt_error_mean[i] = h_ptdist[i]->GetMean();
      pt_error_rms[i] = h_ptdist[i]->GetRMS();
    }
    TGraphErrors gr( num_pt_bins, mean_pt_bins, pt_error_mean, err_pt_bins, pt_error_rms );
    gr.Draw( "ap" );
    gr.SetMarkerStyle( 20 );
    gr.SetTitle( "Single photon p_{T} region\\(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}" );
    gr.SetLineWidth( 2 );
    gr.GetXaxis()->SetRangeUser( pt_bins[0], pt_bins[num_pt_bins] );
    c.Prettify( gr.GetHistogram() );
    c.SetGrid();
    c.Save( "pdf,png", output_dir );
  }
  {
    gStyle->SetOptStat( 111111 );
    Canvas c( "mc_study_match_mass", top_title );
    h_match_mass->Draw( "hist" );
    c.Prettify( h_match_mass );
    c.Save( "pdf,png", output_dir );
    gStyle->SetOptStat( 0 );
  }
  {
    Canvas c( "mc_study_ptpair_ratio", top_title );
    h_diph_ptdist_numer->Divide( h_diph_ptdist_denom );
    h_diph_ptdist_numer->Sumw2();
    h_diph_ptdist_numer->SetMarkerStyle( 20 );
    h_diph_ptdist_numer->Draw( "p" );
    c.Prettify( h_diph_ptdist_numer );
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_match_ptpair", top_title );
    h_match_ptpair->Draw( "hist" );
    c.Prettify( h_match_ptpair );
    c.SetLogy();
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_ptdiff", top_title );
    h_ptdiff->Draw();
    c.Prettify( h_ptdiff );
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_genvtx_dist", top_title );
    h_vtxdist_gen->Draw( "hist" );
    c.Prettify( h_vtxdist_gen );
    c.SetLogy();
    c.Save( "png,pdf", output_dir );
  }
  {
    //Plotter plt( output_dir, "CMS Simulation - Excl. #gamma#gamma #rightarrow #gamma#gamma, #sqrt{s} = 13 TeV, 2016 PU cond." );
    bool ltor = true;
    TH1D* cumul = dynamic_cast<TH1D*>( h_vtxdist_gen->GetCumulative( ltor ) ),
         *cumul_purw = dynamic_cast<TH1D*>( h_vtxdist_gen_purw->GetCumulative( ltor ) );
    Canvas c( "mc_study_genvtx_dist_cumul", top_title );
    THStack st;
    cumul->SetLineColor( kBlack );
    cumul->SetLineWidth( 2 );
    cumul->SetFillColor( kGray );
    cumul->SetFillStyle( 0 );
    st.Add( cumul, "hist" );
    c.AddLegendEntry( cumul, "No PU reweighting", "f" );
    //cumul_purw->SetMarkerStyle( 24 );
    cumul_purw->SetLineColor( kRed+1 );
    cumul_purw->SetLineWidth( 2 );
    st.Add( cumul_purw, "hist" );
    c.AddLegendEntry( cumul_purw, "PU reweighting", "l" );
    st.Draw( "nostack" );
    st.SetTitle( "d(diphoton vtx^{gen}/diphoton vtx^{reco})\\Events fraction (cumulative)\\cm?.2f" );
    //c.SetLegendY1( 0.1 );
    st.SetMaximum( 1.35 );
    c.Prettify( st.GetHistogram() );
    st.SetTitle("");
    c.SetGrid();
    c.Save( "pdf,png", output_dir );
    //Plotter::HistsMap hm;
    //plt.BaseCanvas()->SetLegendY1( 0.1 );
    //hm.push_back( make_pair( "No PU reweighting", cumul ) );
    //hm.push_back( make_pair( "PU reweighting", cumul_purw ) );
    //plt.draw_multiplot( "mc_study_genvtx_dist_cumul", hm );
  }
  Plotter plt( output_dir, "" );
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_ptpair_gen ) );
    hm.push_back( make_pair( "Reco level", h_ptpair_reco ) );
    plt.draw_multiplot( "mc_study_ptpair_comp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_diphvtxz_gen ) );
    hm.push_back( make_pair( "Reco level", h_diphvtxz_reco ) );
    plt.draw_multiplot( "mc_study_diphoton_vtxz", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_ptsingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_ptsingle_reco ) );
    plt.draw_multiplot( "mc_study_ptcomp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_etasingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_etasingle_reco ) );
    plt.draw_multiplot( "mc_study_etacomp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_phisingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_phisingle_reco ) );
    plt.draw_multiplot( "mc_study_phicomp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_energysingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_energysingle_reco ) );
    plt.draw_multiplot( "mc_study_energycomp", hm );
  }
}
