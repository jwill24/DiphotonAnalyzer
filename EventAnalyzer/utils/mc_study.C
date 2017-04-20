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
  const float dr_cut = 0.1;

  //const char* file = "output_0000_GammaGammaToGammaGamma_13TeV_fpmc.root";
  const char* file = "output_GammaGammaToGammaGamma_fpmc_test.root";
  TFile f( file );
  TTree* t = dynamic_cast<TTree*>( f.Get( "ntp" ) );
  unsigned int num_dipho;
  const unsigned short max_dipho = 10;
  float photon1_pt[max_dipho], photon1_eta[max_dipho], photon1_phi[max_dipho],
        photon2_pt[max_dipho], photon2_eta[max_dipho], photon2_phi[max_dipho];
  float dipho_vx[max_dipho], dipho_vy[max_dipho], dipho_vz[max_dipho];
  t->SetBranchAddress( "num_diphoton", &num_dipho );
  t->SetBranchAddress( "diphoton_pt1", photon1_pt );
  t->SetBranchAddress( "diphoton_eta1", photon1_eta );
  t->SetBranchAddress( "diphoton_phi1", photon1_phi );
  t->SetBranchAddress( "diphoton_pt2", photon2_pt );
  t->SetBranchAddress( "diphoton_eta2", photon2_eta );
  t->SetBranchAddress( "diphoton_phi2", photon2_phi );
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
  float genpart_vx[max_genpart], genpart_vy[max_genpart], genpart_vz[max_genpart];
  t->SetBranchAddress( "num_gen_part", &num_genpart );
  t->SetBranchAddress( "gen_part_pdgid", genpart_pdgid );
  t->SetBranchAddress( "gen_part_vertex_x", genpart_vx );
  t->SetBranchAddress( "gen_part_vertex_y", genpart_vy );
  t->SetBranchAddress( "gen_part_vertex_z", genpart_vz );
  float pu_weight;
  t->SetBranchAddress( "pileup_weight", &pu_weight );

  TH1D* h_ptgen = new TH1D( "pt_gen", "Single photon p_{T}\\Events\\GeV", 100, 0., 1000. ),
       *h_ptreco = (TH1D*) h_ptgen->Clone( "pt_reco" );
  TH1D* h_ptpair_gen = new TH1D( "ptpair_gen", "Diphoton p_{T}\\Events\\GeV", 100, 0., 50. ),
       *h_ptpair_reco = (TH1D*) h_ptpair_gen->Clone( "ptpair_reco" );
  TH1D* h_ptdiff = new TH1D( "ptdiff", "p_{T}^{reco} - p_{T}^{gen}\\Events\\GeV?.2f", 100, -25., 25. );
  TH2D* h_ptevol = new TH2D( "ptevol", "", 500, 0., 1000., 100, -25., 25. );
  TH1D* h_vtxdist_gen = new TH1D( "vtxdist_gen", "d(generated/reconstructed diphoton vertex)\\Events fraction\\cm?.2f", 100, 0., 25. ),
       *h_vtxdist_gen_purw = (TH1D*) h_vtxdist_gen->Clone( "vtxdist_gen_purw" );
  TH1D* h_match_mass = new TH1D( "match_mass", "(m_{reco}(#gamma#gamma)-m_{gen}(#gamma#gamma))/m_{gen}(#gamma#gamma)\\Events fraction\\?.4f", 200, -0.25, 0.25 ),
       *h_match_ptpair = new TH1D( "match_pt", "(p_{T}^{reco}(#gamma#gamma)-p_{T}^{gen}(#gamma#gamma))/p_{T}^{gen}(#gamma#gamma)\\Events fraction", 100, -10., 10. );

  const unsigned int num_entries = t->GetEntries();
  for ( unsigned long long i=0; i<num_entries; ++i ) {
    t->GetEntry( i );
    bool has_pho1_match = false,
         has_pho2_match = false;

    TVector3 dipho_gen_vtx;
    for ( unsigned int j=0; j<num_genpart; ++j ) {
      if ( genpart_pdgid[j]!=23 ) continue;
      dipho_gen_vtx = TVector3( genpart_vx[j], genpart_vy[j], genpart_vz[j] );
    }

    for ( unsigned int j=0; j<num_dipho; ++j ) {
      TLorentzVector pho1, pho2,
                     pho1_gen, pho2_gen;
      TVector3 dipho_vtx,
               pho1_gen_vtx, pho2_gen_vtx;

      pho1.SetPtEtaPhiM( photon1_pt[j], photon1_eta[j], photon1_phi[j], 0. );
      pho2.SetPtEtaPhiM( photon2_pt[j], photon2_eta[j], photon2_phi[j], 0. );
      dipho_vtx = TVector3( dipho_vx[j], dipho_vy[j], dipho_vz[j] );
      h_ptreco->Fill( photon1_pt[j] );
      for ( unsigned int k=0; k<num_genpho; ++k ) {
        h_ptdiff->Fill( photon1_pt[j]-genpho_pt[k] );
        h_ptevol->Fill( genpho_pt[k], photon1_pt[j]-genpho_pt[k] );

        TLorentzVector genpho; genpho.SetPtEtaPhiE( genpho_pt[k], genpho_eta[k], genpho_phi[k], genpho_energy[k] );
        if ( genpho.DeltaR( pho1 )<dr_cut ) {
          //cout << "match with pho1" << endl;
          pho1_gen = genpho;
          //pho1_gen_vtx = TVector3( genpho_vx[k], genpho_vy[k], genpho_vz[k] );
          has_pho1_match = true;
        }
        if ( genpho.DeltaR( pho2 )<dr_cut ) {
          //cout << "match with pho2" << endl;
          pho2_gen = genpho;
          //pho2_gen_vtx = TVector3( genpho_vx[k], genpho_vy[k], genpho_vz[k] );
          has_pho2_match = true;
        }
        if ( has_pho1_match || has_pho2_match ) {
          int bin_id = -1;
          for ( unsigned int l=0; l<num_pt_bins; l++ ) {
            if ( genpho_pt[k]<pt_bins[l] || genpho_pt[k]>pt_bins[l+1] ) continue;
            bin_id = l;
          }
          if ( bin_id>=0 ) {
            h_ptdist[bin_id]->Fill( ( photon1_pt[j]-genpho_pt[k] )/genpho_pt[k] );
            //cout << "pt=" << genpho_pt[k] << " between " << pt_bins[bin_id-1] << " and " << pt_bins[bin_id] << endl;
          }
        }
      }
      if ( has_pho1_match && has_pho2_match ) {
        //pho1_gen_vtx.Print(); pho2_gen_vtx.Print();
        //dipho_vtx.Print();
        const double reco_gen_dist = ( dipho_vtx-dipho_gen_vtx ).Mag();
        h_vtxdist_gen->Fill( reco_gen_dist, 1./num_entries );
        h_vtxdist_gen_purw->Fill( reco_gen_dist, pu_weight/num_entries );
        const TLorentzVector pair_gen( pho1_gen+pho2_gen ),
                             pair_reco( pho1+pho2 );

        h_ptpair_gen->Fill( pair_gen.Pt() );
        h_ptpair_reco->Fill( pair_reco.Pt() );
        h_match_mass->Fill( ( pair_reco.M()-pair_gen.M() )/pair_gen.M(), pu_weight/num_entries );
        h_match_ptpair->Fill( ( pair_reco.Pt()-pair_gen.Pt() )/pair_gen.Pt(), pu_weight/num_entries );
      }
    }
    for ( unsigned int k=0; k<num_genpho; ++k ) {
      h_ptgen->Fill( genpho_pt[k] );
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
    Canvas c( "mc_study_match_ptpair", top_title );
    h_match_ptpair->Draw( "hist" );
    c.Prettify( h_match_ptpair );
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
    h_vtxdist_gen->Draw();
    c.Prettify( h_vtxdist_gen );
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
    //cumul->SetLineWidth( 2 );
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
  {
    Plotter plt( output_dir, "" );
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_ptpair_gen ) );
    hm.push_back( make_pair( "Reco level", h_ptpair_reco ) );
    plt.draw_multiplot( "mc_study_ptpair_comp", hm );
  }
  {
    Plotter plt( output_dir, "" );
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_ptgen ) );
    hm.push_back( make_pair( "Reco level", h_ptreco ) );
    plt.draw_multiplot( "mc_study_ptcomp", hm );
  }
}
