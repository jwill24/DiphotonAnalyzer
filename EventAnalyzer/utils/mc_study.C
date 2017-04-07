#include "Canvas.h"
#include "Plotter.h"

void mc_study()
{
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

  TH1D* h_ptgen = new TH1D( "pt_gen", "p_{T} (gen)\\Events", 100, 0., 1000. ),
       *h_ptreco = (TH1D*) h_ptgen->Clone( "pt_reco" );
  TH1D* h_ptdiff = new TH1D( "ptdiff", "p_{T} (reco)-p_{T} (gen)\\Events\\GeV", 100, -25., 25. );
  TH1D* h_vtxdist_gen = new TH1D( "vtxdist_gen", "d(generated/reconstructed diphoton vertex)\\Events fraction\\cm?.2f", 100, 0., 25. ),
       *h_vtxdist_gen_purw = (TH1D*) h_vtxdist_gen->Clone( "vtxdist_gen_purw" );

  for ( unsigned long long i=0; i<t->GetEntries(); ++i ) {
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
        TLorentzVector genpho; genpho.SetPtEtaPhiE( genpho_pt[k], genpho_eta[k], genpho_phi[k], genpho_energy[k] );
        if ( genpho.DeltaR( pho1 )<0.5 ) {
          //cout << "match with pho1" << endl;
          pho1_gen = genpho;
          //pho1_gen_vtx = TVector3( genpho_vx[k], genpho_vy[k], genpho_vz[k] );
          has_pho1_match = true;
        }
        if ( genpho.DeltaR( pho2 )<0.5 ) {
          //cout << "match with pho2" << endl;
          pho2_gen = genpho;
          //pho2_gen_vtx = TVector3( genpho_vx[k], genpho_vy[k], genpho_vz[k] );
          has_pho2_match = true;
        }
      }
      if ( has_pho1_match && has_pho2_match ) {
        //pho1_gen_vtx.Print(); pho2_gen_vtx.Print();
        //dipho_vtx.Print();
        const double reco_gen_dist = ( dipho_vtx-dipho_gen_vtx ).Mag();
        h_vtxdist_gen->Fill( reco_gen_dist, 1./t->GetEntries() );
        h_vtxdist_gen_purw->Fill( reco_gen_dist, pu_weight/t->GetEntries() );
      }
    }
    for ( unsigned int k=0; k<num_genpho; ++k ) {
      h_ptgen->Fill( genpho_pt[k] );
    }
  }
  gStyle->SetOptStat( 0 );
  {
    Canvas c( "mc_study_ptdiff" );
    h_ptdiff->Draw();
    c.Prettify( h_ptdiff );
    c.Save( "png" );
  }
  {
    Canvas c( "mc_study_genvtx_dist", "CMS Simulation - Excl. #gamma#gamma #rightarrow #gamma#gamma, #sqrt{s} = 13 TeV, 2016 PU cond." );
    h_vtxdist_gen->Draw();
    c.Prettify( h_vtxdist_gen );
    c.Save( "png,pdf" );
  }
  {
    Plotter plt( ".", "CMS Simulation - Excl. #gamma#gamma #rightarrow #gamma#gamma, #sqrt{s} = 13 TeV, 2016 PU cond." );
    TH1D* cumul = dynamic_cast<TH1D*>( h_vtxdist_gen->GetCumulative( true ) ),
         *cumul_purw = dynamic_cast<TH1D*>( h_vtxdist_gen_purw->GetCumulative( true ) );
    cumul->SetTitle( "d(generated/reconstructed diphoton vertex)\\Events fraction (cumulative)\\cm?.2f" );
    Plotter::HistsMap hm;
    plt.BaseCanvas()->SetLegendY1( 0.1 );
    hm.push_back( make_pair( "PU reweighting", cumul_purw ) );
    hm.push_back( make_pair( "No PU reweighting", cumul ) );
    plt.draw_multiplot( "mc_study_genvtx_dist_cumul", hm );
  }
  {
    Plotter plt( ".", "" );
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen-level", h_ptgen ) );
    hm.push_back( make_pair( "Reco-level", h_ptreco ) );
    plt.draw_multiplot( "mc_study_ptcomp", hm );
  }
}
