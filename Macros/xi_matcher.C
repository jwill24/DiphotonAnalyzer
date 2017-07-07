#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"

#include <fstream>

#define MAX_PROTONS 30
#define MAX_DIPH 20
#define MAX_ELE 50
#define MAX_MUON 50
#define MAX_JET 200

const float max_xi = 0.25;

bool is_matched( int n_sigma, float xi_rp, float xi_cs, float err_xi_rp, float err_xi_cs )
{
  const double combined_error = sqrt( err_xi_rp*err_xi_rp + err_xi_cs*err_xi_cs );
  const double delta = fabs( xi_cs-xi_rp );

 return ( delta/combined_error<=n_sigma );
}
void fill_matching( float num_sigma, const std::vector< std::pair<float, float> >& tracks, float limit, float xi_gg, float err_xi_gg, TGraphErrors& gr_match, TGraphErrors& gr_unmatch, TGraphErrors& gr_ooa, std::vector< std::pair<float, float> >& candidates );
void plot_matching( const char* name, TGraphErrors& gr_unmatch, TGraphErrors& gr_match, TGraphErrors& gr_ooa, double limits );
void plot_xispectrum( const char* name, TH1D* spec, double limit );
void plot_generic( const char* name, TH1* plot, const char* plot_style="", bool logy=false );
void plot_chisq( const char* name, TH1D* ina, TH1D* ooa );
void plot_xyaccept( const char* name, TH2D* plot, TH2D* plot_ooa );

void xi_matcher()
{
  //TFile f( "Samples/output_Run2016BCG_looseCuts_22jun.root" );
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  //xi_reco::load_file( "TreeProducer/data/ctpps_optics_9mar2017.root" );
  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  const float rel_err_xi_gg = 0.039;
  const float num_sigma = 2.0;

  const float lim_45n = 0.033;
  const float lim_45f = 0.024;
  /*const float lim_56n = 0.050;
  const float lim_56f = 0.037;

  const float lim_45n = 0.034;
  const float lim_45f = 0.023;*/
  const float lim_56n = 0.042;
  const float lim_56f = 0.032;


  /*const float lim_56n = 0.050;
  const float lim_56f = 0.037;*///FIXME!!!

  unsigned int fill_number;
  tr->SetBranchAddress( "fill_number", &fill_number );

  unsigned int num_proton_track;
  float proton_track_x[MAX_PROTONS], proton_track_y[MAX_PROTONS];
  float proton_track_normchi2[MAX_PROTONS];
  unsigned int proton_track_side[MAX_PROTONS], proton_track_pot[MAX_PROTONS];
  tr->SetBranchAddress( "num_proton_track", &num_proton_track );
  tr->SetBranchAddress( "proton_track_x", proton_track_x );
  tr->SetBranchAddress( "proton_track_y", proton_track_y );
  tr->SetBranchAddress( "proton_track_side", proton_track_side );
  tr->SetBranchAddress( "proton_track_pot", proton_track_pot );
  tr->SetBranchAddress( "proton_track_normchi2", proton_track_normchi2 );

  unsigned int num_diphoton;
  float diphoton_pt1[MAX_DIPH], diphoton_pt2[MAX_DIPH];
  float diphoton_eta1[MAX_DIPH], diphoton_eta2[MAX_DIPH];
  float diphoton_r91[MAX_DIPH], diphoton_r92[MAX_DIPH];
  float diphoton_mass[MAX_DIPH], diphoton_pt[MAX_DIPH], diphoton_rapidity[MAX_DIPH], diphoton_dphi[MAX_DIPH];
  float diphoton_vertex_x[MAX_DIPH], diphoton_vertex_y[MAX_DIPH], diphoton_vertex_z[MAX_DIPH];
  tr->SetBranchAddress( "num_diphoton", &num_diphoton );
  tr->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  tr->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  tr->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  tr->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  tr->SetBranchAddress( "diphoton_r91", diphoton_r91 );
  tr->SetBranchAddress( "diphoton_r92", diphoton_r92 );
  tr->SetBranchAddress( "diphoton_mass", diphoton_mass );
  tr->SetBranchAddress( "diphoton_pt", diphoton_pt );
  tr->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
  tr->SetBranchAddress( "diphoton_dphi", diphoton_dphi );
  tr->SetBranchAddress( "diphoton_vertex_x", diphoton_vertex_x );
  tr->SetBranchAddress( "diphoton_vertex_y", diphoton_vertex_y );
  tr->SetBranchAddress( "diphoton_vertex_z", diphoton_vertex_z );

  unsigned int num_electron;
  float electron_vtx_x[MAX_ELE], electron_vtx_y[MAX_ELE], electron_vtx_z[MAX_ELE];
  tr->SetBranchAddress( "num_electron", &num_electron );
  tr->SetBranchAddress( "electron_vtx_x", electron_vtx_x );
  tr->SetBranchAddress( "electron_vtx_y", electron_vtx_y );
  tr->SetBranchAddress( "electron_vtx_z", electron_vtx_z );

  unsigned int num_muon;
  float muon_vtx_x[MAX_ELE], muon_vtx_y[MAX_ELE], muon_vtx_z[MAX_ELE];
  tr->SetBranchAddress( "num_muon", &num_muon );
  tr->SetBranchAddress( "muon_vtx_x", muon_vtx_x );
  tr->SetBranchAddress( "muon_vtx_y", muon_vtx_y );
  tr->SetBranchAddress( "muon_vtx_z", muon_vtx_z );

  unsigned int num_jet;
  float jet_pt[MAX_JET];
  int jet_dipho_match[MAX_JET];
  tr->SetBranchAddress( "num_jet", &num_jet );
  tr->SetBranchAddress( "jet_pt", jet_pt );
  tr->SetBranchAddress( "jet_dipho_match", jet_dipho_match );

  TGraphErrors gr_45n_ooa, gr_45n_unmatch, gr_45n_match;
  TGraphErrors gr_45f_ooa, gr_45f_unmatch, gr_45f_match;
  TGraphErrors gr_56n_ooa, gr_56n_unmatch, gr_56n_match;
  TGraphErrors gr_56f_ooa, gr_56f_unmatch, gr_56f_match;

  TH1D* h_xispectrum_45n = new TH1D( "xisp_45n", "", 100, 0., max_xi ),
       *h_xispectrum_45f = dynamic_cast<TH1D*>( h_xispectrum_45n->Clone( "xisp_45f" ) ),
       *h_xispectrum_56n = dynamic_cast<TH1D*>( h_xispectrum_45n->Clone( "xisp_56n" ) ),
       *h_xispectrum_56f = dynamic_cast<TH1D*>( h_xispectrum_45n->Clone( "xisp_56f" ) );
  TH1D* h_fwdtrk_chisq_45n_ina = new TH1D( "fwdtrk_chisq_45n_ina", "", 100, 0., 60. ),
       *h_fwdtrk_chisq_45n_ooa = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_45n_ooa" ) ),
       *h_fwdtrk_chisq_45f_ina = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_45f_ina" ) ),
       *h_fwdtrk_chisq_45f_ooa = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_45f_ooa" ) ),
       *h_fwdtrk_chisq_56n_ina = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_56n_ina" ) ),
       *h_fwdtrk_chisq_56n_ooa = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_56n_ooa" ) ),
       *h_fwdtrk_chisq_56f_ina = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_56f_ina" ) ),
       *h_fwdtrk_chisq_56f_ooa = dynamic_cast<TH1D*>( h_fwdtrk_chisq_45n_ina->Clone( "fwdtrk_chisq_56f_ooa" ) );
  TH1D* h_nearfar_xi_45 = new TH1D( "nearfar_xi_45", "#xi(45N)-#xi(45F)\\Events", 100, -0.2, 0.2 );
  TH1D* h_nearfar_xi_56 = new TH1D( "nearfar_xi_56", "#xi(56N)-#xi(56F)\\Events", 100, -0.2, 0.2 );
  TH2D* h_nearfar_x_45 = new TH2D( "nearfar_x_45", "x(45N) (cm)\\x(45F) (cm)", 100, 0., 2.5, 100, 0., 2.5 ),
       *h_nearfar_x_45_ooa = dynamic_cast<TH2D*>( h_nearfar_x_45->Clone( "nearfar_x_45_ooa" ) );
  TH2D* h_nearfar_x_56 = new TH2D( "nearfar_x_56", "x(56N) (cm)\\x(56F) (cm)", 100, 0., 2.5, 100, 0., 2.5 ),
       *h_nearfar_x_56_ooa = dynamic_cast<TH2D*>( h_nearfar_x_56->Clone( "nearfar_x_56_ooa" ) );
  TH2D* h_nearfar_y_45 = new TH2D( "nearfar_y_45", "y(45N) (cm)\\y(45F) (cm)", 100, -2.5, 2.5, 100, -2.5, 2.5 ),
       *h_nearfar_y_45_ooa = dynamic_cast<TH2D*>( h_nearfar_y_45->Clone( "nearfar_y_56_ooa" ) );
  TH2D* h_nearfar_y_56 = new TH2D( "nearfar_y_56", "y(56N) (cm)\\y(56F) (cm)", 100, -2.5, 2.5, 100, -2.5, 2.5 ),
       *h_nearfar_y_56_ooa = dynamic_cast<TH2D*>( h_nearfar_y_56->Clone( "nearfar_y_56_ooa" ) );

  TH1D h_mpair( "mpair", "m(#gamma#gamma)\\Events\\GeV", 10, 250., 1500. );
  TH1D h_ypair( "ypair", "Y(#gamma#gamma)\\Events\\?.2f", 10, -2.5, 2.5 );
  TH1D h_ptpair( "ptpair", "p_{T}(#gamma#gamma)\\Events\\GeV", 50, 0., 100. );
  TH1D h_acopl( "acopl", "1-#||{#Delta#phi/#pi}\\Events\\?.3f", 20, 0., 0.008 );

  const unsigned short num_ptcut_bins = 5, num_dist_bins = 4;
  double extrajets_ptcut[num_ptcut_bins] = { 50., 100., 200., 500., 1000. };
  double extraleptons_dist[num_dist_bins] = { 0.1, 0.5, 2.0, 5.0 };
  unsigned short num_extraleptons[num_dist_bins], num_extrajets[num_ptcut_bins];
  TH1D* h_numextrajets[num_ptcut_bins], *h_numextraleptons[5];
  for ( unsigned short i=0; i<num_ptcut_bins; i++ ) h_numextrajets[i] = new TH1D( Form( "num_extrajets_ptgt%.0f", extrajets_ptcut[i] ), Form( "jet p_{T} > %.0f GeV", extrajets_ptcut[i] ), 10, 0., 10. );
  for ( unsigned short i=0; i<num_dist_bins; i++ ) h_numextraleptons[i] = new TH1D( Form( "num_extraleptons_dist%.1f", extraleptons_dist[i] ), Form( "d(e/#mu vertex, #gamma#gamma vertex) < %.1f cm", extraleptons_dist[i] ), 10, 0., 10. );

  ofstream file_cands( "events_list.txt" );

  float min_xi_45n = 1., min_xi_45f = 1., min_xi_56n = 1., min_xi_56f = 1.;
  unsigned short num_ooa_45n = 0, num_ooa_45f = 0, num_ooa_56n = 0, num_ooa_56f = 0;

  const unsigned long long num_events = tr->GetEntriesFast();
  for ( unsigned long long i=0; i<num_events; i++ ) {
    tr->GetEntry( i );
    //cout << "event " << i << ": " << num_proton_track << " proton tracks, " << num_diphoton << " diphoton candidates" << endl;

    // first loop to identify the tracks and their respective pot

    auto align = pot_align::get_alignments( fill_number );

    vector< pair<float, float> > xi_45n, xi_45f, xi_56n, xi_56f;
    vector< pair<float, float> > xy_45n, xy_45f, xy_56n, xy_56f;
    vector< pair<float, float> > xy_45n_ooa, xy_45f_ooa, xy_56n_ooa, xy_56f_ooa;

    for ( unsigned short j=0; j<num_proton_track; j++ ) {
      const unsigned short pot_id = 100*proton_track_side[j]+proton_track_pot[j];
      const auto& al = align[pot_id];

      double xi, xi_err;
      xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
      if ( proton_track_side[j]==0 && proton_track_pot[j]==2 ) {
        if ( xi>=lim_45n ) {
          h_fwdtrk_chisq_45n_ina->Fill( proton_track_normchi2[j] );
        }
        else { // outside pot acceptance
          cout << "45N: " << fill_number << endl;
          h_fwdtrk_chisq_45n_ooa->Fill( proton_track_normchi2[j] );
          xy_45n_ooa.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
          num_ooa_45n++;
        }
        xi_45n.emplace_back( xi, xi_err );
        xy_45n.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
        h_xispectrum_45n->Fill( xi );
//cout << "xi45n: " << proton_track_xi[j] << " +/- " << proton_track_xi_error[j] << " / " << xi << " +/- " << xi_err << endl;
      }
      else if ( proton_track_side[j]==0 && proton_track_pot[j]==3 ) {
        if ( xi>=lim_45f ) {
          h_fwdtrk_chisq_45f_ina->Fill( proton_track_normchi2[j] );
        }
        else { // outside pot acceptance
          cout << "45F: " << fill_number << endl;
          h_fwdtrk_chisq_45f_ooa->Fill( proton_track_normchi2[j] );
          xy_45f_ooa.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
          num_ooa_45f++;
        }
        xi_45f.emplace_back( xi, xi_err );
        xy_45f.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
        h_xispectrum_45f->Fill( xi );
      }
      else if ( proton_track_side[j]==1 && proton_track_pot[j]==2 ) {
        if ( xi>=lim_56n ) {
          h_fwdtrk_chisq_56n_ina->Fill( proton_track_normchi2[j] );
        }
        else { // outside pot acceptance
          cout << "56N: " << fill_number << endl;
          h_fwdtrk_chisq_56n_ooa->Fill( proton_track_normchi2[j] );
          xy_56n_ooa.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
          num_ooa_56n++;
        }
        xi_56n.emplace_back( xi, xi_err );
        xy_56n.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
        h_xispectrum_56n->Fill( xi );
      }
      else if ( proton_track_side[j]==1 && proton_track_pot[j]==3 ) {
        if ( xi>=lim_56f ) {
          h_fwdtrk_chisq_56f_ina->Fill( proton_track_normchi2[j] );
        }
        else { // outside pot acceptance
          cout << "56F: " << fill_number << endl;
          h_fwdtrk_chisq_56f_ooa->Fill( proton_track_normchi2[j] );
          xy_56f_ooa.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
          num_ooa_56f++;
        }
//cout << proton_track_xi[j]-xi << "\txi56f: " << proton_track_xi[j] << " +/- " << proton_track_xi_error[j] << " / " << xi << " +/- " << xi_err << endl;
        xi_56f.emplace_back( xi, xi_err );
        xy_56f.emplace_back( proton_track_x[j]+al.x, proton_track_y[j]-al.y );
        h_xispectrum_56f->Fill( xi );
      }
    }
    for ( const auto& near : xi_45n ) {
      for ( const auto& far : xi_45f ) {
        h_nearfar_xi_45->Fill( near.first-far.first );
      }
    }
    for ( const auto& near : xi_56n ) {
      for ( const auto& far : xi_56f ) {
        h_nearfar_xi_56->Fill( near.first-far.first );
      }
    }
    for ( const auto& near : xy_45n ) {
      for ( const auto& far : xy_45f ) {
        h_nearfar_x_45->Fill( near.first*1.e2, far.first*1.e2 );
        h_nearfar_y_45->Fill( near.second*1.e2, far.second*1.e2 );
      }
    }
    for ( const auto& near : xy_56n ) {
      for ( const auto& far : xy_56f ) {
        h_nearfar_x_56->Fill( near.first*1.e2, far.first*1.e2 );
        h_nearfar_y_56->Fill( near.second*1.e2, far.second*1.e2 );
      }
    }
    for ( const auto& near : xy_45n_ooa ) {
      for ( const auto& far : xy_45f_ooa ) {
        h_nearfar_x_45_ooa->Fill( near.first*1.e2, far.first*1.e2 );
        h_nearfar_y_45_ooa->Fill( near.second*1.e2, far.second*1.e2 );
      }
    }
    for ( const auto& near : xy_56n_ooa ) {
      for ( const auto& far : xy_56f_ooa ) {
        h_nearfar_x_56_ooa->Fill( near.first*1.e2, far.first*1.e2 );
        h_nearfar_y_56_ooa->Fill( near.second*1.e2, far.second*1.e2 );
      }
    }

    // second loop to identify the diphoton candidates and the corresponding xi values
    for ( unsigned short j=0; j<num_diphoton; j++ ) {

      //----- photon quality cuts

      if ( diphoton_pt1[j]<50. ) continue;
      if ( diphoton_pt2[j]<50. ) continue;
      if ( diphoton_eta1[j]>2.5 || ( diphoton_eta1[j]>1.4442 && diphoton_eta1[j]<1.566 ) ) continue;
      if ( diphoton_eta2[j]>2.5 || ( diphoton_eta2[j]>1.4442 && diphoton_eta2[j]<1.566 ) ) continue;
      if ( diphoton_r91[j]<0.94 ) continue;
      if ( diphoton_r92[j]<0.94 ) continue;

      //----- back-to-back photons

      if ( 1.-fabs( diphoton_dphi[j] )/M_PI>0.005 ) continue;

      //----- search for associated leptons

      const double lepton_mindist = 2.0; // in centimeters
      TVector3 diph_vtx( diphoton_vertex_x[j], diphoton_vertex_y[j], diphoton_vertex_z[j] );
      for ( unsigned short k=0; k<num_dist_bins; k++ ) num_extraleptons[k] = 0;
      unsigned short num_close_ele = 0, num_close_muon = 0;
      for ( unsigned short k=0; k<num_electron; k++ ) {
        TVector3 ele_vtx( electron_vtx_x[k], electron_vtx_y[k], electron_vtx_z[k] );
        const float ele_dist = ( ele_vtx-diph_vtx ).Mag();
        for ( unsigned short l=0; l<num_dist_bins; l++ ) if ( ele_dist<extraleptons_dist[l] ) num_extraleptons[l]++;
        if ( ele_dist<lepton_mindist ) num_close_ele++;
      }
      for ( unsigned short k=0; k<num_muon; k++ ) {
        TVector3 mu_vtx( muon_vtx_x[k], muon_vtx_y[k], muon_vtx_z[k] );
        const float mu_dist = ( mu_vtx-diph_vtx ).Mag();
        for ( unsigned short l=0; l<num_dist_bins; l++ ) if ( mu_dist<extraleptons_dist[l] ) num_extraleptons[l]++;
        if ( mu_dist<lepton_mindist ) num_close_muon++;
      }
      for ( unsigned short k=0; k<num_dist_bins; k++ ) h_numextraleptons[k]->Fill( num_extraleptons[k], 1./num_events );

      //----- search for associated jets

      const double min_extrajet_pt = 500.;
      unsigned short num_associated_jet = 0;
      for ( unsigned short k=0; k<num_ptcut_bins; k++ ) num_extrajets[k] = 0;
      for ( unsigned short k=0; k<num_jet; k++ ) {
        if ( jet_dipho_match[k]!=j ) continue;
        for ( unsigned short l=0; l<num_ptcut_bins; l++ ) if ( jet_pt[k]>=extrajets_ptcut[l] ) num_extrajets[l]++;
        if ( jet_pt[k]>min_extrajet_pt ) num_associated_jet++;
      }
      for ( unsigned short k=0; k<num_ptcut_bins; k++ ) h_numextrajets[k]->Fill( num_extrajets[k], 1./num_events );

      //----- another batch of "exclusivity" cuts

      if ( num_close_ele>0 ) continue;
      if ( num_close_muon>0 ) continue;
      if ( num_associated_jet>0 ) continue;

      //----- reconstruct the energy loss from central system

      const float xip_gg = ( diphoton_pt1[j]*exp( +diphoton_eta1[j] ) + diphoton_pt2[j]*exp( +diphoton_eta2[j] ) ) / sqrt_s;
      const float xim_gg = ( diphoton_pt1[j]*exp( -diphoton_eta1[j] ) + diphoton_pt2[j]*exp( -diphoton_eta2[j] ) ) / sqrt_s;

      vector< pair<float, float> > matching_45n, matching_45f, matching_56n, matching_56f;

      fill_matching( num_sigma, xi_45n, lim_45n, xip_gg, xip_gg*rel_err_xi_gg, gr_45n_match, gr_45n_unmatch, gr_45n_ooa, matching_45n );
      fill_matching( num_sigma, xi_45f, lim_45f, xip_gg, xip_gg*rel_err_xi_gg, gr_45f_match, gr_45f_unmatch, gr_45f_ooa, matching_45f );
      fill_matching( num_sigma, xi_56n, lim_56n, xim_gg, xim_gg*rel_err_xi_gg, gr_56n_match, gr_56n_unmatch, gr_56n_ooa, matching_56n );
      fill_matching( num_sigma, xi_56f, lim_56f, xim_gg, xim_gg*rel_err_xi_gg, gr_56f_match, gr_56f_unmatch, gr_56f_ooa, matching_56f );

      // at least one tag (all pots, no distinction on the arm)
      if ( matching_45n.size()>0 || matching_45f.size()>0 || matching_56n.size()>0 || matching_56f.size()>0 ) {
        h_mpair.Fill( diphoton_mass[j] );
        h_ypair.Fill( diphoton_rapidity[j] );
        h_ptpair.Fill( diphoton_pt[j] );
        h_acopl.Fill( 1.0-fabs( diphoton_dphi[j] )/M_PI );
      }

      // at least one tag in each arm (= double tag!)
      if ( ( matching_45n.size()+matching_45f.size()>0 ) && ( matching_56n.size()+matching_56f.size()>0 ) ) {
        cout << "@@@ DOUBLE TAGGING" << endl;
        //cout << "--> " << matching_45n.size() << " / " << matching_45f.size() << " / " << matching_56n.size() << " / " << matching_56f.size() << endl;
        vector<diproton_candidate_t> candidates;
        for ( const auto& p45 : matching_45n ) {
          cout << "45n: " << p45.first << " +/- " << p45.second << endl;
          for ( const auto& p56 : matching_56n ) {
            candidates.emplace_back( p45.first, p45.second, p56.first, p56.second );
            cout << "56n: " << p56.first << " +/- " << p56.second << endl;
          }
          for ( const auto& p56 : matching_56f ) {
            candidates.emplace_back( p45.first, p45.second, p56.first, p56.second );
            cout << "56f: " << p56.first << " +/- " << p56.second << endl;
          }
        }
        for ( const auto& p45 : matching_45f ) {
          cout << "45f: " << p45.first << " +/- " << p45.second << endl;
          for ( const auto& p56 : matching_56n ) {
            candidates.emplace_back( p45.first, p45.second, p56.first, p56.second );
            cout << "56n: " << p56.first << " +/- " << p56.second << endl;
          }
          for ( const auto& p56 : matching_56f ) {
            candidates.emplace_back( p45.first, p45.second, p56.first, p56.second );
            cout << "56f: " << p56.first << " +/- " << p56.second << endl;
          }
        }
        cout << candidates.size() << " candidate(s) in total!" << endl;
        for ( const auto& cand : candidates ) {
          cout << "masses: diphoton: " << diphoton_mass[j] << ", diphoton: " << cand.mass() << " +/- " << cand.mass_error() << endl;
          file_cands << "2\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
        }
      }
      // only one arm tagged
      else if ( ( matching_45n.size()+matching_45f.size()>0 ) || ( matching_56n.size()+matching_56f.size()>0 ) ) {
        file_cands << "1\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
      }

    }
  }
  file_cands.close();

  cout << "outside acceptance: 45N: " << num_ooa_45n << ", 45F: " << num_ooa_45f << ", 56N: " << num_ooa_56n << ", 56F: " << num_ooa_56f << endl;

  cout << "minimum xi expected: 45N: " << lim_45n << ", 45F: " << lim_45f << ", 56N: " << lim_56n << ", 56F: " << lim_56f << endl;
  cout << "minimum xi observed: 45N: " << min_xi_45n << ", 45F: " << min_xi_45f << ", 56N: " << min_xi_56n << ", 56F: " << min_xi_56f << endl;

  //----- plotting part

  gr_45n_match.SetTitle( "45N" );
  gr_45f_match.SetTitle( "45F" );
  gr_56n_match.SetTitle( "56N" );
  gr_56f_match.SetTitle( "56F" );

  plot_matching( "xi_matching_45n", gr_45n_unmatch, gr_45n_match, gr_45n_ooa, lim_45n );
  plot_matching( "xi_matching_45f", gr_45f_unmatch, gr_45f_match, gr_45f_ooa, lim_45f );
  plot_matching( "xi_matching_56n", gr_56n_unmatch, gr_56n_match, gr_56n_ooa, lim_56n );
  plot_matching( "xi_matching_56f", gr_56f_unmatch, gr_56f_match, gr_56f_ooa, lim_56f );

  h_xispectrum_45n->SetTitle( "#xi (45N)\\Events\\?.4f" );
  h_xispectrum_45f->SetTitle( "#xi (45F)\\Events\\?.4f" );
  h_xispectrum_56n->SetTitle( "#xi (56N)\\Events\\?.4f" );
  h_xispectrum_56f->SetTitle( "#xi (56F)\\Events\\?.4f" );

  plot_xispectrum( "xi_spectrum_45n", h_xispectrum_45n, lim_45n );
  plot_xispectrum( "xi_spectrum_45f", h_xispectrum_45f, lim_45f );
  plot_xispectrum( "xi_spectrum_56n", h_xispectrum_56n, lim_56n );
  plot_xispectrum( "xi_spectrum_56f", h_xispectrum_56f, lim_56f );

  plot_generic( "xitagged_mpair", &h_mpair );
  plot_generic( "xitagged_rapidity", &h_ypair );
  plot_generic( "xitagged_ptpair", &h_ptpair );
  plot_generic( "xitagged_acoplanarity", &h_acopl );

  h_fwdtrk_chisq_45n_ina->SetTitle( "45N" );
  h_fwdtrk_chisq_45f_ina->SetTitle( "45F" );
  h_fwdtrk_chisq_56n_ina->SetTitle( "56N" );
  h_fwdtrk_chisq_56f_ina->SetTitle( "56F" );

  plot_chisq( "xi_accept_chisq_45n", h_fwdtrk_chisq_45n_ina, h_fwdtrk_chisq_45n_ooa );
  plot_chisq( "xi_accept_chisq_45f", h_fwdtrk_chisq_45f_ina, h_fwdtrk_chisq_45f_ooa );
  plot_chisq( "xi_accept_chisq_56n", h_fwdtrk_chisq_56n_ina, h_fwdtrk_chisq_56n_ooa );
  plot_chisq( "xi_accept_chisq_56f", h_fwdtrk_chisq_56f_ina, h_fwdtrk_chisq_56f_ooa );

  plot_generic( "xicorr_nearfar_45", h_nearfar_xi_45, "", true );
  plot_generic( "xicorr_nearfar_56", h_nearfar_xi_56, "", true );

  plot_xyaccept( "xi_xcorr_nearfar_45", h_nearfar_x_45, h_nearfar_x_45_ooa );
  plot_xyaccept( "xi_xcorr_nearfar_56", h_nearfar_x_56, h_nearfar_x_56_ooa );
  plot_xyaccept( "xi_ycorr_nearfar_45", h_nearfar_y_45, h_nearfar_y_45_ooa );
  plot_xyaccept( "xi_ycorr_nearfar_56", h_nearfar_y_56, h_nearfar_y_56_ooa );

  int colours[] = { kBlack, kRed+1, kBlue-2, kGreen+2, kOrange };
  {
    Canvas c( "num_extraleptons", "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
    THStack hs;
    c.SetLegendX1( 0.35 );
    double max = -1.;
    for ( unsigned short i=0; i<num_dist_bins; i++ ) {
      h_numextraleptons[i]->Sumw2();
      h_numextraleptons[i]->SetMarkerStyle( 20+i );
      h_numextraleptons[i]->SetMarkerColor( colours[i] );
      h_numextraleptons[i]->SetLineColor( kBlack );
      hs.Add( h_numextraleptons[i] );
      c.AddLegendEntry( h_numextraleptons[i], h_numextraleptons[i]->GetTitle() );
      max = TMath::Max( h_numextraleptons[i]->GetMaximum(), max );
    }
    hs.Draw( "p,nostack" );
    hs.SetTitle( "Number of leptons associated to #gamma#gamma vertex\\Events fraction" );
    c.GetLegend()->SetFillStyle( 0 );
    c.GetLegend()->SetLineWidth( 0 );
    hs.SetMaximum( max*1.1 );
    auto ln = new TLine( 1., 0., 1., max*1.1 );
    ln->SetLineStyle( 2 );
    ln->Draw();
    c.Prettify( hs.GetHistogram() );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    Canvas c( "num_extrajets", "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
    THStack hs;
    double max = -1.;
    for ( unsigned short i=0; i<num_ptcut_bins; i++ ) {
      h_numextrajets[i]->Sumw2();
      h_numextrajets[i]->SetMarkerStyle( 20+i );
      h_numextrajets[i]->SetMarkerColor( colours[i] );
      h_numextrajets[i]->SetLineColor( kBlack );
      hs.Add( h_numextrajets[i] );
      c.AddLegendEntry( h_numextrajets[i], h_numextrajets[i]->GetTitle() );
      max = TMath::Max( h_numextrajets[i]->GetMaximum(), max );
    }
    hs.Draw( "p,nostack" );
    hs.SetTitle( "Number of jets associated to #gamma#gamma vertex\\Events fraction" );
    hs.SetMaximum( max*1.1 );
    auto ln = new TLine( 1., 0., 1., max*1.1 );
    ln->SetLineStyle( 2 );
    ln->Draw();
    c.Prettify( hs.GetHistogram() );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }

}

void
fill_matching( float num_sigma, const std::vector< std::pair<float, float> >& tracks, float limit, float xi_gg, float err_xi_gg, TGraphErrors& gr_match, TGraphErrors& gr_unmatch, TGraphErrors& gr_ooa, std::vector< std::pair<float, float> >& candidates )
{
  for ( const auto& pair : tracks ) {
    if ( xi_gg>limit ) {
      if ( is_matched( num_sigma, pair.first, xi_gg, pair.second, err_xi_gg ) ) {
        unsigned int n = gr_match.GetN();
        gr_match.SetPoint( n, pair.first, xi_gg );
        gr_match.SetPointError( n, pair.second, err_xi_gg );
        candidates.push_back( pair );
      }
      else {
        unsigned int n = gr_unmatch.GetN();
        gr_unmatch.SetPoint( n, pair.first, xi_gg );
        gr_unmatch.SetPointError( n, pair.second, err_xi_gg );
      }
    }
    else {
      unsigned int n = gr_ooa.GetN();
      gr_ooa.SetPoint( n, pair.first, xi_gg );
      gr_ooa.SetPointError( n, pair.second, err_xi_gg );
    }
  }
}

void
plot_matching( const char* name, TGraphErrors& gr_unmatch, TGraphErrors& gr_match, TGraphErrors& gr_ooa, double limits )
{
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );

  TF1 lim( "lim", "x", 0., max_xi );
  lim.SetTitle( "#xi(RP)\\#xi(#gamma#gamma)" );
  lim.Draw("lr+");
  lim.SetLineColor(4);
  lim.SetLineWidth(2);

  int colour = kAzure+10;
  TBox box( 0., 0., limits, max_xi );
  box.SetFillColor( colour );
  box.SetLineColor( kBlack );
  box.SetLineWidth( 1 );
  box.Draw( "l" );

  lim.Draw( "lr,same" );

  gr_unmatch.SetLineWidth(2);
  gr_unmatch.SetMarkerStyle(22);
  gr_unmatch.SetMarkerSize(1.25);
  gr_unmatch.Draw("p same");

  gr_match.SetLineColor(2);
  gr_match.SetLineWidth(2);
  gr_match.SetMarkerColor(2);
  gr_match.SetMarkerStyle(20);
  gr_match.SetMarkerSize(1.25);
  gr_match.Draw( "p same" );

  gr_ooa.SetLineWidth(2);
  gr_ooa.SetMarkerStyle(25);
  gr_ooa.SetMarkerSize(1.25);
  gr_ooa.Draw( "p same" );

  c.SetLegendX1( 0.4 );
  c.AddLegendEntry( &gr_match, "Matching events", "lp" );
  c.AddLegendEntry( &gr_unmatch, "Non-matching events", "lp" );
  c.AddLegendEntry( &gr_ooa, "Out of acceptance events", "lp" );
  c.AddLegendEntry( &box, "No acceptance for RP", "f" );

  PaveText lab( 0.8, 0.2, 0.85, 0.25 );
  lab.SetTextSize( 0.075 );
  lab.SetFillStyle( 0 );
  lab.SetLineWidth( 0 );
  lab.AddText( gr_match.GetTitle() );
  lab.Draw( "same" );

  lim.GetYaxis()->SetRangeUser( 0., max_xi );

  TH2D *axiiis = new TH2D( Form( "axiiis_%s", name ), "", 10, 0, max_xi, 10, 0, max_xi );
  axiiis->Draw( "sameaxis" );
  c.GetLegend()->SetFillStyle( 0 );
  c.GetLegend()->SetLineWidth( 0 );
  c.Prettify( lim.GetHistogram() );

  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

void
plot_xispectrum( const char* name, TH1D* spec, double limit )
{
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
  gStyle->SetStatX( 0.88 );
  gStyle->SetStatY( 0.92 );
  gStyle->SetOptStat( "erm" );
  const double max_y = spec->GetMaximum()*1.1;
  spec->Sumw2();
  spec->Draw( "hist" );
  spec->GetYaxis()->SetRangeUser( 0., max_y );
  spec->SetLineColor( kBlack );
  //spec->SetLineWidth( 2 );
  //spec->SetFillColor( kBlack );
  TBox lim( 0., 0., limit, max_y );
  spec->Draw( "p same" );
  lim.SetFillStyle( 3004 );
  lim.SetFillColor( kBlack );
  lim.SetLineColor( kBlack );
  lim.SetLineWidth( 1 );
  lim.Draw( "l" );
  c.Prettify( spec );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

void
plot_chisq( const char* name, TH1D* ina, TH1D* ooa )
{
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
  THStack st;
  st.Add( ina );
  ina->SetLineColor( kBlack );
  ina->SetLineWidth( 2 );
  ina->SetLineStyle( 2 );
  c.AddLegendEntry( ina, "Within acceptance" );
  st.Add( ooa );
  ooa->SetLineColor( kRed+1 );
  ooa->SetLineWidth( 2 );
  ooa->SetLineStyle( 1 );
  c.AddLegendEntry( ooa, "Outside acceptance" );
  st.Draw( "hist,nostack" );
  st.GetHistogram()->SetTitle( "Local track #chi^{2}/ndf (TOTEM strips)\\Events\\?.3f" );

  PaveText lab( 0.8, 0.45, 0.85, 0.5 );
  lab.SetTextSize( 0.075 );
  lab.SetFillStyle( 0 );
  lab.SetLineWidth( 0 );
  lab.AddText( ina->GetTitle() );
  lab.Draw( "same" );

  c.SetLogy();
  c.Prettify( st.GetHistogram() );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

void
plot_generic( const char* name, TH1* plot, const char* plot_style, bool logy )
{
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
  plot->Sumw2();
  plot->Draw( plot_style );
  plot->SetMarkerStyle( 20 );
  plot->SetLineColor( kBlack );
  if ( logy ) c.SetLogy();
  c.Prettify( plot );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

void
plot_xyaccept( const char* name, TH2D* plot, TH2D* plot_ooa )
{
  gStyle->SetOptStat( 0 );
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
  plot->Draw( "colz" );
  plot_ooa->Draw( "p same" );
  plot_ooa->SetMarkerStyle( 24 );
  plot_ooa->SetMarkerSize( 1.4 );
  plot_ooa->SetMarkerColor( kRed );
  c.SetLegendX1( 0.15 );
  c.SetLegendY1( 0.82 );
  c.AddLegendEntry( plot_ooa, "Out of #xi acceptance", "p" );
  c.GetLegend()->SetFillStyle( 0 );
  c.GetLegend()->SetLineWidth( 0 );
  c.Prettify( plot );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

