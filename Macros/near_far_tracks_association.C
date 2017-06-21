#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStyle.h"

#include "Canvas.h"
#include "Plotter.h"

#include <fstream>
#include <iostream>
#include <map>

#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp/"
#define default_ntp_file "Samples/output_Run2016B_mgg-ov-500.root"

void near_far_tracks_association( TString file=default_ntp_file )
{
  TFile f(file);
  if ( !f.IsOpen() ) return;

  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );
  // general quantities
  unsigned int run_id, lumisection;
  unsigned long long event_number;
  tr->SetBranchAddress( "run_id", &run_id );
  tr->SetBranchAddress( "lumisection", &lumisection );
  tr->SetBranchAddress( "event_number", &event_number );
  // diphoton quantities
  const unsigned short max_diph = 20;
  unsigned int num_diphoton;
  float diphoton_vertex_x[max_diph], diphoton_vertex_y[max_diph], diphoton_vertex_z[max_diph];
  tr->SetBranchAddress( "diphoton_vertex_x", diphoton_vertex_x );
  tr->SetBranchAddress( "diphoton_vertex_y", diphoton_vertex_y );
  tr->SetBranchAddress( "diphoton_vertex_z", diphoton_vertex_z );
  // proton quantities
  unsigned int num_proton;
  const unsigned short max_pr = 20;
  unsigned int proton_side[max_pr], proton_pot[max_pr];
  float proton_xi[max_pr], proton_xi_error[max_pr];
  unsigned int proton_link_id[max_pr];
  float proton_link_dist[max_pr];
  tr->SetBranchAddress( "num_proton_track", &num_proton );
  tr->SetBranchAddress( "proton_track_side", proton_side );
  tr->SetBranchAddress( "proton_track_pot", proton_pot );
  tr->SetBranchAddress( "proton_track_xi", proton_xi );
  tr->SetBranchAddress( "proton_track_xi_error", proton_xi_error );
  tr->SetBranchAddress( "proton_track_link_nearfar", proton_link_id );
  tr->SetBranchAddress( "proton_track_link_mindist", proton_link_dist );
  // vertex quantities
  unsigned int num_vertex;
  tr->SetBranchAddress( "num_vertex", &num_vertex );

  TH1D* h_num_proton = new TH1D( "num_proton", "Forward track multiplicity\\Events", 6, 0., 6. ),
       *h_num_vtx = new TH1D( "num_vertex", "Vertex multiplicity\\Events", 40, 0., 40. );
  TGraphErrors h_ximatch_45nf_thr0p1, h_ximatch_45nf_thr0p2, h_ximatch_45nf_thr0p5, h_ximatch_45nf_thr1p0,
               h_ximatch_56nf_thr0p1, h_ximatch_56nf_thr0p2, h_ximatch_56nf_thr0p5, h_ximatch_56nf_thr1p0;

  unsigned int num_evts_notag = 0, num_evts_with_tag = 0;
  TLorentzVector pho1, pho2, electron, muon, jet;
  // tree readout stage
  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );
    // dump the list of events in a text file

    //----- proton tracks retrieval part -----

    h_num_proton->Fill( num_proton );

    // xi / error
    typedef std::pair<float, float> trackparam_t;
    typedef std::map<unsigned int, trackparam_t> tracks_map;
    tracks_map tracks_45_near, tracks_45_far, tracks_56_near, tracks_56_far;

    // first loop to identify the RP and unflatten the collection
    for ( unsigned int j=0; j<num_proton; j++ ) {
      switch ( proton_side[j] ) {
        case 0: { // 4-5
          switch ( proton_pot[j] ) {
            case 2: { tracks_45_near.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // near pot
            case 3: { tracks_45_far.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // far pot
          }
        } break;
        case 1: { // 5-6
          switch ( proton_pot[j] ) {
            case 2: { tracks_56_near.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // near pot
            case 3: { tracks_56_far.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // far pot
          }
        } break;
      }
    }

    const float sqs = 13.e3;
    float xi_56 = -1., err_xi_56 = -1.,
          xi_45 = -1., err_xi_45 = -1.;
    float max_diproton_mass = -1., max_diproton_mass_error = -1.,
          max_diproton_mass_rap = -999.;

    if ( tracks_45_near.size()>0 && tracks_56_near.size()>0 ) {
      for ( tracks_map::const_iterator it_45=tracks_45_near.begin(); it_45!=tracks_45_near.end(); it_45++ ) {
        const trackparam_t xi45 = it_45->second;
        for ( tracks_map::const_iterator it_56=tracks_56_near.begin(); it_56!=tracks_56_near.end(); it_56++ ) {
          const trackparam_t xi56 = it_56->second;
          const float diproton_mass = sqs * sqrt( xi45.first * xi56.first ),
                      err_diproton_mass = sqrt( pow( xi45.second, 2 ) + pow( xi56.second, 2 ) ) * diproton_mass/2.;
          if ( diproton_mass>max_diproton_mass ) {
            max_diproton_mass = diproton_mass; max_diproton_mass_error = err_diproton_mass;
            max_diproton_mass_rap = log( xi56.first/xi45.first );
            xi_45 = xi45.first; err_xi_45 = xi45.second;
            xi_56 = xi56.first; err_xi_56 = xi56.second;
          }
        }
      }
    }

    const bool has_singletag = ( tracks_45_near.size()>0 || tracks_45_far.size()>0 || tracks_56_near.size()>0 || tracks_56_far.size()>0 ),
               has_doubletag = ( ( tracks_45_near.size()>0 || tracks_45_far.size()>0 ) && ( tracks_56_near.size()>0 || tracks_56_far.size()>0 ) );


    // near-far matching test
    for ( tracks_map::const_iterator it=tracks_45_near.begin(); it!=tracks_45_near.end(); it++ ) {
      const unsigned int linked_id = proton_link_id[it->first];
      if ( linked_id>100 ) continue;
      if ( proton_link_dist[it->first]<0.1 ) { int id = h_ximatch_45nf_thr0p1.GetN(); h_ximatch_45nf_thr0p1.SetPoint( id, it->second.first, tracks_45_far[linked_id].first ); h_ximatch_45nf_thr0p1.SetPointError( id, it->second.second, tracks_45_far[linked_id].second ); }
      if ( proton_link_dist[it->first]<0.2 ) { int id = h_ximatch_45nf_thr0p2.GetN(); h_ximatch_45nf_thr0p2.SetPoint( id, it->second.first, tracks_45_far[linked_id].first ); h_ximatch_45nf_thr0p2.SetPointError( id, it->second.second, tracks_45_far[linked_id].second ); }
      if ( proton_link_dist[it->first]<0.5 ) { int id = h_ximatch_45nf_thr0p5.GetN(); h_ximatch_45nf_thr0p5.SetPoint( id, it->second.first, tracks_45_far[linked_id].first ); h_ximatch_45nf_thr0p5.SetPointError( id, it->second.second, tracks_45_far[linked_id].second ); }
      if ( proton_link_dist[it->first]<1.0 ) { int id = h_ximatch_45nf_thr1p0.GetN(); h_ximatch_45nf_thr1p0.SetPoint( id, it->second.first, tracks_45_far[linked_id].first ); h_ximatch_45nf_thr1p0.SetPointError( id, it->second.second, tracks_45_far[linked_id].second ); }
    }
    for ( tracks_map::const_iterator it=tracks_56_near.begin(); it!=tracks_56_near.end(); it++ ) {
      const unsigned int linked_id = proton_link_id[it->first];
      if ( linked_id>100 ) continue;
      if ( proton_link_dist[it->first]<0.1 ) { int id = h_ximatch_56nf_thr0p1.GetN(); h_ximatch_56nf_thr0p1.SetPoint( id, it->second.first, tracks_56_far[linked_id].first ); h_ximatch_56nf_thr0p1.SetPointError( id, it->second.second, tracks_56_far[linked_id].second ); }
      if ( proton_link_dist[it->first]<0.2 ) { int id = h_ximatch_56nf_thr0p2.GetN(); h_ximatch_56nf_thr0p2.SetPoint( id, it->second.first, tracks_56_far[linked_id].first ); h_ximatch_56nf_thr0p2.SetPointError( id, it->second.second, tracks_56_far[linked_id].second ); }
      if ( proton_link_dist[it->first]<0.5 ) { int id = h_ximatch_56nf_thr0p5.GetN(); h_ximatch_56nf_thr0p5.SetPoint( id, it->second.first, tracks_56_far[linked_id].first ); h_ximatch_56nf_thr0p5.SetPointError( id, it->second.second, tracks_56_far[linked_id].second ); }
      if ( proton_link_dist[it->first]<1.0 ) { int id = h_ximatch_56nf_thr1p0.GetN(); h_ximatch_56nf_thr1p0.SetPoint( id, it->second.first, tracks_56_far[linked_id].first ); h_ximatch_56nf_thr1p0.SetPointError( id, it->second.second, tracks_56_far[linked_id].second ); }
    }
    h_num_vtx->Fill( num_vertex );

  }

  // plotting stage

  gStyle->SetOptStat( 0 );

  const string top_label_str = "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV";
  const char* top_label = top_label_str.c_str();

  {
    Plotter plt( out_path, top_label );

    // sector 45
    Plotter::GraphsMap gm_45;
    gm_45.insert( std::make_pair( "d #leq 0.1 cm", &h_ximatch_45nf_thr0p1 ) );
    gm_45.insert( std::make_pair( "d #leq 0.2 cm", &h_ximatch_45nf_thr0p2 ) );
    gm_45.insert( std::make_pair( "d #leq 0.5 cm", &h_ximatch_45nf_thr0p5 ) );
    gm_45.insert( std::make_pair( "d #leq 1.0 cm", &h_ximatch_45nf_thr1p0 ) );
    plt.plot_xi_correlations( "45", gm_45 );
    // sector 56
    Plotter::GraphsMap gm_56;
    gm_56.insert( std::make_pair( "d #leq 0.1 cm", &h_ximatch_56nf_thr0p1 ) );
    gm_56.insert( std::make_pair( "d #leq 0.2 cm", &h_ximatch_56nf_thr0p2 ) );
    gm_56.insert( std::make_pair( "d #leq 0.5 cm", &h_ximatch_56nf_thr0p5 ) );
    gm_56.insert( std::make_pair( "d #leq 1.0 cm", &h_ximatch_56nf_thr1p0 ) );
    plt.plot_xi_correlations( "56", gm_56 );
  }

  {
    Canvas c( "num_proton", top_label );
    h_num_proton->Sumw2();
    h_num_proton->Draw();
    h_num_proton->SetMarkerStyle( 20 );
    h_num_proton->SetLineColor( kBlack );
    c.Prettify( h_num_proton );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
}
