#include "Plotter.h"

typedef struct track_t
{
  track_t( float x, float y, float xi, float xi_error ) :
    x( x ), y( y ), xi( xi ), xi_error( xi_error ) {}
  track_t() :
    x( 0. ), y( 0. ), xi( 0. ), xi_error( 0. ) {}
  float x, y;
  float xi, xi_error;
} track_t;
typedef pair<track_t,track_t> candidate_t;

void
xi_cutter( TString file="Samples/output_Run2016BCG_looseCuts_10mar_xifix.root" )
{
  const char* out_path = "/afs/cern.ch/user/l/lforthom/www/private/twophoton/xi_study";

  const float xidiff_cuts[] = { 0.6, 0.2, 0.1 };
  const size_t num_xi_cuts = sizeof( xidiff_cuts )/sizeof( xidiff_cuts[0] );

  TFile in_file( file );
  TTree* tr = dynamic_cast<TTree*>( in_file.Get( "ntp" ) );
  const unsigned short max_num_tracks = 20;
  unsigned int num_tracks, track_side[max_num_tracks], track_pot[max_num_tracks];
  float track_x[max_num_tracks], track_y[max_num_tracks],
        track_xi[max_num_tracks], err_track_xi[max_num_tracks];
  tr->SetBranchAddress( "num_proton_track", &num_tracks );
  tr->SetBranchAddress( "proton_track_side", track_side );
  tr->SetBranchAddress( "proton_track_pot", track_pot );
  tr->SetBranchAddress( "proton_track_x", track_x );
  tr->SetBranchAddress( "proton_track_y", track_y );
  tr->SetBranchAddress( "proton_track_xi", track_xi );
  tr->SetBranchAddress( "proton_track_xi_error", err_track_xi );

  //enum pot { near45 = 2, far45 = 3, near56 = 102, far56 = 103 };
  enum sector { sector45 = 0, sector56 = 1 };
  enum pot { nearPot = 2, farPot = 3 };
  TString potName[4] = { "45n", "45f", "56n", "56f" };

  TH1D* h_per_pot[4],
       *h_per_pot_cut[num_xi_cuts][4];
  TH2D* h_per_pot_hitmap[4],
       *h_per_pot_hitmap_cut[num_xi_cuts][4];
  for ( unsigned int i=0; i<4; i++ ) {
    h_per_pot[i] = new TH1D( Form( "per_pot_%s", potName[i].Data() ), Form( "Tracks' #xi (sector %s)\\Events\\?.3f", potName[i].Data() ), 70, 0.03, 0.38 ),
    h_per_pot_hitmap[i] = new TH2D( Form( "per_pot_hm_%s", potName[i].Data() ), "Track x (cm)\\Track y (cm)", 100, 0., 5., 100, -2.5, 2.5 );
    for ( unsigned int j=0; j<num_xi_cuts; j++ ) {
      h_per_pot_cut[j][i] = dynamic_cast<TH1D*>( h_per_pot[i]->Clone( Form( "per_pot_cut%d_%s", j, potName[i].Data() ) ) );
      h_per_pot_hitmap_cut[j][i] = dynamic_cast<TH2D*>( h_per_pot_hitmap[i]->Clone( Form( "per_pot_hm_cut%d_%s", j, potName[i].Data() ) ) );
    }
  }
  TH2D* h_corr_x_45 = new TH2D( "corr_nf_45", "Track x position (near pot)\\Track x position (far pot)", 100, 0., 2.5, 100, 0., 2.5 ),
       *h_corr_x_56 = dynamic_cast<TH2D*>( h_corr_x_45->Clone( "corr_nf_56" ) );
  TH1D* h_diff_x_45 = new TH1D( "diff_nf_45", "Near track x - far track x\\Events\\cm", 100, -1., 1. ),
       *h_diff_x_56 = dynamic_cast<TH1D*>( h_diff_x_45->Clone( "diff_nf_56" ) );

  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );

    vector<track_t> track_45n, track_45f, track_56n, track_56f;
    for ( unsigned short j=0; j<num_tracks; j++ ) {
      if ( track_side[j]==sector45 ) {
        if ( track_pot[j]==nearPot ) {
          track_t trk( track_x[j], track_y[j], track_xi[j], err_track_xi[j] );
          track_45n.push_back( trk );
          h_per_pot[0]->Fill( trk.xi );
          h_per_pot_hitmap[0]->Fill( trk.x*100, trk.y*100 );
        }
        if ( track_pot[j]==farPot ) {
          track_t trk( track_x[j], track_y[j], track_xi[j], err_track_xi[j] );
          track_45f.push_back( trk );
          h_per_pot[1]->Fill( trk.xi );
          h_per_pot_hitmap[1]->Fill( trk.x*100, trk.y*100 );
        }
      }
      if ( track_side[j]==sector56 ) {
        if ( track_pot[j]==nearPot ) {
          track_t trk( track_x[j], track_y[j], track_xi[j], err_track_xi[j] );
          track_56n.push_back( trk );
          h_per_pot[2]->Fill( trk.xi );
          h_per_pot_hitmap[2]->Fill( trk.x*100, trk.y*100 );
        }
        if ( track_pot[j]==farPot ) {
          track_t trk( track_x[j], track_y[j], track_xi[j], err_track_xi[j] );
          track_56f.push_back( trk );
          h_per_pot[3]->Fill( trk.xi );
          h_per_pot_hitmap[3]->Fill( trk.x*100, trk.y*100 );
        }
      }
    }

    vector<candidate_t> cand_45_cut[num_xi_cuts], cand_56_cut[num_xi_cuts];

    //----- sector 45

    for ( vector<track_t>::const_iterator itn=track_45n.begin(); itn!=track_45n.end(); itn++ ) {
      for ( vector<track_t>::const_iterator itf=track_45f.begin(); itf!=track_45f.end(); itf++ ) {
        h_corr_x_45->Fill( itn->x*100, itf->x*100 );
        h_diff_x_45->Fill( ( itn->x-itf->x )*100 );
        for ( unsigned short j=0; j<num_xi_cuts; j++ ) {
          if ( fabs( itn->x-itf->x )<xidiff_cuts[j] ) { cand_45_cut[j].push_back( make_pair( *itn, *itf ) ); }
        }
      }
    }
    for ( unsigned short j=0; j<num_xi_cuts; j++ ) {
      for ( vector<candidate_t>::const_iterator it=cand_45_cut[j].begin(); it!=cand_45_cut[j].end(); it++ ) {
        h_per_pot_hitmap_cut[j][0]->Fill( it->first.x*100, it->first.y*100 );
        h_per_pot_cut[j][0]->Fill( it->first.xi );
        h_per_pot_hitmap_cut[j][1]->Fill( it->second.x*100, it->second.y*100 );
        h_per_pot_cut[j][1]->Fill( it->second.xi );
      }
    }

    //----- sector 56

    for ( vector<track_t>::const_iterator itn=track_56n.begin(); itn!=track_56n.end(); itn++ ) {
      for ( vector<track_t>::const_iterator itf=track_56f.begin(); itf!=track_56f.end(); itf++ ) {
        h_corr_x_56->Fill( itn->x*100, itf->x*100 );
        h_diff_x_56->Fill( ( itn->x-itf->x )*100 );
        for ( unsigned short j=0; j<num_xi_cuts; j++ ) {
          if ( fabs( itn->x-itf->x )<xidiff_cuts[j] ) { cand_56_cut[j].push_back( make_pair( *itn, *itf ) ); }
        }
        //if ( itn->x-itf->x>-0.1 ) {}
      }
    }
    for ( unsigned short j=0; j<num_xi_cuts; j++ ) {
      for ( vector<candidate_t>::const_iterator it=cand_56_cut[j].begin(); it!=cand_56_cut[j].end(); it++ ) {
        h_per_pot_hitmap_cut[j][2]->Fill( it->first.x*100, it->first.y*100 );
        h_per_pot_cut[j][2]->Fill( it->first.xi );
        h_per_pot_hitmap_cut[j][3]->Fill( it->second.x*100, it->second.y*100 );
        h_per_pot_cut[j][3]->Fill( it->second.xi );
      }
    }

  } // loop over events

  //----- drawing part

  gStyle->SetOptStat( 0 );

  Plotter plt( out_path, "" );
  /*{
    Plotter::HistsMap hm;
    for ( unsigned short i=0; i<4; i++ ) {
      hm.push_back( make_pair( potName[i], h_per_pot[i] ) );
    }
    plt.draw_4plots( "xistudy_per_pot", hm );
  }*/
  for ( unsigned short i=0; i<4; i++ ) {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "All tracks", h_per_pot[i] ) );
    for ( unsigned short j=0; j<num_xi_cuts; j++ ) {
      hm.push_back( make_pair( Form( "|#xi_{N}-#xi_{F}|<%.1f", xidiff_cuts[j] ), h_per_pot_cut[j][i] ) );
    }
    plt.plot_multihists( Form( "xistudy_per_pot_%s", potName[i].Data() ), hm, 0., 1.24 );
  }
  {
    Plotter::HistsMap hm;
    for ( unsigned short i=0; i<4; i++ ) {
      hm.push_back( make_pair( potName[i], h_per_pot_hitmap[i] ) );
    }
    plt.draw_4plots( "xistudy_per_pot_hitmap", hm );
  }
  for ( unsigned short i=0; i<num_xi_cuts; i++ ) {
    Plotter::HistsMap hm;
    for ( unsigned short j=0; j<4; j++ ) {
      hm.push_back( make_pair( potName[j], h_per_pot_hitmap_cut[i][j] ) );
    }
    plt.draw_4plots( Form( "xistudy_per_pot_hitmap_cut%d", i+1 ), hm );
  }
  {
    Canvas c( "xistudy_corr_nf_45" );
    h_corr_x_45->Draw( "colz" );
    c.Prettify( h_corr_x_45 );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "xistudy_corr_nf_56" );
    h_corr_x_56->Draw( "colz" );
    c.Prettify( h_corr_x_56 );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "xistudy_diff_nf_45" );
    h_diff_x_45->Draw();
    c.Prettify( h_diff_x_45 );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "xistudy_diff_nf_56" );
    h_diff_x_56->Draw();
    c.Prettify( h_diff_x_56 );
    c.Save( "pdf,png", out_path );
  }
}
