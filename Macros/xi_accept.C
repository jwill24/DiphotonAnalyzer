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

const float max_xi = 0.25;

void plot_xispectrum( const char* name, TH1D* spec, double limit );

void xi_accept()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  //xi_reco::load_file( "TreeProducer/data/ctpps_optics_9mar2017.root" );
  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  const float rel_err_xi_gg = 0.028;
  const float num_sigma = 2.0;

  const float lim_45n = 0.033;
  const float lim_45f = 0.024;
  const float lim_56n = 0.050;
  const float lim_56f = 0.037;

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

  for ( unsigned long long i=0; i<tr->GetEntriesFast(); i++ ) {
    tr->GetEntry( i );

    auto align = pot_align::get_alignments( fill_number );

    vector< pair<float, float> > xi_45n, xi_45f, xi_56n, xi_56f;

    for ( unsigned short j=0; j<num_proton_track; j++ ) {
      double xi, xi_err;
      const unsigned short pot_id = proton_track_side[j]*100+proton_track_pot[j];
      const auto& al = align[pot_id];
      xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );

      if ( proton_track_side[j]==0 && proton_track_pot[j]==2 ) {
        if ( xi>=lim_45n ) {
        }
        else { // outside pot acceptance
        }
      }
      else if ( proton_track_side[j]==0 && proton_track_pot[j]==3 ) {
        if ( xi>=lim_45f ) {
        }
        else { // outside pot acceptance
        }
      }
      else if ( proton_track_side[j]==1 && proton_track_pot[j]==2 ) {
        if ( xi>=lim_56n ) {
        }
        else { // outside pot acceptance
        }
      }
      else if ( proton_track_side[j]==1 && proton_track_pot[j]==3 ) {
        if ( xi>=lim_56f ) {
        }
        else { // outside pot acceptance
        }
      }
    }

  }

  //----- plotting part

  h_xispectrum_45n->SetTitle( "#xi (45N)\\Events\\?.4f" );
  h_xispectrum_45f->SetTitle( "#xi (45F)\\Events\\?.4f" );
  h_xispectrum_56n->SetTitle( "#xi (56N)\\Events\\?.4f" );
  h_xispectrum_56f->SetTitle( "#xi (56F)\\Events\\?.4f" );

  plot_xispectrum( "xi_spectrum_45n", h_xispectrum_45n, lim_45n );
  plot_xispectrum( "xi_spectrum_45f", h_xispectrum_45f, lim_45f );
  plot_xispectrum( "xi_spectrum_56n", h_xispectrum_56n, lim_56n );
  plot_xispectrum( "xi_spectrum_56f", h_xispectrum_56f, lim_56f );

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
  lim.SetFillStyle( 3003 );
  lim.SetFillColor( kBlack );
  lim.SetLineColor( kBlack );
  lim.SetLineWidth( 1 );
  lim.Draw( "l" );
  c.Prettify( spec );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
