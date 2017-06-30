#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

#define MAX_PROTONS 30

void xi_accept()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  const float num_sigma = 2.0;

  const map<unsigned short,string> pot_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };
  const map<unsigned short,float> xi_accept = { { 2, 0.033 }, { 3, 0.024 }, { 102, 0.050 }, { 103, 0.037 } };

  const vector<unsigned short> fills_to_study = { 4947, 4985, 5017, 5030 };

  unsigned int fill_number;
  tr->SetBranchAddress( "fill_number", &fill_number );

  unsigned int num_proton_track;
  float proton_track_x[MAX_PROTONS], proton_track_y[MAX_PROTONS];
  unsigned int proton_track_side[MAX_PROTONS], proton_track_pot[MAX_PROTONS];
  float proton_track_normchi2[MAX_PROTONS];
  tr->SetBranchAddress( "num_proton_track", &num_proton_track );
  tr->SetBranchAddress( "proton_track_x", proton_track_x );
  tr->SetBranchAddress( "proton_track_y", proton_track_y );
  tr->SetBranchAddress( "proton_track_side", proton_track_side );
  tr->SetBranchAddress( "proton_track_pot", proton_track_pot );
  tr->SetBranchAddress( "proton_track_normchi2", proton_track_normchi2 );

  typedef map<unsigned short,TH1D*> plots_t;
  map<unsigned short,plots_t> plots_vs_fill;

  for ( unsigned long long i=0; i<tr->GetEntriesFast(); i++ ) {
    tr->GetEntry( i );

    auto align = pot_align::get_alignments( fill_number );
    if ( plots_vs_fill.find( fill_number )==plots_vs_fill.end() ) { // fill info not yet created
      // create the plots for this fill
      plots_t plots;
      for ( const auto& pot : pot_names ) {
        plots[pot.first] = new TH1D( Form( "xi_spectrum_%s_fill%d", pot.second.c_str(), fill_number ), Form( "#xi(%s)\\Events", pot.second.c_str() ), 100, 0., 0.25 );
      }
      plots_vs_fill.insert( make_pair( fill_number, plots ) );
    }
    // retrieve the plots for this fill
    auto& plots = plots_vs_fill[fill_number];

    for ( unsigned short j=0; j<num_proton_track; j++ ) {
      double xi, xi_err;
      const unsigned short pot_id = proton_track_side[j]*100+proton_track_pot[j];
      const auto& al = align[pot_id];
      xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );

      //float xi_cut = xi_accept[pot_id];
      plots[pot_id]->Fill( xi );
      /*if ( xi>=xi_cut ) {
      }
      else { // outside pot acceptance
      }*/
    }
  }

  //----- plotting part

  for ( const auto& pot : xi_accept ) {
    Canvas c( Form( "xiaccept_%d", pot.first ), "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
    THStack st;
    unsigned short i = 0;
    //for ( auto& fill : plots_vs_fill ) {
    for ( auto& fill : fills_to_study ) {
      if ( i>10 ) break; //FIXME
      //auto plot = fill.second[pot.first];
      auto plot = plots_vs_fill[fill][pot.first];
      plot->SetLineWidth( 2 );
      plot->SetLineColor( i+1 );
      st.Add( plot );
      //c.AddLegendEntry( plot, Form( "Fill %d", fill.first ) );
      c.AddLegendEntry( plot, Form( "Fill %d", fill ) );
      if ( i==0 ) st.SetTitle( plot->GetTitle() );
      i++;
    }
    st.Draw( "nostack" );
    c.Prettify( st.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}

