#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"

#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

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
struct track_t {
  track_t() : xi( 0. ), err_xi( 0. ), x( 0. ), y( 0. ) {}
  track_t( float xi, float err_xi, float x, float y ) : xi( xi ), err_xi( err_xi ), x( x ), y( y ) {}
  float xi, err_xi, x, y;
};

bool is_matched( int n_sigma, float x1, float x2, float err_x1, float err_x2 )
{
  const double combined_error = sqrt( err_x1*err_x1 + err_x2*err_x2 );
  const double delta = fabs( x1-x2 );

 return ( delta/combined_error<=n_sigma );
}
void plot_matching( const char* name, TGraphErrors&, TGraphErrors&, TGraphErrors&, TGraphErrors&, double min, double max, double xleg=0.5 );
vector<pair<float,float> > merge_nearfar( const vector<track_t>& near_tracks, const vector<track_t>& far_tracks, float xdiff_cut=0.01 );

map<string,float> pots_accept = { { "45N", 0.033 }, { "45F", 0.024 }, { "56N", 0.050 }, { "56F", 0.037 } };

void massrap_matcher()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  //const float rel_err_xi_gg = 0.039;
  const float num_sigma = 2.0;

  TreeEvent ev;
  ev.attach( tr, true );

  TGraphErrors gr_mass_massrapmatch, gr_mass_massmatch, gr_mass_rapmatch, gr_mass_nomatch;
  TGraphErrors gr_rap_massrapmatch, gr_rap_massmatch, gr_rap_rapmatch, gr_rap_nomatch;

  unsigned int num_massmatch = 0, num_rapmatch = 0, num_massrapmatch = 0, num_nomatch = 0;

  TH1D* h_mass_all = new TH1D( "mass_all", "Diproton missing mass@@Events@@GeV", 12, 200., 2000. );
  TH1D* h_rap_all = new TH1D( "rap_all", "Diproton rapidity@@Events", 20, -1., 1. );

  const unsigned long long num_events = tr->GetEntriesFast();
  for ( unsigned long long i = 0; i < num_events; ++i ) {
    tr->GetEntry( i );
    //cout << "event " << i << ": " << ev.num_proton_track << " proton tracks, " << ev.num_diphoton << " diphoton candidates" << endl;

    // first loop to identify the tracks and their respective pot

    auto align = pot_align::get_alignments( ev.fill_number );

    vector<track_t> xi_45n, xi_45f, xi_56n, xi_56f;

    for ( unsigned short j = 0; j < ev.num_proton_track; ++j ) {
      const unsigned short pot_id = 100*ev.proton_track_side[j]+ev.proton_track_pot[j];
      const auto& al = align[pot_id];

      //----- reconstruct the kinematics
      double xi, xi_err;
      xi_reco::reconstruct( ev.proton_track_x[j]+al.x, ev.proton_track_side[j], ev.proton_track_pot[j], xi, xi_err );

      //----- associate each track to a RP
      if      ( ev.proton_track_side[j] == 0 && ev.proton_track_pot[j] == 2 ) xi_45n.emplace_back( xi, xi_err, ev.proton_track_x[j]+al.x, ev.proton_track_y[j]-al.y );
      else if ( ev.proton_track_side[j] == 0 && ev.proton_track_pot[j] == 3 ) xi_45f.emplace_back( xi, xi_err, ev.proton_track_x[j]+al.x, ev.proton_track_y[j]-al.y );
      else if ( ev.proton_track_side[j] == 1 && ev.proton_track_pot[j] == 2 ) xi_56n.emplace_back( xi, xi_err, ev.proton_track_x[j]+al.x, ev.proton_track_y[j]-al.y );
      else if ( ev.proton_track_side[j] == 1 && ev.proton_track_pot[j] == 3 ) xi_56f.emplace_back( xi, xi_err, ev.proton_track_x[j]+al.x, ev.proton_track_y[j]-al.y );
    }

    //----- merge 2 tracks in one if N-F pot content is similar
    vector<pair<float,float> > xi_45, xi_56;

    const float xdiff_cut = 0.01;

    //--- sector 45

    xi_45 = merge_nearfar( xi_45n, xi_45f, xdiff_cut );
    xi_56 = merge_nearfar( xi_56n, xi_56f, xdiff_cut );

    /*//FIXME FIXME
    for ( const auto& trk : xi_45n ) xi_45.emplace_back( trk.xi, trk.err_xi );
    for ( const auto& trk : xi_45f ) xi_45.emplace_back( trk.xi, trk.err_xi );
    for ( const auto& trk : xi_56n ) xi_56.emplace_back( trk.xi, trk.err_xi );
    for ( const auto& trk : xi_56f ) xi_56.emplace_back( trk.xi, trk.err_xi );
    //FIXME FIXME*/

    //---- identify the diproton candidates

    vector<diproton_candidate_t> candidates;
    for ( const auto trk45 : xi_45 ) {
      for ( const auto trk56 : xi_56 ) {
        candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
      }
    }
    //cout << candidates.size() << " diproton candidate(s) in total!" << endl;

    //----- identify the diphoton candidates

    for ( unsigned short j = 0; j < ev.num_diphoton; ++j ) {

      //----- photon quality cuts

      if ( ev.diphoton_pt1[j] < 50. ) continue;
      if ( ev.diphoton_pt2[j] < 50. ) continue;
      if ( ev.diphoton_eta1[j] > 2.5 || ( ev.diphoton_eta1[j] > 1.4442 && ev.diphoton_eta1[j] < 1.566 ) ) continue;
      if ( ev.diphoton_eta2[j] > 2.5 || ( ev.diphoton_eta2[j] > 1.4442 && ev.diphoton_eta2[j] < 1.566 ) ) continue;
      if ( ev.diphoton_r91[j] < 0.94 ) continue;
      if ( ev.diphoton_r92[j] < 0.94 ) continue;
      if ( ev.diphoton_mass[j] < 350. ) continue;

      //----- back-to-back photons

      if ( 1.-fabs( ev.diphoton_dphi[j] )/M_PI > 0.005 ) continue;

      const float xip = ( ev.diphoton_pt1[j]*exp( +ev.diphoton_eta1[j] ) + ev.diphoton_pt2[j]*exp( +ev.diphoton_eta2[j] ) ) / sqrt_s,
                  xim = ( ev.diphoton_pt1[j]*exp( -ev.diphoton_eta1[j] ) + ev.diphoton_pt2[j]*exp( -ev.diphoton_eta2[j] ) ) / sqrt_s;

      //if ( ( xim < pots_accept["56N"] && xim < pots_accept["56F"] ) || xim > 0.15 ) continue;
      //if ( ( xip < pots_accept["45N"] && xip < pots_accept["45F"] ) || xip > 0.15 ) continue;

      //----- search for associated leptons

      /*const double lepton_mindist = 5.0; // in centimeters
      TVector3 diph_vtx( ev.diphoton_vertex_x[j], ev.diphoton_vertex_y[j], ev.diphoton_vertex_z[j] );
      unsigned short num_close_ele = 0, num_close_muon = 0;
      for ( unsigned short k = 0; k < ev.num_electron; ++k ) {
        TVector3 ele_vtx( electron_vtx_x[k], electron_vtx_y[k], electron_vtx_z[k] );
        const float ele_dist = ( ele_vtx-diph_vtx ).Mag();
        if ( ele_dist < lepton_mindist ) num_close_ele++;
      }
      for ( unsigned short k = 0; k < ev.num_muon; ++k ) {
        TVector3 mu_vtx( muon_vtx_x[k], muon_vtx_y[k], muon_vtx_z[k] );
        const float mu_dist = ( mu_vtx-diph_vtx ).Mag();
        if ( mu_dist < lepton_mindist ) num_close_muon++;
      }*/

      //----- search for associated jets

      TLorentzVector pho1, pho2;
      pho1.SetPtEtaPhiM( ev.diphoton_pt1[j], ev.diphoton_eta1[j], ev.diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( ev.diphoton_pt2[j], ev.diphoton_eta2[j], ev.diphoton_phi2[j], 0. );
      TLorentzVector cms = pho1+pho2;

      //const double min_extrajet_pt = 500.;
      /*unsigned short num_associated_jet = 0;
      for ( unsigned short k = 0; k < ev.num_jet; ++k ) {
        if ( jet_dipho_match[k] != j ) continue;
        //if ( jet_pt[k] > min_extrajet_pt ) num_associated_jet++;
        TLorentzVector jet;
        jet.SetPtEtaPhiE( jet_pt[k], jet_eta[k], jet_phi[k], jet_energy[k] );
        cms += jet;
      }*/

      //----- another batch of "exclusivity" cuts

      /*if ( num_close_ele > 0 ) continue;
      if ( num_close_muon > 0 ) continue;*/
      //if ( num_associated_jet > 0 ) continue;

      //----- reconstruct the energy loss from central system

      float diphoton_mass_error = ev.diphoton_mass[j]*0.02;
      float diphoton_rapidity_error = fabs( ev.diphoton_rapidity[j] )*0.061;

      for ( const auto cand : candidates ) {
        h_mass_all->Fill( cand.mass() );
        h_rap_all->Fill( cand.rapidity() );
        bool mass_match = is_matched( num_sigma, cms.M(), cand.mass(), diphoton_mass_error, cand.mass_error() );
        bool rap_match = is_matched( num_sigma, cms.Rapidity(), cand.rapidity(), diphoton_rapidity_error, cand.rapidity_error() );
        if ( mass_match && rap_match ) {
          cout << "@@@ DOUBLE TAGGING" << endl;
          cout << "masses: central system: " << cms.M() << ", diphoton: " << ev.diphoton_mass[j] << " +/- " << diphoton_mass_error << ", diphoton: " << cand.mass() << " +/- " << cand.mass_error() << endl;
          cout << "rapidities: central system: " << cms.Rapidity() << ", diphoton: " << ev.diphoton_rapidity[j] << " +/- " << diphoton_rapidity_error << ", diphoton: " << cand.rapidity() << " +/- " << cand.rapidity_error() << endl;
          gr_mass_massrapmatch.SetPoint( num_massrapmatch, cand.mass(), cms.M() );
          gr_mass_massrapmatch.SetPointError( num_massrapmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_massrapmatch.SetPoint( num_massrapmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_massrapmatch.SetPointError( num_massrapmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_massrapmatch++;
        }
        else if ( mass_match ) { // only match in mass
          gr_mass_massmatch.SetPoint( num_massmatch, cand.mass(), cms.M() );
          gr_mass_massmatch.SetPointError( num_massmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_massmatch.SetPoint( num_massmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_massmatch.SetPointError( num_massmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_massmatch++;
        }
        else if ( rap_match ) { // only match in rapidity
          gr_mass_rapmatch.SetPoint( num_rapmatch, cand.mass(), cms.M() );
          gr_mass_rapmatch.SetPointError( num_rapmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_rapmatch.SetPoint( num_rapmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_rapmatch.SetPointError( num_rapmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_rapmatch++;
        }
        else { // no matching at all
          gr_mass_nomatch.SetPoint( num_nomatch, cand.mass(), cms.M() );
          gr_mass_nomatch.SetPointError( num_nomatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_nomatch.SetPoint( num_nomatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_nomatch.SetPointError( num_nomatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_nomatch++;
        }
        cout << "matching: " << mass_match << "\t" << rap_match << endl;
      }

    }
  }
cout << "in plot:\n\t" << "not matching: " << num_nomatch << "\n\tmass match: " << num_massmatch << "\n\trap match: " << num_rapmatch << "\n\tboth match: " << num_massrapmatch << endl;

  //----- plotting part

  gr_mass_massrapmatch.SetTitle( "Diproton system missing mass (GeV)@@Diphoton mass (GeV)" );
  gr_rap_massrapmatch.SetTitle( "Diproton system rapidity@@Diphoton rapidity" );

  plot_matching( "2d_massmatch", gr_mass_nomatch, gr_mass_rapmatch, gr_mass_massmatch, gr_mass_massrapmatch, 300., 1800. );
  plot_matching( "2d_rapmatch", gr_rap_nomatch, gr_rap_rapmatch, gr_rap_massmatch, gr_rap_massrapmatch, -3., 3., 0.15 );

  {
    Canvas c( "1d_massmatch", "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
    h_mass_all->Sumw2();
    h_mass_all->Draw();
    c.Prettify( h_mass_all );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    Canvas c( "1d_rapmatch", "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
    h_rap_all->Sumw2();
    h_rap_all->Draw();
    c.Prettify( h_rap_all );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }

}

void plot_matching( const char* name, TGraphErrors& gr_nomatch, TGraphErrors& gr_rapmatch, TGraphErrors& gr_massmatch, TGraphErrors& gr_massrapmatch, double min, double max, double xleg )
{
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
  /*auto tmp = new TH2D( Form( "tmp_%s", name ), gr_massrapmatch.GetTitle(), 2, min, max, 2, min, max );
  tmp->Draw();*/
  auto diag = new TF1( "diag", "x", min, max );
  diag->SetLineColor( kGray );

  c.SetLegendX1( xleg );
  c.SetLegendY1( 0.7 );
  TMultiGraph mg;
  gr_nomatch.SetMarkerStyle( 24 );
  gr_rapmatch.SetMarkerStyle( 22 );
  gr_rapmatch.SetMarkerSize( 1.2 );
  gr_massmatch.SetMarkerStyle( 20 );
  gr_massmatch.SetMarkerSize( 1.0 );
  gr_massrapmatch.SetMarkerStyle( 21 );
  gr_massrapmatch.SetMarkerColor( kBlue+1 );
  gr_massrapmatch.SetMarkerSize( 1.3 );
  gr_massmatch.SetMarkerColor( kRed+1 );
  gr_rapmatch.SetMarkerColor( kGreen+2 );
  if ( gr_massrapmatch.GetN() > 0 ) c.AddLegendEntry( &gr_massrapmatch, "2D matching", "p" );
  c.AddLegendEntry( &gr_massmatch, "Mass matching", "lp" );
  c.AddLegendEntry( &gr_rapmatch, "Rapidity matching", "lp" );
  c.AddLegendEntry( &gr_nomatch, "No matching", "lp" );
  mg.Add( &gr_nomatch );
  mg.Add( &gr_massmatch );
  mg.Add( &gr_rapmatch );
  mg.Add( &gr_massrapmatch );

  diag->Draw( "l" );
  mg.Draw( "p" );

  //PaveText pt( 0.15, 0.8, 0.4, 0.95 );
  PaveText pt( xleg, 0.8, xleg+0.25, 0.95 );
  pt.SetTextAlign( kHAlignLeft+kVAlignBottom );
  pt.AddText( "1-|#Delta#phi(#gamma,#gamma)/#pi| < 0.005" );
  pt.Draw();

  diag->GetHistogram()->SetTitle( gr_massrapmatch.GetTitle() );
  c.Prettify( diag->GetHistogram() );
  //c.Prettify( tmp );
  //mg.GetXaxis()->SetLimits( min, max );
  diag->GetHistogram()->GetYaxis()->SetRangeUser( min, max );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}


/*void
plot_matching( const char* name, TGraphErrors& gr_unmatch, TGraphErrors& gr_match, TGraphErrors& gr_ooa, double limits )
{
  Canvas c( name, "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );

  TF1 lim( "lim", "x", 0., max_xi );
  lim.SetTitle( "#xi(RP)@@#xi(#gamma#gamma)" );
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
  lim.SetFillStyle( 3003 );
  lim.SetFillColor( kBlack );
  lim.SetLineColor( kBlack );
  lim.SetLineWidth( 1 );
  lim.Draw( "l" );
  c.Prettify( spec );
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
*/
  
vector<pair<float,float> > merge_nearfar( const vector<track_t>& near_tracks, const vector<track_t>& far_tracks, float xdiff_cut ) {
  vector<pair<float,float> > out;
  set<unsigned short> matched_far_ids;

  //--- first loop to extract near tracks with matching
  for ( const auto& near : near_tracks ) {
    float min_xidiff = 999.;
    int matched_far_id = -1;
    pair<track_t,track_t> assoc;
    unsigned short id_far = 0;
    for ( const auto& far : far_tracks ) {
      if ( fabs( near.x-far.x ) < xdiff_cut ) {
        //--- near-far track association candidate
        float xidiff = fabs( near.xi-far.xi );
        if ( xidiff < min_xidiff ) {
          assoc = make_pair( near, far );
          min_xidiff = xidiff;
          matched_far_id = id_far;
        }
      }
      id_far++;
    }
    if ( min_xidiff < 999. ) {
      out.emplace_back( assoc.second.xi, assoc.second.err_xi ); // store the far tracks info
      matched_far_ids.insert( matched_far_id );
    }
    //--- store the near track if no matching found with far track
    else out.emplace_back( near.xi, near.err_xi );
  }
  //--- second loop to add the remaining, unmatched far tracks
  unsigned short id_far = 0;
  for ( const auto& far : far_tracks ) {
    if ( matched_far_ids.count( id_far ) == 0 ) {
      //--- discard far track if a mapping is already found
      out.emplace_back( far.xi, far.err_xi );
    }
  }
  return out;
}
