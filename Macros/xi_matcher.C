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

#define MAX_PROTONS 10
#define MAX_DIPH 20

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

void xi_matcher()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_22jun.root" );
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
  float proton_track_x[MAX_PROTONS];
  float proton_track_xi[MAX_PROTONS], proton_track_xi_error[MAX_PROTONS];
  unsigned int proton_track_side[MAX_PROTONS], proton_track_pot[MAX_PROTONS];
  tr->SetBranchAddress( "num_proton_track", &num_proton_track );
  tr->SetBranchAddress( "proton_track_x", proton_track_x );
  tr->SetBranchAddress( "proton_track_xi", proton_track_xi );
  tr->SetBranchAddress( "proton_track_xi_error", proton_track_xi_error );
  tr->SetBranchAddress( "proton_track_side", proton_track_side );
  tr->SetBranchAddress( "proton_track_pot", proton_track_pot );

  unsigned int num_diphoton;
  float diphoton_pt1[MAX_DIPH], diphoton_pt2[MAX_DIPH];
  float diphoton_eta1[MAX_DIPH], diphoton_eta2[MAX_DIPH];
  float diphoton_r91[MAX_DIPH], diphoton_r92[MAX_DIPH];
  float diphoton_mass[MAX_DIPH], diphoton_rapidity[MAX_DIPH], diphoton_dphi[MAX_DIPH];
  tr->SetBranchAddress( "num_diphoton", &num_diphoton );
  tr->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  tr->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  tr->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  tr->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  tr->SetBranchAddress( "diphoton_r91", diphoton_r91 );
  tr->SetBranchAddress( "diphoton_r92", diphoton_r92 );
  tr->SetBranchAddress( "diphoton_mass", diphoton_mass );
  tr->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
  tr->SetBranchAddress( "diphoton_dphi", diphoton_dphi );

  TGraphErrors gr_45n_ooa, gr_45n_unmatch, gr_45n_match;
  TGraphErrors gr_45f_ooa, gr_45f_unmatch, gr_45f_match;
  TGraphErrors gr_56n_ooa, gr_56n_unmatch, gr_56n_match;
  TGraphErrors gr_56f_ooa, gr_56f_unmatch, gr_56f_match;

  TH1D* h_xispectrum_45n = new TH1D( "xisp_45n", "", 100, 0., max_xi ),
       *h_xispectrum_45f = dynamic_cast<TH1D*>( h_xispectrum_45n->Clone( "xisp_45f" ) ),
       *h_xispectrum_56n = dynamic_cast<TH1D*>( h_xispectrum_45n->Clone( "xisp_56n" ) ),
       *h_xispectrum_56f = dynamic_cast<TH1D*>( h_xispectrum_45n->Clone( "xisp_56f" ) );

  ofstream file_cands( "events_list.txt" );

  float min_xi_45n = 1., min_xi_45f = 1., min_xi_56n = 1., min_xi_56f = 1.;

  for ( unsigned long long i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );
    //cout << "event " << i << ": " << num_proton_track << " proton tracks, " << num_diphoton << " diphoton candidates" << endl;

    // first loop to identify the tracks and their respective pot

    auto align = pot_align::get_alignments( fill_number );

    vector< pair<float, float> > xi_45n, xi_45f, xi_56n, xi_56f;

    for ( unsigned short j=0; j<num_proton_track; j++ ) {
      double xi, xi_err;
      if ( proton_track_side[j]==0 && proton_track_pot[j]==2 ) {
        const auto& al = align[2];
        xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
        if ( proton_track_xi[j]<min_xi_45n ) { cout << "45N: " << fill_number << endl; min_xi_45n = proton_track_xi[j]; }
        xi_45n.emplace_back( xi, xi_err );
        h_xispectrum_45n->Fill( xi );
//cout << "xi45n: " << proton_track_xi[j] << " +/- " << proton_track_xi_error[j] << " / " << xi << " +/- " << xi_err << endl;
      }
      if ( proton_track_side[j]==0 && proton_track_pot[j]==3 ) {
        const auto& al = align[3];
        xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
        if ( proton_track_xi[j]<min_xi_45f ) { cout << "45F: " << fill_number << endl; min_xi_45f = proton_track_xi[j]; }
        xi_45f.emplace_back( xi, xi_err );
        h_xispectrum_45f->Fill( xi );
      }
      if ( proton_track_side[j]==1 && proton_track_pot[j]==2 ) {
        const auto& al = align[102];
        xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
        if ( proton_track_xi[j]<min_xi_56n ) { cout << "56N: " << fill_number << endl; min_xi_56n = proton_track_xi[j]; }
        xi_56n.emplace_back( xi, xi_err );
        h_xispectrum_56n->Fill( xi );
      }
      if ( proton_track_side[j]==1 && proton_track_pot[j]==3 ) {
        const auto& al = align[103];
        xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
        if ( proton_track_xi[j]<min_xi_56f ) { cout << "56F: " << fill_number << endl; min_xi_56f = proton_track_xi[j]; }
//cout << proton_track_xi[j]-xi << "\txi56f: " << proton_track_xi[j] << " +/- " << proton_track_xi_error[j] << " / " << xi << " +/- " << xi_err << endl;
        xi_56f.emplace_back( xi, xi_err );
        h_xispectrum_56f->Fill( xi );
      }
    }

    // second loop to identify the diphoton candidates and the corresponding xi values
    for ( unsigned short j=0; j<num_diphoton; j++ ) {

      if ( diphoton_pt1[j]<50. ) continue;
      if ( diphoton_pt2[j]<50. ) continue;
      if ( diphoton_eta1[j]>2.5 || ( diphoton_eta1[j]>1.4442 && diphoton_eta1[j]<1.566 ) ) continue;
      if ( diphoton_eta2[j]>2.5 || ( diphoton_eta2[j]>1.4442 && diphoton_eta2[j]<1.566 ) ) continue;
      if ( diphoton_r91[j]<0.94 ) continue;
      if ( diphoton_r92[j]<0.94 ) continue;

      if ( 1.-fabs( diphoton_dphi[j] )/M_PI>0.005 ) continue;

      const float xip_gg = ( diphoton_pt1[j]*exp( +diphoton_eta1[j] ) + diphoton_pt2[j]*exp( +diphoton_eta2[j] ) ) / sqrt_s;
      const float xim_gg = ( diphoton_pt1[j]*exp( -diphoton_eta1[j] ) + diphoton_pt2[j]*exp( -diphoton_eta2[j] ) ) / sqrt_s;

      vector< pair<float, float> > matching_45n, matching_45f, matching_56n, matching_56f;

      fill_matching( num_sigma, xi_45n, lim_45n, xip_gg, xip_gg*rel_err_xi_gg, gr_45n_match, gr_45n_unmatch, gr_45n_ooa, matching_45n );
      fill_matching( num_sigma, xi_45f, lim_45f, xip_gg, xip_gg*rel_err_xi_gg, gr_45f_match, gr_45f_unmatch, gr_45f_ooa, matching_45f );
      fill_matching( num_sigma, xi_56n, lim_56n, xim_gg, xim_gg*rel_err_xi_gg, gr_56n_match, gr_56n_unmatch, gr_56n_ooa, matching_56n );
      fill_matching( num_sigma, xi_56f, lim_56f, xim_gg, xim_gg*rel_err_xi_gg, gr_56f_match, gr_56f_unmatch, gr_56f_ooa, matching_56f );

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
      else if ( ( matching_45n.size()+matching_45f.size()>0 ) || ( matching_56n.size()+matching_56f.size()>0 ) ) {
        file_cands << "1\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
      }

    }
  }
  file_cands.close();

  cout << "minimum xi expected: 45N: " << lim_45n << ", 45F: " << lim_45f << ", 56N: " << lim_56n << ", 56F: " << lim_56f << endl;
  cout << "minimum xi observed: 45N: " << min_xi_45n << ", 45F: " << min_xi_45f << ", 56N: " << min_xi_56n << ", 56F: " << min_xi_56f << endl;

  //----- plotting part

  plot_matching( "xi_matching_45n", gr_45n_unmatch, gr_45n_match, gr_45n_ooa, lim_45n );
  plot_matching( "xi_matching_45f", gr_45f_unmatch, gr_45f_match, gr_45f_ooa, lim_45f );
  plot_matching( "xi_matching_56n", gr_56n_unmatch, gr_56n_match, gr_56n_ooa, lim_56n );
  plot_matching( "xi_matching_56f", gr_56f_unmatch, gr_56f_match, gr_56f_ooa, lim_56f );

  plot_xispectrum( "xi_spectrum_45n", h_xispectrum_45n, lim_45n );
  plot_xispectrum( "xi_spectrum_45f", h_xispectrum_45f, lim_45f );
  plot_xispectrum( "xi_spectrum_56n", h_xispectrum_56n, lim_56n );
  plot_xispectrum( "xi_spectrum_56f", h_xispectrum_56f, lim_56f );

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
  box.SetLineColor( colour );
  box.SetLineWidth( 0 );
  box.Draw();

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

  c.AddLegendEntry( &gr_match, "Matching events", "lp" );
  c.AddLegendEntry( &gr_unmatch, "Non-matching events", "lp" );
  c.AddLegendEntry( &gr_ooa, "Out of acceptance events", "lp" );
  c.AddLegendEntry( &box, "No acceptance for RP", "f" );

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
  Canvas c( name );
  spec->Draw();
  TLine lim( limit, 0., limit, spec->GetMaximum() );
  lim.SetLineColor( kRed );
  lim.SetLineWidth( 2 );
  lim.Draw();
  c.Prettify( spec );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
