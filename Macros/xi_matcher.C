#include "Canvas.h"

#define MAX_PROTONS 10
#define MAX_DIPH 20

bool is_matched( int n_sigma, float xi_rp, float xi_cs, float err_xi_rp, float err_xi_cs )
{
  const double combined_error = sqrt( err_xi_rp*err_xi_rp + err_xi_cs*err_xi_cs );
  const double delta = fabs( xi_cs-xi_rp );

 return ( delta/combined_error<=n_sigma );
}
void fill_matching( float num_sigma, const std::vector< std::pair<float, float> >& tracks, float limit, float xi_gg, float err_xi_gg, TGraphErrors& gr_match, TGraphErrors& gr_unmatch, TGraphErrors& gr_ooa );
void plot_matching( const char* name, TGraphErrors& gr_unmatch, TGraphErrors& gr_match, TGraphErrors& gr_ooa, double limits );

void xi_matcher()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_9may.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  const float sqrt_s = 13.0e3;
  const float rel_err_xi_gg = 0.028;
  const float num_sigma = 2.0;

  const float lim_45n = 0.033;
  const float lim_45f = 0.024;
  const float lim_56n = 0.0496;
  const float lim_56f = 0.037;

  unsigned int num_proton_track;
  float proton_track_xi[MAX_PROTONS], proton_track_xi_error[MAX_PROTONS];
  unsigned int proton_track_side[MAX_PROTONS], proton_track_pot[MAX_PROTONS];
  tr->SetBranchAddress( "num_proton_track", &num_proton_track );
  tr->SetBranchAddress( "proton_track_xi", proton_track_xi );
  tr->SetBranchAddress( "proton_track_xi_error", proton_track_xi_error );
  tr->SetBranchAddress( "proton_track_side", proton_track_side );
  tr->SetBranchAddress( "proton_track_pot", proton_track_pot );

  unsigned int num_diphoton;
  float diphoton_pt1[MAX_DIPH], diphoton_pt2[MAX_DIPH];
  float diphoton_eta1[MAX_DIPH], diphoton_eta2[MAX_DIPH];
  float diphoton_r91[MAX_DIPH], diphoton_r92[MAX_DIPH];
  float diphoton_dphi[MAX_DIPH];
  tr->SetBranchAddress( "num_diphoton", &num_diphoton );
  tr->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  tr->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  tr->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  tr->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  tr->SetBranchAddress( "diphoton_r91", diphoton_r91 );
  tr->SetBranchAddress( "diphoton_r92", diphoton_r92 );
  tr->SetBranchAddress( "diphoton_dphi", diphoton_dphi );

  TGraphErrors gr_45n_ooa, gr_45n_unmatch, gr_45n_match;
  TGraphErrors gr_45f_ooa, gr_45f_unmatch, gr_45f_match;
  TGraphErrors gr_56n_ooa, gr_56n_unmatch, gr_56n_match;
  TGraphErrors gr_56f_ooa, gr_56f_unmatch, gr_56f_match;

  for ( unsigned long long i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );
    //cout << "event " << i << ": " << num_proton_track << " proton tracks, " << num_diphoton << " diphoton candidates" << endl;

    // first loop to identify the tracks and their respective pot

    vector< pair<float, float> > xi_45n, xi_45f, xi_56n, xi_56f;

    for ( unsigned short j=0; j<num_proton_track; j++ ) {
      if ( proton_track_side[j]==0 && proton_track_pot[j]==2 ) { xi_45n.emplace_back( proton_track_xi[j], proton_track_xi_error[j] ); }
      if ( proton_track_side[j]==0 && proton_track_pot[j]==3 ) { xi_45f.emplace_back( proton_track_xi[j], proton_track_xi_error[j] ); }
      if ( proton_track_side[j]==1 && proton_track_pot[j]==2 ) { xi_56n.emplace_back( proton_track_xi[j], proton_track_xi_error[j] ); }
      if ( proton_track_side[j]==1 && proton_track_pot[j]==3 ) { xi_56f.emplace_back( proton_track_xi[j], proton_track_xi_error[j] ); }
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

      fill_matching( num_sigma, xi_45n, lim_45n, xip_gg, xip_gg*rel_err_xi_gg, gr_45n_match, gr_45n_unmatch, gr_45n_ooa );
      fill_matching( num_sigma, xi_45f, lim_45f, xip_gg, xip_gg*rel_err_xi_gg, gr_45f_match, gr_45f_unmatch, gr_45f_ooa );
      fill_matching( num_sigma, xi_56n, lim_56n, xim_gg, xim_gg*rel_err_xi_gg, gr_56n_match, gr_56n_unmatch, gr_56n_ooa );
      fill_matching( num_sigma, xi_56f, lim_56f, xim_gg, xim_gg*rel_err_xi_gg, gr_56f_match, gr_56f_unmatch, gr_56f_ooa );

      

      cout << "--> " << gr_45n_match.GetN() << " / " << gr_45f_match.GetN() << " / " << gr_56n_match.GetN() << " / " << gr_56f_match.GetN() << endl;

    }
  }

  //----- plotting part

  plot_matching( "xi_matching_45n", gr_45n_unmatch, gr_45n_match, gr_45n_ooa, lim_45n );
  plot_matching( "xi_matching_45f", gr_45f_unmatch, gr_45f_match, gr_45f_ooa, lim_45f );
  plot_matching( "xi_matching_56n", gr_56n_unmatch, gr_56n_match, gr_56n_ooa, lim_56n );
  plot_matching( "xi_matching_56f", gr_56f_unmatch, gr_56f_match, gr_56f_ooa, lim_56f );

}

void
fill_matching( float num_sigma, const std::vector< std::pair<float, float> >& tracks, float limit, float xi_gg, float err_xi_gg, TGraphErrors& gr_match, TGraphErrors& gr_unmatch, TGraphErrors& gr_ooa )
{
  for ( const auto& pair : tracks ) {
    if ( xi_gg>limit ) {
      if ( is_matched( num_sigma, pair.first, xi_gg, pair.second, err_xi_gg ) ) {
        unsigned int n = gr_match.GetN();
        gr_match.SetPoint( n, pair.first, xi_gg );
        gr_match.SetPointError( n, pair.second, err_xi_gg );
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

  TF1 lim( "lim", "x", 0., 0.2 );
  lim.SetTitle( "#xi(RP)\\#xi(#gamma#gamma)" );
  lim.Draw("lr+");
  lim.SetLineColor(4);
  lim.SetLineWidth(2);

  int colour = kAzure+10;
  TBox box( 0., 0., limits, 0.2 );
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

  lim.GetYaxis()->SetRangeUser( 0., 0.2 );

  TH2D *axiiis = new TH2D( Form( "axiiis_%s", name ), "", 10, 0, 0.2, 10, 0, 0.2 );
  axiiis->Draw( "sameaxis" );
  c.GetLegend()->SetFillStyle( 0 );
  c.GetLegend()->SetLineWidth( 0 );
  c.Prettify( lim.GetHistogram() );

  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

