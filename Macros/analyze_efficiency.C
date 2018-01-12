#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"
#include "TStyle.h"
#include <map>
#include <iostream>

#include "Canvas.h"
#include "xi_reconstruction.h"
#include "pot_alignment.h"

#define loc_www "/afs/cern.ch/user/l/lforthom/www/private/ctpps/efficiency"

void analyze_efficiency()
{
  const unsigned int ref_fill = 4985;

  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );
  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );

  auto f = TFile::Open( "/eos/cms/store/user/lforthom/ctpps/efficiency_study/merged_eff_outputBCG.root" );
  auto tree = dynamic_cast<TTree*>( f->Get( "ntp" ) );
  tree->SetBranchStatus( "*", 0 );
  unsigned int fill_number;
  tree->SetBranchStatus( "fill_number", 1 ); tree->SetBranchAddress( "fill_number", &fill_number );
  const unsigned short max_tracks = 20;
  unsigned int num_proton_track, proton_track_side[max_tracks], proton_track_pot[max_tracks];
  float proton_track_x[max_tracks], proton_track_y[max_tracks];
  tree->SetBranchStatus( "num_proton_track", 1 ); tree->SetBranchAddress( "num_proton_track", &num_proton_track );
  tree->SetBranchStatus( "proton_track_x", 1 ); tree->SetBranchAddress( "proton_track_x", proton_track_x );
  tree->SetBranchStatus( "proton_track_y", 1 ); tree->SetBranchAddress( "proton_track_y", proton_track_y );
  tree->SetBranchStatus( "proton_track_side", 1 ); tree->SetBranchAddress( "proton_track_side", proton_track_side );
  tree->SetBranchStatus( "proton_track_pot", 1 ); tree->SetBranchAddress( "proton_track_pot", proton_track_pot );

  map<unsigned short,const char*> pot_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };
  map<unsigned short,pair<double,double> > pot_fit_limits = { { 2, { 9., 15. } }, { 3, { 7., 15. } }, { 102, { 8., 13. }  }, { 103, { 6.5, 13. } } };
  map<unsigned short,TH1D*> h_num_x, h_denom_x, h_num_y, h_denom_y, h_num_xi, h_denom_xi;
  //map<unsigned short,TEfficiency*> eff_x, eff_y, eff_xi;
  double y_bins[] = { -10., -9., -8., -7., -6.,
                      -5., -4.5, -4., -3.5, -3., -2.5, -2.25,
                      -2., -1.75, -1.5, -1.375, -1.25, -1.125, -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125,
                       0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875,
                       1., 1.125, 1.25, 1.375, 1.5, 1.75,
                       2., 2.25, 2.5, 3., 3.5, 4., 4.5,
                       5., 6., 7., 8., 9., 10. };
  for ( const auto& p : pot_names ) {
    h_num_x[p.first] = new TH1D( Form( "h_num_x_%d", p.first ), Form( "Track x (%s)@@Entries@@mm?.2f", p.second ), 40, 0., 20. );
    h_denom_x[p.first] = dynamic_cast<TH1D*>( h_num_x[p.first]->Clone( Form( "h_denom_x_%d", p.first ) ) );
    //h_num_y[p.first] = new TH1D( Form( "h_num_y_%d", p.first ), Form( "Track y (%s)@@Entries@@mm?.2f", p.second ), 100, -12.5, 12.5 );
    h_num_y[p.first] = new TH1D( Form( "h_num_y_%d", p.first ), Form( "Track y (%s)@@Entries@@mm?.2f", p.second ), sizeof( y_bins )/sizeof( double )-1, y_bins );
    h_denom_y[p.first] = dynamic_cast<TH1D*>( h_num_y[p.first]->Clone( Form( "h_denom_y_%d", p.first ) ) );
    h_num_xi[p.first] = new TH1D( Form( "h_num_xi_%d", p.first ), Form( "Track #xi (%s)@@Entries@@?.3f", p.second ), 40, 0., 0.2 );
    h_denom_xi[p.first] = dynamic_cast<TH1D*>( h_num_xi[p.first]->Clone( Form( "h_denom_xi_%d", p.first ) ) );
    /*eff_x[p.first] = new TEfficiency( Form( "eff_x_%d", p.first ), Form( ";Track x (%s);Efficiency", p.second ), 40, 0., 20. );
    eff_y[p.first] = new TEfficiency( Form( "eff_x_%d", p.first ), Form( ";Track y (%s);Efficiency", p.second ), sizeof( y_bins )/sizeof( double )-1, y_bins );
    eff_xi[p.first] = new TEfficiency( Form( "eff_x_%d", p.first ), Form( ";Track #xi (%s);Efficiency", p.second ), 50, 0., 0.2 );*/
  }

  const unsigned long long num_entries = tree->GetEntriesFast();
  for ( unsigned long long i = 0; i < num_entries; ++i ) {
    tree->GetEntry( i );
    if ( i % 1000000 == 0 ) cout << "-- event " << i << "/" << num_entries << endl;
    const bool is_ref_fill = ( fill_number == ref_fill );

    auto align = pot_align::get_alignments( fill_number );

    for ( unsigned int j = 0; j < num_proton_track; ++j ) {
      const unsigned short pid = 100 * proton_track_side[j] + proton_track_pot[j];
      double xi = 0., xi_err = 0.;
      xi_reco::reconstruct( proton_track_x[j]+align[pid].x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
      const double trk_x = ( proton_track_x[j]+align[pid].x )*1.e3;
      const double trk_y = ( proton_track_y[j]-align[pid].y )*1.e3;
      if ( !is_ref_fill ) {
        h_num_x[pid]->Fill( trk_x );
        h_num_y[pid]->Fill( trk_y );
        h_num_xi[pid]->Fill( xi );
      }
      else {
        h_denom_x[pid]->Fill( trk_x );
        h_denom_y[pid]->Fill( trk_y );
        h_denom_xi[pid]->Fill( xi );
      }
      /*eff_x[pid]->Fill( !is_ref_fill, trk_x );
      eff_y[pid]->Fill( !is_ref_fill, trk_y );
      eff_xi[pid]->Fill( !is_ref_fill, xi );*/
    }
  }

  // rescaling for a given x range

  for ( const auto& p : pot_names ) {
    // first start by scaling the two histograms to 1
    h_num_x[p.first]->Sumw2();
    h_denom_x[p.first]->Sumw2();
    h_num_y[p.first]->Sumw2();
    h_denom_y[p.first]->Sumw2();
    h_num_xi[p.first]->Sumw2();
    h_denom_xi[p.first]->Sumw2();
    h_num_x[p.first]->Scale( 1./h_num_x[p.first]->Integral( "width" ) );
    h_denom_x[p.first]->Scale( 1./h_denom_x[p.first]->Integral( "width" ) );
    const auto limits = pot_fit_limits[p.first];
    const double weight_num = h_num_x[p.first]->Integral( h_num_x[p.first]->GetXaxis()->FindBin( limits.first ), h_num_x[p.first]->GetXaxis()->FindBin( limits.second ) );
    const double weight_denom = h_denom_x[p.first]->Integral( h_denom_x[p.first]->GetXaxis()->FindBin( limits.first ), h_denom_x[p.first]->GetXaxis()->FindBin( limits.second ) );
    const double norm = weight_denom/weight_num;
    cout << p.second << "::" << weight_num << "|" << weight_denom << "|ratio=" << norm << endl;
    h_num_x[p.first]->Scale( norm );
    h_num_y[p.first]->Scale( norm/h_num_y[p.first]->Integral( "width" ) );
    h_denom_y[p.first]->Scale( 1./h_denom_y[p.first]->Integral( "width" ) );
    h_num_xi[p.first]->Scale( norm/h_num_xi[p.first]->Integral( "width" ) );
    h_denom_xi[p.first]->Scale( 1./h_denom_xi[p.first]->Integral( "width" ) );
  }

  // plotting part

  const string top_title = "CMS-TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}";

  auto text = new TText();
  text->SetTextColor( kGray+3 );

  const vector<string> distrib = { "x", "y", "xi" };
  vector<pair<map<unsigned short,TH1D*>,map<unsigned short,TH1D*> > > hists = { { h_num_x, h_denom_x }, { h_num_y, h_denom_y }, { h_num_xi, h_denom_xi } };
  unsigned short i = 0;
  for ( auto& hist : hists ) {
    for ( const auto& p : pot_names ) {
      { // plot both the numerator and denominator
        Canvas c( Form( "dist_%s_%s", distrib[i].c_str(), p.second ), top_title.c_str() );
        c.SetLegendX1( 0.475 );
        THStack hs;
        hs.Add( hist.first[p.first] );
        c.AddLegendEntry( hist.first[p.first], "All runs in 2016B/C/G", "f" );
        hist.first[p.first]->SetLineColor( kBlack );
        hist.first[p.first]->SetFillStyle( 3004 );
        hist.first[p.first]->SetFillColor( kBlack );
        hs.Add( hist.second[p.first] );
        c.AddLegendEntry( hist.second[p.first], Form( "Reference fill (%d)", ref_fill ), "l" );
        hist.second[p.first]->SetLineColor( kRed+1 );
        hist.second[p.first]->SetLineWidth( 2 );
        hs.Draw( "hist,nostack" );
        hs.GetHistogram()->SetTitle( hist.first[p.first]->GetTitle() );
        c.Prettify( hs.GetHistogram() );
        if ( distrib[i] == "x" ) {
          auto limits = pot_fit_limits[p.first];
          auto range = new TArrow( limits.first, 0.015, limits.second, 0.015, 0.015, "<>" );
          //range->SetNDC( true );
          range->SetLineColor( kGray+3 );
          range->SetLineWidth( 3 );
          range->Draw();
          auto range_leg = text->DrawText( 0.5*( limits.first+limits.second ), 0.0175, "Norm. range" );
          range_leg->SetTextFont( 42 );
          range_leg->SetTextAlign( 21 );
          range_leg->SetTextSize( 0.035 );
          range_leg->Draw();
        }
        c.Save( "pdf,png", loc_www );
      }
      { // plot the efficiencies
        Canvas c( Form( "ratio_%s_%s", distrib[i].c_str(), p.second ), top_title.c_str() );
        gStyle->SetOptStat( 0 );
        //Option_t* estim_mode = "cp v"; // Clopper-Pearson estimator
        //auto gr = new TGraphAsymmErrors( hist.first[p.first], hist.second[p.first], estim_mode );
        /*auto gr = new TGraphAsymmErrors( hist.first[p.first], hist.second[p.first] );
        gr->Draw( "alp" );
        gr->GetHistogram()->SetTitle( hist.first[p.first]->GetTitle() );
        c.Prettify( gr->GetHistogram() );
        gr->GetYaxis()->SetTitle( "Efficiency" );
        gr->GetYaxis()->SetRangeUser( 0., 1.1 );
        gr->SetLineColor( kBlack );
        gr->SetMarkerStyle( 24 );
        gr->SetMarkerColor( kRed+1 );*/
        auto ratio = dynamic_cast<TH1D*>( hist.first[p.first]->Clone() );
        ratio->SetTitle( TString( ratio->GetTitle() ).ReplaceAll( "Entries", "Efficiency" ) );
        auto den = dynamic_cast<TH1D*>( hist.second[p.first]->Clone() );
        /*ratio->Sumw2();
        den->Sumw2();*/
        ratio->Divide( den );
        //ratio->Draw( "p" );
        ratio->Draw( "e0" );
        //ratio->Draw( "hist" );
        //ratio->SetLineStyle( 2 );
        c.Prettify( ratio );
        ratio->SetMarkerStyle( 24 );
        ratio->SetLineColor( kBlack );
        ratio->GetYaxis()->SetRangeUser( 0., 1.1 );
        //ratio->GetYaxis()->SetTitle( "Efficiency" );
        c.SetGrid( 1, 1 );
        c.Save( "pdf,png", loc_www );
      }
      /*{ // plot the efficiencies
        Canvas c( Form( "eff_%s_%s", distrib[i].c_str(), p.second ), top_title.c_str() );
        auto eff = new TEfficiency( *hist.first[p.first], *hist.second[p.first] );
        eff->Draw( "alp" );
        eff->SetMarkerStyle( 24 );
        c.Save( "pdf,png", loc_www );
      }*/
    }
    i++;
  }
  
  /*vector<map<unsigned short,TEfficiency*> > effs = { { eff_x, eff_y, eff_xi } };
  i = 0;
  for ( auto& eff : effs ) {
    for ( const auto& p : pot_names ) {
      { // plot the efficiencies
        Canvas c( Form( "eff_%s_%s", distrib[i].c_str(), p.second ), top_title.c_str() );
        eff[p.first]->Draw();
        eff[p.first]->SetMarkerStyle( 24 );
        c.Save( "pdf,png", loc_www );
      }
    }
    i++;
  }*/
}

