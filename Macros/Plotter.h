#ifndef Plotter_h
#define Plotter_h

#include "Canvas.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "THStack.h"
#include "TStyle.h"

#include <map>
#include <vector>

class Plotter
{
  public:
    typedef std::vector< std::pair<std::string, TGraphErrors*> > GraphsMap;
    typedef std::vector< std::pair<std::string, TH1*> > HistsMap;

  public:
    Plotter( const char* out_path, const char* top_label ) : out_path_( out_path ), top_label_( top_label ) {}
    ~Plotter() {}

    void plot_multihists( const char* name, HistsMap hm, float min_ratio_y=0., float max_ratio_y=0.55, bool draw_overflow=true ) const {
      if ( hm.size()==0 ) return;

      Canvas c( name, top_label_, hm.size()>1 );
      TH1* hist = 0;

      unsigned short i = 0;
      double min = 0.001,
             max = -1.;
      HistsMap hm2;
      for ( HistsMap::iterator it=hm.begin(); it!=hm.end(); it++ ) {
        hist = dynamic_cast<TH1*>( ( draw_overflow ) ? WithOverflow( it->second ) : it->second );
        hm2.push_back( std::make_pair( it->first, hist ) );
        if ( i>0 ) {
          if ( !draw_overflow ) hist->Sumw2();
          hist->SetMarkerStyle( marker_pool_[i] );
          hist->SetMarkerColor( colour_pool_[i] );
        }
        //hist->Draw( ( i==0 ) ? "" : "e1 same" );
        hist->Draw( ( i==0 ) ? "hist" : "same" );
        hist->SetLineColor( kBlack );
	//hist->SetFillStyle( 3004 );
	//hist->SetFillColor( kBlack );
        //h->SetMarkerStyle( 20 );
        min = TMath::Min( min, hist->GetMinimum() );
        max = TMath::Max( max, hist->GetMaximum() );
        if ( i>2 ) hist->SetFillStyle( 3005 );
        if ( hm.size()>1 ) {
          c.AddLegendEntry( hist, it->first.c_str(), ( i==0 ) ? "f" : "ep" );
        }
        i++;
      }
      hist = dynamic_cast<TH1*>( hm.begin()->second );
      hist->SetMinimum( min );
      hist->SetMaximum( max*1.1 );
      c.Prettify( hist );
      /*if ( hm.size()==2 ) c.RatioPlot( hm[0].second, hm[1].second, 0, min_ratio_y, max_ratio_y );
      if ( hm.size()>2 ) c.RatioPlot( hm[0].second, hm[1].second, hm[2].second, min_ratio_y, max_ratio_y );*/
      if ( hm.size()>1 ) c.RatioPlot( hm2, min_ratio_y, max_ratio_y );
      c.Save( "pdf,png", out_path_ );
    }

    void plot_3hists( const char* name, TH1D* h, TH1D* h_1tag=0, TH1D* h_2tag=0, bool draw_overflow=true ) const {
      HistsMap hm;
      hm.push_back( std::make_pair( "All diphotons", h ) );
      hm.push_back( std::make_pair( "track(s) in 1 arm", h_1tag ) );
      hm.push_back( std::make_pair( "track(s) in 2 arms", h_2tag ) );
      plot_multihists( name, hm, 0., 0.64, draw_overflow );
    }

    void plot_balances( const char* name, const char* title, GraphsMap gr_map,
                        const float& min_xy=-999., const float& max_xy=-999., const float& rp_acc=-1. ) const {
      Canvas c( name, top_label_ );
      //c.DrawFrame( min_xy, min_xy, max_xy, max_xy );

      TMultiGraph mg;
      TGraphErrors* gr = 0;

      unsigned short i = 0;
      for ( GraphsMap::iterator it=gr_map.begin(); it!=gr_map.end(); it++ ) {
        gr = ( TGraphErrors* )it->second;
        mg.Add( gr, "ep" );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        //c.SetLegendY1( 0.18 );
        if ( gr_map.size()>0 && strcmp( it->first.c_str(), "" )!=0 ) c.AddLegendEntry( gr, it->first.c_str(), "p" );
        i++;
      }


      mg.Draw( "ap" );
      mg.GetHistogram()->SetTitle( title );
      c.Prettify( mg.GetHistogram() );
      draw_diagonal( min_xy, max_xy, 0. );
      if ( min_xy!=-999. and max_xy!=-999. ) {
        mg.GetXaxis()->SetLimits( min_xy, max_xy );
        mg.GetYaxis()->SetRangeUser( min_xy, max_xy );
      }

      if ( rp_acc>0. ) {
        TPave* acc_h = new TPave( min_xy, min_xy, max_xy, rp_acc, 1 );
        acc_h->SetLineWidth( 1 );
        acc_h->SetLineColor( kBlack );
        acc_h->SetFillColorAlpha( kGray, 0.4 );
        //acc_h->SetFillColor( kGray );
        acc_h->Draw( "same" );
        TPave* acc_v = new TPave( min_xy, min_xy, rp_acc, max_xy, 1 );
        acc_v->SetLineWidth( 1 );
        acc_v->SetLineColor( kBlack );
        acc_v->SetFillColorAlpha( kGray, 0.4 );
        //acc_v->SetFillColor( kGray );
        acc_v->Draw( "same" );
      }
      c.Save( "pdf,png", out_path_ );
    }

    void plot_balances( const char* name, const char* title, TGraphErrors* gr, const float& min_xy=-999., const float& max_xy=-999., const float& rp_acc=-1. ) const {
      GraphsMap gm; gm.push_back( std::pair<const char*, TGraphErrors*>( gr->GetTitle(), gr ) );
      return plot_balances( name, title, gm, min_xy, max_xy, rp_acc );
    }

    void plot_xi_correlations( const char* sector, GraphsMap graphs_map ) const {
      Canvas c( Form( "xi_nearfar_corr_%s", sector ), top_label_ );
      TMultiGraph mg;

      TGraphErrors* gr = 0;

      unsigned short i = 0;
      for ( GraphsMap::iterator it=graphs_map.begin(); it!=graphs_map.end(); it++ ) {
        gr = ( TGraphErrors* )it->second;
        mg.Add( gr );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        if ( strcmp( it->first.c_str(), "" )!=0 ) c.AddLegendEntry( gr, it->first.c_str(), "p" );
        i++;
      }

      mg.Draw( "ap" );
      mg.GetHistogram()->SetTitle( Form( "#xi_{%s} (near pot)\\#xi_{%s} (far pot)", sector, sector ) );

      c.Prettify( mg.GetHistogram() );
      draw_diagonal( 0., 0.2 );

      mg.GetXaxis()->SetLimits( 0., 0.2 );
      mg.GetYaxis()->SetRangeUser( 0., 0.2 );

      c.Save( "pdf,png", out_path_ );
    }

    void draw_multigraph( const char* filename, GraphsMap g_map, float lim_min=-999., float lim_max=-999., bool draw_diag=false, float y_limit=-1. ) const {
      Canvas c( filename, top_label_ );
      //if ( lim_min!=-999. && lim_max!=-999. )
      //  c.DrawFrame( lim_min, lim_min, lim_max, lim_max );

      TGraphErrors* gr = 0;

      TMultiGraph mg;
      unsigned short i = 0;
      for ( GraphsMap::iterator it=g_map.begin(); it!=g_map.end(); it++ ) {
        gr = ( TGraphErrors* )it->second;
        mg.Add( gr );
        if ( strcmp( gr->GetTitle(), "" )!=0 ) mg.SetTitle( gr->GetTitle() );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        if ( strcmp( it->first.c_str(), "" )!=0 ) c.AddLegendEntry( gr, it->first.c_str(), "p" );
        i++;
      }

      mg.Draw( "ap" );
      c.Prettify( mg.GetHistogram() );
      if ( draw_diag ) draw_diagonal( lim_min, lim_max );

      if ( lim_min>0 || lim_max>0. ) {
        mg.GetXaxis()->SetLimits( ( lim_min>0 ? lim_min : 0. ), ( lim_max>0 ? lim_max : 500. ) );
        mg.GetYaxis()->SetRangeUser( ( lim_min>0 ? lim_min : 0. ), ( lim_max>0 ? lim_max : 500. ) );
      }

      if ( y_limit>0. && y_limit<lim_max ) {
        TPave* acc_h = new TPave( lim_min, lim_min, lim_max, y_limit, 1 );
        acc_h->SetLineWidth( 1 );
        acc_h->SetLineColor( kBlack );
        acc_h->SetFillColorAlpha( kGray, 0.4 );
        acc_h->Draw( "same" );
        TPave* acc_v = new TPave( lim_min, lim_min, y_limit, lim_max, 1 );
        acc_v->SetLineWidth( 1 );
        acc_v->SetLineColor( kBlack );
        acc_v->SetFillColorAlpha( kGray, 0.2 );
        acc_v->Draw( "same" );
      }
      c.Save( "pdf,png", out_path_ );
    }

    void draw_multiplot( const char* filename, HistsMap h_map_data, HistsMap h_map_mc, HistsMap h_map_sig, TString label="", bool colours=true, bool logy=false ) const {
      Canvas c( filename, top_label_, true );
      TH1D* hist = 0;
      TH1D* h_data, *h_mc = 0;
      double max_bin = -1.;
      unsigned short i = 0;
      THStack hs_mc, hs_data, hs_sig;
      const float leg_size_y = TMath::Max( 0.15, 0.03*( h_map_data.size()+h_map_mc.size()+h_map_sig.size() ) );
      c.SetLegendSizeY( leg_size_y );
      c.SetLegendY1( 0.75-leg_size_y+0.15 );
      for ( HistsMap::iterator it=h_map_data.begin(); it!=h_map_data.end(); it++ ) {
        hist = ( TH1D* )it->second;
        if ( i==0 ) h_data = dynamic_cast<TH1D*>( hist );
        else h_data->Add( hist );
        // draw the data distributions unstacked
        hist->Sumw2();
        hist->SetMarkerStyle( 20+i );
        hist->SetMarkerColor( kBlack );
        hist->SetLineColor( kBlack );
        if ( strcmp( it->first.c_str(), "" )!=0 ) c.AddLegendEntry( hist, it->first.c_str(), "lp" );
        max_bin = TMath::Max( max_bin, hist->GetMaximum() );
        hs_data.Add( hist );
        i++;
      }
      i = 0;
      for ( HistsMap::iterator it=h_map_mc.begin(); it!=h_map_mc.end(); it++ ) {
        hist = dynamic_cast<TH1D*>( it->second );
        if ( i==0 ) h_mc = dynamic_cast<TH1D*>( hist->Clone() );
        else { h_mc->Add( hist ); }
        if ( colours ) hist->SetFillColorAlpha( colour_pool_[i+1], 0.5 );
        //hist->SetFillStyle( 3002 );
        hist->SetLineColor( kBlack );
        if ( strcmp( it->first.c_str(), "" )!=0 ) { c.AddLegendEntry( hist, it->first.c_str(), "f" ); }
        hs_mc.Add( hist );
        if ( i==0 ) hs_mc.SetTitle( hist->GetTitle() );
        i++;
      }
      i = 0;
      for ( HistsMap::iterator it=h_map_sig.begin(); it!=h_map_sig.end(); it++ ) {
        hist = dynamic_cast<TH1D*>( it->second );
        TH1D* hist_stacked = dynamic_cast<TH1D*>( hist->Clone() );
        hist->SetLineColor( kGreen );
        hist->SetLineWidth( 3 );
        hist->SetLineStyle( i );
        if ( strcmp( it->first.c_str(), "" )!=0 ) { c.AddLegendEntry( hist, it->first.c_str(), "l" ); }
        hs_sig.Add( hist );
        if ( i==0 ) hs_sig.SetTitle( hist->GetTitle() );
        hist_stacked->SetLineColor( kBlack );
        hist_stacked->SetLineWidth( 1 );
        hist_stacked->SetFillColorAlpha( kGreen-9, 0.5 );
        hs_mc.Add( hist_stacked );
        //h_mc->Add( hist_stacked );
        i++;
      }
      hs_mc.Draw( "hist" );
      hs_data.Draw( "same nostack" );
      hs_sig.Draw( "hist same nostack" );
      if ( h_mc ) {
        h_mc->Draw( "e2 same" );
        h_mc->SetFillColor( kBlack );
        h_mc->SetFillStyle( 3004 );
      }
      hist = ( TH1D* )hs_mc.GetHistogram();
      max_bin = TMath::Max( hist->GetMaximum(), max_bin );
      if ( logy ) {
        hs_mc.SetMaximum( max_bin*5. );
        if ( hs_mc.GetMinimum()==0. ) hs_mc.SetMinimum( 0.1 );
      }
      else { hs_mc.SetMaximum( max_bin*1.4 ); }
      c.Prettify( hist );
      HistsMap hm;
      hm.emplace_back( "mc", h_mc );
      hm.emplace_back( "data", h_data );
      hs_mc.SetTitle( "" );
      c.RatioPlot( hm, -0.4, 2.4, 1.0 );
      c.cd( 1 );
      if ( !label.IsNull() ) {
        PaveText* lab = new PaveText( 0.135, 0.96, 0.2, 0.97 );
        lab->SetTextSize( 0.05 );
        lab->SetTextAlign( kVAlignTop+kHAlignLeft );
        lab->AddText( label );
        lab->Draw( "same" );
      }
      if ( logy ) dynamic_cast<TPad*>( c.GetPad( 1 ) )->SetLogy();
      c.Save( "pdf,png", out_path_ );
    }

    void draw_multiplot( const char* filename, HistsMap h_map, bool compute_w2=true, bool logy=false, float min_y=0.5 ) const {
      Canvas c( filename, top_label_ );
      THStack hs;
      TH1D* hist = 0;
      unsigned short i = 0;
      double max_bin = 0.;
      for ( HistsMap::iterator it=h_map.begin(); it!=h_map.end(); it++ ) {
        hist = ( TH1D* )it->second;
        if ( i==0 ) hs.SetTitle( hist->GetTitle() );
        if ( compute_w2 ) hist->Sumw2();
        hist->SetMarkerStyle( marker_pool_[i] );
        hist->SetMarkerColor( colour_pool_[i] );
        hist->SetLineColor( kBlack );
        if ( strcmp( it->first.c_str(), "" )!=0 ) c.AddLegendEntry( hist, it->first.c_str(), ( compute_w2 ) ? "elp" : "lp" );
        hs.Add( hist );
        i++;
      }
      hs.Draw( "nostack" );
      c.Prettify( hs.GetHistogram() );
      if ( logy ) {
        c.SetLogy();
        hs.GetYaxis()->SetRangeUser( min_y, max_bin*3 );
        hs.SetMinimum( min_y );
      }
      else {
        hs.GetYaxis()->SetRangeUser( 0., max_bin*1.1 );
      }
      hs.SetTitle( "" );
      c.Save( "pdf,png", out_path_ );
    }

    void draw_4plots( const char* filename, HistsMap h_map, Option_t* draw_style="colz" ) const {
      Canvas c( filename, top_label_ );
      c.Divide( 2, 2 );
      TH1D* hist = 0;
      unsigned int i = 0;
      for ( HistsMap::iterator it=h_map.begin(); it!=h_map.end(); it++ ) {
        if ( i>3 ) break;
        hist = ( TH1D* )it->second;
        TPad* pad = (TPad*)c.GetPad( i+1 );
        pad->SetTicks();
        pad->SetLeftMargin( 0.15 );
        pad->SetRightMargin( 0.15 );
        pad->SetTopMargin( 0.1 );
        pad->SetBottomMargin( 0.13 );
        c.cd( i+1 );
        hist->Draw( draw_style );
        c.Prettify( hist );

        PaveText* title = new PaveText( 0.82, 0.85, 0.82, 0.85 );
        title->SetTextSize( 0.05 );
        title->SetTextAlign( kHAlignRight+kVAlignTop );
        title->AddText( it->first.c_str() );
        title->Draw( "same" );

        hist->GetXaxis()->SetLabelSize( 12 );
        hist->GetYaxis()->SetLabelSize( 12 );
        hist->GetZaxis()->SetLabelSize( 0.045 );
        hist->GetXaxis()->SetTitleOffset( 2.0 );
        hist->GetYaxis()->SetTitleOffset( 2.0 );
        hist->GetXaxis()->SetTitleSize( 16 );
        hist->GetYaxis()->SetTitleSize( 16 );
        i++;
      }
      c.Save( "pdf,png", out_path_ );
    }
    void draw_stonratio( const char* filename, TH1D* h_sig, TH1D* h_bck ) const {
      Canvas c( filename, top_label_ );
      TH1D* ratio = dynamic_cast<TH1D*>( h_sig->Clone( "ratio" ) );
      ratio->Clear();
      for ( int i=0; i<=ratio->GetNbinsX(); i++ ) {
        const double sig = h_sig->GetBinContent( i ), bck = h_bck->GetBinContent( i );
        const double r = ( sig+bck!=0. ) ? sig/sqrt( sig+bck ) : 0.;
        ratio->SetBinContent( i, r );
        ratio->SetBinError( i, 0. ); //FIXME
      }
      ratio->Draw();
      c.Prettify( ratio );
      c.SetLogy();
      c.Save( "pdf,png", out_path_ );
    }

  private:

    void draw_diagonal( const float& min, const float& max, const float& x_resol=-1., const float& y_resol=-1., bool abs_unc=false) const {
      //FIXME to do: implement x resolution
      if (y_resol>0.) {
        TGraph* sigmay_pm = new TGraph(4);
        if (!abs_unc) {
          sigmay_pm->SetPoint(0, min, min+y_resol*min);
          sigmay_pm->SetPoint(1, max, max+y_resol*max);
          sigmay_pm->SetPoint(2, max, max-y_resol*max);
          sigmay_pm->SetPoint(3, min, min-y_resol*min);
        }
        else {
          sigmay_pm->SetPoint(0, min, min+y_resol);
          sigmay_pm->SetPoint(1, max, max+y_resol);
          sigmay_pm->SetPoint(2, max, max-y_resol);
          sigmay_pm->SetPoint(3, min, min-y_resol);
        }
        sigmay_pm->SetFillColorAlpha(kBlack, 0.1);
        sigmay_pm->Draw("f");
      }
      TLine l;
      l.SetLineWidth( 3 );
      l.SetLineColor( kGray );
      l.SetLineStyle( 2 );
      l.DrawLine( min, min, max, max );
    }

    TString out_path_;
    TString top_label_;
    static const int* marker_pool_;
    static const int* colour_pool_;

};

static const int markers[] = { 24, 20, 25, 21, 26, 22, 27, 23, 28, 24 };
static const int colours[] = { kBlack, kRed+1, kGreen+2, kBlue+1, kMagenta+1, kOrange+1, kGray };

const int* Plotter::marker_pool_ = markers;
const int* Plotter::colour_pool_ = colours;

#endif
