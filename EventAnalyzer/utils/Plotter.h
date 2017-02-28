#ifndef Plotter_h
#define Plotter_h

#include "Canvas.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TStyle.h"
#include <map>
#include <vector>

class Plotter
{
  public:
    typedef std::map<TString, TGraph*> GraphsMap;
    typedef std::vector< std::pair<TString, TH1*> > HistsMap;

  public:
    Plotter( const char* out_path, const char* top_label ) : out_path_( out_path ), top_label_( top_label ) {}
    ~Plotter() {}

    void plot_3hists( const char* name, TH1* h, TH1* h_1tag=0, TH1* h_2tag=0 ) const {
      Canvas c( name, top_label_, true );
      h->Draw();
      h->SetLineColor( kBlack );
      double min = TMath::Min( 0.001, h->GetMinimum() );
      //h->SetMarkerStyle( 20 );
      if ( h_1tag ) {
        h_1tag->Sumw2();
        h_1tag->Draw("e1 same");
        h_1tag->SetLineColor( kBlack );
        h_1tag->SetMarkerStyle( 20 );
        h_1tag->SetMarkerColor( kRed+1 );
        min = TMath::Min( min, h_1tag->GetMinimum() );
      }
      if ( h_2tag ) {
        h_2tag->Sumw2();
        h_2tag->Draw("e1 same");
        h_2tag->SetLineColor( kBlack );
        h_2tag->SetMarkerStyle( 24 );
        h_2tag->SetMarkerColor( kGreen+2 );
        h_2tag->SetFillStyle( 3005 );
        h_2tag->SetFillColor( kGreen+2 );
        min = TMath::Min( min, h_2tag->GetMinimum() );
      }
      if ( h_1tag && h_2tag ) {
        c.AddLegendEntry( h, "All diphotons", "f" );
        if ( h_1tag ) c.AddLegendEntry( h_1tag, "track(s) in 1 arm", "ep" );
        if ( h_2tag ) c.AddLegendEntry( h_2tag, "track(s) in 2 arms", "ep" );
      }
      h->SetMinimum( min );
      c.Prettify( h );
      if ( h_1tag && h_2tag )c .RatioPlot( h, h_1tag, h_2tag, 0., 0.55 );
      c.Save( "png,pdf", out_path_ );
    }

    void plot_balances( const char* name, const char* title,
                        TGraphErrors* h2, TGraphErrors* h2_mtag=0, TGraphErrors* h2_ytag=0,
                        const float& min_xy=-1., const float& max_xy=-1., const float& rp_acc=-1. ) const {
      Canvas c( name, top_label_ );
      c.DrawFrame( min_xy, min_xy, max_xy, max_xy );

      TMultiGraph mg;
      h2->SetMarkerStyle( 24 );
      //c.SetLegendY1( 0.18 );
      if ( h2_mtag || h2_ytag ) c.AddLegendEntry( h2, "All candidates", "p" );

      mg.Add( h2, "ep" );
      if ( h2_mtag ) {
        h2_mtag->SetMarkerStyle( 20 );
        h2_mtag->SetMarkerSize( .95 );
        h2_mtag->SetMarkerColor( kRed+1 );
        c.AddLegendEntry( h2_mtag, "Mass matching", "p" );
        mg.Add( h2_mtag, "p" );
      }
      if ( h2_ytag ) {
        h2_ytag->SetMarkerStyle( 31 );
        h2_ytag->SetMarkerSize( .7 );
        h2_ytag->SetMarkerColor( kGreen+2 );
        c.AddLegendEntry( h2_ytag, "Rapidity matching", "p" );
        mg.Add( h2_ytag, "p" );
      }

      mg.Draw( "same" );
      if ( min_xy!=-1. and max_xy!=-1. ) {
        mg.GetXaxis()->SetLimits( min_xy, max_xy );
        mg.GetYaxis()->SetRangeUser( min_xy, max_xy );
      }

      mg.GetHistogram()->SetTitle( title );
      c.Prettify( mg.GetHistogram() );

      draw_diagonal( min_xy, max_xy, 0. );
      if ( rp_acc>0. ) {
        TPave* acc = new TPave( min_xy, min_xy, max_xy, rp_acc, 1 );
        acc->SetLineWidth( 1 );
        acc->SetLineColor( kBlack );
        acc->SetFillColorAlpha( kGray, 0.4 );
        //acc->SetFillColor( kGray );
        acc->Draw( "same" );
      }


      c.Save( "pdf,png", out_path_ );
    }

    void plot_balances_old( const char* name, const char* title,
                            TGraphErrors* h2, TGraphErrors* h2_mtag=0, TGraphErrors* h2_ytag=0,
                            const float& min_xy=-1., const float& max_xy=-1., const float& mass_resol=-1., bool abs_unc=false ) const {
      Canvas c( name, top_label_ );
      TMultiGraph mg;
      h2->SetMarkerStyle( 24 );
      //c.SetLegendY1( 0.18 );
      if ( h2_mtag or h2_ytag ) c.AddLegendEntry( h2, "All candidates", "p" );
      mg.Add( h2, "ep" );
      if ( h2_mtag ) {
        h2_mtag->Draw( "p same" );
        h2_mtag->SetMarkerStyle( 20 );
        h2_mtag->SetMarkerSize( .95 );
        h2_mtag->SetMarkerColor( kRed+1 );
        c.AddLegendEntry( h2_mtag, "Mass matching", "p" );
        mg.Add( h2_mtag );
      }
      if ( h2_ytag ) {
        h2_ytag->Draw( "p same" );
        h2_ytag->SetMarkerStyle( 31 );
        h2_ytag->SetMarkerSize( .7 );
        h2_ytag->SetMarkerColor( kGreen+2 );
        c.AddLegendEntry( h2_ytag, "Rapidity matching", "p" );
        mg.Add( h2_ytag );
      }
      mg.Draw( "ap" );
      if ( min_xy!=-1. and max_xy!=-1. ) {
        mg.GetXaxis()->SetLimits( min_xy, max_xy );
        mg.GetYaxis()->SetRangeUser( min_xy, max_xy );
      }
      if ( mass_resol>0. ) {
        for ( unsigned int i=0; i<h2->GetN(); i++ ) {
          double x, y; h2->GetPoint( i, x, y );
          if ( !abs_unc ) h2->SetPointError( i, 0., y*mass_resol );
          else            h2->SetPointError( i, 0., mass_resol );
        }
      }

      draw_diagonal( min_xy, max_xy, 0., mass_resol, abs_unc );
      mg.GetHistogram()->SetTitle( title );
      c.Prettify( mg.GetHistogram() );
      c.Save( "pdf,png", out_path_ );
    }

    void plot_xi_correlations( const char* sector, GraphsMap graphs_map ) const {
      Canvas c( Form( "xi_nearfar_corr_%s", sector ), top_label_ );
      TMultiGraph mg;

      TGraph* gr = 0;

      unsigned short i = 0;
      for ( GraphsMap::iterator it=graphs_map.begin(); it!=graphs_map.end(); it++ ) {
        gr = ( TGraph* )it->second->Clone();
        mg.Add( gr );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        c.AddLegendEntry( gr, it->first, "p" );
        i++;
      }

      mg.Draw( "ap" );
      mg.GetXaxis()->SetLimits( 0., 0.2 );
      mg.GetYaxis()->SetRangeUser( 0., 0.2 );
      mg.GetHistogram()->SetTitle( Form( "#xi_{%s} (near pot)\\#xi_{%s} (far pot)", sector, sector ) );
      c.Prettify( mg.GetHistogram() );

      draw_diagonal( 0., 0.2 );

      c.Save( "pdf,png", out_path_ );
      delete gr;
    }

    void draw_multiplot( const char* filename, HistsMap h_map, bool compute_w2=true ) const {
      Canvas c( filename, top_label_ );
      TH1* hist = 0;
      unsigned short i = 0;
      double max_bin = 0.;
      for ( HistsMap::iterator it=h_map.begin(); it!=h_map.end(); it++ ) {
        hist = ( TH1* )it->second;
        if ( compute_w2 ) hist->Sumw2();
        hist->Draw( ( i>0 ) ? "same" : "" );
        hist->SetMarkerStyle( marker_pool_[i] );
        hist->SetMarkerColor( colour_pool_[i] );
        hist->SetLineColor( kBlack );
        c.AddLegendEntry( hist, it->first, ( compute_w2 ) ? "elp" : "lp" );
        max_bin = TMath::Max( max_bin, hist->GetMaximum() );
        i++;
      }
      c.Prettify( h_map.begin()->second );
      h_map.begin()->second->GetYaxis()->SetRangeUser( 0., max_bin*1.1 );
      c.Save( "pdf,png", out_path_ );

    }

    void draw_4hitmaps( const char* filename, HistsMap h_map ) const {
      Canvas c( filename, top_label_ );
      c.Divide( 2, 2 );
      TH2* hist = 0;
      unsigned int i = 0;
      for ( HistsMap::iterator it=h_map.begin(); it!=h_map.end(); it++ ) {
        if ( i>3 ) break;
        hist = ( TH2* )it->second;
        TPad* pad = (TPad*)c.GetPad( i+1 );
        pad->SetTicks();
        pad->SetLeftMargin( 0.15 );
        pad->SetBottomMargin( 0.15 );
        c.cd( i+1 );
        hist->Draw( "colz" );
        c.Prettify( hist );
        hist->GetXaxis()->SetLabelSize( 12 );
        hist->GetYaxis()->SetLabelSize( 12 );
        hist->GetXaxis()->SetTitleOffset( 2.0 );
        hist->GetYaxis()->SetTitleOffset( 2.0 );
        hist->GetXaxis()->SetTitleSize( 16 );
        hist->GetYaxis()->SetTitleSize( 16 );
        i++;
      }
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

static const int markers[] = { 20, 24, 21, 25, 22, 26, 23, 27, 24, 28 };
static const int colours[] = { kBlack, kBlue+1, kRed+1, kGreen+1, kMagenta, kOrange, kGray };

const int* Plotter::marker_pool_ = markers;
const int* Plotter::colour_pool_ = colours;

#endif
