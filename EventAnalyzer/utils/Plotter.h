#include "Canvas.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TStyle.h"

class Plotter
{
  public:
    Plotter( const char* out_path, const char* top_label ) : out_path_( out_path ), top_label_( top_label ) {}
    ~Plotter() {}

    void plot_3hists( const char* name, TH1* h, TH1* h_1tag, TH1* h_2tag ) const {
      Canvas c( name, top_label_, true );
      h_1tag->Sumw2();
      h_2tag->Sumw2();
      h->Draw();
      h->SetLineColor( kBlack );
      double min = TMath::Min( 0.001, h->GetMinimum() );
      //h->SetMarkerStyle( 20 );
      h_1tag->Draw("e1 same");
      h_1tag->SetLineColor( kBlack );
      h_1tag->SetMarkerStyle( 20 );
      h_1tag->SetMarkerColor( kRed+1 );
      min = TMath::Min( min, h_1tag->GetMinimum() );
      h_2tag->Draw("e1 same");
      h_2tag->SetLineColor( kBlack );
      h_2tag->SetMarkerStyle( 24 );
      h_2tag->SetMarkerColor( kGreen+2 );
      h_2tag->SetFillStyle( 3005 );
      h_2tag->SetFillColor( kGreen+2 );
      min = TMath::Min( min, h_2tag->GetMinimum() );
      c.AddLegendEntry( h, "All diphotons", "f" );
      c.AddLegendEntry( h_1tag, "#geq 1 proton tag", "ep" );
      c.AddLegendEntry( h_2tag, "#geq 2 proton tags", "ep" );
      h->SetMinimum( min );
      c.Prettify( h );
      c.RatioPlot( h, h_1tag, h_2tag, 0., 0.55 );
      c.Save( "png", out_path_ );
      c.Save( "pdf", out_path_ );
    }

    void plot_balances( const char* name, const char* title,
                        TGraphErrors* h2, TGraphErrors* h2_mtag=0, TGraphErrors* h2_ytag=0,
                        const float& min_xy=-1., const float& max_xy=-1., const float& rp_acc=-1. ) const {
      Canvas c( name, top_label_ );
      c.DrawFrame( min_xy, min_xy, max_xy, max_xy );
      draw_diagonal( min_xy, max_xy, 0. );

      if ( rp_acc>0. ) {
        TPave* acc = new TPave( min_xy, min_xy, max_xy, rp_acc, 1 );
        acc->SetLineWidth( 1 );
        acc->SetLineColor( kBlack );
        acc->SetFillColorAlpha( kGray, 0.4 );
        //acc->SetFillColor( kGray );
        acc->Draw( "same" );
      }

      TMultiGraph mg;
      h2->SetMarkerStyle( 24 );
      //c.SetLegendY1( 0.18 );
      if ( h2_mtag!=0 or h2_ytag!=0 ) c.AddLegendEntry( h2, "All candidates", "p" );

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
      c.Save( "pdf", out_path_ );
      c.Save( "png", out_path_ );
    }

    void plot_balances_old( const char* name, const char* title,
                            TGraphErrors* h2, TGraphErrors* h2_mtag=0, TGraphErrors* h2_ytag=0,
                            const float& min_xy=-1., const float& max_xy=-1., const float& mass_resol=-1., bool abs_unc=false ) const {
      Canvas c( name, top_label_ );
      TMultiGraph mg;
      h2->SetMarkerStyle( 24 );
      //c.SetLegendY1( 0.18 );
      if ( h2_mtag!=0 or h2_ytag!=0 ) c.AddLegendEntry( h2, "All candidates", "p" );
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
      c.Save( "pdf", out_path_ );
      c.Save( "png", out_path_ );
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
};
