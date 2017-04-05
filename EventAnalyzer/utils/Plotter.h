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
    typedef std::vector< std::pair<const char*, TGraphErrors*> > GraphsMap;
    typedef std::vector< std::pair<const char*, TH1*> > HistsMap;

  public:
    Plotter( const char* out_path, const char* top_label ) :
      out_path_( out_path ),
      canv_( 0 ),
      canv_template_( new Canvas( "canv_template", top_label ) )
    {
      reset_canvas();
    }
    ~Plotter() {
      cout << __PRETTY_FUNCTION__ << endl;
      if ( canv_ ) delete canv_;
cout << "canv del" << endl;
      if ( canv_template_ ) delete canv_template_;
cout << "canv tmpl del" << endl;
    }

    Canvas* BaseCanvas() { return canv_; }

    void plot_multihists( const char* name, HistsMap hm, float min_ratio_y=0., float max_ratio_y=0.55, bool draw_overflow=true )
    {
      if ( hm.size()==0 ) return;

      Canvas* c = new Canvas( *canv_, name );
      if ( hm.size()>1 ) c->SetRatioPlot();

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
        //h->SetMarkerStyle( 20 );
        min = TMath::Min( min, hist->GetMinimum() );
        max = TMath::Max( max, hist->GetMaximum() );
        if ( i>2 ) hist->SetFillStyle( 3005 );
        if ( hm.size()>1 ) {
          c->AddLegendEntry( hist, it->first, ( i==0 ) ? "f" : "ep" );
        }
        i++;
      }
      hist = dynamic_cast<TH1*>( hm.begin()->second );
      hist->SetMinimum( min );
      hist->SetMaximum( max*1.1 );
      c->Prettify( hist );
      /*if ( hm.size()==2 ) c->RatioPlot( hm[0].second, hm[1].second, 0, min_ratio_y, max_ratio_y );
      if ( hm.size()>2 ) c->RatioPlot( hm[0].second, hm[1].second, hm[2].second, min_ratio_y, max_ratio_y );*/
      c->RatioPlot( hm2, min_ratio_y, max_ratio_y );
      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

    void plot_3hists( const char* name, TH1D* h, TH1D* h_1tag=0, TH1D* h_2tag=0, bool draw_overflow=true )
    {
      HistsMap hm;
      hm.push_back( std::make_pair( "All diphotons", h ) );
      hm.push_back( std::make_pair( "track(s) in 1 arm", h_1tag ) );
      hm.push_back( std::make_pair( "track(s) in 2 arms", h_2tag ) );
      plot_multihists( name, hm, 0., 0.55, draw_overflow );
    }

    void plot_balances( const char* name, const char* title, GraphsMap gr_map,
                        const float& min_xy=-999., const float& max_xy=-999., const float& rp_acc=-1. )
    {
      Canvas* c = new Canvas( *canv_, name );
      //c->DrawFrame( min_xy, min_xy, max_xy, max_xy );

      TMultiGraph mg;
      TGraphErrors* gr = 0;

      unsigned short i = 0;
      for ( GraphsMap::iterator it=gr_map.begin(); it!=gr_map.end(); it++ ) {
        gr = ( TGraphErrors* )it->second;
        mg.Add( gr, "ep" );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        //c->SetLegendY1( 0.18 );
        if ( gr_map.size()>0 && strcmp( it->first, "" )!=0 ) c->AddLegendEntry( gr, it->first, "p" );
        i++;
      }


      mg.Draw( "ap" );
      mg.GetHistogram()->SetTitle( title );
      c->Prettify( mg.GetHistogram() );
      draw_diagonal( min_xy, max_xy, 0. );
      if ( min_xy!=-999. and max_xy!=-999. ) {
        mg.GetXaxis()->SetLimits( min_xy, max_xy );
        mg.GetYaxis()->SetRangeUser( min_xy, max_xy );
      }

      if ( rp_acc>0. ) {
        TPave* acc = new TPave( min_xy, min_xy, max_xy, rp_acc, 1 );
        acc->SetLineWidth( 1 );
        acc->SetLineColor( kBlack );
        acc->SetFillColorAlpha( kGray, 0.4 );
        //acc->SetFillColor( kGray );
        acc->Draw( "same" );
      }
      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

    void plot_balances( const char* name, const char* title, TGraphErrors* gr, const float& min_xy=-999., const float& max_xy=-999., const float& rp_acc=-1. )
    {
      GraphsMap gm; gm.push_back( std::pair<const char*, TGraphErrors*>( gr->GetTitle(), gr ) );
      plot_balances( name, title, gm, min_xy, max_xy, rp_acc );
    }

    void plot_xi_correlations( const char* sector, GraphsMap graphs_map )
    {
      Canvas* c = new Canvas( *canv_, Form( "xi_nearfar_corr_%s", sector ) );
      TMultiGraph mg;

      TGraphErrors* gr = 0;

      unsigned short i = 0;
      for ( GraphsMap::iterator it=graphs_map.begin(); it!=graphs_map.end(); it++ ) {
        gr = ( TGraphErrors* )it->second;
        mg.Add( gr );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        if ( strcmp( it->first, "" )!=0 ) c->AddLegendEntry( gr, it->first, "p" );
        i++;
      }

      mg.Draw( "ap" );
      mg.GetHistogram()->SetTitle( Form( "#xi_{%s} (near pot)\\#xi_{%s} (far pot)", sector, sector ) );

      c->Prettify( mg.GetHistogram() );
      draw_diagonal( 0., 0.2 );

      mg.GetXaxis()->SetLimits( 0., 0.2 );
      mg.GetYaxis()->SetRangeUser( 0., 0.2 );

      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

    void draw_multigraph( const char* name, GraphsMap g_map, float lim_min=-999., float lim_max=-999., bool draw_diag=false, float y_limit=-1. )
    {
      Canvas* c = dynamic_cast<Canvas*>( canv_->Clone( name ) );
      //if ( lim_min!=-999. && lim_max!=-999. )
      //  c->DrawFrame( lim_min, lim_min, lim_max, lim_max );

      TGraphErrors* gr = 0;

      TMultiGraph mg;
      unsigned short i = 0;
      for ( GraphsMap::iterator it=g_map.begin(); it!=g_map.end(); it++ ) {
        gr = ( TGraphErrors* )it->second;
        mg.Add( gr );
        if ( strcmp( gr->GetTitle(), "" )!=0 ) mg.SetTitle( gr->GetTitle() );
        gr->SetMarkerStyle( marker_pool_[i] );
        gr->SetMarkerColor( colour_pool_[i] );
        if ( strcmp( it->first, "" )!=0 ) c->AddLegendEntry( gr, it->first, "p" );
        i++;
      }

      mg.Draw( "ap" );
      c->Prettify( mg.GetHistogram() );
      if ( draw_diag ) draw_diagonal( lim_min, lim_max );

      if ( lim_min>0 || lim_max>0. ) {
        mg.GetXaxis()->SetLimits( ( lim_min>0 ? lim_min : 0. ), ( lim_max>0 ? lim_max : 500. ) );
        mg.GetYaxis()->SetRangeUser( ( lim_min>0 ? lim_min : 0. ), ( lim_max>0 ? lim_max : 500. ) );
      }

      if ( y_limit>0. && y_limit<lim_max ) {
        TPave* acc = new TPave( lim_min, lim_min, lim_max, y_limit, 1 );
        acc->SetLineWidth( 1 );
        acc->SetLineColor( kBlack );
        acc->SetFillColorAlpha( kGray, 0.4 );
        acc->Draw( "same" );
      }
      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

    void draw_multiplot( const char* name, HistsMap h_map_data, HistsMap h_map_mc, HistsMap h_map_sig, bool logy=false )
    {
      Canvas* c = new Canvas( *canv_, name );
      c->SetRatioPlot();
      TH1D* hist = 0;
      TH1D* h_data, *h_mc = 0;
      double max_bin = -1.;
      unsigned short i = 0;
      THStack hs_mc, hs_data, hs_sig;
      for ( HistsMap::iterator it=h_map_data.begin(); it!=h_map_data.end(); it++ ) {
        hist = ( TH1D* )it->second;
        if ( i==0 ) h_data = dynamic_cast<TH1D*>( hist );
        else h_data->Add( hist );
        // draw the data distributions unstacked
        hist->Sumw2();
        hist->SetMarkerStyle( 20+i );
        hist->SetMarkerColor( kBlack );
        hist->SetLineColor( kBlack );
        if ( strcmp( it->first, "" )!=0 ) c->AddLegendEntry( hist, it->first, "lp" );
        max_bin = TMath::Max( max_bin, hist->GetMaximum() );
        hs_data.Add( hist );
        i++;
      }
      i = 0;
      for ( HistsMap::iterator it=h_map_mc.begin(); it!=h_map_mc.end(); it++ ) {
        hist = dynamic_cast<TH1D*>( it->second );
        if ( i==0 ) h_mc = dynamic_cast<TH1D*>( hist->Clone() );
        else { h_mc->Add( hist ); }
        hist->SetFillColorAlpha( colour_pool_[i+1], 0.5 );
        //hist->SetFillStyle( 3002 );
        hist->SetLineColor( kBlack );
        if ( strcmp( it->first, "" )!=0 ) { c->AddLegendEntry( hist, it->first, "f" ); }
        hs_mc.Add( hist );
        if ( i==0 ) hs_mc.SetTitle( hist->GetTitle() );
        i++;
      }
      i = 0;
      for ( HistsMap::iterator it=h_map_sig.begin(); it!=h_map_sig.end(); it++ ) {
        hist = dynamic_cast<TH1D*>( it->second );
        hist->SetLineColor( kGreen );
        hist->SetLineWidth( 3 );
        hist->SetLineStyle( i );
        if ( strcmp( it->first, "" )!=0 ) { c->AddLegendEntry( hist, it->first, "l" ); }
        hs_sig.Add( hist );
        if ( i==0 ) hs_sig.SetTitle( hist->GetTitle() );
        i++;
      }
      hs_mc.Draw( "hist" );
      hs_sig.Draw( "hist same nostack" );
      hs_data.Draw( "same nostack" );
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
      c->Prettify( hist );

      HistsMap hm;
      hm.push_back( std::make_pair( "mc", h_mc ) );
      hm.push_back( std::make_pair( "data", h_data ) );
      hs_mc.SetTitle( "" );
      c->RatioPlot( hm, -0.4, 2.4, 1.0 );
      c->cd( 1 );
      if ( logy ) dynamic_cast<TPad*>( c->GetPad( 1 ) )->SetLogy();
      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

    void draw_multiplot( const char* name, HistsMap h_map, bool compute_w2=true, bool logy=false )
    {
      Canvas* c = new Canvas( *canv_, name );
      TH1D* hist = 0;
      unsigned short i = 0;
      double max_bin = 0.;
      for ( HistsMap::iterator it=h_map.begin(); it!=h_map.end(); it++ ) {
        hist = ( TH1D* )it->second;
        if ( compute_w2 ) {
          hist->Sumw2();
          hist->Draw( ( i>0 ) ? "p same" : "p" );
        }
        else hist->Draw( ( i>0 ) ? "same" : "" );
        hist->SetMarkerStyle( marker_pool_[i] );
        hist->SetMarkerColor( colour_pool_[i] );
        hist->SetLineColor( kBlack );
        if ( strcmp( it->first, "" )!=0 ) c->AddLegendEntry( hist, it->first, ( compute_w2 ) ? "elp" : "lp" );
        max_bin = TMath::Max( max_bin, hist->GetMaximum() );
        i++;
      }
      c->Prettify( h_map.begin()->second );
      if ( logy ) {
        h_map.begin()->second->GetYaxis()->SetRangeUser( 0.5, max_bin*3 );
        c->SetLogy();
      }
      else {
        h_map.begin()->second->GetYaxis()->SetRangeUser( 0., max_bin*1.1 );
      }
      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

    void draw_4plots( const char* name, HistsMap h_map, Option_t* draw_style="colz" )
    {
      Canvas* c = new Canvas( *canv_, name );
      c->Divide( 2, 2 );
      TH1D* hist = 0;
      unsigned int i = 0;
      for ( HistsMap::iterator it=h_map.begin(); it!=h_map.end(); it++ ) {
        if ( i>3 ) break;
        hist = ( TH1D* )it->second;
        TPad* pad = dynamic_cast<TPad*>( c->GetPad( i+1 ) );
        pad->SetTicks();
        pad->SetLeftMargin( 0.15 );
        pad->SetRightMargin( 0.15 );
        pad->SetTopMargin( 0.1 );
        pad->SetBottomMargin( 0.13 );
        c->cd( i+1 );
        hist->Draw( draw_style );
        c->Prettify( hist );

        TPaveText* title = new TPaveText( 0.82, 0.85, 0.82, 0.85, "NDC" );
        title->SetFillStyle( 0 );
        title->SetLineColor( 0 );
        title->SetTextFont( 42 );
        title->SetTextSize( 0.05 );
        title->SetTextAlign( kHAlignRight+kVAlignTop );
        title->AddText( it->first );
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
      c->Save( "pdf,png", out_path_ );
      delete c;
      reset_canvas();
    }

  private:
    void reset_canvas()
    {
cout << __PRETTY_FUNCTION__ << endl;
return; //FIXME
      if ( !canv_template_ ) return;
cout << "mid0>>>" << endl;
      if ( canv_ ) delete canv_;

cout << "mid>>>" << endl;
      canv_ = new Canvas( *canv_template_, "" );
      //*canv_ = *canv_template_;
cout << "end>>>" << endl;
    }

    void draw_diagonal( const float& min, const float& max, const float& x_resol=-1., const float& y_resol=-1., bool abs_unc=false) const
    {
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

    Canvas* canv_;
    Canvas* canv_template_;

    TString out_path_;
    static const int* marker_pool_;
    static const int* colour_pool_;

};

static const int markers[] = { 24, 20, 25, 21, 26, 22, 27, 23, 28, 24 };
static const int colours[] = { kBlack, kRed+1, kGreen+2, kBlue+1, kMagenta+1, kOrange+1, kGray };

const int* Plotter::marker_pool_ = markers;
const int* Plotter::colour_pool_ = colours;

#endif
