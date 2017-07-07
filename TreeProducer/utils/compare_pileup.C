#include "DiphotonAnalyzer/Macros/Canvas.h"

void compare_pileup()
{
  TFile f1( "data/pileup_data16BCG_PPSruns_v1.root" ), f2( "data/pileup_data16BCG_PPSruns_v2.root" ), f_mc( "data/pileup_mc.root" );
  auto pu_1 = dynamic_cast<TH1D*>( f1.Get( "pileup" ) ), pu_2 = dynamic_cast<TH1D*>( f2.Get( "pileup" ) ), pu_mc = dynamic_cast<TH1D*>( f_mc.Get( "pileup" ) );

  pu_1->Scale( 1./pu_1->Integral() );
  pu_2->Scale( 1./pu_2->Integral() );
  pu_mc->Scale( 1./pu_mc->Integral() );

  Canvas c( "pileup_dist_comp" );
  THStack st;
  c.SetLegendX1( 0.42 );
  st.Add( pu_mc, "hist" );
  pu_mc->SetLineStyle( 2 );
  pu_mc->SetLineWidth( 2 );
  pu_mc->SetLineColor( kBlack );
  c.AddLegendEntry( pu_mc, "MC (25ns_Moriond17MC)" );
  st.Add( pu_1, "hist" );
  pu_1->SetLineWidth( 2 );
  pu_1->SetLineColor( kBlack );
  c.AddLegendEntry( pu_1, "Old" );
  st.Add( pu_2, "p" );
  pu_2->SetMarkerStyle( 20 );
  pu_2->SetMarkerColor( kRed+1 );
  c.AddLegendEntry( pu_2, "Corrected", "p" );
  st.Draw( "nostack" );
  st.GetHistogram()->SetTitle( "Number of primary vertices\\Events fraction" );
  c.Prettify( st.GetHistogram() );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
