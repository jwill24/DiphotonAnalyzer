void
optics_comparator()
{
  TFile* f[2];
  f[0] = new TFile( "data/ctpps_optics_v1.root" );
  f[1] = new TFile( "data/ctpps_optics_9mar2017.root" );
  TGraph* g_x_to_xi_R_1_N[2],
         *g_x_to_xi_R_1_F[2],
         *g_x_to_xi_L_1_N[2],
         *g_x_to_xi_L_1_F[2];

  TMultiGraph mg_45n, mg_45f, mg_56n, mg_56f;

  for ( unsigned short i=0; i<2; i++ ) {
    g_x_to_xi_R_1_N[i] = ( TGraph* )f[i]->Get( "g_x_to_xi_R_1_N" ),
    g_x_to_xi_R_1_F[i] = ( TGraph* )f[i]->Get( "g_x_to_xi_R_1_F" ),
    g_x_to_xi_L_1_N[i] = ( TGraph* )f[i]->Get( "g_x_to_xi_L_1_N" ),
    g_x_to_xi_L_1_F[i] = ( TGraph* )f[i]->Get( "g_x_to_xi_L_1_F" );
    g_x_to_xi_L_1_N[i]->SetLineColor( ( i==0 ) ? kBlack : kRed+1 );
    g_x_to_xi_L_1_F[i]->SetLineColor( ( i==0 ) ? kBlack : kRed+1 );
    g_x_to_xi_R_1_N[i]->SetLineColor( ( i==0 ) ? kBlack : kRed+1 );
    g_x_to_xi_R_1_F[i]->SetLineColor( ( i==0 ) ? kBlack : kRed+1 );
    mg_45n.Add( g_x_to_xi_L_1_N[i] );
    mg_45f.Add( g_x_to_xi_L_1_F[i] );
    mg_56n.Add( g_x_to_xi_R_1_N[i] );
    mg_56f.Add( g_x_to_xi_R_1_F[i] );
  }
  TLegend leg( 0.5, 0.7, 0.7, 0.8 );
  leg.AddEntry( g_x_to_xi_L_1_N[0], "old optics file", "l" );
  leg.AddEntry( g_x_to_xi_L_1_N[1], "new optics file", "l" );

  TCanvas c;
  c.Divide( 2, 2 );
  c.cd( 1 );
  mg_45n.Draw( "al" );
  mg_45n.GetXaxis()->SetTitle( "45N" );
  c.cd( 2 );
  mg_45f.Draw( "al" );
  mg_45f.GetXaxis()->SetTitle( "45F" );
  c.cd( 3 );
  mg_56n.Draw( "al" );
  mg_56n.GetXaxis()->SetTitle( "56N" );
  c.cd( 4 );
  mg_56f.Draw( "al" );
  mg_56f.GetXaxis()->SetTitle( "56F" );
  c.cd();
  leg.Draw();
  c.SaveAs( "comparison_optics_curves.pdf" );
}
