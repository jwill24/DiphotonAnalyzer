#include "Canvas.h"

void test_rap()
{
  //TFile f( "Samples/output_DiPhotonJetsBox_MGG-80toInf_Sherpa_v3.root" );
  TFile f( "Samples/output_GammaGammaToGammaGamma_fpmc_r9gt0.root" );
  TTree* t = ( TTree* )f.Get( "ntp" );
  unsigned int num_diphoton;
  float diphoton_mass[10], diphoton_rap[10];
  float diphoton_pt1[10], diphoton_eta1[10], diphoton_phi1[10], diphoton_energy1[10];
  float diphoton_pt2[10], diphoton_eta2[10], diphoton_phi2[10], diphoton_energy2[10];
  t->SetBranchAddress( "num_diphoton", &num_diphoton );
  t->SetBranchAddress( "diphoton_mass", diphoton_mass );
  t->SetBranchAddress( "diphoton_rapidity", diphoton_rap );
  t->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  t->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  t->SetBranchAddress( "diphoton_phi1", diphoton_phi1 );
  t->SetBranchAddress( "diphoton_energy1", diphoton_energy1 );
  t->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  t->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  t->SetBranchAddress( "diphoton_phi2", diphoton_phi2 );
  t->SetBranchAddress( "diphoton_energy2", diphoton_energy2 );

  TH1D* h_mass_meth1 = new TH1D( "mass_meth1", "", 200, 200., 2200. ),
       *h_mass_meth2 = ( TH1D* )h_mass_meth1->Clone( "mass_meth2" );
  TH1D* h_rap_meth1 = new TH1D( "rap_meth1", "", 100, -5., 5. ),
       *h_rap_meth2 = ( TH1D* )h_rap_meth1->Clone( "rap_meth2" );
  TLorentzVector lep1, lep2;
  for ( unsigned long long i=0; i<t->GetEntries(); i++ ) {
    t->GetEntry( i );
    for ( unsigned int j=0; j<num_diphoton; j++ ) {
      /*lep1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      lep2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );*/
      lep1.SetPtEtaPhiE( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], diphoton_energy1[j] );
      lep2.SetPtEtaPhiE( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], diphoton_energy2[j] );
      cout << lep1.M() << "\t" << lep2.M() << endl;

      h_mass_meth1->Fill( ( lep1+lep2 ).M() );
      h_mass_meth2->Fill( diphoton_mass[j] );
      h_rap_meth1->Fill( ( lep1+lep2 ).Rapidity() );
      h_rap_meth2->Fill( diphoton_rap[j] );
    }
  }

  {
    Canvas c( "rap_comparison", "" );
    //h_rap_meth2->Draw( );
    //h_rap_meth2->SetLineColor( kRed );
    h_rap_meth1->Divide( h_rap_meth2 );
    h_rap_meth1->Draw();
    c.Prettify( h_rap_meth1 );
    c.Save( "png,pdf" );
  }
  {
    Canvas c( "mass_comparison", "" );
    //h_mass_meth2->Draw( );
    //h_mass_meth2->SetLineColor( kRed );
    h_mass_meth1->Divide( h_mass_meth2 );
    h_mass_meth1->Draw();
    c.Prettify( h_mass_meth1 );
    c.Save( "png,pdf" );
  }
}
