#include "Canvas.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>

//#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/"
#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp/"
#define default_ntp_file "run2016BCG_21oct.root"

//#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/runBC"
//#define default_ntp_file "run2016BC_28sep.root"

//#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/matching_with_met"
//#define default_ntp_file "run2016BCG_21oct.root"

void plot_3hists( const char* name, const char* top_label, TH1* h, TH1* h_1tag, TH1* h_2tag )
{
  Canvas c( name, top_label, true );
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
  c.Save( "png", out_path );
  c.Save( "pdf", out_path );
}

void plot_balances_new( const char* name, const char* top_label, const char* title, TGraphErrors* h2, TGraphErrors* h2_mtag=0, TGraphErrors* h2_ytag=0,
                    const float& min_xy=-1., const float& max_xy=-1. )
{
  Canvas c( name, top_label );
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

  c.DrawDiagonal( min_xy, max_xy, 0. );
  mg.GetHistogram()->SetTitle( title );
  c.Prettify( mg.GetHistogram() );
  c.Save( "pdf", out_path );
  c.Save( "png", out_path );
}

void plot_balances( const char* name, const char* top_label, const char* title, TGraphErrors* h2, TGraphErrors* h2_mtag=0, TGraphErrors* h2_ytag=0,
                    const float& min_xy=-1., const float& max_xy=-1., const float& mass_resol=-1., bool abs_unc=false )
{
  Canvas c( name, top_label );
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

  c.DrawDiagonal( min_xy, max_xy, 0., mass_resol, abs_unc );
  mg.GetHistogram()->SetTitle( title );
  c.Prettify( mg.GetHistogram() );
  c.Save( "pdf", out_path );
  c.Save( "png", out_path );
}

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 )
{
  if ( r9<0.94 ) return 0.004;
  if ( fabs( eta )<1.4442 ) return 0.002;
  if ( fabs( eta )>1.566 ) return 0.005;
  return 1.;
}

void tree_reader_fallback( TString file=default_ntp_file )
{
  TFile f(file);
  if ( !f.IsOpen() ) return;

  const float sqrt_s = 13.e3;
  const bool compute_with_met = true;
  const double rel_err_xi = 0.05; // 15% error on xi determination

  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );
  // general quantities
  unsigned int run_id, lumisection;
  unsigned long long event_number;
  tr->SetBranchAddress( "run_id", &run_id );
  tr->SetBranchAddress( "lumisection", &lumisection );
  tr->SetBranchAddress( "event_number", &event_number );
  // diphoton quantities
  unsigned int num_diphoton;
  float diphoton_pt1[5], diphoton_eta1[5], diphoton_phi1[5], diphoton_r9_1[5],
        diphoton_pt2[5], diphoton_eta2[5], diphoton_phi2[5], diphoton_r9_2[5];
  float diphoton_pt[5], diphoton_mass[5], diphoton_rapidity[5], diphoton_dphi[5];
  tr->SetBranchAddress( "num_diphoton", &num_diphoton );
  tr->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  tr->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  tr->SetBranchAddress( "diphoton_phi1", diphoton_phi1 );
  tr->SetBranchAddress( "diphoton_r91", diphoton_r9_1 );
  tr->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  tr->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  tr->SetBranchAddress( "diphoton_phi2", diphoton_phi2 );
  tr->SetBranchAddress( "diphoton_r92", diphoton_r9_2 );
  tr->SetBranchAddress( "diphoton_pt", diphoton_pt );
  tr->SetBranchAddress( "diphoton_dphi", diphoton_dphi );
  tr->SetBranchAddress( "diphoton_mass", diphoton_mass );
  tr->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
  unsigned int diphoton_vertex_vtx1mmdist[5], diphoton_vertex_vtx2mmdist[5], diphoton_vertex_vtx5mmdist[5], diphoton_vertex_vtx1cmdist[5];
  float diphoton_vertex_nearestvtxdist[5];
  unsigned int diphoton_vertex_tracks[5];
  tr->SetBranchAddress( "diphoton_vertex_vtx1mmdist", diphoton_vertex_vtx1mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx2mmdist", diphoton_vertex_vtx2mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx5mmdist", diphoton_vertex_vtx5mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx1cmdist", diphoton_vertex_vtx1cmdist );
  tr->SetBranchAddress( "diphoton_vertex_nearestvtxdist", diphoton_vertex_nearestvtxdist );
  tr->SetBranchAddress( "diphoton_vertex_tracks", diphoton_vertex_tracks );
  // proton quantities
  unsigned int num_proton;
  unsigned int proton_side[10], proton_pot[10];
  float proton_xi[10], proton_xi_error[10];
  tr->SetBranchAddress( "num_proton_track", &num_proton );
  tr->SetBranchAddress( "proton_track_side", proton_side );
  tr->SetBranchAddress( "proton_track_pot", proton_pot );
  tr->SetBranchAddress( "proton_track_xi", proton_xi );
  tr->SetBranchAddress( "proton_track_xi_error", proton_xi_error );
  // vertex quantities
  unsigned int num_vertex;
  tr->SetBranchAddress( "num_vertex", &num_vertex );
  // other quantities
  float met, met_phi;
  tr->SetBranchAddress( "met", &met );
  tr->SetBranchAddress( "met_phi", &met_phi );

  //                JW
  // Electron Quantities
  unsigned int num_electron;
  float electron_pt[10], electron_eta[10], electron_phi[10], electron_energy[10];
  tr->SetBranchAddress( "num_electron", &num_electron );
  tr->SetBranchAddress( "electron_pt", electron_pt );
  tr->SetBranchAddress( "electron_eta", electron_eta );
  tr->SetBranchAddress( "electron_phi", electron_phi );
  tr->SetBranchAddress( "electron_energy", electron_energy );
  // Muon Quantities
  unsigned int num_muon;
  float muon_pt[10], muon_eta[10], muon_phi[10], muon_energy[10];
  tr->SetBranchAddress("num_muon",&num_muon );
  tr->SetBranchAddress("muon_pt", muon_pt );
  tr->SetBranchAddress("muon_eta",muon_eta );
  tr->SetBranchAddress("muon_phi",muon_phi );
  tr->SetBranchAddress("muon_energy", muon_energy );
  //Jet Quantities
  unsigned int num_jet;
  float jet_pt[10], jet_eta[10], jet_phi[10], jet_energy[10], jet_mass[10];
  tr->SetBranchAddress("num_jet",&num_jet );
  tr->SetBranchAddress("jet_pt", jet_pt );
  tr->SetBranchAddress("jet_eta",jet_eta );
  tr->SetBranchAddress("jet_phi",jet_phi );
  tr->SetBranchAddress("jet_energy", jet_energy );
  tr->SetBranchAddress("jet_mass", jet_mass );
  //

  TH1D* h_num_proton = new TH1D( "num_proton", "Number of protons reconstructed in event\\Events", 6, 0., 6. );
  TH1D* h_mpp_over_mgg = new TH1D( "mpp_over_mgg", "m_{pp}^{missing} / m_{#gamma#gamma} for double-tag events\\Events\\?.2f", 30, -2., 4. ),
       *h_ypp_minus_ygg = new TH1D( "ypp_minus_ygg", "y_{pp}^{missing} - y_{#gamma#gamma} for double-tag events\\Events\\?.2f", 50, -2.5, 2.5 );
  TH1D* h_met = new TH1D( "met", "Missing E_{T}\\Events\\GeV?.0f", 48, 0., 240. ),
       *h_met_1tag = (TH1D*)h_met->Clone( "met_1tag" ),
       *h_met_2tag = (TH1D*)h_met->Clone( "met_2tag" );
  TH1D* h_diphoton_pt = new TH1D( "diphoton_pt", "Diphoton p_{T}\\Events\\GeV?.0f", 40, 0., 400. ),
       *h_diphoton_pt_1tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_1tag" ),
       *h_diphoton_pt_2tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_2tag" );
  TH1D* h_diphoton_pt_zoom = new TH1D( "diphoton_pt_zoom", "Diphoton p_{T}\\Events\\GeV?.0f", 10, 0., 10. ),
       *h_diphoton_pt_zoom_1tag = (TH1D*)h_diphoton_pt_zoom->Clone( "diphoton_pt_zoom_1tag" ),
       *h_diphoton_pt_zoom_2tag = (TH1D*)h_diphoton_pt_zoom->Clone( "diphoton_pt_zoom_2tag" );
  TH1D* h_diphoton_leadpt = new TH1D( "leadphoton_pt", "Leading photon p_{T}\\Events\\GeV?.0f", 35, 50., 750. ),
       *h_diphoton_leadpt_1tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_1tag" ),
       *h_diphoton_leadpt_2tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_2tag" );
  TH1D* h_diphoton_subleadpt = new TH1D( "subleadphoton_pt", "Subleading photon p_{T}\\Events\\GeV?.0f", 35, 50., 750. ),
       *h_diphoton_subleadpt_1tag = (TH1D*)h_diphoton_subleadpt->Clone( "subleadphoton_pt_1tag" ),
       *h_diphoton_subleadpt_2tag = (TH1D*)h_diphoton_subleadpt->Clone( "subleadphoton_pt_2tag" );
  TH1D* h_diphoton_leadeta = new TH1D( "leadphoton_eta", "Leading photon #eta\\Events\\?.3f", 40, -2.5, 2.5 ),
       *h_diphoton_leadeta_1tag = (TH1D*)h_diphoton_leadeta->Clone( "leadphoton_eta_1tag" ),
       *h_diphoton_leadeta_2tag = (TH1D*)h_diphoton_leadeta->Clone( "leadphoton_eta_2tag" );
  TH1D* h_diphoton_subleadeta = new TH1D( "subleadphoton_eta", "Subleading photon #eta\\Events\\?.3f", 40, -2.5, 2.5 ),
       *h_diphoton_subleadeta_1tag = (TH1D*)h_diphoton_subleadeta->Clone( "subleadphoton_eta_1tag" ),
       *h_diphoton_subleadeta_2tag = (TH1D*)h_diphoton_subleadeta->Clone( "subleadphoton_eta_2tag" );
  TH1D* h_diphoton_dphi = new TH1D( "diphoton_dphi", "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.2f", 50, 0., 1. ),
       *h_diphoton_dphi_1tag = (TH1D*)h_diphoton_dphi->Clone( "diphoton_dphi_1tag" ),
       *h_diphoton_dphi_2tag = (TH1D*)h_diphoton_dphi->Clone( "diphoton_dphi_2tag" );
  TH1D* h_diphoton_dphi_zoom = new TH1D( "diphoton_dphi_zoom", "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.2f", 20, 0., 0.1 ),
       *h_diphoton_dphi_zoom_1tag = (TH1D*)h_diphoton_dphi_zoom->Clone( "diphoton_dphi_zoom_1tag" ),
       *h_diphoton_dphi_zoom_2tag = (TH1D*)h_diphoton_dphi_zoom->Clone( "diphoton_dphi_zoom_2tag" );
  TH1D* h_diphoton_mass = new TH1D( "diphoton_mass", "Diphoton mass\\Events\\GeV?.1f", 50, 500., 2000. ),
       *h_diphoton_mass_1tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_1tag" ),
       *h_diphoton_mass_2tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_2tag" );
  TH1D* h_diphoton_mass_withmet = new TH1D( "diphotonmet_mass", "Diphoton + missing E_{T} mass\\Events\\GeV?.0f", 50, 500., 2000. ),
       *h_diphoton_mass_withmet_1tag = (TH1D*)h_diphoton_mass_withmet->Clone( "diphotonmet_mass_1tag" ),
       *h_diphoton_mass_withmet_2tag = (TH1D*)h_diphoton_mass_withmet->Clone( "diphotonmet_mass_2tag" );
  TH1D* h_diphoton_mass_incl = new TH1D( "diphotonincl_mass", "Diphoton + other objects mass\\Events\\GeV?.0f", 50, 500., 2000. ),
       *h_diphoton_mass_incl_1tag = (TH1D*)h_diphoton_mass_incl->Clone( "diphotonmet_mass_1tag" ),
       *h_diphoton_mass_incl_2tag = (TH1D*)h_diphoton_mass_incl->Clone( "diphotonmet_mass_2tag" );
  TH1D* h_diphoton_rap = new TH1D( "diphoton_rap", "Diphoton rapidity\\Events\\?.1f", 40, -2.5, 2.5 ),
       *h_diphoton_rap_1tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_1tag" ),
       *h_diphoton_rap_2tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_2tag" );
  TH1D* h_diphoton_ntrk = new TH1D( "diphoton_ntrk", "Number of tracks on diphoton vertex\\Events\\?.1f", 20, 0., 20. ),
       *h_diphoton_ntrk_1tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_1tag" ),
       *h_diphoton_ntrk_2tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_2tag" );
  TH1D* h_num_vtx = new TH1D( "num_vtx", "Number of primary vertices in event\\Events", 40, 0., 40. ),
       *h_num_vtx_1tag = (TH1D*)h_num_vtx->Clone( "num_vtx_1tag" ),
       *h_num_vtx_2tag = (TH1D*)h_num_vtx->Clone( "num_vtx_2tag" );
  TH1D* h_num_vtx_1mm = new TH1D( "num_vtx_1mm", "Number of primary vertices near diphoton vertex\\Events", 12, 0., 12. ),
       *h_num_vtx_2mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_2mm" ),
       *h_num_vtx_5mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_5mm" ),
       *h_num_vtx_1cm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_1cm" );
  TH1D* h_diphoton_closestvtx = new TH1D( "diphoton_closestvtx", "Distance diphoton/nearest vertex\\Events\\mm?.1f", 25, 0., 2.5 ),
       *h_diphoton_closestvtx_1tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_1tag" ),
       *h_diphoton_closestvtx_2tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_2tag" );
  //  JW
  TH1D* h_lep_pt = new TH1D("lep_pt", "Lepton pT\\Events", 35, 0., 750.),
       *h_lep_pt_1tag = (TH1D*)h_lep_pt->Clone( "lep_pt_1tag" ),
       *h_lep_pt_2tag = (TH1D*)h_lep_pt->Clone( "lep_pt_2tag" );
  TH1D* h_jet_pt = new TH1D("jet_pt", "Jet pT\\Events", 35, 0., 750.),
       *h_jet_pt_1tag = (TH1D*)h_jet_pt->Clone( "jet_pt_1tag" ),
       *h_jet_pt_1tag = (TH1D*)h_jet_pt->Clone( "jet_pt_2tag" );
  //
  TH2D* h_met_vs_pt = new TH2D( "met_vs_pt", "Missing E_{T} (GeV)\\Diphoton p_{T} (GeV)", 40, 0., 400., 40, 0., 400. ),
       *h_met_vs_pt_2tag = (TH2D*)h_met_vs_pt->Clone( "met_vs_pt_2tag" ),
       *h_metx_vs_mety = new TH2D( "metx_vs_mety", "#slash{E}_{T,x} (GeV)\\#slash{E}_{T,y} (GeV)", 50, -125., 125., 50, -125., 125. ),
       *h_metx_vs_mety_2tag = (TH2D*)h_metx_vs_mety->Clone( "metx_vs_mety_2tag" );
  TH2D* h_mggmet_vs_mgg = new TH2D( "mggmet_vs_mgg", "Diphoton + #slash{E}_{T} mass (GeV)\\Diphoton mass (GeV)", 100, 500., 2000., 100, 500., 2000. );
  TGraphErrors* h_ptgg_vs_mpp = new TGraphErrors,
               *h_ptgg_vs_mgg = new TGraphErrors;
  TGraphErrors* h_mgg_vs_mpp = new TGraphErrors,
               *h_mgg_vs_mpp_candm = new TGraphErrors,
               *h_mgg_vs_mpp_candy = new TGraphErrors,
               *h_mggmet_vs_mpp = new TGraphErrors,
               *h_mggmet_vs_mpp_candm = new TGraphErrors,
               *h_mggmet_vs_mpp_candy = new TGraphErrors;
  TGraphErrors* h_ygg_vs_ypp = new TGraphErrors,
               *h_ygg_vs_ypp_candm = new TGraphErrors,
               *h_ygg_vs_ypp_candy = new TGraphErrors,
               *h_yggmet_vs_ypp = new TGraphErrors,
               *h_yggmet_vs_ypp_candm = new TGraphErrors,
               *h_yggmet_vs_ypp_candy = new TGraphErrors;
  TGraphErrors* h_xi1gg_vs_xi1pp = new TGraphErrors,
               *h_xi1gg_vs_xi1pp_candm = new TGraphErrors,
               *h_xi1gg_vs_xi1pp_candy = new TGraphErrors,
               *h_xi1ggmet_vs_xi1pp = new TGraphErrors,
               *h_xi2gg_vs_xi2pp = new TGraphErrors,
               *h_xi2gg_vs_xi2pp_candm = new TGraphErrors,
               *h_xi2gg_vs_xi2pp_candy = new TGraphErrors,
               *h_xi2ggmet_vs_xi2pp = new TGraphErrors;
  TGraphErrors* h_ximatch_45n = new TGraphErrors,
               *h_ximatch_45f = new TGraphErrors,
               *h_ximatch_56n = new TGraphErrors,
               *h_ximatch_56f = new TGraphErrors;

  ofstream events_list( "events_list.txt" );

  unsigned int num_evts_notag = 0, num_evts_with_tag = 0;
  TLorentzVector pho1, pho2, electron, muon, jet;
  // tree readout stage
  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );
    // dump the list of events in a text file

    //----- proton tracks retrieval part -----

    h_num_proton->Fill( num_proton );
    unsigned int num_1tag = 0, num_2tag = 0;
    float xi_prot1 = -1., err_xi_prot1 = -1.,
          xi_prot2 = -1., err_xi_prot2 = -1.;
    float max_diproton_mass = -1., max_diproton_mass_error = -1.,
          max_diproton_mass_rap = -999.;

    unsigned int num_diproton = 0;
    bool has_45 = false, has_56 = false;

    // xi / error
    typedef std::pair<float, float> trackparam_t;
    std::vector<trackparam_t> tracks_45_near, tracks_45_far,
                              tracks_56_near, tracks_56_far;

    // first loop to identify the RP and unflatten the collection
    for ( unsigned int j=0; j<num_proton; j++ ) {
      switch ( proton_side[j] ) {
        case 0: { // 4-5
          switch ( proton_pot[j] ) {
            case 2: { tracks_45_near.push_back( trackparam_t( proton_xi[j], proton_xi_error[j] ) ); } break; // near pot
            case 3: { tracks_45_far.push_back( trackparam_t( proton_xi[j], proton_xi_error[j] ) ); } break; // far pot
          }
          has_45 = true;
        } break;
        case 1: { // 5-6
          switch ( proton_pot[j] ) {
            case 2: { tracks_56_near.push_back( trackparam_t( proton_xi[j], proton_xi_error[j] ) ); } break; // near pot
            case 3: { tracks_56_far.push_back( trackparam_t( proton_xi[j], proton_xi_error[j] ) ); } break; // far pot
          }
          has_56 = true;
        } break;
      }
    }

    if ( has_45 ) num_1tag++;
    if ( has_56 ) num_1tag++;

    const float sqs = 13.e3;

    if ( has_45 && has_56 ) {
      for ( unsigned short j=0; j<tracks_45_near.size(); j++ ) {
        for ( unsigned short k=0; k<tracks_56_near.size(); k++ ) {
          const float diproton_mass = sqs*sqrt( tracks_45_near[j].first*tracks_56_near[k].first );
          if ( diproton_mass>max_diproton_mass ) {
            max_diproton_mass = diproton_mass;
            max_diproton_mass_error = sqs*sqs/4./diproton_mass * tracks_45_near[j].second * tracks_56_near[k].second;
            max_diproton_mass_rap = log( tracks_56_near[k].first/tracks_45_near[j].first )/2.;
          }
        }
      }
      num_2tag++;
    }

    //----- leptons+jets retrieval part -----

    //          JW
    for ( unsigned int j=0; j<num_electron; j++ ) {
      electron.SetPtEtaPhiM( electron_pt[j], electron_eta[j], electron_phi[j], 0.000510998928 );
    }
    for ( unsigned int j=0; j<num_muon; j++ ) {
      muon.SetPtEtaPhiM( muon_pt[j], muon_eta[j], muon_phi[j], .1056583715 );
    }
    for ( unsigned int j=0; j<num_jet; j++ ) {
      jet.SetPtEtaPhiM( jet_pt[j], jet_eta[j], jet_phi[j], jet_mass[j] );
    }
    //

    //----- diphotons retrieval part -----

    for ( unsigned int j=0; j<num_diphoton; j++ ) {

      const float energy_corr_pho1 = photon_rel_energy_scale( diphoton_pt1[j], diphoton_eta1[j], diphoton_r9_1[j] ) * diphoton_pt1[j],
                  energy_corr_pho2 = photon_rel_energy_scale( diphoton_pt2[j], diphoton_eta2[j], diphoton_r9_2[j] ) * diphoton_pt2[j];

      const float xi_reco1 = ( diphoton_pt1[j] * exp( -diphoton_eta1[j] ) + diphoton_pt2[j] * exp( -diphoton_eta2[j] ) )/sqrt_s,
                  xi_reco2 = ( diphoton_pt1[j] * exp(  diphoton_eta1[j] ) + diphoton_pt2[j] * exp(  diphoton_eta2[j] ) )/sqrt_s,
                  err_xi_reco1 = sqrt( pow( energy_corr_pho1, 2 ) * exp( -2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp( -2.*diphoton_eta2[j] ) )/sqrt_s,
                  err_xi_reco2 = sqrt( pow( energy_corr_pho1, 2 ) * exp(  2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp(  2.*diphoton_eta2[j] ) )/sqrt_s; //FIXME
cout << energy_corr_pho1 << "\t" << energy_corr_pho2 << "\t" << err_xi_reco1 << "\t" << err_xi_reco2 << endl;
      const float xi_reco1_withmet = xi_reco1 + met/sqrt_s,
                  xi_reco2_withmet = xi_reco2 + met/sqrt_s;

      //const float met_x = met*cos( met_phi ),
      //            met_y = met*sin( met_phi );
      //pho1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      //pho2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );
      //const TLorentzVector lv_met( Met_x, met_y, 0., met ),
      //                     dipho_met = pho1+pho2+lv_met;
      ////cout << dipho_met.M() << " <---> " << diphoton_mass[j] << endl;
      //const float diphoton_plus_met_mass = dipho_met.M(),
      //            diphoton_plus_met_rap = dipho_met.Rapidity();

      const float met_x = met*cos( met_phi ),
	met_y = met*sin( met_phi );
      pho1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );

      //                  JW                                                                                                  
      const TLorentzVector lv_met( met_x, met_y, 0., met ),
                           dipho_met = pho1+pho2+lv_met,
                           dipho_incl = pho1+pho2+electron+muon+jet+lv_met;
      //cout << dipho_met.M() << " <---> " << diphoton_mass[j] << endl;
      const float diphoton_plus_met_mass = dipho_met.M(),
                  diphoton_plus_met_rap = dipho_met.Rapidity(),
                  diphoton_incl_mass = dipho_incl.M(),
                  diphoton_incl_rap = dipho_incl.Rapidity();

      //cout << xi_reco1 << ", " << xi_reco2 << endl;
      h_diphoton_pt->Fill( diphoton_pt[j] );
      h_diphoton_pt_zoom->Fill( diphoton_pt[j] );
      h_diphoton_mass->Fill( diphoton_mass[j] );
      h_diphoton_mass_withmet->Fill( diphoton_plus_met_mass );
      h_diphoton_mass_incl->Fill( diphoton_incl_mass );
      h_diphoton_rap->Fill( diphoton_rapidity[j] );
      h_diphoton_closestvtx->Fill( diphoton_vertex_nearestvtxdist[j] );
      h_diphoton_ntrk->Fill( diphoton_vertex_tracks[j] );
      h_diphoton_dphi->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
      h_diphoton_dphi_zoom->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
      h_diphoton_leadpt->Fill( diphoton_pt1[j] );
      h_diphoton_subleadpt->Fill( diphoton_pt2[j] );
      h_diphoton_leadeta->Fill( diphoton_eta1[j] );
      h_diphoton_subleadeta->Fill( diphoton_eta2[j] );

      h_met_vs_pt->Fill( met, diphoton_pt[j] );
      h_metx_vs_mety->Fill( met_x, met_y );
      h_mggmet_vs_mgg->Fill( diphoton_plus_met_mass, diphoton_mass[j] );

      for ( unsigned short k=0; k<tracks_45_near.size(); k++ ) {
        const unsigned short id = h_ximatch_45n->GetN();
        h_ximatch_45n->SetPoint( id, tracks_45_near[k].first, xi_reco1 );
        h_ximatch_45n->SetPointError( id, tracks_45_near[k].second, err_xi_reco1 );
      }
      for ( unsigned short k=0; k<tracks_45_far.size(); k++ ) {
        const unsigned short id = h_ximatch_45f->GetN();
        h_ximatch_45f->SetPoint( id, tracks_45_far[k].first, xi_reco1 );
        h_ximatch_45f->SetPointError( id, tracks_45_far[k].second, err_xi_reco1 );
      }
      for ( unsigned short k=0; k<tracks_56_near.size(); k++ ) {
        const unsigned short id = h_ximatch_56n->GetN();
        h_ximatch_56n->SetPoint( id, tracks_56_near[k].first, xi_reco2 );
        h_ximatch_56n->SetPointError( id, tracks_56_near[k].second, err_xi_reco2 );
      }
      for ( unsigned short k=0; k<tracks_56_far.size(); k++ ) {
        const unsigned short id = h_ximatch_56f->GetN();
        h_ximatch_56f->SetPoint( id, tracks_56_far[k].first, xi_reco2 );
        h_ximatch_56f->SetPointError( id, tracks_56_far[k].second, err_xi_reco2 );
      }

      if ( num_1tag>0 ) {
        h_diphoton_pt_1tag->Fill( diphoton_pt[j] );
        h_diphoton_pt_zoom_1tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_1tag->Fill( diphoton_mass[j] );
        h_diphoton_mass_withmet_1tag->Fill( diphoton_plus_met_mass );
        h_diphoton_mass_incl_1tag->Fill( diphoton_incl_mass );
        h_diphoton_rap_1tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_1tag->Fill( diphoton_vertex_nearestvtxdist[j] );
        h_diphoton_ntrk_1tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_1tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
        h_diphoton_dphi_zoom_1tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
        h_diphoton_leadpt_1tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_1tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_1tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_1tag->Fill( diphoton_eta2[j] );
	
	//              JW
	h_lep_pt_1tag->Fill( electron_pt[j] + muon_pt[j] );
	h_jet_pt_1tag->Fill( jet_pt[j] )
	//
      }
      if ( num_2tag>0 ) {
        h_diphoton_pt_2tag->Fill( diphoton_pt[j] );
        h_diphoton_pt_zoom_2tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_2tag->Fill( diphoton_mass[j] );
        h_diphoton_mass_withmet_2tag->Fill( diphoton_plus_met_mass );
        h_diphoton_mass_incl_2tag->Fill( diphoton_incl_mass );
        h_diphoton_rap_2tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_2tag->Fill( diphoton_vertex_nearestvtxdist[j] );
        h_diphoton_ntrk_2tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_2tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
        h_diphoton_dphi_zoom_2tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
        h_diphoton_leadpt_2tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_2tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_2tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_2tag->Fill( diphoton_eta2[j] );

	//              JW                                                                                                                                      
	h_lep_pt_2tag->Fill( electron_pt[j] + muon_pt[j] );
	h_jet_pt_2tag->Fill( jet_pt[j] )
	//  

        h_xi1gg_vs_xi1pp->SetPoint( h_xi1gg_vs_xi1pp->GetN(), xi_reco1, xi_prot1 );
        h_xi2gg_vs_xi2pp->SetPoint( h_xi2gg_vs_xi2pp->GetN(), xi_reco2, xi_prot2 );

        h_xi1ggmet_vs_xi1pp->SetPoint( h_xi1ggmet_vs_xi1pp->GetN(), xi_reco1_withmet, xi_prot1 );
        h_xi2ggmet_vs_xi2pp->SetPoint( h_xi2ggmet_vs_xi2pp->GetN(), xi_reco2_withmet, xi_prot2 );

        h_met_vs_pt_2tag->Fill( met, diphoton_pt[j] );
        h_metx_vs_mety_2tag->Fill( met_x, met_y );

        events_list << run_id << ":" << lumisection << ":" << event_number << endl;
      }

      if ( num_diproton>0 ) {
        h_mgg_vs_mpp->SetPoint( h_mgg_vs_mpp->GetN(), diphoton_mass[j], max_diproton_mass );
        h_mgg_vs_mpp->SetPointError( h_mgg_vs_mpp->GetN(), 0., max_diproton_mass_error );
        h_mggmet_vs_mpp->SetPoint( h_mggmet_vs_mpp->GetN(), diphoton_plus_met_mass, max_diproton_mass );
        h_ygg_vs_ypp->SetPoint( h_ygg_vs_ypp->GetN(), diphoton_rapidity[j], max_diproton_mass_rap );
	h_yggmet_vs_ypp->SetPoint( h_yggmet_vs_ypp->GetN(), diphoton_plus_met_rap, max_diproton_mass_rap );
        h_mpp_over_mgg->Fill( max_diproton_mass/diphoton_mass[j] );
        h_ypp_minus_ygg->Fill( max_diproton_mass_rap - diphoton_rapidity[j] );
        h_ptgg_vs_mpp->SetPoint( h_ptgg_vs_mpp->GetN(), diphoton_pt[j], max_diproton_mass );
        h_ptgg_vs_mgg->SetPoint( h_ptgg_vs_mgg->GetN(), diphoton_pt[j], diphoton_mass[j] );

        if ( ( !compute_with_met && fabs( diphoton_mass[j]-max_diproton_mass )<max_diproton_mass*rel_err_xi )
          || ( compute_with_met && fabs( diphoton_plus_met_mass-max_diproton_mass )<max_diproton_mass*rel_err_xi ) ) {

          h_mgg_vs_mpp_candm->SetPoint( h_mgg_vs_mpp_candm->GetN(), diphoton_mass[j], max_diproton_mass );
          h_mggmet_vs_mpp_candm->SetPoint( h_mggmet_vs_mpp_candm->GetN(), diphoton_plus_met_mass, max_diproton_mass );
          h_ygg_vs_ypp_candm->SetPoint( h_ygg_vs_ypp_candm->GetN(), diphoton_rapidity[j], max_diproton_mass_rap );
	  h_yggmet_vs_ypp_candm->SetPoint( h_yggmet_vs_ypp_candm->GetN(), diphoton_plus_met_rap, max_diproton_mass_rap );
          h_xi1gg_vs_xi1pp_candm->SetPoint( h_xi1gg_vs_xi1pp_candm->GetN(), xi_reco1, xi_prot1 );
          h_xi2gg_vs_xi2pp_candm->SetPoint( h_xi2gg_vs_xi2pp_candm->GetN(), xi_reco2, xi_prot2 );

          if ( ( !compute_with_met && fabs( diphoton_rapidity[j]-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) )
            || ( compute_with_met && fabs( diphoton_plus_met_rap-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) ) ) {
            /*TLorentzVector p1, p2;
            p1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
            p2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );
            cout << "---------> " << (p1+p2).M() << endl;*/

            cout << " ---> event: " << run_id << ":" << lumisection << ":" << event_number << endl;
            cout << "      diphoton: pt=" << diphoton_pt[j] << ", mass=" << diphoton_mass[j] << ", rapidity=" << diphoton_rapidity[j] << ", dphi=" << diphoton_dphi[j] << endl
                 << "      diproton: mass=" << max_diproton_mass << ", rapidity=" << max_diproton_mass_rap << endl
                 << "      MET: " << met << endl;
            cout << "      single photon: pt=" << diphoton_pt1[j] << ", eta=" << diphoton_eta1[j] << ", phi=" << diphoton_phi1[j] << endl
                 << "                     pt=" << diphoton_pt2[j] << ", eta=" << diphoton_eta2[j] << ", phi=" << diphoton_phi2[j] << endl;
          }
        }
        if ( ( !compute_with_met && fabs( diphoton_rapidity[j]-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) )
          || ( compute_with_met && fabs( diphoton_plus_met_rap-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) ) ) {
          h_mgg_vs_mpp_candy->SetPoint( h_mgg_vs_mpp_candy->GetN(), diphoton_mass[j], max_diproton_mass );
          h_mggmet_vs_mpp_candy->SetPoint( h_mggmet_vs_mpp_candy->GetN(), diphoton_plus_met_mass, max_diproton_mass );
          h_ygg_vs_ypp_candy->SetPoint( h_ygg_vs_ypp_candy->GetN(), diphoton_rapidity[j], max_diproton_mass_rap );
	  h_yggmet_vs_ypp_candy->SetPoint( h_yggmet_vs_ypp_candy->GetN(), diphoton_plus_met_rap, max_diproton_mass_rap );
          h_xi1gg_vs_xi1pp_candy->SetPoint( h_xi1gg_vs_xi1pp_candy->GetN(), xi_reco1, xi_prot1 );
          h_xi2gg_vs_xi2pp_candy->SetPoint( h_xi2gg_vs_xi2pp_candy->GetN(), xi_reco2, xi_prot2 );
        }
        num_evts_with_tag++;
      }

      h_num_vtx_1mm->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm->Fill( diphoton_vertex_vtx1cmdist[j] );

      num_evts_notag++;
    }
    h_num_vtx->Fill( num_vertex );
    h_met->Fill( met );
    if ( num_1tag>0 ) {
      h_num_vtx_1tag->Fill( num_vertex );
      h_met_1tag->Fill( met );
    }
    if ( num_2tag>0 ) {
      h_num_vtx_2tag->Fill( num_vertex );
      h_met_2tag->Fill( met );
    }

  }

  // plotting stage
  cout << "events: " << num_evts_notag << ", with tag: " << num_evts_with_tag << endl;

  const float lumi_b = 5.060924481910, // fb-1
              lumi_c = 1.490748474431, // fb-1
              lumi_g = 3.742171002882; // fb-1

  const float lumi = ( file.Contains( "2016BCG" ) ) ? lumi_b+lumi_c+lumi_g : lumi_b+lumi_c;
  /*float lumi = 0.; string run_name;
  switch ( run ) {
    case 'B': lumi = lumi_b; run_name = "B"; break;
    case 'C': lumi = lumi_c; run_name = "C"; break;
    case '0': default: lumi = lumi_b+lumi_c; run_name = "BC"; break;
  }*/

  gStyle->SetOptStat( 0 );

  const string top_label_str = Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.2f fb^{-1}", lumi );
  const char* top_label = top_label_str.c_str();

  {
    plot_3hists( "diphoton_mass", top_label, h_diphoton_mass, h_diphoton_mass_1tag, h_diphoton_mass_2tag );
    plot_3hists( "diphoton_mass_withmet", top_label, h_diphoton_mass_withmet, h_diphoton_mass_withmet_1tag, h_diphoton_mass_withmet_2tag );
    plot_3hists( "diphoton_mass_incl", top_label, h_diphoton_mass_incl, h_diphoton_mass_incl_1tag, h_diphoton_mass_incl_2tag );
    plot_3hists( "diphoton_pt", top_label, h_diphoton_pt, h_diphoton_pt_1tag, h_diphoton_pt_2tag );
    plot_3hists( "diphoton_pt_zoom", top_label, h_diphoton_pt_zoom, h_diphoton_pt_zoom_1tag, h_diphoton_pt_zoom_2tag );
    plot_3hists( "diphoton_lead_pt", top_label, h_diphoton_leadpt, h_diphoton_leadpt_1tag, h_diphoton_leadpt_2tag );
    plot_3hists( "diphoton_sublead_pt", top_label, h_diphoton_subleadpt, h_diphoton_subleadpt_1tag, h_diphoton_subleadpt_2tag );

    /*h_diphoton_leadeta->GetYaxis()->SetRangeUser(0., 55.);
    h_diphoton_subleadeta->GetYaxis()->SetRangeUser(0., 55.);*/
    h_diphoton_leadeta->SetMinimum( 0. );
    h_diphoton_subleadeta->SetMinimum( 0. );

    plot_3hists( "diphoton_lead_eta", top_label, h_diphoton_leadeta, h_diphoton_leadeta_1tag, h_diphoton_leadeta_2tag );
    plot_3hists( "diphoton_sublead_eta", top_label, h_diphoton_subleadeta, h_diphoton_subleadeta_1tag, h_diphoton_subleadeta_2tag );
    plot_3hists( "diphoton_dphi", top_label, h_diphoton_dphi, h_diphoton_dphi_1tag, h_diphoton_dphi_2tag );
    plot_3hists( "diphoton_dphi_zoom", top_label, h_diphoton_dphi_zoom, h_diphoton_dphi_zoom_1tag, h_diphoton_dphi_zoom_2tag );
    plot_3hists( "diphoton_rapidity", top_label, h_diphoton_rap, h_diphoton_rap_1tag, h_diphoton_rap_2tag );
    plot_3hists( "diphoton_closest_vtx", top_label, h_diphoton_closestvtx, h_diphoton_closestvtx_1tag, h_diphoton_closestvtx_2tag );
    plot_3hists( "diphoton_vtx_numtracks", top_label, h_diphoton_ntrk, h_diphoton_ntrk_1tag, h_diphoton_ntrk_2tag );
    plot_3hists( "num_vertex", top_label, h_num_vtx, h_num_vtx_1tag, h_num_vtx_2tag );
    plot_3hists( "event_met", top_label, h_met, h_met_1tag, h_met_2tag );
    //              JW
    plot_3hists( "lep_pt", top_label, h_lep_pt, h_lep_pt_1tag, h_lep_pt_2tag );
    plot_3hists( "jet_pt", top_label, h_jet_pt, h_jet_pt_1tag, h_lep_pt_2tag );
    //

    cout << "total candidates: " << h_mgg_vs_mpp->Integral() << endl
         << " -> with mass matching: " << h_mgg_vs_mpp_candm->Integral() << endl
         << " -> with rapiditiy matching: " << h_mgg_vs_mpp_candy->Integral() << endl;

    plot_balances( "mass_balance", top_label, "Diphoton mass (GeV)\\Diproton missing mass (GeV)", h_mgg_vs_mpp, h_mgg_vs_mpp_candm, h_mgg_vs_mpp_candy, 500., 2000., rel_err_xi );
    plot_balances( "mass_balance_withmet", top_label, "Diphoton + #slash{E}_{T} mass (GeV)\\Diproton missing mass (GeV)", h_mggmet_vs_mpp, h_mggmet_vs_mpp_candm, h_mggmet_vs_mpp_candy, 500., 2000., rel_err_xi );
    plot_balances( "rapidity_balance", top_label, "Diphoton rapidity\\Diproton rapidity", h_ygg_vs_ypp, h_ygg_vs_ypp_candm, h_ygg_vs_ypp_candy, -3., 3., rel_err_xi/sqrt( 2. ), true );
    plot_balances( "rapidity_balance_withmet", top_label, "Diphoton + #slash{E}_{T} rapidity\\Diproton rapidity", h_yggmet_vs_ypp, h_yggmet_vs_ypp_candm, h_yggmet_vs_ypp_candy, -3., 3., rel_err_xi/sqrt( 2. ), true );
    plot_balances( "xi1_balance", top_label, "#xi_{1} from diphoton system\\Proton #xi_{1}", h_xi1gg_vs_xi1pp, h_xi1gg_vs_xi1pp_candm, h_xi1gg_vs_xi1pp_candy, 0., 0.45, rel_err_xi );
    plot_balances( "xi2_balance", top_label, "#xi_{2} from diphoton system\\Proton #xi_{2}", h_xi2gg_vs_xi2pp, h_xi2gg_vs_xi2pp_candm, h_xi2gg_vs_xi2pp_candy, 0., 0.45, rel_err_xi );
    plot_balances( "xi1_balance_withmet", top_label, "#xi_{1} from diphoton + #slash{E}_{T} system\\Proton #xi_{1}", h_xi1ggmet_vs_xi1pp, 0, 0, 0., 0.45, rel_err_xi );
    plot_balances( "xi2_balance_withmet", top_label, "#xi_{2} from diphoton + #slash{E}_{T} system\\Proton #xi_{2}", h_xi2ggmet_vs_xi2pp, 0, 0, 0., 0.45, rel_err_xi );

    plot_balances_new( "xi1_balance_45n", top_label, "Proton #xi (sector 45 - near pot)\\Diphoton #xi", h_ximatch_45n, 0, 0, 0., 0.45 );
    plot_balances_new( "xi1_balance_45f", top_label, "Proton #xi (sector 45 - far pot)\\Diphoton #xi", h_ximatch_45f, 0, 0, 0., 0.45 );
    plot_balances_new( "xi1_balance_56n", top_label, "Proton #xi (sector 56 - near pot)\\Diphoton #xi", h_ximatch_56n, 0, 0, 0., 0.45 );
    plot_balances_new( "xi1_balance_56f", top_label, "Proton #xi (sector 56 - far pot)\\Diphoton #xi", h_ximatch_56f, 0, 0, 0., 0.45 );

    plot_balances( "diphoton_pt_vs_diproton_mass", top_label, "Diphoton p_{T} (GeV)\\Diproton mass (GeV)", h_ptgg_vs_mpp );
    plot_balances( "diphoton_pt_vs_diphoton_mass", top_label, "Diphoton p_{T} (GeV)\\Diphoton mass (GeV)", h_ptgg_vs_mgg );
  }

  {
    Canvas c( "met_x_vs_y", top_label );
    h_metx_vs_mety->Draw( "colz" );
    h_metx_vs_mety_2tag->Draw( "p same" );
    h_metx_vs_mety_2tag->SetMarkerStyle( 24 );
    h_metx_vs_mety_2tag->SetMarkerColor( kRed );
    c.SetLegendX1( 0.15 );
    c.SetLegendY1( 0.15 );
    c.AddLegendEntry( h_metx_vs_mety_2tag, "#geq 2 proton tags", "p" );
    c.Prettify( h_metx_vs_mety );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "diphoton_mass_vs diphotonmet_mass", top_label );
    h_mggmet_vs_mgg->Draw( "colz" );
    c.DrawDiagonal( h_mggmet_vs_mgg->GetXaxis()->GetXmin(), h_mggmet_vs_mgg->GetXaxis()->GetXmax() );
    h_mggmet_vs_mgg->Draw( "colz same" );
    c.Prettify( h_mggmet_vs_mgg );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "diphoton_pt_vs_met", top_label );
    h_met_vs_pt->Draw( "colz" );
    h_met_vs_pt_2tag->Draw( "p same" );
    h_met_vs_pt_2tag->SetMarkerStyle( 24 );
    h_met_vs_pt_2tag->SetMarkerColor( kRed );
    c.AddLegendEntry( h_met_vs_pt_2tag, "#geq 2 proton tags", "p" );
    c.DrawDiagonal( h_met_vs_pt->GetXaxis()->GetXmin(), h_met_vs_pt->GetXaxis()->GetXmax() );
    c.Prettify( h_met_vs_pt );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  /*{
    Canvas c( "diphoton_pt_vs_diproton_mass", top_label );
    h_ptgg_vs_mpp->Draw( "colz" );
    c.DrawDiagonal( h_ptgg_vs_mpp->GetXaxis()->GetXmin(), h_ptgg_vs_mpp->GetXaxis()->GetXmax() );
    c.Prettify( h_ptgg_vs_mpp );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }*/
  {
    Canvas c( "num_close_vertex", top_label );
    h_num_vtx_1mm->Sumw2();
    h_num_vtx_2mm->Sumw2();
    h_num_vtx_5mm->Sumw2();
    h_num_vtx_1cm->Sumw2();
    h_num_vtx_1mm->Draw();
    h_num_vtx_1mm->SetLineColor( kBlack );
    h_num_vtx_2mm->Draw( "same" );
    h_num_vtx_2mm->SetLineColor( kBlack );
    h_num_vtx_2mm->SetMarkerColor( kRed+1 );
    h_num_vtx_5mm->Draw( "same" );
    h_num_vtx_5mm->SetMarkerColor( kGreen+2 );
    h_num_vtx_5mm->SetLineColor( kBlack );
    h_num_vtx_1cm->Draw( "same" );
    h_num_vtx_1cm->SetMarkerColor( kBlue+1 );
    h_num_vtx_1cm->SetLineColor( kBlack );
    h_num_vtx_1mm->SetMarkerStyle( 20 );
    h_num_vtx_2mm->SetMarkerStyle( 21 );
    h_num_vtx_5mm->SetMarkerStyle( 22 );
    h_num_vtx_1cm->SetMarkerStyle( 23 );
    c.AddLegendEntry( h_num_vtx_1mm, "at 1 mm distance" );
    c.AddLegendEntry( h_num_vtx_2mm, "at 2 mm distance" );
    c.AddLegendEntry( h_num_vtx_5mm, "at 5 mm distance" );
    c.AddLegendEntry( h_num_vtx_1cm, "at 1 cm distance" );
    c.Prettify( h_num_vtx_1mm );
    c.SetLogy();
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "num_proton", top_label );
    h_num_proton->Sumw2();
    h_num_proton->Draw();
    h_num_proton->SetMarkerStyle( 20 );
    h_num_proton->SetLineColor( kBlack );
    c.Prettify( h_num_proton );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "mass_ratio", top_label );
    h_mpp_over_mgg->Sumw2();
    h_mpp_over_mgg->Draw();
    h_mpp_over_mgg->SetMarkerStyle( 20 );
    h_mpp_over_mgg->SetLineColor( kBlack );
    c.Prettify( h_mpp_over_mgg );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "rapidity_difference", top_label );
    h_ypp_minus_ygg->Sumw2();
    h_ypp_minus_ygg->Draw();
    h_ypp_minus_ygg->SetMarkerStyle( 20 );
    h_ypp_minus_ygg->SetLineColor( kBlack );
    c.Prettify( h_ypp_minus_ygg );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }

}

//  LocalWords:  SetMarkerStyle
