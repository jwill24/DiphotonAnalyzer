#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStyle.h"

#include "Canvas.h"
#include "Plotter.h"

#include <fstream>
#include <iostream>

#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton"
//#define default_ntp_file "output_Run2016B_mgg-ov-500.root"
//#define default_ntp_file "output_Run2016BC_looseCuts.root"
#define default_ntp_file "output_Run2016BCG_looseCuts_prelim.root"

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 );

enum Sector { sector45 = 0, sector56 = 1 };
enum Pot { nearPot = 2, farPot = 3 };

void tree_reader( TString file=default_ntp_file )
{
  TFile f( file );
  if ( !f.IsOpen() ) return;

  const float sqrt_s = 13.e3;
  const bool compute_with_met = true;

  //----- cuts -----

  const float max_dist_vtx_inclobjects = 0.2, // in cm
              min_pt_photon = 50.,
              min_r9_photon = 0.94,
              max_eta_photon = 2.5,
              min_etaveto_photon = 1.4442,
              max_etaveto_photon = 1.566,
              min_mass_diphoton = 300.;
  /*const float max_pt_diphoton = 500.,
              max_acopl = 0.01;*/
  const float max_pt_diphoton = 10000.,
              max_acopl = 0.1;
  const float oth_obj_minpt = 150.;
  const float n_sigma = 1.0; // xi matching sigma number

  /*210-N (beam 1): 0.0334862
    210-N (beam 2): 0.0496449
    210-F (beam 1): 0.023913
    210-F (beam 2): 0.0368263*/
  const float lim_45n = 0.034,
              lim_45f = 0.024,
              lim_56n = 0.050,
              lim_56f = 0.037;

  TH1D* h_pt_incl = new TH1D( "pt_incl", "", 100, 0., 2000. );

  //----- tree loading -----

  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );
  // general quantities
  unsigned int run_id, lumisection;
  unsigned long long event_number;
  tr->SetBranchAddress( "run_id", &run_id );
  tr->SetBranchAddress( "lumisection", &lumisection );
  tr->SetBranchAddress( "event_number", &event_number );
  // diphoton quantities
  unsigned int num_diphoton;
  const unsigned short max_diph = 10;
  float diphoton_pt1[max_diph], diphoton_eta1[max_diph], diphoton_phi1[max_diph], diphoton_r9_1[max_diph],
        diphoton_pt2[max_diph], diphoton_eta2[max_diph], diphoton_phi2[max_diph], diphoton_r9_2[max_diph];
  float diphoton_pt[max_diph], diphoton_mass[max_diph], diphoton_rapidity[max_diph], diphoton_dphi[max_diph];
  float diphoton_vertex_x[max_diph], diphoton_vertex_y[max_diph], diphoton_vertex_z[max_diph];
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
  tr->SetBranchAddress( "diphoton_vertex_x", diphoton_vertex_x );
  tr->SetBranchAddress( "diphoton_vertex_y", diphoton_vertex_y );
  tr->SetBranchAddress( "diphoton_vertex_z", diphoton_vertex_z );
  tr->SetBranchAddress( "diphoton_vertex_vtx1mmdist", diphoton_vertex_vtx1mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx2mmdist", diphoton_vertex_vtx2mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx5mmdist", diphoton_vertex_vtx5mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx1cmdist", diphoton_vertex_vtx1cmdist );
  tr->SetBranchAddress( "diphoton_vertex_nearestvtxdist", diphoton_vertex_nearestvtxdist );
  tr->SetBranchAddress( "diphoton_vertex_tracks", diphoton_vertex_tracks );
  // proton quantities
  unsigned int num_proton;
  const unsigned short max_pr = 20;
  unsigned int proton_side[max_pr], proton_pot[max_pr];
  float proton_x[max_pr], proton_y[max_pr];
  float proton_xi[max_pr], proton_xi_error[max_pr];
  unsigned int proton_link_id[max_pr];
  float proton_link_dist[max_pr];
  tr->SetBranchAddress( "num_proton_track", &num_proton );
  tr->SetBranchAddress( "proton_track_side", proton_side );
  tr->SetBranchAddress( "proton_track_pot", proton_pot );
  tr->SetBranchAddress( "proton_track_x", proton_x );
  tr->SetBranchAddress( "proton_track_y", proton_y );
  tr->SetBranchAddress( "proton_track_xi", proton_xi );
  tr->SetBranchAddress( "proton_track_xi_error", proton_xi_error );
  tr->SetBranchAddress( "proton_track_link_nearfar", proton_link_id );
  tr->SetBranchAddress( "proton_track_link_mindist", proton_link_dist );
  // vertex quantities
  unsigned int num_vertex;
  const unsigned short max_vtx = 80;
  float vertex_x[max_vtx], vertex_y[max_vtx], vertex_z[max_vtx];
  tr->SetBranchAddress( "num_vertex", &num_vertex );
  tr->SetBranchAddress( "vertex_x", vertex_x );
  tr->SetBranchAddress( "vertex_y", vertex_y );
  tr->SetBranchAddress( "vertex_z", vertex_z );

  // other quantities
  float met, met_phi;
  tr->SetBranchAddress( "met", &met );
  tr->SetBranchAddress( "met_phi", &met_phi );

  //                JW
  // Electron Quantities
  unsigned int num_electron;
  const unsigned short max_ele = 20;
  float electron_pt[max_ele], electron_eta[max_ele], electron_phi[max_ele], electron_energy[max_ele];
  float electron_vtx_x[max_ele], electron_vtx_y[max_ele], electron_vtx_z[max_ele];
  tr->SetBranchAddress( "num_electron", &num_electron );
  tr->SetBranchAddress( "electron_pt", electron_pt );
  tr->SetBranchAddress( "electron_eta", electron_eta );
  tr->SetBranchAddress( "electron_phi", electron_phi );
  tr->SetBranchAddress( "electron_energy", electron_energy );
  tr->SetBranchAddress( "electron_vtx_x", electron_vtx_x );
  tr->SetBranchAddress( "electron_vtx_y", electron_vtx_y );
  tr->SetBranchAddress( "electron_vtx_z", electron_vtx_z );
  // Muon Quantities
  unsigned int num_muon;
  const unsigned short max_mu = 20;
  float muon_pt[max_mu], muon_eta[max_mu], muon_phi[max_mu], muon_energy[max_mu];
  float muon_vtx_x[max_mu], muon_vtx_y[max_mu], muon_vtx_z[max_mu];
  tr->SetBranchAddress( "num_muon",&num_muon );
  tr->SetBranchAddress( "muon_pt", muon_pt );
  tr->SetBranchAddress( "muon_eta",muon_eta );
  tr->SetBranchAddress( "muon_phi",muon_phi );
  tr->SetBranchAddress( "muon_energy", muon_energy );
  tr->SetBranchAddress( "muon_vtx_x", muon_vtx_x );
  tr->SetBranchAddress( "muon_vtx_y", muon_vtx_y );
  tr->SetBranchAddress( "muon_vtx_z", muon_vtx_z );
  //Jet Quantities
  unsigned int num_jet;
  const unsigned short max_jet = 50;
  float jet_pt[max_jet], jet_eta[max_jet], jet_phi[max_jet], jet_energy[max_jet];
  float jet_vtx_x[max_jet], jet_vtx_y[max_jet], jet_vtx_z[max_jet];
  tr->SetBranchAddress( "num_jet", &num_jet );
  tr->SetBranchAddress( "jet_pt",  jet_pt );
  tr->SetBranchAddress( "jet_eta", jet_eta );
  tr->SetBranchAddress( "jet_phi", jet_phi );
  tr->SetBranchAddress( "jet_energy", jet_energy );
  tr->SetBranchAddress( "jet_vtx_x", jet_vtx_x );
  tr->SetBranchAddress( "jet_vtx_y", jet_vtx_y );
  tr->SetBranchAddress( "jet_vtx_z", jet_vtx_z );
  //

  TH1D* h_mpp_over_mgg = new TH1D( "mpp_over_mgg", "m_{pp}^{missing} / m_{#gamma#gamma} for double-tag events\\Events\\?.2f", 30, -2., 4. ),
       *h_ypp_minus_ygg = new TH1D( "ypp_minus_ygg", "y_{pp}^{missing} - y_{#gamma#gamma} for double-tag events\\Events\\?.2f", 50, -2.5, 2.5 );
  TH1D* h_met = new TH1D( "met", "Missing E_{T}\\Events\\GeV?.0f", 35, 0., 3500. ),
       *h_met_1tag = (TH1D*)h_met->Clone( "met_1tag" ),
       *h_met_2tag = (TH1D*)h_met->Clone( "met_2tag" );
  // diphoton only
  TH1D* h_num_diphoton = new TH1D( "num_diphoton", "Diphoton multiplicity\\Event", 10, 0., 10. ),
       *h_num_diphoton_1tag = (TH1D*)h_num_diphoton->Clone( "num_diphoton_1tag" ),
       *h_num_diphoton_2tag = (TH1D*)h_num_diphoton->Clone( "num_diphoton_2tag" );
  TH1D* h_diphoton_pt = new TH1D( "diphoton_pt", "Diphoton p_{T}\\Events\\GeV?.0f", 40, 0., 400. ),
       *h_diphoton_pt_1tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_1tag" ),
       *h_diphoton_pt_2tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_2tag" );
  TH1D* h_diphoton_pt_zoom = new TH1D( "diphoton_pt_zoom", "Diphoton p_{T}\\Events\\GeV?.0f", 10, 0., 10. ),
       *h_diphoton_pt_zoom_1tag = (TH1D*)h_diphoton_pt_zoom->Clone( "diphoton_pt_zoom_1tag" ),
       *h_diphoton_pt_zoom_2tag = (TH1D*)h_diphoton_pt_zoom->Clone( "diphoton_pt_zoom_2tag" );
  TH1D* h_diphoton_leadpt = new TH1D( "leadphoton_pt", "Leading photon p_{T}\\Events\\GeV?.0f", ( 750.-min_pt_photon )/20, min_pt_photon, 750. ),
       *h_diphoton_leadpt_1tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_1tag" ),
       *h_diphoton_leadpt_2tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_2tag" );
  TH1D* h_diphoton_subleadpt = new TH1D( "subleadphoton_pt", "Subleading photon p_{T}\\Events\\GeV?.0f", ( 750.-min_pt_photon )/20, min_pt_photon, 750. ),
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
  TH1D* h_diphoton_mass = new TH1D( "diphoton_mass", "Diphoton mass\\Events\\GeV?.0f", ( int )( ( 2000.-min_mass_diphoton )/25 ), min_mass_diphoton, 2000. ),
       *h_diphoton_mass_1tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_1tag" ),
       *h_diphoton_mass_2tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_2tag" );
  TH1D* h_diphoton_mass_withmet = new TH1D( "diphotonmet_mass", "Diphoton + missing E_{T} mass\\Events\\GeV?.0f", ( int )( ( 3000.-min_mass_diphoton )/40 ), min_mass_diphoton, 3000. ),
       *h_diphoton_mass_withmet_1tag = (TH1D*)h_diphoton_mass_withmet->Clone( "diphotonmet_mass_1tag" ),
       *h_diphoton_mass_withmet_2tag = (TH1D*)h_diphoton_mass_withmet->Clone( "diphotonmet_mass_2tag" );
  TH1D* h_diphoton_mass_incl = new TH1D( "diphotonincl_mass", "Diphoton + other objects mass\\Events\\GeV?.0f", ( int )( ( 3500.-min_mass_diphoton )/50 ), min_mass_diphoton, 3500. ),
       *h_diphoton_mass_incl_1tag = (TH1D*)h_diphoton_mass_incl->Clone( "diphotonmet_mass_1tag" ),
       *h_diphoton_mass_incl_2tag = (TH1D*)h_diphoton_mass_incl->Clone( "diphotonmet_mass_2tag" );
  TH1D* h_diphoton_rap = new TH1D( "diphoton_rap", "Diphoton rapidity\\Events\\?.1f", 40, -2.5, 2.5 ),
       *h_diphoton_rap_1tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_1tag" ),
       *h_diphoton_rap_2tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_2tag" );
  TH1D* h_diphoton_ntrk = new TH1D( "diphoton_ntrk", "Number of tracks on diphoton vertex\\Events\\?.1f", 20, 0., 20. ),
       *h_diphoton_ntrk_1tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_1tag" ),
       *h_diphoton_ntrk_2tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_2tag" );
  TH1D* h_diphoton_vtxz = new TH1D( "diphoton_vtxz", "Diphoton longitudinal vertex position\\Events\\cm?.1f", 40, -10., 10. ),
       *h_diphoton_vtxz_1tag = (TH1D*)h_diphoton_vtxz->Clone( "diphoton_vtxz_1tag" ),
       *h_diphoton_vtxz_2tag = (TH1D*)h_diphoton_vtxz->Clone( "diphoton_vtxz_2tag" );
  // vertexing study
  TH1D* h_num_vtx = new TH1D( "num_vtx", "Number of primary vertices in event\\Events", 40, 0., 40. ),
       *h_num_vtx_1tag = (TH1D*)h_num_vtx->Clone( "num_vtx_1tag" ),
       *h_num_vtx_2tag = (TH1D*)h_num_vtx->Clone( "num_vtx_2tag" );
  TH1D* h_vtxz = new TH1D( "vtxz", "z position, all primary vertices in event\\Events\\cm", 40, -20., 20. ),
       *h_vtxz_1tag = (TH1D*)h_vtxz->Clone( "vtxz_1tag" ),
       *h_vtxz_2tag = (TH1D*)h_vtxz->Clone( "vtxz_2tag" );
  TH1D* h_num_vtx_1mm = new TH1D( "num_vtx_1mm", "Number of primary vertices near diphoton vertex\\Events", 12, 0., 12. ),
       *h_num_vtx_2mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_2mm" ),
       *h_num_vtx_5mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_5mm" ),
       *h_num_vtx_1cm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_1cm" ),
       *h_num_vtx_1mm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_1mm" ),
       *h_num_vtx_2mm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_2mm" ),
       *h_num_vtx_5mm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_5mm" ),
       *h_num_vtx_1cm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_1cm" );
  TH1D* h_diphoton_closestvtx = new TH1D( "diphoton_closestvtx", "Distance diphoton/nearest vertex\\Events\\cm?.1f", 25, 0., 5.0 ),
       *h_diphoton_closestvtx_1tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_1tag" ),
       *h_diphoton_closestvtx_2tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_2tag" );
  TH1D* h_electron_pt_1mm = new TH1D( "electron_pt_1mm", "Close electron E_{T}\\Events\\GeV", 25, 0., 500. ),
       *h_electron_pt_2mm = (TH1D*)h_electron_pt_1mm->Clone( "electron_pt_2mm" ),
       *h_electron_pt_5mm = (TH1D*)h_electron_pt_1mm->Clone( "electron_pt_5mm" ),
       *h_electron_pt_1cm = (TH1D*)h_electron_pt_1mm->Clone( "electron_pt_1cm" );
  TH1D* h_muon_pt_1mm = new TH1D( "muon_pt_1mm", "Close muon p_{T}\\Events\\GeV", 25, 0., 250. ),
       *h_muon_pt_2mm = (TH1D*)h_muon_pt_1mm->Clone( "muon_pt_2mm" ),
       *h_muon_pt_5mm = (TH1D*)h_muon_pt_1mm->Clone( "muon_pt_5mm" ),
       *h_muon_pt_1cm = (TH1D*)h_muon_pt_1mm->Clone( "muon_pt_1cm" );
  TH1D* h_jet_pt_1mm = new TH1D( "jet_pt_1mm", "Close jet p_{T}\\Events\\GeV", 50, 0., 500. ),
       *h_jet_pt_2mm = (TH1D*)h_jet_pt_1mm->Clone( "jet_pt_2mm" ),
       *h_jet_pt_5mm = (TH1D*)h_jet_pt_1mm->Clone( "jet_pt_5mm" ),
       *h_jet_pt_1cm = (TH1D*)h_jet_pt_1mm->Clone( "jet_pt_1cm" );
  //  JW
  TH1D* h_lep_pt = new TH1D("lep_pt", "Lepton pT\\Events", 35, 0., 750.),
       *h_lep_pt_1tag = (TH1D*)h_lep_pt->Clone( "lep_pt_1tag" ),
       *h_lep_pt_2tag = (TH1D*)h_lep_pt->Clone( "lep_pt_2tag" );
  TH1D* h_jet_pt = new TH1D("jet_pt", "Jet pT\\Events", 35, 0., 750.),
       *h_jet_pt_1tag = (TH1D*)h_jet_pt->Clone( "jet_pt_1tag" ),
       *h_jet_pt_2tag = (TH1D*)h_jet_pt->Clone( "jet_pt_2tag" );
  //
  TH2D* h_met_vs_pt = new TH2D( "met_vs_pt", "Missing E_{T} (GeV)\\Diphoton p_{T} (GeV)", 40, 0., 400., 40, 0., 400. ),
       *h_met_vs_pt_2tag = (TH2D*)h_met_vs_pt->Clone( "met_vs_pt_2tag" ),
       *h_metx_vs_mety = new TH2D( "metx_vs_mety", "#slash{E}_{T,x} (GeV)\\#slash{E}_{T,y} (GeV)", 50, -4500., 4500., 50, -4500., 4500. ),
       *h_metx_vs_mety_1tag = (TH2D*)h_metx_vs_mety->Clone( "metx_vs_mety_1tag" ),
       *h_metx_vs_mety_2tag = (TH2D*)h_metx_vs_mety->Clone( "metx_vs_mety_2tag" );
  TH2D* h_mggmet_vs_mgg = new TH2D( "mggmet_vs_mgg", "Diphoton + #slash{E}_{T} mass (GeV)\\Diphoton mass (GeV)", ( int )( ( 2000.-min_mass_diphoton )/50 ), min_mass_diphoton, 2000., ( int )( ( 2000.-min_mass_diphoton )/50 ), min_mass_diphoton, 2000. );
  TH2D* h_hitmap_45n = new TH2D( "hitmap_45n", "RP track x (cm)\\RP track y (cm)", 50, 0., 5., 50, -2.5, 2.5 ),
       *h_hitmap_45f = (TH2D*)h_hitmap_45n->Clone( "hitmap_45f" ),
       *h_hitmap_56n = (TH2D*)h_hitmap_45n->Clone( "hitmap_56n" ),
       *h_hitmap_56f = (TH2D*)h_hitmap_45n->Clone( "hitmap_56f" ),
       *h_hitmap_45n_2tag = (TH2D*)h_hitmap_45n->Clone( "hitmap_45n_2tag" ),
       *h_hitmap_45f_2tag = (TH2D*)h_hitmap_45n->Clone( "hitmap_45f_2tag" ),
       *h_hitmap_56n_2tag = (TH2D*)h_hitmap_45n->Clone( "hitmap_56n_2tag" ),
       *h_hitmap_56f_2tag = (TH2D*)h_hitmap_45n->Clone( "hitmap_56f_2tag" );
  TH1D* h_xidist_45n_all = new TH1D( "xidist_45n_all", "#xi^{RP}\\Events", 100, 0., 0.3 ),
       *h_xidist_45f_all = (TH1D*)h_xidist_45n_all->Clone( "xidist_45f_all" ),
       *h_xidist_56n_all = (TH1D*)h_xidist_45n_all->Clone( "xidist_56n_all" ),
       *h_xidist_56f_all = (TH1D*)h_xidist_45n_all->Clone( "xidist_56f_all" ),
       *h_xidist_45n_2tag = (TH1D*)h_xidist_45n_all->Clone( "xidist_45n_2tag" ),
       *h_xidist_45f_2tag = (TH1D*)h_xidist_45n_all->Clone( "xidist_45f_2tag" ),
       *h_xidist_56n_2tag = (TH1D*)h_xidist_45n_all->Clone( "xidist_56n_2tag" ),
       *h_xidist_56f_2tag = (TH1D*)h_xidist_45n_all->Clone( "xidist_56f_2tag" );
  // balances
  TGraphErrors* h_ptgg_vs_mpp = new TGraphErrors, *h_ptgg_vs_mgg = new TGraphErrors;
  TGraphErrors* h_mgg_vs_mpp = new TGraphErrors, *h_mgg_vs_mpp_candm = new TGraphErrors, *h_mgg_vs_mpp_candy = new TGraphErrors,
               *h_mggmet_vs_mpp = new TGraphErrors, *h_mggmet_vs_mpp_candm = new TGraphErrors, *h_mggmet_vs_mpp_candy = new TGraphErrors,
               *h_mggincl_vs_mpp = new TGraphErrors, *h_mggincl_vs_mpp_candm = new TGraphErrors, *h_mggincl_vs_mpp_candy = new TGraphErrors;
  TGraphErrors* h_ygg_vs_ypp = new TGraphErrors, *h_ygg_vs_ypp_candm = new TGraphErrors, *h_ygg_vs_ypp_candy = new TGraphErrors,
               *h_yggmet_vs_ypp = new TGraphErrors, *h_yggmet_vs_ypp_candm = new TGraphErrors, *h_yggmet_vs_ypp_candy = new TGraphErrors,
               *h_yggincl_vs_ypp = new TGraphErrors, *h_yggincl_vs_ypp_candm = new TGraphErrors, *h_yggincl_vs_ypp_candy = new TGraphErrors;
  TGraphErrors* h_ximatch_45n = new TGraphErrors, *h_ximatch_45f = new TGraphErrors, *h_ximatch_56n = new TGraphErrors, *h_ximatch_56f = new TGraphErrors,
               *h_ximatch_45n_matched = new TGraphErrors, *h_ximatch_45f_matched = new TGraphErrors, *h_ximatch_56n_matched = new TGraphErrors, *h_ximatch_56f_matched = new TGraphErrors,
               *h_ximatch_45n_withmet = new TGraphErrors, *h_ximatch_45f_withmet = new TGraphErrors, *h_ximatch_56n_withmet = new TGraphErrors, *h_ximatch_56f_withmet = new TGraphErrors;
  // proton reco study
  TH1D* h_num_proton = new TH1D( "num_proton", "Forward track multiplicity\\Events", 6, 0., 6. ),
       *h_num_proton_45 = (TH1D*)h_num_proton->Clone( "num_proton_45" ),
       *h_num_proton_56 = (TH1D*)h_num_proton->Clone( "num_proton_56" );

  const unsigned short num_bins = 50;
  const float max_bin_dist = 5.;
  float dist_bin[num_bins];
  unsigned short num_ele_bin[num_bins], num_ele_ptcut_bin[num_bins],
                 num_mu_bin[num_bins], num_mu_ptcut_bin[num_bins],
                 num_jet_bin[num_bins], num_jet_ptcut_bin[num_bins];
  for ( unsigned short i=0; i<num_bins; i++ ) {
    dist_bin[i] = max_bin_dist*i/( num_bins+1 );
  }

  TH1D* gr_veto_ele = new TH1D( "gr_veto_ele", "Other objects multiplicity at distance (cm)\\Diphotons", num_bins-1, dist_bin ),
       *gr_veto_mu = ( TH1D* )gr_veto_ele->Clone( "gr_veto_mu" ),
       *gr_veto_jet = ( TH1D* )gr_veto_ele->Clone( "gr_veto_jet" ),
       *gr_veto_ele_ptcut = ( TH1D* )gr_veto_ele->Clone( "gr_veto_ele_ptcut" ),
       *gr_veto_mu_ptcut = ( TH1D* )gr_veto_ele->Clone( "gr_veto_mu_ptcut" ),
       *gr_veto_jet_ptcut = ( TH1D* )gr_veto_ele->Clone( "gr_veto_jet_ptcut" );

  ofstream events_list( "events_list.txt" );

  unsigned int num_evts_notag = 0, num_evts_with_tag = 0;
  TLorentzVector pho1, pho2, electron, muon, jet;


  unsigned short cand_1tag_id = 0, cand_2tag_id = 0;

  // tree readout stage
  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );

    if ( run_id>280385 ) { continue; } // skip post-TS2 runs

    //----- diphotons retrieval part -----

    unsigned short num_diphoton_cand = 0, num_diphoton_cand_1tag = 0, num_diphoton_cand_2tag = 0;
    unsigned short num_proton_45 = 0, num_proton_56 = 0;

    // first fill the proton hitmaps outside the diphoton loop
    for ( unsigned int j=0; j<num_proton; j++ ) {
      if ( proton_side[j]==0 && proton_pot[j]==2 ) {
        h_hitmap_45n->Fill( proton_x[j]*100., proton_y[j]*100. );
        h_xidist_45n_all->Fill( proton_xi[j] );
        num_proton_45++;
      }
      if ( proton_side[j]==0 && proton_pot[j]==3 ) {
        h_hitmap_45f->Fill( proton_x[j]*100., proton_y[j]*100. );
        h_xidist_45f_all->Fill( proton_xi[j] );
        num_proton_45++;
      }
      if ( proton_side[j]==1 && proton_pot[j]==2 ) {
        h_hitmap_56n->Fill( proton_x[j]*100., proton_y[j]*100. );
        h_xidist_56n_all->Fill( proton_xi[j] );
        num_proton_56++;
      }
      if ( proton_side[j]==1 && proton_pot[j]==3 ) {
        h_hitmap_56f->Fill( proton_x[j]*100., proton_y[j]*100. );
        h_xidist_56f_all->Fill( proton_xi[j] );
        num_proton_56++;
      }
    }

    for ( unsigned short i=0; i<num_bins; i++ ) {
      num_ele_bin[i] = num_mu_bin[i] = num_jet_bin[i] = 0;
      num_ele_ptcut_bin[i] = num_mu_ptcut_bin[i] = num_jet_ptcut_bin[i] = 0;
    }

    //----- BEFORE THE LOOP ON DIPHOTON -----

    for ( unsigned int j=0; j<num_diphoton; j++ ) {

      const float acopl = 1-fabs( diphoton_dphi[j]/TMath::Pi() );

      const float energy_corr_pho1 = photon_rel_energy_scale( diphoton_pt1[j], diphoton_eta1[j], diphoton_r9_1[j] ) * diphoton_pt1[j],
                  energy_corr_pho2 = photon_rel_energy_scale( diphoton_pt2[j], diphoton_eta2[j], diphoton_r9_2[j] ) * diphoton_pt2[j];

      const float xim_reco = ( diphoton_pt1[j] * exp( -diphoton_eta1[j] ) + diphoton_pt2[j] * exp( -diphoton_eta2[j] ) )/sqrt_s,
                  xip_reco = ( diphoton_pt1[j] * exp(  diphoton_eta1[j] ) + diphoton_pt2[j] * exp(  diphoton_eta2[j] ) )/sqrt_s,
                  err_xim_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp( -2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp( -2.*diphoton_eta2[j] ) )/sqrt_s,
                  err_xip_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp(  2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp(  2.*diphoton_eta2[j] ) )/sqrt_s; //FIXME

      const float xim_reco_withmet = xim_reco + met/sqrt_s, err_xim_reco_withmet = 0.,
                  xip_reco_withmet = xip_reco + met/sqrt_s, err_xip_reco_withmet = 0.;

      //----- PRESELECTION ------

      //----- quality cuts for the photons -----

      if ( diphoton_pt1[j]<min_pt_photon || diphoton_pt2[j]<min_pt_photon ) continue;
      if ( fabs( diphoton_eta1[j] )>max_eta_photon || fabs( diphoton_eta2[j] )>max_eta_photon ) continue;
      if ( ( fabs( diphoton_eta1[j] )>=min_etaveto_photon && fabs( diphoton_eta1[j] )<=max_etaveto_photon ) ||
           ( fabs( diphoton_eta2[j] )>=min_etaveto_photon && fabs( diphoton_eta2[j] )<=max_etaveto_photon ) ) continue;
      if ( diphoton_r9_1[j]<min_r9_photon || diphoton_r9_2[j]<min_r9_photon ) continue;
      if ( diphoton_mass[j]<min_mass_diphoton ) continue;

      //----- build the event's kinematics -----

      pho1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );

      const TVector3 diph_vtx( diphoton_vertex_x[j], diphoton_vertex_y[j], diphoton_vertex_z[j] );
      float min_dist_vtx = 999.999;
      for ( unsigned int k=0; k<num_vertex; k++ ) {
        const TVector3 oth_vtx( vertex_x[k], vertex_y[k], vertex_z[k] );
        float dist = ( oth_vtx-diph_vtx ).Mag();
        if ( dist>0. && dist<min_dist_vtx ) {
          min_dist_vtx = dist;
        }
      }
      h_diphoton_closestvtx->Fill( min_dist_vtx );

      //----- leptons+jets retrieval part -----

      bool has_ele = false, has_muon = false, has_jet = false;
      TLorentzVector electrons, muons, jets;

      unsigned short num_ele_included = 0,
                     num_mu_included = 0,
                     num_jet_included = 0;

      for ( unsigned int k=0; k<num_electron; k++ ) {
        TLorentzVector ele; ele.SetPtEtaPhiE( electron_pt[k], electron_eta[k], electron_phi[k], electron_energy[k] );
        const TVector3 ele_vtx( electron_vtx_x[k], electron_vtx_y[k], electron_vtx_z[k] );
        const float dist = ( ele_vtx-diph_vtx ).Mag();
        if ( dist<0.1 ) { h_electron_pt_1mm->Fill( ele.Pt() ); }
        if ( dist<0.2 ) { h_electron_pt_2mm->Fill( ele.Pt() ); }
        if ( dist<0.5 ) { h_electron_pt_5mm->Fill( ele.Pt() ); }
        if ( dist<1.0 ) { h_electron_pt_1cm->Fill( ele.Pt() ); }
        for ( unsigned short l=0; l<num_bins; l++ ) {
          if ( dist<dist_bin[l] ) num_ele_bin[l]++;
          if ( ele.Pt()>oth_obj_minpt && dist<dist_bin[l] ) num_ele_ptcut_bin[l]++;
        }
        if ( dist<max_dist_vtx_inclobjects ) {
          electrons += ele;
          has_ele = true;
          num_ele_included++;
        }
      }

      for ( unsigned int k=0; k<num_muon; k++ ) {
        TLorentzVector mu; mu.SetPtEtaPhiE( muon_pt[k], muon_eta[k], muon_phi[k], muon_energy[k] );
        const TVector3 mu_vtx( muon_vtx_x[k], muon_vtx_y[k], muon_vtx_z[k] );
        const float dist = ( mu_vtx-diph_vtx ).Mag();
        if ( dist<0.1 ) { h_muon_pt_1mm->Fill( mu.Pt() ); }
        if ( dist<0.2 ) { h_muon_pt_2mm->Fill( mu.Pt() ); }
        if ( dist<0.5 ) { h_muon_pt_5mm->Fill( mu.Pt() ); }
        if ( dist<1.0 ) { h_muon_pt_1cm->Fill( mu.Pt() ); }
        for ( unsigned short l=0; l<num_bins; l++ ) {
          if ( dist<dist_bin[l] ) num_mu_bin[l]++;
          if ( mu.Pt()>oth_obj_minpt && dist<dist_bin[l] ) num_mu_ptcut_bin[l]++;
        }
        if ( dist<max_dist_vtx_inclobjects ) {
          muons += mu;
          has_muon = true;
          num_mu_included++;
        }
      }

      for ( unsigned int k=0; k<num_jet; k++ ) {
        TLorentzVector jet; jet.SetPtEtaPhiE( jet_pt[k], jet_eta[k], jet_phi[k], jet_energy[k] );
        const TVector3 jet_vtx( jet_vtx_x[k], jet_vtx_y[k], jet_vtx_z[k] );
        const float dist = ( jet_vtx-diph_vtx ).Mag();
        if ( dist<0.1 ) { h_jet_pt_1mm->Fill( jet.Pt() ); }
        if ( dist<0.2 ) { h_jet_pt_2mm->Fill( jet.Pt() ); }
        if ( dist<0.5 ) { h_jet_pt_5mm->Fill( jet.Pt() ); }
        if ( dist<1.0 ) { h_jet_pt_1cm->Fill( jet.Pt() ); }
        for ( unsigned short l=0; l<num_bins; l++ ) {
          if ( dist<dist_bin[l] ) num_jet_bin[l]++;
          if ( jet.Pt()>oth_obj_minpt && dist<dist_bin[l] ) num_jet_ptcut_bin[l]++;
        }
        if ( dist<max_dist_vtx_inclobjects ) {
          jets += jet;
          has_jet = true;
          num_jet_included++;
        }
      }

      for ( unsigned short k=0; k<num_bins; k++ ) {
        gr_veto_ele->Fill( dist_bin[k], num_ele_bin[k]/num_diphoton );
        gr_veto_mu->Fill( dist_bin[k], num_mu_bin[k]/num_diphoton );
        gr_veto_jet->Fill( dist_bin[k], num_jet_bin[k]/num_diphoton );
        gr_veto_ele_ptcut->Fill( dist_bin[k], num_ele_ptcut_bin[k]/num_diphoton );
        gr_veto_mu_ptcut->Fill( dist_bin[k], num_mu_ptcut_bin[k]/num_diphoton );
        gr_veto_jet_ptcut->Fill( dist_bin[k], num_jet_ptcut_bin[k]/num_diphoton );
      }

      const float met_x = met*cos( met_phi ),
                  met_y = met*sin( met_phi );

      const TLorentzVector lv_met( met_x, met_y, 0., met ),
                           dipho_met = pho1+pho2+lv_met,
                           //dipho_incl = pho1+pho2+electrons+muons+jets+lv_met;
                           dipho_incl = pho1+pho2+electrons+muons+jets;
      //cout << diphoton_mass[j] << "\t" << dipho_met.M() << "\t" << dipho_incl.M() << endl;
      const float diphoton_plus_met_mass = dipho_met.M(),
                  diphoton_plus_met_rap = dipho_met.Rapidity(),
                  diphoton_incl_mass = dipho_incl.M(),
                  diphoton_incl_rap = dipho_incl.Rapidity();

      //----- N-1 plots -----

      h_diphoton_pt->Fill( diphoton_pt[j] );
      h_diphoton_pt_zoom->Fill( diphoton_pt[j] );
      h_diphoton_mass->Fill( diphoton_mass[j] );
      h_diphoton_mass_withmet->Fill( diphoton_plus_met_mass );
      h_diphoton_mass_incl->Fill( diphoton_incl_mass );
      h_diphoton_rap->Fill( diphoton_rapidity[j] );
      h_diphoton_ntrk->Fill( diphoton_vertex_tracks[j] );
      h_diphoton_dphi->Fill( acopl );
      h_diphoton_dphi_zoom->Fill( acopl );
      h_diphoton_leadpt->Fill( diphoton_pt1[j] );
      h_diphoton_subleadpt->Fill( diphoton_pt2[j] );
      h_diphoton_leadeta->Fill( diphoton_eta1[j] );
      h_diphoton_subleadeta->Fill( diphoton_eta2[j] );
      h_diphoton_vtxz->Fill( diphoton_vertex_z[j] );

      h_lep_pt->Fill( ( electrons+muons ).Pt() );
      h_jet_pt->Fill( jets.Pt() );

      h_met_vs_pt->Fill( met, diphoton_pt[j] );
      h_metx_vs_mety->Fill( met_x, met_y );
      h_mggmet_vs_mgg->Fill( diphoton_plus_met_mass, diphoton_mass[j] );

      // rough check of the proton tagging information for N-1 plots
      // (to extract the tagged arms only ; not yet at a Roman pot granularity)

      bool has_singletag_45 = false, has_singletag_56 = false;

      for ( unsigned short k=0; k<num_proton; k++ ) {
        if ( proton_side[k]==sector45 ) {
          if ( proton_pot[k]==nearPot ) { has_singletag_45 = true; }
          if ( proton_pot[k]==farPot ) { has_singletag_45 = true; }
        }
        if ( proton_side[k]==sector56 ) {
          if ( proton_pot[k]==nearPot ) { has_singletag_56 = true; }
          if ( proton_pot[k]==farPot ) { has_singletag_56 = true; }
        }
      }

      bool has_singletag = ( has_singletag_45 || has_singletag_56 ),
           has_doubletag = ( has_singletag_45 && has_singletag_56 );

      if ( has_singletag ) {
        h_diphoton_pt_1tag->Fill( diphoton_pt[j] );
        h_diphoton_pt_zoom_1tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_1tag->Fill( diphoton_mass[j] );
        h_diphoton_mass_withmet_1tag->Fill( diphoton_plus_met_mass );
        h_diphoton_mass_incl_1tag->Fill( diphoton_incl_mass );
        h_diphoton_rap_1tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_1tag->Fill( min_dist_vtx );
        h_diphoton_ntrk_1tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_1tag->Fill( acopl );
        h_diphoton_dphi_zoom_1tag->Fill( acopl );
        h_diphoton_leadpt_1tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_1tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_1tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_1tag->Fill( diphoton_eta2[j] );
        h_diphoton_vtxz_1tag->Fill( diphoton_vertex_z[j] );
	
        h_lep_pt_1tag->Fill( ( electrons+muons ).Pt() );
        h_jet_pt_1tag->Fill( jets.Pt() );

        h_metx_vs_mety_1tag->Fill( met_x, met_y );

        num_diphoton_cand_1tag++;
      }
      if ( has_doubletag ) {
        //cout << "maximal diproton mass: " << max_diproton_mass << " +- " << max_diproton_mass_error << " (rapidity=" << max_diproton_mass_rap << ")" << endl;
        h_diphoton_pt_2tag->Fill( diphoton_pt[j] );
        h_diphoton_pt_zoom_2tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_2tag->Fill( diphoton_mass[j] );
        h_diphoton_mass_withmet_2tag->Fill( diphoton_plus_met_mass );
        h_diphoton_mass_incl_2tag->Fill( diphoton_incl_mass );
        h_diphoton_rap_2tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_2tag->Fill( min_dist_vtx );
        h_diphoton_ntrk_2tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_2tag->Fill( acopl );
        h_diphoton_dphi_zoom_2tag->Fill( acopl );
        h_diphoton_leadpt_2tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_2tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_2tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_2tag->Fill( diphoton_eta2[j] );
        h_diphoton_vtxz_2tag->Fill( diphoton_vertex_z[j] );

        h_ptgg_vs_mgg->SetPoint( h_ptgg_vs_mgg->GetN(), diphoton_pt[j], diphoton_mass[j] );

        h_lep_pt_2tag->Fill( ( electrons+muons ).Pt() );
        h_jet_pt_2tag->Fill( jets.Pt() );

        h_met_vs_pt_2tag->Fill( met, diphoton_pt[j] );
        h_metx_vs_mety_2tag->Fill( met_x, met_y );

        num_evts_with_tag++;
        num_diphoton_cand_2tag++;
      }

      h_num_vtx_1mm->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm->Fill( diphoton_vertex_vtx1cmdist[j] );

      num_evts_notag++;
      num_diphoton_cand++;

      //----- exclusivity cuts -----

      if ( has_ele || has_muon || has_jet ) continue; //FIXME FIXME FIXME FIXME FIXME
      if ( diphoton_pt[j]>max_pt_diphoton ) continue;
      if ( acopl>max_acopl ) continue;

      //----- FROM THAT POINT ON, EXCLUSIVE DIPHOTON CANDIDATE -----

      h_num_vtx_1mm_excl->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm_excl->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm_excl->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm_excl->Fill( diphoton_vertex_vtx1cmdist[j] );

      //----- forward tracks retrieval part -----

      float xi_45 = -1., err_xi_45 = -1.,
            xi_56 = -1., err_xi_56 = -1.;
      int id_45 = -1, id_56 = -1;

      bool compat_45f = false, compat_45n = false,
           compat_56f = false, compat_56n = false;

      for ( unsigned short k=0; k<num_proton; k++ ) {
        if ( proton_side[k]==sector45 ) {
          if ( proton_pot[k]==farPot ) {
            const unsigned short id = h_ximatch_45f->GetN();

            if ( fabs( proton_xi[k]-xip_reco )<n_sigma*proton_xi_error[k] && xip_reco>lim_45f ) {
              compat_45f = true;
              h_ximatch_45f_matched->SetPoint( id, proton_xi[k], xip_reco );
            }

            h_ximatch_45f->SetPoint( id, proton_xi[k], xip_reco );
            h_ximatch_45f->SetPointError( id, proton_xi_error[k], err_xip_reco );
            h_ximatch_45f_withmet->SetPoint( id, proton_xi[k], xip_reco_withmet );
            h_ximatch_45f_withmet->SetPointError( id, proton_xi_error[k], err_xip_reco_withmet );
            {
              xi_45 = proton_xi[k];
              err_xi_45 = proton_xi_error[k];
              id_45 = k;
            }
          }
          if ( proton_pot[k]==nearPot ) {
            const unsigned short id = h_ximatch_45n->GetN();

            if ( fabs( proton_xi[k]-xip_reco )<n_sigma*proton_xi_error[k] && xip_reco>lim_45n ) {
              compat_45n = true;
              h_ximatch_45n_matched->SetPoint( id, proton_xi[k], xip_reco );
            }

            h_ximatch_45n->SetPoint( id, proton_xi[k], xip_reco );
            h_ximatch_45n->SetPointError( id, proton_xi_error[k], err_xip_reco );
            h_ximatch_45n_withmet->SetPoint( id, proton_xi[k], xip_reco_withmet );
            h_ximatch_45n_withmet->SetPointError( id, proton_xi_error[k], err_xip_reco_withmet );
            if ( xi_45<0. ) {
              xi_45 = proton_xi[k];
              err_xi_45 = proton_xi_error[k];
              id_45 = k;
            }
          }
        }
        if ( proton_side[k]==sector56 ) {
          if ( proton_pot[k]==farPot ) {
            const unsigned short id = h_ximatch_56f->GetN();

            if ( fabs( proton_xi[k]-xim_reco )<n_sigma*proton_xi_error[k] && xim_reco>lim_56f ) {
              compat_56f = true;
              h_ximatch_56f_matched->SetPoint( id, proton_xi[k], xim_reco );
            }

            h_ximatch_56f->SetPoint( id, proton_xi[k], xim_reco );
            h_ximatch_56f->SetPointError( id, proton_xi_error[k], err_xim_reco );
            h_ximatch_56f_withmet->SetPoint( id, proton_xi[k], xim_reco_withmet );
            h_ximatch_56f_withmet->SetPointError( id, proton_xi_error[k], err_xim_reco_withmet );
            {
              xi_56 = proton_xi[k];
              err_xi_56 = proton_xi_error[k];
              id_56 = k;
            }
          }
          if ( proton_pot[k]==nearPot ) {
            const unsigned short id = h_ximatch_56n->GetN();

            if ( fabs( proton_xi[k]-xim_reco )<n_sigma*proton_xi_error[k] && xim_reco>lim_56n ) {
              compat_56n = true;
              h_ximatch_56n_matched->SetPoint( id, proton_xi[k], xim_reco );
            }

            h_ximatch_56n->SetPoint( id, proton_xi[k], xim_reco );
            h_ximatch_56n->SetPointError( id, proton_xi_error[k], err_xim_reco );
            h_ximatch_56n_withmet->SetPoint( id, proton_xi[k], xim_reco_withmet );
            h_ximatch_56n_withmet->SetPointError( id, proton_xi_error[k], err_xim_reco_withmet );
            if ( xi_56<0. ) {
              xi_56 = proton_xi[k];
              err_xi_56 = proton_xi_error[k];
              id_56 = k;
            }
          }
        }
      } // loop over proton tracks

      // dump the list of events in a text file

      if ( compat_45n || compat_45f ) {
        events_list << "1\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
      }
      if ( compat_56n || compat_56f ) {
        events_list << "1\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
      }
      if ( ( compat_45n || compat_45f ) && ( compat_56n || compat_56f ) ) {
        events_list << "2\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
      }

      if ( xi_45>=0. && xi_56>=0. ) { // double tag

        const float miss_mass = sqrt_s * sqrt( xi_45*xi_56 ),
                    err_miss_mass = ( sqrt_s/2. ) * ( 1./xi_45/xi_56 ) * sqrt( pow( err_xi_45/xi_45, 2 ) + pow( err_xi_56/xi_56, 2 ) ) / miss_mass,
                    dipr_rapidity = log( xi_45/xi_56 ) / 2.,
                    err_dipr_rapidity = sqrt( pow( err_xi_45/xi_45, 2 ) + pow( err_xi_56/xi_56, 2 ) ) / 2.;
        //cout << xi_45 << "\t" << xi_56 << "\t" << miss_mass << "\t" << err_miss_mass << "\t" << diphoton_mass[j] << endl;

        h_mgg_vs_mpp->SetPoint( cand_2tag_id, diphoton_mass[j], miss_mass );
        h_mgg_vs_mpp->SetPointError( cand_2tag_id, 0., err_miss_mass );


        h_ygg_vs_ypp->SetPoint( cand_2tag_id, diphoton_rapidity[j], dipr_rapidity );
        h_ygg_vs_ypp->SetPointError( cand_2tag_id, 0., err_dipr_rapidity );
        if ( fabs( diphoton_mass[j]-miss_mass )<1.0*err_miss_mass ) {
          h_mgg_vs_mpp_candm->SetPoint( cand_2tag_id, diphoton_mass[j], miss_mass );
          h_mgg_vs_mpp_candm->SetPointError( cand_2tag_id, 0., err_miss_mass );
          h_ygg_vs_ypp_candm->SetPoint( cand_2tag_id, diphoton_rapidity[j], dipr_rapidity );
          h_ygg_vs_ypp_candm->SetPointError( cand_2tag_id, 0., err_dipr_rapidity );
        }
        if ( fabs( diphoton_rapidity[j]-dipr_rapidity )<1.0*err_dipr_rapidity ) {
          h_mgg_vs_mpp_candy->SetPoint( cand_2tag_id, diphoton_mass[j], miss_mass );
          h_mgg_vs_mpp_candy->SetPointError( cand_2tag_id, 0., err_miss_mass );
          h_ygg_vs_ypp_candy->SetPoint( cand_2tag_id, diphoton_rapidity[j], dipr_rapidity );
          h_ygg_vs_ypp_candy->SetPointError( cand_2tag_id, 0., err_dipr_rapidity );
        }

        h_mggmet_vs_mpp->SetPoint( cand_2tag_id, diphoton_plus_met_mass, miss_mass );
        h_yggmet_vs_ypp->SetPoint( cand_2tag_id, diphoton_plus_met_rap, dipr_rapidity );
        if ( fabs( diphoton_plus_met_mass-miss_mass )<1.0*err_miss_mass ) {
          h_mggmet_vs_mpp_candm->SetPoint( cand_2tag_id, diphoton_plus_met_mass, miss_mass );
          h_yggmet_vs_ypp_candm->SetPoint( cand_2tag_id, diphoton_plus_met_rap, dipr_rapidity );
        }
        if ( fabs( diphoton_plus_met_rap-dipr_rapidity )<1.0*err_dipr_rapidity ) {
          h_mggmet_vs_mpp_candy->SetPoint( cand_2tag_id, diphoton_plus_met_mass, miss_mass );
          h_yggmet_vs_ypp_candy->SetPoint( cand_2tag_id, diphoton_plus_met_rap, dipr_rapidity );
        }

        h_mggincl_vs_mpp->SetPoint( cand_2tag_id, diphoton_incl_mass, miss_mass );
        h_yggincl_vs_ypp->SetPoint( cand_2tag_id, diphoton_incl_rap, dipr_rapidity );
        if ( fabs( diphoton_incl_mass-miss_mass )<1.0*err_miss_mass ) {
          h_mggincl_vs_mpp_candm->SetPoint( cand_2tag_id, diphoton_incl_mass, miss_mass );
          h_yggincl_vs_ypp_candm->SetPoint( cand_2tag_id, diphoton_incl_rap, dipr_rapidity );
          //cout << run_id << ":" << lumisection << ":" << event_number << "\t" << "pt===" << dipho_incl.Pt() << "\t" << diphoton_incl_mass << "\t" << diphoton_mass[j] << "\t" << num_ele_included << "\t" << num_mu_included << "\t" << num_jet_included << endl;
          h_pt_incl->Fill( dipho_incl.Pt() );
        }
        if ( fabs( diphoton_incl_rap-dipr_rapidity )<1.0*err_dipr_rapidity ) {
          h_mggincl_vs_mpp_candy->SetPoint( cand_2tag_id, diphoton_incl_mass, miss_mass );
          h_yggincl_vs_ypp_candy->SetPoint( cand_2tag_id, diphoton_incl_rap, dipr_rapidity );
        }

        h_mpp_over_mgg->Fill( miss_mass/diphoton_mass[j] );
        h_ypp_minus_ygg->Fill( dipr_rapidity - diphoton_rapidity[j] );

        if ( proton_pot[id_45]==2 ) {
          h_hitmap_45n_2tag->Fill( proton_x[id_45]*100., proton_y[id_45]*100. );
          h_xidist_45n_2tag->Fill( proton_xi[id_45] );
        }
        if ( proton_pot[id_45]==3 ) {
          h_hitmap_45f_2tag->Fill( proton_x[id_45]*100., proton_y[id_45]*100. );
          h_xidist_45f_2tag->Fill( proton_xi[id_45] );
        }
        if ( proton_pot[id_56]==2 ) {
          h_hitmap_56n_2tag->Fill( proton_x[id_56]*100., proton_y[id_56]*100. );
          h_xidist_56n_2tag->Fill( proton_xi[id_56] );
        }
        if ( proton_pot[id_56]==3 ) {
          h_hitmap_56f_2tag->Fill( proton_x[id_56]*100., proton_y[id_56]*100. );
          h_xidist_56f_2tag->Fill( proton_xi[id_56] );
        }

        cand_2tag_id++;
      }

      else if ( xi_45>=0. || xi_56>=0. ) { // single tag
        // dump the list of events in a text file
        cand_1tag_id++;
      }

    } // loop over diphotons

    //----- AFTER THE LOOP ON DIPHOTON -----

    h_num_vtx->Fill( num_vertex );
    h_met->Fill( met );
    if ( num_proton_45>0 || num_proton_56>0 ) {
      h_num_vtx_1tag->Fill( num_vertex );
      h_met_1tag->Fill( met );
    }
    if ( num_proton_45>0 && num_proton_56>0 ) {
      h_num_vtx_2tag->Fill( num_vertex );
      h_met_2tag->Fill( met );
    }

    for ( unsigned short j=0; j<num_vertex; j++ ) {
      h_vtxz->Fill( vertex_z[j] );
      if ( num_proton_45>0 || num_proton_56>0 ) { // double tag
        h_vtxz_1tag->Fill( vertex_z[j] );
      }
      if ( num_proton_45>0 && num_proton_56>0 ) { // single tag
        h_vtxz_2tag->Fill( vertex_z[j] );
      }
    }

    h_num_diphoton->Fill( num_diphoton_cand );
    //h_num_diphoton->Fill( num_diphoton );
    h_num_diphoton_1tag->Fill( num_diphoton_cand_1tag );
    h_num_diphoton_2tag->Fill( num_diphoton_cand_2tag );

    h_num_proton->Fill( num_proton );
    h_num_proton_45->Fill( num_proton_45 );
    h_num_proton_56->Fill( num_proton_56 );

  }

  // plotting stage
  cout << "events: " << num_evts_notag << ", with tag: " << num_evts_with_tag << endl;

  /*const float lumi_b = 5.060924481910, // fb-1
              lumi_c = 1.490748474431, // fb-1
              lumi_g = 3.742171002882; // fb-1*/
  const float lumi_b = 5.055026851041,
              lumi_c = 1.470474265456,
              lumi_g = 6.757328202581; //FIXME update me!!

  const float lumi = ( file.Contains( "2016BCG" ) ) ? lumi_b+lumi_c+lumi_g : lumi_b+lumi_c;

  gStyle->SetOptStat( 0 );

  const string top_label_str = ( lumi<10. ) 
    ? Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.2f fb^{-1}", lumi )
    : Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", lumi );
  const char* top_label = top_label_str.c_str();

  const Plotter plt( out_path, top_label );

  {
    //----- preselection level -----

    plt.plot_3hists( "presel_diphoton_mult", h_num_diphoton, h_num_diphoton_1tag, h_num_diphoton_2tag );

    plt.plot_3hists( "presel_diphoton_mass", h_diphoton_mass, h_diphoton_mass_1tag, h_diphoton_mass_2tag );
    plt.plot_3hists( "presel_diphoton_mass_withmet", h_diphoton_mass_withmet, h_diphoton_mass_withmet_1tag, h_diphoton_mass_withmet_2tag );
    plt.plot_3hists( "presel_diphoton_mass_incl", h_diphoton_mass_incl, h_diphoton_mass_incl_1tag, h_diphoton_mass_incl_2tag );
    plt.plot_3hists( "presel_diphoton_pt", h_diphoton_pt, h_diphoton_pt_1tag, h_diphoton_pt_2tag );
    plt.plot_3hists( "presel_diphoton_pt_zoom", h_diphoton_pt_zoom, h_diphoton_pt_zoom_1tag, h_diphoton_pt_zoom_2tag );
    plt.plot_3hists( "presel_diphoton_lead_pt", h_diphoton_leadpt, h_diphoton_leadpt_1tag, h_diphoton_leadpt_2tag );
    plt.plot_3hists( "presel_diphoton_sublead_pt", h_diphoton_subleadpt, h_diphoton_subleadpt_1tag, h_diphoton_subleadpt_2tag );
    plt.plot_3hists( "presel_diphoton_vtxz", h_diphoton_vtxz, h_diphoton_vtxz_1tag, h_diphoton_vtxz_2tag );

    /*h_diphoton_leadeta->GetYaxis()->SetRangeUser(0., 55.);
    h_diphoton_subleadeta->GetYaxis()->SetRangeUser(0., 55.);*/
    h_diphoton_leadeta->SetMinimum( 0. );
    h_diphoton_subleadeta->SetMinimum( 0. );

    plt.plot_3hists( "presel_diphoton_lead_eta", h_diphoton_leadeta, h_diphoton_leadeta_1tag, h_diphoton_leadeta_2tag );
    plt.plot_3hists( "presel_diphoton_sublead_eta", h_diphoton_subleadeta, h_diphoton_subleadeta_1tag, h_diphoton_subleadeta_2tag );
    plt.plot_3hists( "presel_diphoton_dphi", h_diphoton_dphi, h_diphoton_dphi_1tag, h_diphoton_dphi_2tag );
    plt.plot_3hists( "presel_diphoton_dphi_zoom", h_diphoton_dphi_zoom, h_diphoton_dphi_zoom_1tag, h_diphoton_dphi_zoom_2tag );
    plt.plot_3hists( "presel_diphoton_rapidity", h_diphoton_rap, h_diphoton_rap_1tag, h_diphoton_rap_2tag );
    plt.plot_3hists( "presel_diphoton_closest_vtx", h_diphoton_closestvtx, h_diphoton_closestvtx_1tag, h_diphoton_closestvtx_2tag );
    plt.plot_3hists( "presel_diphoton_vtx_numtracks", h_diphoton_ntrk, h_diphoton_ntrk_1tag, h_diphoton_ntrk_2tag );
    plt.plot_3hists( "presel_num_vertex", h_num_vtx, h_num_vtx_1tag, h_num_vtx_2tag );
    plt.plot_3hists( "presel_vertexz_allvtx", h_vtxz, h_vtxz_1tag, h_vtxz_2tag );
    plt.plot_3hists( "presel_event_met", h_met, h_met_1tag, h_met_2tag );
    //              JW
    plt.plot_3hists( "lep_pt", h_lep_pt, h_lep_pt_1tag, h_lep_pt_2tag );
    plt.plot_3hists( "jet_pt", h_jet_pt, h_jet_pt_1tag, h_jet_pt_2tag );
    //

    //plt.plot_balances( "presel_diphoton_pt_vs_diproton_mass", "Diphoton p_{T} (GeV)\\Diproton mass (GeV)", &h_ptgg_vs_mpp );

    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "Sector 45", h_num_proton_45 ) );
      hm.push_back( std::make_pair( "Sector 56", h_num_proton_56 ) );
      hm.push_back( std::make_pair( "Both sides", h_num_proton ) );
      plt.draw_multiplot( "presel_num_proton", hm );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "45-N", h_hitmap_45n ) );
      hm.push_back( std::make_pair( "45-F", h_hitmap_45f ) );
      hm.push_back( std::make_pair( "56-N", h_hitmap_56n ) );
      hm.push_back( std::make_pair( "56-F", h_hitmap_56f ) );
      plt.draw_4plots( "presel_hitmaps_allevents", hm );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "45-N", h_xidist_45n_all ) );
      hm.push_back( std::make_pair( "45-F", h_xidist_45f_all ) );
      hm.push_back( std::make_pair( "56-N", h_xidist_56n_all ) );
      hm.push_back( std::make_pair( "56-F", h_xidist_56f_all ) );
      plt.draw_4plots( "presel_xidists_allevents", hm, "hist,e" );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_electron_pt_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_electron_pt_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_electron_pt_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_electron_pt_1cm ) );
      plt.draw_multiplot( "presel_electron_pt_close_vertex", hm, true, true );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_muon_pt_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_muon_pt_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_muon_pt_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_muon_pt_1cm ) );
      plt.draw_multiplot( "presel_muon_pt_close_vertex", hm, true, true );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_jet_pt_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_jet_pt_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_jet_pt_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_jet_pt_1cm ) );
      plt.draw_multiplot( "presel_jet_pt_close_vertex", hm, true, true );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_num_vtx_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_num_vtx_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_num_vtx_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_num_vtx_1cm ) );
      plt.draw_multiplot( "presel_num_close_vertex", hm, true, true );
    }

    {
      gr_veto_ele->Scale( 1./tr->GetEntries() );
      gr_veto_mu->Scale( 1./tr->GetEntries() );
      gr_veto_jet->Scale( 1./tr->GetEntries() );
      Plotter::HistsMap gm;
      gm.push_back( std::make_pair( "all electrons", gr_veto_ele ) );
      gm.push_back( std::make_pair( "all muons", gr_veto_mu ) );
      gm.push_back( std::make_pair( "all jets", gr_veto_jet ) );
      plt.draw_multiplot( "presel_veto_eff_dist", gm );
    }
    {
      gr_veto_ele_ptcut->Scale( 1./tr->GetEntries() );
      gr_veto_mu_ptcut->Scale( 1./tr->GetEntries() );
      gr_veto_jet_ptcut->Scale( 1./tr->GetEntries() );
      Plotter::HistsMap gm;
      gm.push_back( std::make_pair( Form( "electrons, p_{T} > %.f GeV", oth_obj_minpt ), gr_veto_ele_ptcut ) );
      gm.push_back( std::make_pair( Form( "muons, p_{T} > %.f GeV", oth_obj_minpt ), gr_veto_mu_ptcut ) );
      gm.push_back( std::make_pair( Form( "jets, p_{T} > %.f GeV", oth_obj_minpt ), gr_veto_jet_ptcut ) );
      plt.draw_multiplot( "presel_veto-ptcut_eff_dist", gm );
    }

    //----- exclusive selection -----

    cout << "total candidates: " << h_mgg_vs_mpp->Integral() << endl
         << " -> with mass matching: " << h_mgg_vs_mpp_candm->Integral() << endl
         << " -> with rapiditiy matching: " << h_mgg_vs_mpp_candy->Integral() << endl;

    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "Exclusive diphotons", h_mgg_vs_mpp ) );
      gm.push_back( std::make_pair( "Mass matching", h_mgg_vs_mpp_candm ) );
      gm.push_back( std::make_pair( "Rapidity matching", h_mgg_vs_mpp_candy ) );
      plt.plot_balances( "excl_mass_balance", "Diphoton mass (GeV)\\Diproton missing mass (GeV)", gm, min_mass_diphoton, 2500. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "Exclusive diphotons", h_mggmet_vs_mpp ) );
      gm.push_back( std::make_pair( "Mass matching", h_mggmet_vs_mpp_candm ) );
      gm.push_back( std::make_pair( "Rapidity matching", h_mggmet_vs_mpp_candy ) );
      plt.plot_balances( "excl_mass_balance_withmet", "Diphoton + #slash{E}_{T} mass (GeV)\\Diproton missing mass (GeV)", gm, min_mass_diphoton, 3000. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "Exclusive diphotons", h_mggincl_vs_mpp ) );
      gm.push_back( std::make_pair( "Mass matching", h_mggincl_vs_mpp_candm ) );
      gm.push_back( std::make_pair( "Rapidity matching", h_mggincl_vs_mpp_candy ) );
      plt.plot_balances( "excl_mass_balance_incl", "Diphoton + surrounding objects mass (GeV)\\Diproton missing mass (GeV)", gm, min_mass_diphoton, 3000. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "Exclusive diphotons", h_ygg_vs_ypp ) );
      gm.push_back( std::make_pair( "Mass matching", h_ygg_vs_ypp_candm ) );
      gm.push_back( std::make_pair( "Rapidity matching", h_ygg_vs_ypp_candy ) );
      plt.plot_balances( "excl_rapidity_balance", "Diphoton rapidity\\Diproton rapidity", gm, -3., 3. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "Exclusive diphotons", h_yggmet_vs_ypp ) );
      gm.push_back( std::make_pair( "Mass matching", h_yggmet_vs_ypp_candm ) );
      gm.push_back( std::make_pair( "Rapidity matching", h_yggmet_vs_ypp_candy ) );
      plt.plot_balances( "excl_rapidity_balance_withmet", "Diphoton + #slash{E}_{T} rapidity\\Diproton rapidity", gm, -1., 1. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "Exclusive diphotons", h_yggincl_vs_ypp ) );
      gm.push_back( std::make_pair( "Mass matching", h_yggincl_vs_ypp_candm ) );
      gm.push_back( std::make_pair( "Rapidity matching", h_yggincl_vs_ypp_candy ) );
      plt.plot_balances( "excl_rapidity_balance_incl", "Diphoton + surrounding objects rapidity\\Diproton rapidity", gm, -1.5, 1.5 );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "All candidates", h_ximatch_45n ) );
      gm.push_back( std::make_pair( Form( "#xi within %.0f #sigma", n_sigma ), h_ximatch_45n_matched ) );
      //gm.push_back( std::make_pair( "Including #slash{E}_{T}", h_ximatch_45n_withmet ) );
      plt.draw_multigraph( "excl_xi1_balance_45n", gm, 0., 0.3, true, lim_45n );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "All candidates", h_ximatch_45f ) );
      gm.push_back( std::make_pair( Form( "#xi within %.0f #sigma", n_sigma ), h_ximatch_45f_matched ) );
      plt.draw_multigraph( "excl_xi1_balance_45f", gm, 0., 0.3, true, lim_45f );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "All candidates", h_ximatch_56n ) );
      gm.push_back( std::make_pair( Form( "#xi within %.0f #sigma", n_sigma ), h_ximatch_56n_matched ) );
      plt.draw_multigraph( "excl_xi1_balance_56n", gm, 0., 0.3, true, lim_56n );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( std::make_pair( "All candidates", h_ximatch_56f ) );
      gm.push_back( std::make_pair( Form( "#xi within %.0f #sigma", n_sigma ), h_ximatch_56f_matched ) );
      plt.draw_multigraph( "excl_xi1_balance_56f", gm, 0., 0.3, true, lim_56f );
    }
    {
      plt.plot_balances( "excl_xi1_balance_45n_withmet", "Proton #xi (sector 45 - near pot)\\#xi from diphoton + #slash{E}_{T}", h_ximatch_45n_withmet, 0., 0.3, 0.037 );
      plt.plot_balances( "excl_xi1_balance_45f_withmet", "Proton #xi (sector 45 - far pot)\\#xi from diphoton + #slash{E}_{T}", h_ximatch_45f_withmet, 0., 0.3, 0.026 );
      plt.plot_balances( "excl_xi1_balance_56n_withmet", "Proton #xi (sector 56 - near pot)\\#xi from diphoton + #slash{E}_{T}", h_ximatch_56n_withmet, 0., 0.3, 0.048 );
      plt.plot_balances( "excl_xi1_balance_56f_withmet", "Proton #xi (sector 56 - far pot)\\#xi from diphoton + #slash{E}_{T}", h_ximatch_56f_withmet, 0., 0.3, 0.037 );
    }

    /*{
      h_ptgg_vs_mgg.GetXaxis()->SetRangeUser( 0., 100. );
      plt.plot_balances( "excl_diphoton_pt_vs_diphoton_mass", "Diphoton p_{T} (GeV)\\Diphoton mass (GeV)", h_ptgg_vs_mgg );
    }*/
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "45-N", h_hitmap_45n_2tag ) );
      hm.push_back( std::make_pair( "45-F", h_hitmap_45f_2tag ) );
      hm.push_back( std::make_pair( "56-N", h_hitmap_56n_2tag ) );
      hm.push_back( std::make_pair( "56-F", h_hitmap_56f_2tag ) );
      plt.draw_4plots( "excl_hitmaps", hm );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "45-N", h_xidist_45n_2tag ) );
      hm.push_back( std::make_pair( "45-F", h_xidist_45f_2tag ) );
      hm.push_back( std::make_pair( "56-N", h_xidist_56n_2tag ) );
      hm.push_back( std::make_pair( "56-F", h_xidist_56f_2tag ) );
      plt.draw_4plots( "excl_xidists_2tags", hm, "hist,e" );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_num_vtx_1mm_excl ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_num_vtx_2mm_excl ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_num_vtx_5mm_excl ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_num_vtx_1cm_excl ) );
      plt.draw_multiplot( "excl_num_close_vertex", hm, true );
    }
  }

  {
    Canvas c( "met_x_vs_y", top_label );
    h_metx_vs_mety->Draw( "colz" );
    h_metx_vs_mety_1tag->Draw( "p same" );
    h_metx_vs_mety_1tag->SetMarkerStyle( 26 );
    h_metx_vs_mety_1tag->SetMarkerColor( kBlack );
    h_metx_vs_mety_2tag->Draw( "p same" );
    h_metx_vs_mety_2tag->SetMarkerStyle( 24 );
    h_metx_vs_mety_2tag->SetMarkerColor( kRed );
    c.SetLegendX1( 0.15 );
    c.SetLegendY1( 0.15 );
    c.AddLegendEntry( h_metx_vs_mety_1tag, "#geq 1 proton tags", "p" );
    c.AddLegendEntry( h_metx_vs_mety_2tag, "#geq 2 proton tags", "p" );
    c.Prettify( h_metx_vs_mety );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "diphoton_mass_vs diphotonmet_mass", top_label );
    h_mggmet_vs_mgg->Draw( "colz" );
    //c.DrawDiagonal( h_mggmet_vs_mgg->GetXaxis()->GetXmin(), h_mggmet_vs_mgg->GetXaxis()->GetXmax() );
    h_mggmet_vs_mgg->Draw( "colz same" );
    c.Prettify( h_mggmet_vs_mgg );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "diphoton_pt_vs_met", top_label );
    h_met_vs_pt->Draw( "colz" );
    h_met_vs_pt_2tag->Draw( "p same" );
    h_met_vs_pt_2tag->SetMarkerStyle( 24 );
    h_met_vs_pt_2tag->SetMarkerColor( kRed );
    c.AddLegendEntry( h_met_vs_pt_2tag, "#geq 2 proton tags", "p" );
    //c.DrawDiagonal( h_met_vs_pt->GetXaxis()->GetXmin(), h_met_vs_pt->GetXaxis()->GetXmax() );
    c.Prettify( h_met_vs_pt );
    c.Save( "pdf,png", out_path );
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
    Canvas c( "mass_ratio", top_label );
    h_mpp_over_mgg->Sumw2();
    h_mpp_over_mgg->Draw();
    h_mpp_over_mgg->SetMarkerStyle( 20 );
    h_mpp_over_mgg->SetLineColor( kBlack );
    c.Prettify( h_mpp_over_mgg );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "rapidity_difference", top_label );
    h_ypp_minus_ygg->Sumw2();
    h_ypp_minus_ygg->Draw();
    h_ypp_minus_ygg->SetMarkerStyle( 20 );
    h_ypp_minus_ygg->SetLineColor( kBlack );
    c.Prettify( h_ypp_minus_ygg );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "pt_incl_evts", top_label );
    h_pt_incl->Draw();
    c.Prettify( h_pt_incl );
    c.Save( "pdf,png", out_path );
  }
}

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 )
{
  if ( r9<0.94 ) return 0.004;
  if ( fabs( eta )<1.4442 ) return 0.002;
  if ( fabs( eta )>1.566 ) return 0.005;
  return 1.;
}
