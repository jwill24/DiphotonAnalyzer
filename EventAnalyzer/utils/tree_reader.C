#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStyle.h"

#include "Canvas.h"
#include "Plotter.h"

#include <fstream>
#include <iostream>

//#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/"
#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp/"
#define default_ntp_file "run2016BCG_21oct.root"

//#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/runBC"
//#define default_ntp_file "run2016BC_28sep.root"

//#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/matching_with_met"
//#define default_ntp_file "run2016BCG_21oct.root"

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 );

void tree_reader( TString file=default_ntp_file )
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
  float proton_xi[max_pr], proton_xi_error[max_pr];
  unsigned int proton_link_id[max_pr];
  float proton_link_dist[max_pr];
  tr->SetBranchAddress( "num_proton_track", &num_proton );
  tr->SetBranchAddress( "proton_track_side", proton_side );
  tr->SetBranchAddress( "proton_track_pot", proton_pot );
  tr->SetBranchAddress( "proton_track_xi", proton_xi );
  tr->SetBranchAddress( "proton_track_xi_error", proton_xi_error );
  tr->SetBranchAddress( "proton_track_link_nearfar", proton_link_id );
  tr->SetBranchAddress( "proton_track_link_mindist", proton_link_dist );
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
  const unsigned short max_jet = 20;
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

  TH1D* h_num_proton = new TH1D( "num_proton", "Forward track multiplicity\\Events", 6, 0., 6. );
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
       *h_jet_pt_2tag = (TH1D*)h_jet_pt->Clone( "jet_pt_2tag" );
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

    // xi / error
    typedef std::pair<float, float> trackparam_t;
    typedef std::map<unsigned int, trackparam_t> tracks_map;
    tracks_map tracks_45_near, tracks_45_far, tracks_56_near, tracks_56_far;

    // first loop to identify the RP and unflatten the collection
    for ( unsigned int j=0; j<num_proton; j++ ) {
      switch ( proton_side[j] ) {
        case 0: { // 4-5
          switch ( proton_pot[j] ) {
            case 2: { tracks_45_near.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // near pot
            case 3: { tracks_45_far.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // far pot
          }
        } break;
        case 1: { // 5-6
          switch ( proton_pot[j] ) {
            case 2: { tracks_56_near.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // near pot
            case 3: { tracks_56_far.insert( std::make_pair( j, trackparam_t( proton_xi[j], proton_xi_error[j] ) ) ); } break; // far pot
          }
        } break;
      }
    }

    const float sqs = 13.e3;
    float xi_56 = -1., err_xi_56 = -1.,
          xi_45 = -1., err_xi_45 = -1.;
    float max_diproton_mass = -1., max_diproton_mass_error = -1.,
          max_diproton_mass_rap = -999.;

    for ( tracks_map::const_iterator it_45=tracks_45_near.begin(); it_45!=tracks_45_near.end(); it_45++ ) {
      const unsigned int linked_id = proton_link_id[it_45->first];
      if ( linked_id<999 && proton_link_dist[it_45->first]<0.1 ) {
        cout << "--> " << it_45->second.first << " <> " << tracks_45_far[linked_id].first << endl;
      }
      /*if ( track )unsigned short j=0; j<tracks_45_near.size(); j++ ) {
      for ( unsigned short k=0; k<tracks_56_near.size(); k++ ) {
        const float diproton_mass = sqs*sqrt( tracks_45_near[j].first*tracks_56_near[k].first );
        if ( diproton_mass>max_diproton_mass ) {
          max_diproton_mass = diproton_mass;
          max_diproton_mass_error = sqs*sqs/4./diproton_mass * tracks_45_near[j].second * tracks_56_near[k].second;
          max_diproton_mass_rap = log( tracks_56_near[k].first/tracks_45_near[j].first )/2.;
        }
      }*/
    }

    if ( tracks_45_near.size()>0 && tracks_56_near.size()>0 ) {
      for ( tracks_map::const_iterator it_45=tracks_45_near.begin(); it_45!=tracks_45_near.end(); it_45++ ) {
        const trackparam_t xi45 = it_45->second;
        for ( tracks_map::const_iterator it_56=tracks_56_near.begin(); it_56!=tracks_56_near.end(); it_56++ ) {
          const trackparam_t xi56 = it_56->second;
          const float diproton_mass = sqs * sqrt( xi45.first * xi56.first ),
                      err_diproton_mass = sqrt( pow( xi45.second, 2 ) + pow( xi56.second, 2 ) ) * diproton_mass/2.;
          if ( diproton_mass>max_diproton_mass ) {
            max_diproton_mass = diproton_mass; max_diproton_mass_error = err_diproton_mass;
            max_diproton_mass_rap = log( xi56.first/xi45.first );
            xi_45 = xi45.first; err_xi_45 = xi45.second;
            xi_56 = xi56.first; err_xi_56 = xi56.second;
          }
        }
      }
    }

    const bool has_singletag = ( tracks_45_near.size()>0 || tracks_45_far.size() || tracks_56_near.size() || tracks_56_far.size() ),
               has_doubletag = ( ( tracks_45_near.size()>0 || tracks_45_far.size() ) && ( tracks_56_near.size() || tracks_56_far.size() ) );


    //          JW
    //

    //----- diphotons retrieval part -----

    for ( unsigned int j=0; j<num_diphoton; j++ ) {

      const float energy_corr_pho1 = photon_rel_energy_scale( diphoton_pt1[j], diphoton_eta1[j], diphoton_r9_1[j] ) * diphoton_pt1[j],
                  energy_corr_pho2 = photon_rel_energy_scale( diphoton_pt2[j], diphoton_eta2[j], diphoton_r9_2[j] ) * diphoton_pt2[j];

      const float xim_reco = ( diphoton_pt1[j] * exp( -diphoton_eta1[j] ) + diphoton_pt2[j] * exp( -diphoton_eta2[j] ) )/sqrt_s,
                  xip_reco = ( diphoton_pt1[j] * exp(  diphoton_eta1[j] ) + diphoton_pt2[j] * exp(  diphoton_eta2[j] ) )/sqrt_s,
                  err_xim_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp( -2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp( -2.*diphoton_eta2[j] ) )/sqrt_s,
                  err_xip_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp(  2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp(  2.*diphoton_eta2[j] ) )/sqrt_s; //FIXME
//cout << energy_corr_pho1 << "\t" << energy_corr_pho2 << "\t" << err_xim_reco << "\t" << err_xip_reco << endl;
      const float xim_reco_withmet = xim_reco + met/sqrt_s,
                  xip_reco_withmet = xip_reco + met/sqrt_s;

      const TVector3 diph_vtx( diphoton_vertex_x[j], diphoton_vertex_y[j], diphoton_vertex_z[j] );

      const float max_dist_vtx = 0.1; // in cm

      //----- leptons+jets retrieval part -----

      bool has_ele = false, has_muon = false, has_jet = false;
      TLorentzVector electrons, muons, jets;

      for ( unsigned int k=0; k<num_electron; k++ ) {
        TLorentzVector ele; ele.SetPtEtaPhiE( electron_pt[k], electron_eta[k], electron_phi[k], electron_energy[k] );
        const TVector3 ele_vtx( electron_vtx_x[k], electron_vtx_y[k], electron_vtx_z[k] );
        if ( ( ele_vtx-diph_vtx ).Mag()<max_dist_vtx ) {
          electrons += ele;
          has_ele = true;
        }
      }

      for ( unsigned int k=0; k<num_muon; k++ ) {
        TLorentzVector mu; mu.SetPtEtaPhiE( muon_pt[k], muon_eta[k], muon_phi[k], muon_energy[k] );
        const TVector3 mu_vtx( muon_vtx_x[k], muon_vtx_y[k], muon_vtx_z[k] );
        if ( ( mu_vtx-diph_vtx ).Mag()<max_dist_vtx ) {
          muons += mu;
          has_muon = true;
        }
      }

      for ( unsigned int k=0; k<num_jet; k++ ) {
        TLorentzVector jet; jet.SetPtEtaPhiE( jet_pt[k], jet_eta[k], jet_phi[k], jet_energy[k] );
        const TVector3 jet_vtx( jet_vtx_x[k], jet_vtx_y[k], jet_vtx_z[k] );
        if ( ( jet_vtx-diph_vtx ).Mag()<max_dist_vtx ) {
          jets += jet;
          has_jet = true;
        }
      }

///FIXME FIXME FIXME FIXME FIXME

// """exclusivity""" cut
//if ( has_ele || has_muon || has_jet ) continue;

///FIXME FIXME FIXME FIXME FIXME


      const float met_x = met*cos( met_phi ),
                  met_y = met*sin( met_phi );
      pho1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );

      //                  JW                                                                                                  
      const TLorentzVector lv_met( met_x, met_y, 0., met ),
                           dipho_met = pho1+pho2+lv_met,
                           dipho_incl = pho1+pho2+electrons+muons+jets+lv_met;
      //cout << diphoton_mass[j] << "\t" << dipho_met.M() << "\t" << dipho_incl.M() << endl;
      const float diphoton_plus_met_mass = dipho_met.M(),
                  diphoton_plus_met_rap = dipho_met.Rapidity(),
                  diphoton_incl_mass = dipho_incl.M(),
                  diphoton_incl_rap = dipho_incl.Rapidity();


      //cout << ( pho1+pho2 ).Rapidity() << "\t" << xip_reco << ", " << xim_reco << endl;
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

      /*for ( unsigned short k=0; k<tracks_45_near.size(); k++ ) {
        const unsigned short id = h_ximatch_45n->GetN();
        h_ximatch_45n->SetPoint( id, tracks_45_near[k].first, xip_reco );
        h_ximatch_45n->SetPointError( id, tracks_45_near[k].second, err_xip_reco );
      }
      for ( unsigned short k=0; k<tracks_45_far.size(); k++ ) {
        const unsigned short id = h_ximatch_45f->GetN();
        h_ximatch_45f->SetPoint( id, tracks_45_far[k].first, xip_reco );
        h_ximatch_45f->SetPointError( id, tracks_45_far[k].second, err_xip_reco );
      }
      for ( unsigned short k=0; k<tracks_56_near.size(); k++ ) {
        const unsigned short id = h_ximatch_56n->GetN();
        h_ximatch_56n->SetPoint( id, tracks_56_near[k].first, xim_reco );
        h_ximatch_56n->SetPointError( id, tracks_56_near[k].second, err_xim_reco );
      }
      for ( unsigned short k=0; k<tracks_56_far.size(); k++ ) {
        const unsigned short id = h_ximatch_56f->GetN();
        h_ximatch_56f->SetPoint( id, tracks_56_far[k].first, xim_reco );
        h_ximatch_56f->SetPointError( id, tracks_56_far[k].second, err_xim_reco );
      }*/

      if ( has_singletag ) {
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
	h_jet_pt_1tag->Fill( jet_pt[j] );
	//
      }
      if ( has_doubletag ) {
cout << "maximal diproton mass: " << max_diproton_mass << " +- " << max_diproton_mass_error << " (rapidity=" << max_diproton_mass_rap << ")" << endl;
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
	h_jet_pt_2tag->Fill( jet_pt[j] );
	//  

        h_xi1gg_vs_xi1pp->SetPoint( h_xi1gg_vs_xi1pp->GetN(), xim_reco, xi_56 );
        h_xi2gg_vs_xi2pp->SetPoint( h_xi2gg_vs_xi2pp->GetN(), xip_reco, xi_45 );

        h_xi1ggmet_vs_xi1pp->SetPoint( h_xi1ggmet_vs_xi1pp->GetN(), xim_reco_withmet, xi_56 );
        h_xi2ggmet_vs_xi2pp->SetPoint( h_xi2ggmet_vs_xi2pp->GetN(), xip_reco_withmet, xi_45 );

        h_met_vs_pt_2tag->Fill( met, diphoton_pt[j] );
        h_metx_vs_mety_2tag->Fill( met_x, met_y );

        events_list << run_id << ":" << lumisection << ":" << event_number << endl;

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
          h_xi1gg_vs_xi1pp_candm->SetPoint( h_xi1gg_vs_xi1pp_candm->GetN(), xim_reco, xi_56 );
          h_xi2gg_vs_xi2pp_candm->SetPoint( h_xi2gg_vs_xi2pp_candm->GetN(), xip_reco, xi_45 );

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
          h_xi1gg_vs_xi1pp_candy->SetPoint( h_xi1gg_vs_xi1pp_candy->GetN(), xim_reco, xi_56 );
          h_xi2gg_vs_xi2pp_candy->SetPoint( h_xi2gg_vs_xi2pp_candy->GetN(), xip_reco, xi_45 );
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
    if ( has_singletag ) {
      h_num_vtx_1tag->Fill( num_vertex );
      h_met_1tag->Fill( met );
    }
    if ( has_doubletag ) {
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
    Plotter plt( out_path, top_label );
    plt.plot_3hists( "diphoton_mass", h_diphoton_mass, h_diphoton_mass_1tag, h_diphoton_mass_2tag );
    plt.plot_3hists( "diphoton_mass_withmet", h_diphoton_mass_withmet, h_diphoton_mass_withmet_1tag, h_diphoton_mass_withmet_2tag );
    plt.plot_3hists( "diphoton_mass_incl", h_diphoton_mass_incl, h_diphoton_mass_incl_1tag, h_diphoton_mass_incl_2tag );
    plt.plot_3hists( "diphoton_pt", h_diphoton_pt, h_diphoton_pt_1tag, h_diphoton_pt_2tag );
    plt.plot_3hists( "diphoton_pt_zoom", h_diphoton_pt_zoom, h_diphoton_pt_zoom_1tag, h_diphoton_pt_zoom_2tag );
    plt.plot_3hists( "diphoton_lead_pt", h_diphoton_leadpt, h_diphoton_leadpt_1tag, h_diphoton_leadpt_2tag );
    plt.plot_3hists( "diphoton_sublead_pt", h_diphoton_subleadpt, h_diphoton_subleadpt_1tag, h_diphoton_subleadpt_2tag );

    /*h_diphoton_leadeta->GetYaxis()->SetRangeUser(0., 55.);
    h_diphoton_subleadeta->GetYaxis()->SetRangeUser(0., 55.);*/
    h_diphoton_leadeta->SetMinimum( 0. );
    h_diphoton_subleadeta->SetMinimum( 0. );

    plt.plot_3hists( "diphoton_lead_eta", h_diphoton_leadeta, h_diphoton_leadeta_1tag, h_diphoton_leadeta_2tag );
    plt.plot_3hists( "diphoton_sublead_eta", h_diphoton_subleadeta, h_diphoton_subleadeta_1tag, h_diphoton_subleadeta_2tag );
    plt.plot_3hists( "diphoton_dphi", h_diphoton_dphi, h_diphoton_dphi_1tag, h_diphoton_dphi_2tag );
    plt.plot_3hists( "diphoton_dphi_zoom", h_diphoton_dphi_zoom, h_diphoton_dphi_zoom_1tag, h_diphoton_dphi_zoom_2tag );
    plt.plot_3hists( "diphoton_rapidity", h_diphoton_rap, h_diphoton_rap_1tag, h_diphoton_rap_2tag );
    plt.plot_3hists( "diphoton_closest_vtx", h_diphoton_closestvtx, h_diphoton_closestvtx_1tag, h_diphoton_closestvtx_2tag );
    plt.plot_3hists( "diphoton_vtx_numtracks", h_diphoton_ntrk, h_diphoton_ntrk_1tag, h_diphoton_ntrk_2tag );
    plt.plot_3hists( "num_vertex", h_num_vtx, h_num_vtx_1tag, h_num_vtx_2tag );
    plt.plot_3hists( "event_met", h_met, h_met_1tag, h_met_2tag );
    //              JW
    plt.plot_3hists( "lep_pt", h_lep_pt, h_lep_pt_1tag, h_lep_pt_2tag );
    plt.plot_3hists( "jet_pt", h_jet_pt, h_jet_pt_1tag, h_lep_pt_2tag );
    //

    cout << "total candidates: " << h_mgg_vs_mpp->Integral() << endl
         << " -> with mass matching: " << h_mgg_vs_mpp_candm->Integral() << endl
         << " -> with rapiditiy matching: " << h_mgg_vs_mpp_candy->Integral() << endl;

    plt.plot_balances( "mass_balance", "Diphoton mass (GeV)\\Diproton missing mass (GeV)", h_mgg_vs_mpp, h_mgg_vs_mpp_candm, h_mgg_vs_mpp_candy, 500., 2000., rel_err_xi );
    plt.plot_balances( "mass_balance_withmet", "Diphoton + #slash{E}_{T} mass (GeV)\\Diproton missing mass (GeV)", h_mggmet_vs_mpp, h_mggmet_vs_mpp_candm, h_mggmet_vs_mpp_candy, 500., 2000., rel_err_xi );
    plt.plot_balances_old( "rapidity_balance", "Diphoton rapidity\\Diproton rapidity", h_ygg_vs_ypp, h_ygg_vs_ypp_candm, h_ygg_vs_ypp_candy, -3., 3., rel_err_xi/sqrt( 2. ), true );
    plt.plot_balances_old( "rapidity_balance_withmet", "Diphoton + #slash{E}_{T} rapidity\\Diproton rapidity", h_yggmet_vs_ypp, h_yggmet_vs_ypp_candm, h_yggmet_vs_ypp_candy, -3., 3., rel_err_xi/sqrt( 2. ), true );
    plt.plot_balances( "xi1_balance", "#xi_{1} from diphoton system\\Proton #xi_{1}", h_xi1gg_vs_xi1pp, h_xi1gg_vs_xi1pp_candm, h_xi1gg_vs_xi1pp_candy, 0., 0.45, rel_err_xi );
    plt.plot_balances( "xi2_balance", "#xi_{2} from diphoton system\\Proton #xi_{2}", h_xi2gg_vs_xi2pp, h_xi2gg_vs_xi2pp_candm, h_xi2gg_vs_xi2pp_candy, 0., 0.45, rel_err_xi );
    plt.plot_balances( "xi1_balance_withmet", "#xi_{1} from diphoton + #slash{E}_{T} system\\Proton #xi_{1}", h_xi1ggmet_vs_xi1pp, 0, 0, 0., 0.45, rel_err_xi );
    plt.plot_balances( "xi2_balance_withmet", "#xi_{2} from diphoton + #slash{E}_{T} system\\Proton #xi_{2}", h_xi2ggmet_vs_xi2pp, 0, 0, 0., 0.45, rel_err_xi );

    plt.plot_balances( "xi1_balance_45n", "Proton #xi (sector 45 - near pot)\\#xi from diphoton", h_ximatch_45n, 0, 0, 0., 0.25, 0.037 );
    plt.plot_balances( "xi1_balance_45f", "Proton #xi (sector 45 - far pot)\\#xi from diphoton", h_ximatch_45f, 0, 0, 0., 0.25, 0.026 );
    plt.plot_balances( "xi1_balance_56n", "Proton #xi (sector 56 - near pot)\\#xi from diphoton", h_ximatch_56n, 0, 0, 0., 0.25, 0.048 );
    plt.plot_balances( "xi1_balance_56f", "Proton #xi (sector 56 - far pot)\\#xi from diphoton", h_ximatch_56f, 0, 0, 0., 0.25, 0.037 );

    plt.plot_balances( "diphoton_pt_vs_diproton_mass", "Diphoton p_{T} (GeV)\\Diproton mass (GeV)", h_ptgg_vs_mpp );
    plt.plot_balances( "diphoton_pt_vs_diphoton_mass", "Diphoton p_{T} (GeV)\\Diphoton mass (GeV)", h_ptgg_vs_mgg );
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
    //c.DrawDiagonal( h_mggmet_vs_mgg->GetXaxis()->GetXmin(), h_mggmet_vs_mgg->GetXaxis()->GetXmax() );
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
    //c.DrawDiagonal( h_met_vs_pt->GetXaxis()->GetXmin(), h_met_vs_pt->GetXaxis()->GetXmax() );
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

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 )
{
  if ( r9<0.94 ) return 0.004;
  if ( fabs( eta )<1.4442 ) return 0.002;
  if ( fabs( eta )>1.566 ) return 0.005;
  return 1.;
}
