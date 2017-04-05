#include "DatasetHandler.h"
#include "Plotter.h"
#include "Canvas.h"

#include "THStack.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

void
plotter()
{
  gSystem->Load( "libEventFilterUtilities.so" );
  DatasetHandler dsh( "utils/datasets_list.json" );

  const unsigned short num_samples = dsh.size();

  const string kin_regions[4] = { "nosel", "presel", "elastic", "inelastic" };
  typedef enum {
    nosel = 0,
    presel = 1,
    elastic = 2,
    inelastic = 3
  } regions;
  const unsigned short num_regions = sizeof( kin_regions )/sizeof( kin_regions[0] );

  TH1D* h_mass[num_regions][num_samples], *h_ptpair[num_regions][num_samples], *h_dphi[num_regions][num_samples],
       *h_met[num_regions][num_samples],
       *h_ptlead[num_regions][num_samples], *h_ptsublead[num_regions][num_samples], *h_pt[num_regions][num_samples],
       *h_r9lead[num_regions][num_samples], *h_r9sublead[num_regions][num_samples], *h_r9[num_regions][num_samples],
       *h_xip[num_regions][num_samples], *h_xim[num_regions][num_samples],
       *h_ndiph[num_regions][num_samples], *h_nvtx[num_regions][num_samples],
       *h_diph_vtxz[num_regions][num_samples];
  Plotter::HistsMap hm_ndiph[num_regions][3],
                    hm_mass[num_regions][3], hm_ptpair[num_regions][3], hm_dphi[num_regions][3],
                    hm_met[num_regions][3],
                    hm_ptlead[num_regions][3], hm_ptsublead[num_regions][3], hm_pt[num_regions][3],
                    hm_r9lead[num_regions][3], hm_r9sublead[num_regions][3], hm_r9[num_regions][3],
                    hm_xip[num_regions][3], hm_xim[num_regions][3],
                    hm_nvtx[num_regions][3],
                    hm_diph_vtxz[num_regions][3];

  unsigned int num_diph;
  const unsigned short max_diph = 100;
  float diph_mass[max_diph], diph_ptpair[max_diph], diph_dphi[max_diph],
        diph_pt1[max_diph], diph_eta1[max_diph],
        diph_pt2[max_diph], diph_eta2[max_diph],
        diph_r9_1[max_diph], diph_r9_2[max_diph];
  int diph_vtx_id[max_diph];

  unsigned int num_electron;
  const unsigned short max_ele = 50;
  float ele_pt[max_ele], ele_eta[max_ele], ele_phi[max_ele], ele_energy[max_ele];

  unsigned int num_muon;
  const unsigned short max_mu = 50;
  float muon_pt[max_mu], muon_eta[max_mu], muon_phi[max_mu], muon_energy[max_mu];

  unsigned int num_jet;
  const unsigned short max_jet = 50;
  float jet_pt[max_jet], jet_eta[max_jet], jet_phi[max_jet], jet_energy[max_jet];

  unsigned int num_vtx;
  const unsigned short max_vtx = 75;
  float vtx_x[max_vtx], vtx_y[max_vtx], vtx_z[max_vtx];
  //const unsigned short max_vtx = 60;
  float met;
  float pu_weight;

  const float sqrt_s = 13.e3;
  const float lumi = ( 5.055026851041 + 1.470474265456 + 7.487318307770 ) * 1.e3;

  for ( unsigned short i=0; i<num_samples; i++ ) {
    Sample s = dsh.sample( i );
    cout << ">>> Processing " << s.name() << endl;
    const float sample_weight = ( s.type()!=Sample::kData )
      ? s.cross_section()/s.num_events()*lumi
      : 1.;

    for ( unsigned short j=0; j<num_regions; j++ ) {
      h_ndiph[j][i] = new TH1D( Form( "h_%s_ndiph_%d", kin_regions[j].c_str(), i ), "Number of diphoton candidates\\Events", 10, 0., 10. );
      h_mass[j][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[j].c_str(), i ), "Diphoton mass\\Events\\GeV", 66, 350., 2000. );
      h_ptlead[j][i] = new TH1D( Form( "h_%s_ptlead_%d", kin_regions[j].c_str(), i ), "Leading photon p_{T}\\Events\\GeV", 50, 0., 500. );
      h_ptsublead[j][i] = new TH1D( Form( "h_%s_ptsublead_%d", kin_regions[j].c_str(), i ), "Subleading photon p_{T}\\Events\\GeV", 50, 0., 500. );
      h_pt[j][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[j].c_str(), i ), "Single photon p_{T}\\Events\\GeV", 50, 0., 500. );
      h_xip[j][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[j].c_str(), i ), "Diphoton #xi^{+}\\Events", 50, 0., 0.5 );
      h_xim[j][i] = new TH1D( Form( "h_%s_xim_%d", kin_regions[j].c_str(), i ), "Diphoton #xi^{-}\\Events", 50, 0., 0.5 );
      h_diph_vtxz[j][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[j].c_str(), i ), "Diphoton vertex z\\Events\\cm", 40, -20., 20. );
      h_nvtx[j][i] = new TH1D( Form( "h_%s_nvtx_%d", kin_regions[j].c_str(), i ), "Number of primary vertices\\Events", 50, 0., 50. );
      h_met[j][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[j].c_str(), i ), "Event #slash{E}_{T}\\Events\\GeV", 50, 0., 150. );
    }
    h_ptpair[nosel][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[nosel].c_str(), i ), "Diphoton p_{T}\\Events\\GeV", 25, 0., 500. );
    h_ptpair[presel][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[presel].c_str(), i ), "Diphoton p_{T}\\Events\\GeV", 50, 0., 500. );
    h_ptpair[elastic][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[elastic].c_str(), i ), "Diphoton p_{T}\\Events\\GeV", 50, 0., 50. );
    h_ptpair[inelastic][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[inelastic].c_str(), i ), "Diphoton p_{T}\\Events\\GeV", 25, 0., 500. );

    h_r9lead[nosel][i] = new TH1D( Form( "h_%s_r9lead_%d", kin_regions[nosel].c_str(), i ), "Leading photon r_{9}\\Events", 30, 0.85, 1. );
    h_r9lead[presel][i] = new TH1D( Form( "h_%s_r9lead_%d", kin_regions[presel].c_str(), i ), "Leading photon r_{9}\\Events", 30, 0.94, 1. );
    h_r9lead[elastic][i] = new TH1D( Form( "h_%s_r9lead_%d", kin_regions[elastic].c_str(), i ), "Leading photon r_{9}\\Events", 30, 0.94, 1. );
    h_r9lead[inelastic][i] = new TH1D( Form( "h_%s_r9lead_%d", kin_regions[inelastic].c_str(), i ), "Leading photon r_{9}\\Events", 30, 0.94, 1. );

    h_r9sublead[nosel][i] = new TH1D( Form( "h_%s_r9sublead_%d", kin_regions[nosel].c_str(), i ), "Subleading photon r_{9}\\Events", 30, 0.85, 1. );
    h_r9sublead[presel][i] = new TH1D( Form( "h_%s_r9sublead_%d", kin_regions[presel].c_str(), i ), "Subleading photon r_{9}\\Events", 30, 0.94, 1. );
    h_r9sublead[elastic][i] = new TH1D( Form( "h_%s_r9sublead_%d", kin_regions[elastic].c_str(), i ), "Subleading photon r_{9}\\Events", 30, 0.94, 1. );
    h_r9sublead[inelastic][i] = new TH1D( Form( "h_%s_r9sublead_%d", kin_regions[inelastic].c_str(), i ), "Subleading photon r_{9}\\Events", 30, 0.94, 1. );

    h_r9[nosel][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[nosel].c_str(), i ), "Single photon r_{9}\\Events", 30, 0.85, 1. );
    h_r9[presel][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[presel].c_str(), i ), "Single photon r_{9}\\Events", 30, 0.94, 1. );
    h_r9[elastic][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[elastic].c_str(), i ), "Single photon r_{9}\\Events", 30, 0.94, 1. );
    h_r9[inelastic][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[inelastic].c_str(), i ), "Single photon r_{9}\\Events", 30, 0.94, 1. );

    h_dphi[nosel][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[nosel].c_str(), i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 20, 0., 1. );
    h_dphi[presel][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[presel].c_str(), i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 40, 0., 1. );
    h_dphi[elastic][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[elastic].c_str(), i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 20, 0., 0.01 );
    h_dphi[inelastic][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[inelastic].c_str(), i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 40, 0., 1. );

    TTree* tree = s.tree();
    tree->SetBranchAddress( "num_diphoton", &num_diph );
    tree->SetBranchAddress( "diphoton_mass", diph_mass );
    tree->SetBranchAddress( "diphoton_pt", diph_ptpair );
    tree->SetBranchAddress( "diphoton_dphi", diph_dphi );
    tree->SetBranchAddress( "diphoton_pt1", diph_pt1 );
    tree->SetBranchAddress( "diphoton_eta1", diph_eta1 );
    tree->SetBranchAddress( "diphoton_pt2", diph_pt2 );
    tree->SetBranchAddress( "diphoton_eta2", diph_eta2 );
    tree->SetBranchAddress( "diphoton_r91", diph_r9_1 );
    tree->SetBranchAddress( "diphoton_r92", diph_r9_2 );
    tree->SetBranchAddress( "diphoton_vertex_id", diph_vtx_id );
    tree->SetBranchAddress( "num_electron", &num_electron );
    tree->SetBranchAddress( "electron_pt", ele_pt );
    tree->SetBranchAddress( "electron_eta", ele_eta );
    tree->SetBranchAddress( "electron_phi", ele_phi );
    tree->SetBranchAddress( "electron_energy", ele_energy );
    tree->SetBranchAddress( "num_muon", &num_muon );
    tree->SetBranchAddress( "muon_pt", muon_pt );
    tree->SetBranchAddress( "muon_eta", muon_eta );
    tree->SetBranchAddress( "muon_phi", muon_phi );
    tree->SetBranchAddress( "muon_energy", muon_energy );
    tree->SetBranchAddress( "num_jet", &num_jet );
    tree->SetBranchAddress( "jet_pt", jet_pt );
    tree->SetBranchAddress( "jet_eta", jet_eta );
    tree->SetBranchAddress( "jet_phi", jet_phi );
    tree->SetBranchAddress( "jet_energy", jet_energy );
    tree->SetBranchAddress( "num_vertex", &num_vtx );
    tree->SetBranchAddress( "vertex_x", vtx_x );
    tree->SetBranchAddress( "vertex_y", vtx_y );
    tree->SetBranchAddress( "vertex_z", vtx_z );
    tree->SetBranchAddress( "met", &met );
    if ( tree->SetBranchAddress( "pileup_weight", &pu_weight )!=0 ) pu_weight = 1.;
    const unsigned long long num_entries = tree->GetEntries();

    const float pt_cut = 50.;
    const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;
    const float r9_cut = 0.94;
    const float mass_cut = 350.;

    // loop on events
    for ( unsigned long long j=0; j<num_entries; j++ ) {
      if ( fmod( j*1./num_entries, 0.1 )==0 ) {
        cout << "   event " << j << " / " << num_entries << endl;
      }
      tree->GetEntry( j );
      bool has_diphoton_cand = false,
           has_elastic_diphoton_cand = false,
           has_inelastic_diphoton_cand = false;

      //const float s_weight = ( s.type()==Sample::kData ) ? 1. : ( sample_weight * pu_weight );
      const double s_weight = pu_weight;

      // loop on diphotons
      for ( unsigned short k=0; k<num_diph; k++ ) {

        const float acop = 1.-fabs( diph_dphi[k]/TMath::Pi() );
        const float xip = ( diph_pt1[k]*exp( +diph_eta1[k] ) + diph_pt2[k]*exp( +diph_eta2[k] ) ) / sqrt_s,
                    xim = ( diph_pt1[k]*exp( -diph_eta1[k] ) + diph_pt2[k]*exp( -diph_eta2[k] ) ) / sqrt_s;
        const TVector3 diph_vtx( vtx_x[diph_vtx_id[k]], vtx_y[diph_vtx_id[k]], vtx_z[diph_vtx_id[k]] );

        h_ptlead[nosel][i]->Fill( diph_pt1[k], s_weight );
        h_ptsublead[nosel][i]->Fill( diph_pt2[k], s_weight );
        h_pt[nosel][i]->Fill( diph_pt1[k], s_weight );
        h_pt[nosel][i]->Fill( diph_pt2[k], s_weight );

        h_r9lead[nosel][i]->Fill( diph_r9_1[k], s_weight );
        h_r9sublead[nosel][i]->Fill( diph_r9_2[k], s_weight );
        h_r9[nosel][i]->Fill( diph_r9_1[k], s_weight );
        h_r9[nosel][i]->Fill( diph_r9_2[k], s_weight );

        h_mass[nosel][i]->Fill( diph_mass[k], s_weight );
        h_ptpair[nosel][i]->Fill( diph_ptpair[k], s_weight );
        h_dphi[nosel][i]->Fill( acop, s_weight );
        h_met[nosel][i]->Fill( met, s_weight );

        h_xip[nosel][i]->Fill( xip, s_weight );
        h_xim[nosel][i]->Fill( xim, s_weight );

        h_diph_vtxz[nosel][i]->Fill( diph_vtx.z(), s_weight );

        if ( diph_pt1[k]<pt_cut || diph_pt2[k]<pt_cut ) continue;
        if ( ( fabs( diph_eta1[k] )>min_etaveto && fabs( diph_eta1[k] )<max_etaveto ) ||
             ( fabs( diph_eta2[k] )>min_etaveto && fabs( diph_eta2[k] )<max_etaveto ) ||
             ( fabs( diph_eta1[k] )>eta_cut ) ||
             ( fabs( diph_eta2[k] )>eta_cut ) ) continue;
        if ( diph_r9_1[k]<r9_cut || diph_r9_2[k]<r9_cut ) continue;
        if ( diph_mass[k]<mass_cut ) continue;

        h_ptlead[presel][i]->Fill( diph_pt1[k], s_weight );
        h_ptsublead[presel][i]->Fill( diph_pt2[k], s_weight );
        h_pt[presel][i]->Fill( diph_pt1[k], s_weight );
        h_pt[presel][i]->Fill( diph_pt2[k], s_weight );

        h_r9lead[presel][i]->Fill( diph_r9_1[k], s_weight );
        h_r9sublead[presel][i]->Fill( diph_r9_2[k], s_weight );
        h_r9[presel][i]->Fill( diph_r9_1[k], s_weight );
        h_r9[presel][i]->Fill( diph_r9_2[k], s_weight );

        h_mass[presel][i]->Fill( diph_mass[k], s_weight );
        h_ptpair[presel][i]->Fill( diph_ptpair[k], s_weight );
        h_dphi[presel][i]->Fill( acop, s_weight );
        h_met[presel][i]->Fill( met, s_weight );

        h_xip[presel][i]->Fill( xip, s_weight );
        h_xim[presel][i]->Fill( xim, s_weight );

        h_diph_vtxz[presel][i]->Fill( diph_vtx.z(), s_weight );

        for ( unsigned short l=0; l<num_electron; l++ ) {
          TLorentzVector ele; ele.SetPtEtaPhiE( ele_pt[l], ele_eta[l], ele_phi[l], ele_energy[l] );
        }
        for ( unsigned short l=0; l<num_muon; l++ ) {
          TLorentzVector mu; mu.SetPtEtaPhiE( muon_pt[l], muon_eta[l], muon_phi[l], muon_energy[l] );
        }
        for ( unsigned short l=0; l<num_jet; l++ ) {
          TLorentzVector jet; jet.SetPtEtaPhiE( jet_pt[l], jet_eta[l], jet_phi[l], jet_energy[l] );
        }

        has_diphoton_cand = true;

        bool is_elastic = false;
        if ( diph_ptpair[k]<20. ) {
          if ( acop<0.005 ) {
            is_elastic = true;
          }
        }

        if ( is_elastic ) {
          h_ptlead[elastic][i]->Fill( diph_pt1[k], s_weight );
          h_ptsublead[elastic][i]->Fill( diph_pt2[k], s_weight );
          h_pt[elastic][i]->Fill( diph_pt1[k], s_weight );
          h_pt[elastic][i]->Fill( diph_pt2[k], s_weight );

          h_r9lead[elastic][i]->Fill( diph_r9_1[k], s_weight );
          h_r9sublead[elastic][i]->Fill( diph_r9_2[k], s_weight );
          h_r9[elastic][i]->Fill( diph_r9_1[k], s_weight );
          h_r9[elastic][i]->Fill( diph_r9_2[k], s_weight );

          h_mass[elastic][i]->Fill( diph_mass[k], s_weight );
          h_ptpair[elastic][i]->Fill( diph_ptpair[k], s_weight );
          h_dphi[elastic][i]->Fill( acop, s_weight );
          h_met[elastic][i]->Fill( met, s_weight );

          h_xip[elastic][i]->Fill( xip, s_weight );
          h_xim[elastic][i]->Fill( xim, s_weight );

          h_diph_vtxz[elastic][i]->Fill( diph_vtx.z(), s_weight );

          has_elastic_diphoton_cand = true;
        }
        if ( !is_elastic ) {
          h_ptlead[inelastic][i]->Fill( diph_pt1[k], s_weight );
          h_ptsublead[inelastic][i]->Fill( diph_pt2[k], s_weight );
          h_pt[inelastic][i]->Fill( diph_pt1[k], s_weight );
          h_pt[inelastic][i]->Fill( diph_pt2[k], s_weight );

          h_r9lead[inelastic][i]->Fill( diph_r9_1[k], s_weight );
          h_r9sublead[inelastic][i]->Fill( diph_r9_2[k], s_weight );
          h_r9[inelastic][i]->Fill( diph_r9_1[k], s_weight );
          h_r9[inelastic][i]->Fill( diph_r9_2[k], s_weight );

          h_mass[inelastic][i]->Fill( diph_mass[k], s_weight );
          h_ptpair[inelastic][i]->Fill( diph_ptpair[k], s_weight );
          h_dphi[inelastic][i]->Fill( acop, s_weight );
          h_met[inelastic][i]->Fill( met, s_weight );

          h_xip[inelastic][i]->Fill( xip, s_weight );
          h_xim[inelastic][i]->Fill( xim, s_weight );

          h_diph_vtxz[inelastic][i]->Fill( diph_vtx.z(), s_weight );

          has_inelastic_diphoton_cand = true;
        }

        if ( s.type()==Sample::kSignal ) {
          //cout << ">> MC signal candidate!" << endl;
        }

      }
      h_ndiph[nosel][i]->Fill( num_diph, s_weight );
      h_nvtx[nosel][i]->Fill( num_vtx, s_weight );
      if ( has_diphoton_cand ) {
        h_ndiph[presel][i]->Fill( num_diph, s_weight );
        h_nvtx[presel][i]->Fill( num_vtx, s_weight );
      }
      if ( has_elastic_diphoton_cand ) {
        h_ndiph[elastic][i]->Fill( num_diph, s_weight );
        h_nvtx[elastic][i]->Fill( num_vtx, s_weight );
      }
      if ( has_inelastic_diphoton_cand ) {
        h_ndiph[inelastic][i]->Fill( num_diph, s_weight );
        h_nvtx[inelastic][i]->Fill( num_vtx, s_weight );
      }
    } // loop on diphotons

    unsigned short type = 0;
    if ( s.type()==Sample::kData ) { type = 0; }
    else if ( s.type()==Sample::kBackground ) { type = 1; }
    else if ( s.type()==Sample::kSignal ) { type = 2; }

    for ( unsigned short j=0; j<num_regions; j++ ) {
      h_ndiph[j][i]->Scale( sample_weight );
      h_mass[j][i]->Scale( sample_weight );
      h_ptpair[j][i]->Scale( sample_weight );
      h_met[j][i]->Scale( sample_weight );
      h_ptlead[j][i]->Scale( sample_weight );
      h_ptsublead[j][i]->Scale( sample_weight );
      h_pt[j][i]->Scale( sample_weight );
      h_r9lead[j][i]->Scale( sample_weight );
      h_r9sublead[j][i]->Scale( sample_weight );
      h_r9[j][i]->Scale( sample_weight );
      h_dphi[j][i]->Scale( sample_weight );
      h_xip[j][i]->Scale( sample_weight );
      h_xim[j][i]->Scale( sample_weight );
      h_nvtx[j][i]->Scale( sample_weight );
      h_diph_vtxz[j][i]->Scale( sample_weight );

      hm_ndiph[j][type].push_back( make_pair( s.name(), h_ndiph[j][i] ) );
      hm_mass[j][type].push_back( make_pair( s.name(), h_mass[j][i] ) );
      hm_ptpair[j][type].push_back( make_pair( s.name(), h_ptpair[j][i] ) );
      hm_met[j][type].push_back( make_pair( s.name(), h_met[j][i] ) );
      hm_ptlead[j][type].push_back( make_pair( s.name(), h_ptlead[j][i] ) );
      hm_ptsublead[j][type].push_back( make_pair( s.name(), h_ptsublead[j][i] ) );
      hm_pt[j][type].push_back( make_pair( s.name(), h_pt[j][i] ) );
      hm_r9lead[j][type].push_back( make_pair( s.name(), h_r9lead[j][i] ) );
      hm_r9sublead[j][type].push_back( make_pair( s.name(), h_r9sublead[j][i] ) );
      hm_r9[j][type].push_back( make_pair( s.name(), h_r9[j][i] ) );
      hm_dphi[j][type].push_back( make_pair( s.name(), h_dphi[j][i] ) );
      hm_xip[j][type].push_back( make_pair( s.name(), h_xip[j][i] ) );
      hm_xim[j][type].push_back( make_pair( s.name(), h_xim[j][i] ) );
      hm_nvtx[j][type].push_back( make_pair( s.name(), h_nvtx[j][i] ) );
      hm_diph_vtxz[j][type].push_back( make_pair( s.name(), h_diph_vtxz[j][i] ) );
    }

  } // loop on events

  gStyle->SetOptStat( 0 );
  Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/mc_comparison", "CMS Preliminary 2016, #sqrt{s} = 13 TeV" );
  for ( unsigned short i=0; i<num_regions; i++ ) {
    plt.draw_multiplot( Form( "%s_diphoton_mult", kin_regions[i].c_str() ), hm_ndiph[i][0], hm_ndiph[i][1], hm_ndiph[i][2], true );
    plt.draw_multiplot( Form( "%s_diphoton_mass", kin_regions[i].c_str() ), hm_mass[i][0], hm_mass[i][1], hm_mass[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_mass_logscale", kin_regions[i].c_str() ), hm_mass[i][0], hm_mass[i][1], hm_mass[i][2], true );
    plt.draw_multiplot( Form( "%s_diphoton_ptpair", kin_regions[i].c_str() ), hm_ptpair[i][0], hm_ptpair[i][1], hm_ptpair[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_ptpair_logscale", kin_regions[i].c_str() ), hm_ptpair[i][0], hm_ptpair[i][1], hm_ptpair[i][2], true );
    plt.draw_multiplot( Form( "%s_diphoton_met", kin_regions[i].c_str() ), hm_met[i][0], hm_met[i][1], hm_met[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_pt_singlepho", kin_regions[i].c_str() ), hm_pt[i][0], hm_pt[i][1], hm_pt[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_pt_leadpho", kin_regions[i].c_str() ), hm_ptlead[i][0], hm_ptlead[i][1], hm_ptlead[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_pt_subleadpho", kin_regions[i].c_str() ), hm_ptsublead[i][0], hm_ptsublead[i][1], hm_ptsublead[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_r9_singlepho", kin_regions[i].c_str() ), hm_r9[i][0], hm_r9[i][1], hm_r9[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_r9_leadpho", kin_regions[i].c_str() ), hm_r9lead[i][0], hm_r9lead[i][1], hm_r9lead[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_r9_subleadpho", kin_regions[i].c_str() ), hm_r9sublead[i][0], hm_r9sublead[i][1], hm_r9sublead[i][2], false );
    plt.draw_multiplot( Form( "%s_diphoton_dphi", kin_regions[i].c_str() ), hm_dphi[i][0], hm_dphi[i][1], hm_dphi[i][2] );
    plt.draw_multiplot( Form( "%s_diphoton_xip", kin_regions[i].c_str() ), hm_xip[i][0], hm_xip[i][1], hm_xip[i][2] );
    plt.draw_multiplot( Form( "%s_diphoton_xim", kin_regions[i].c_str() ), hm_xim[i][0], hm_xim[i][1], hm_xim[i][2] );
    plt.draw_multiplot( Form( "%s_diphoton_nvtx", kin_regions[i].c_str() ), hm_nvtx[i][0], hm_nvtx[i][1], hm_nvtx[i][2] );
    plt.draw_multiplot( Form( "%s_diphoton_vtx_z", kin_regions[i].c_str() ), hm_diph_vtxz[i][0], hm_diph_vtxz[i][1], hm_diph_vtxz[i][2] );
  }

}
