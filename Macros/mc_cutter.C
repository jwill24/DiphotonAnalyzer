#include "DatasetHandler.h"
#include "EventsSelector.h"
#include "Plotter.h"
#include "Canvas.h"
#include "PhotonScalesParser.h"

#include "THStack.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

void
mc_cutter()
{
  gSystem->Load( "libEventFilterUtilities.so" );
  DatasetHandler dsh( "datasets_list.json" );
  EventsSelector ev_selector( "/afs/cern.ch/user/l/lforthom/public/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns_preTS2.txt" );

  const unsigned short num_samples = dsh.size();

  const PhotonScalesParser pho_scales( "egammaEffi.txt_EGM2D.root" );

  const float pt_cut = 50.;
  const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;
  const float r9_cut = 0.94;
  const float mass_cut = 350.;

  unsigned int run_id, fill_number, lumisection;
  unsigned long long event_number;

  unsigned int num_diph;
  const unsigned short max_diph = 10;
  float diph_mass[max_diph], diph_ptpair[max_diph], diph_dphi[max_diph],
        diph_pt1[max_diph], diph_eta1[max_diph],
        diph_pt2[max_diph], diph_eta2[max_diph],
        diph_r9_1[max_diph], diph_r9_2[max_diph];
  float diph_sc_x_1[max_diph], diph_sc_y_1[max_diph], diph_sc_z_1[max_diph],
        diph_sc_x_2[max_diph], diph_sc_y_2[max_diph], diph_sc_z_2[max_diph];
  int diph_vtx_id[max_diph];

  unsigned int num_electron;
  const unsigned short max_ele = 20;
  float ele_pt[max_ele], ele_eta[max_ele], ele_phi[max_ele], ele_energy[max_ele];
  float ele_vtxx[max_ele], ele_vtxy[max_ele], ele_vtxz[max_ele];

  unsigned int num_muon;
  const unsigned short max_mu = 20;
  float muon_pt[max_mu], muon_eta[max_mu], muon_phi[max_mu], muon_energy[max_mu];
  float muon_vtxx[max_mu], muon_vtxy[max_mu], muon_vtxz[max_mu];

  unsigned int num_jet;
  const unsigned short max_jet = 50;
  float jet_pt[max_jet], jet_eta[max_jet], jet_phi[max_jet], jet_energy[max_jet];

  unsigned int num_vtx;
  const unsigned short max_vtx = 100;
  float vtx_x[max_vtx], vtx_y[max_vtx], vtx_z[max_vtx];
  float met;
  float pu_weight;

  const float sqrt_s = 13.e3;
  const float lumi = 9.412003739742 * 1.e3; // pre-TS2 runs with pots inserted

  TH1D* h_r9lead[2], *h_r9sublead[2];
  TH1D* h_ptlead[2], *h_ptsublead[2];
  TH1D* h_ptpair[2], *h_mpair[2], *h_acopl[2];

  for ( unsigned short i=0; i<2; i++ ) {
    h_r9lead[i] = new TH1D( Form( "r9_lead_%i", i ), "Leading photon r_{9}@@Events@@", 100, 0.5, 1.0 );
    h_r9sublead[i] = new TH1D( Form( "r9_sublead_%i", i ), "Subleading photon r_{9}@@Events@@", 100, 0.5, 1.0 );
    h_ptlead[i] = new TH1D( Form( "pt_lead_%i", i ), "Leading photon p_{T}@@Events@@GeV", 100, 50., 250. );
    h_ptsublead[i] = new TH1D( Form( "pt_sublead_%i", i ), "Subleading photon p_{T}@@Events@@GeV", 100, 0.5, 250. );
    h_ptpair[i] = new TH1D( Form( "ptpair_%i", i ), "Diphoton p_{T}@@Events@@GeV", 100, 0., 50. );
    h_mpair[i] = new TH1D( Form( "mpair_%i", i ), "Diphoton mass@@Events@@GeV", 100, 300., 1300. );
    h_acopl[i] = new TH1D( Form( "acopl_%i", i ), "Acoplanarity@@Events@@GeV", 50, 0., 0.25 );
  }

  for ( unsigned short i=0; i<num_samples; i++ ) {
    Sample s = dsh.sample( i );
    if ( s.type()==Sample::kData ) continue;
    cout << ">>> Processing " << s.name() << endl;
    const float sample_weight = s.cross_section()/s.num_events()*lumi;
    const unsigned short s_type = ( s.type()==Sample::kSignal ) ? 0 : 1;

    TTree* tree = s.tree();
    tree->SetBranchAddress( "run_id", &run_id );
    tree->SetBranchAddress( "fill_number", &fill_number );
    tree->SetBranchAddress( "lumisection", &lumisection );
    tree->SetBranchAddress( "event_number", &event_number );
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
    tree->SetBranchAddress( "diphoton_supercluster_x1", diph_sc_x_1 );
    tree->SetBranchAddress( "diphoton_supercluster_y1", diph_sc_y_1 );
    tree->SetBranchAddress( "diphoton_supercluster_z1", diph_sc_z_1 );
    tree->SetBranchAddress( "diphoton_supercluster_x2", diph_sc_x_2 );
    tree->SetBranchAddress( "diphoton_supercluster_y2", diph_sc_y_2 );
    tree->SetBranchAddress( "diphoton_supercluster_z2", diph_sc_z_2 );
    tree->SetBranchAddress( "num_electron", &num_electron );
    tree->SetBranchAddress( "electron_pt", ele_pt );
    tree->SetBranchAddress( "electron_eta", ele_eta );
    tree->SetBranchAddress( "electron_phi", ele_phi );
    tree->SetBranchAddress( "electron_energy", ele_energy );
    tree->SetBranchAddress( "electron_vtx_x", ele_vtxx );
    tree->SetBranchAddress( "electron_vtx_y", ele_vtxy );
    tree->SetBranchAddress( "electron_vtx_z", ele_vtxz );
    tree->SetBranchAddress( "num_muon", &num_muon );
    tree->SetBranchAddress( "muon_pt", muon_pt );
    tree->SetBranchAddress( "muon_eta", muon_eta );
    tree->SetBranchAddress( "muon_phi", muon_phi );
    tree->SetBranchAddress( "muon_energy", muon_energy );
    tree->SetBranchAddress( "muon_vtx_x", muon_vtxx );
    tree->SetBranchAddress( "muon_vtx_y", muon_vtxy );
    tree->SetBranchAddress( "muon_vtx_z", muon_vtxz );
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

    // loop on events
    for ( unsigned long long j=0; j<num_entries; j++ ) {

      tree->GetEntry( j );
      if ( fmod( j*1./num_entries, 0.1 )==0 ) cout << "   event " << j << " / " << num_entries << endl;

      const double sp_weight = pu_weight*sample_weight;

      // loop on diphotons
      for ( unsigned short k=0; k<num_diph; k++ ) {

        const float acop = 1.-fabs( diph_dphi[k]/TMath::Pi() );
        const float xip = ( diph_pt1[k]*exp( +diph_eta1[k] ) + diph_pt2[k]*exp( +diph_eta2[k] ) ) / sqrt_s,
                    xim = ( diph_pt1[k]*exp( -diph_eta1[k] ) + diph_pt2[k]*exp( -diph_eta2[k] ) ) / sqrt_s;
        const TVector3 diph_vtx( vtx_x[diph_vtx_id[k]], vtx_y[diph_vtx_id[k]], vtx_z[diph_vtx_id[k]] );
        const TVector3 sc_pho1( diph_sc_x_1[k], diph_sc_y_1[k], diph_sc_z_1[k] ),
                       sc_pho2( diph_sc_x_2[k], diph_sc_y_2[k], diph_sc_z_2[k] );

        // EB: 0 < |eta| < 1.4442
        // EE: |eta| > 1.566
        /*unsigned short ev_class = invalid;
        if ( diph_eta1[k]<min_etaveto && diph_eta2[k]>max_etaveto ) ev_class = ebee;
        if ( diph_eta1[k]>max_etaveto && diph_eta2[k]<min_etaveto ) ev_class = ebee;
        else if ( diph_eta1[k]<min_etaveto && diph_eta2[k]<min_etaveto ) ev_class = ebeb;
        else if ( diph_eta1[k]>max_etaveto && diph_eta2[k]>max_etaveto ) ev_class = eeee;
        if ( ev_class==invalid || ev_class==eeee ) continue;*/
        //cout << ">> evt class: " << classes[ev_class] << "  " << diph_eta1[k] << "\t" << diph_eta2[k] << endl;

        float s_weight = 1.;
        const float eff_pho1 = pho_scales.efficiency( diph_pt1[k], sc_pho1.Eta() ),
                    eff_pho2 = pho_scales.efficiency( diph_pt2[k], sc_pho2.Eta() );

        s_weight = sp_weight * ( eff_pho1*eff_pho2 );

	//--- preselection

        if ( diph_pt1[k]<pt_cut || diph_pt2[k]<pt_cut ) continue;
        if ( ( fabs( diph_eta1[k] )>min_etaveto && fabs( diph_eta1[k] )<max_etaveto ) || ( fabs( diph_eta1[k] )>eta_cut ) ||
             ( fabs( diph_eta2[k] )>min_etaveto && fabs( diph_eta2[k] )<max_etaveto ) || ( fabs( diph_eta2[k] )>eta_cut ) ) continue;
        if ( diph_r9_1[k]<r9_cut || diph_r9_2[k]<r9_cut ) continue;
        if ( diph_mass[k]<mass_cut ) continue;

        h_r9lead[s_type]->Fill( diph_r9_1[k], s_weight );
        h_r9sublead[s_type]->Fill( diph_r9_2[k], s_weight );
        h_ptlead[s_type]->Fill( diph_pt1[k], s_weight );
        h_ptsublead[s_type]->Fill( diph_pt2[k], s_weight );
        h_ptpair[s_type]->Fill( diph_ptpair[k], s_weight );
        h_mpair[s_type]->Fill( diph_mass[k], s_weight );
        h_acopl[s_type]->Fill( acop, s_weight );

        /*if ( diph_ptpair[k]<20. && acop<0.005 ) {
          is_elastic = true;
        }*/

      } // loop on diphotons
    } // loop on events
  } // loop on samples

  gStyle->SetOptStat( 0 );
  //Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/mc_comparison", Form( "CMS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", lumi/1.e3 ) );
  Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp", Form( "CMS Simulation 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", lumi/1.e3 ) );
  plt.draw_stonratio( "ston_r9lead", h_r9lead[0], h_r9lead[1] );
  plt.draw_stonratio( "ston_ptpair", h_ptpair[0], h_ptpair[1] );
  plt.draw_stonratio( "ston_mpair", h_mpair[0], h_mpair[1] );
  plt.draw_stonratio( "ston_acopl", h_acopl[0], h_acopl[1] );
  //const string class_name = Form( "#font[62]{%s}", classes[j].c_str() );
}

