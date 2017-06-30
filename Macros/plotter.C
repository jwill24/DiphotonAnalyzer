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

void logarithmicBins( TAxis* axis );
void plot_2ddiscrim( const char* name, TH2D* h2[], bool logx );

typedef enum {
  the_data = 0,
  mc_inclusive = 1,
  mc_signal = 2,
  num_types
} sample_types;

const float sqrt_s = 13.e3;
//const float lumi = ( 5.055026851041 + 1.470474265456 + 7.487318307770 ) * 1.e3;
//const float lumi = ( 5.060574552184 + 1.485056762163 + 7.500305757473 ) * 1.e3; // computed from processedLumis.json
const float the_lumi = 9.412003739742 * 1.e3; // pre-TS2 runs with pots inserted

map<string,float> pots_accept = { { "45N", 0.033 }, { "45F", 0.024 }, { "56N", 0.050 }, { "56F", 0.037 } };

void
plotter()
{
  const float scaling_signal = 100.;

  gSystem->Load( "libEventFilterUtilities.so" );
  DatasetHandler dsh( "datasets_list.json" );
  EventsSelector ev_selector( "/afs/cern.ch/user/l/lforthom/public/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns_preTS2.txt" );

  const unsigned short num_samples = dsh.size();

  typedef enum {
    nosel = 0,
    presel = 1,
    elastic = 2,
    inelastic = 3,
    num_regions
  } regions;
  const string kin_regions[num_regions] = { "nosel", "presel", "elastic", "inelastic" };
  typedef enum {
    ebeb = 0,
    ebee = 1,
    eeee = 2,
    invalid,
    num_classes
  } classes_;
  const string classes[num_classes] = { "EBEB", "EBEE", "EEEE", "invalid" };

  const PhotonScalesParser pho_scales( "egammaEffi.txt_EGM2D.root" );

  TH1D* h_mass[num_regions][num_classes-2][num_samples], *h_ptpair[num_regions][num_classes-2][num_samples], *h_dphi[num_regions][num_classes-2][num_samples],
       *h_met[num_regions][num_classes-2][num_samples],
       *h_ptlead[num_regions][num_classes-2][num_samples], *h_ptsublead[num_regions][num_classes-2][num_samples], *h_pt[num_regions][num_classes-2][num_samples],
       *h_r9lead[num_regions][num_classes-2][num_samples], *h_r9sublead[num_regions][num_classes-2][num_samples], *h_r9[num_regions][num_classes-2][num_samples],
       *h_xip[num_regions][num_classes-2][num_samples], *h_xim[num_regions][num_classes-2][num_samples],
       *h_diph_vtxz[num_regions][num_classes-2][num_samples];
  TH1D* h_ndiph[num_regions][num_classes-2][num_samples],
       *h_nvtx[num_regions][num_classes-2][num_samples];
  TH2D* h2_excl_sel[num_classes], *h2_excl_acop_dpt[num_classes], *h2_excl_jet[num_classes];

  for ( unsigned short i=0; i<num_types; i++ ) {
    h2_excl_sel[i] = new TH2D( Form( "exclusivity_cuts_%d", i ), "Distance to the closest lepton vertex (cm)\\Acoplanarity 1-#||{#Delta#phi/#pi}", 10, -2.5, 0.8, 10, -5., 0. ); logarithmicBins( h2_excl_sel[i]->GetXaxis() ); logarithmicBins( h2_excl_sel[i]->GetYaxis() );
    h2_excl_acop_dpt[i] = new TH2D( Form( "excl_cuts_acop_dpt_%d", i ), "Diphoton #Deltap_{T} (GeV)\\Acoplanarity 1-#||{#Delta#phi/#pi}", 10, -2., 2.3, 10, -5., 0. ); logarithmicBins( h2_excl_acop_dpt[i]->GetXaxis() ); logarithmicBins( h2_excl_acop_dpt[i]->GetYaxis() );
    h2_excl_jet[i] = new TH2D( Form( "excl_cuts_jet_%d", i ), "High-p_{T} jets associated to the #gamma#gamma vertex\\Acoplanarity 1-#||{#Delta#phi/#pi}", 5, 0., 5., 10, -5., 0. ); logarithmicBins( h2_excl_jet[i]->GetYaxis() );
  }

  //Plotter::HistsMap hm_nvtx[num_regions][3];
  //vector< pair<const char*, TH1*> > hm_nvtx[num_regions][3];
  Plotter::HistsMap hm_ndiph[num_regions][num_classes-2][3], hm_nvtx[num_regions][num_classes-2][3],
                    hm_mass[num_regions][num_classes-2][3], hm_ptpair[num_regions][num_classes-2][3], hm_dphi[num_regions][num_classes-2][3],
                    hm_met[num_regions][num_classes-2][3],
                    hm_ptlead[num_regions][num_classes-2][3], hm_ptsublead[num_regions][num_classes-2][3], hm_pt[num_regions][num_classes-2][3],
                    hm_r9lead[num_regions][num_classes-2][3], hm_r9sublead[num_regions][num_classes-2][3], hm_r9[num_regions][num_classes-2][3],
                    hm_xip[num_regions][num_classes-2][3], hm_xim[num_regions][num_classes-2][3],
                    hm_diph_vtxz[num_regions][num_classes-2][3];

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
  int jet_dipho_match[max_jet];

  unsigned int num_vtx;
  const unsigned short max_vtx = 100;
  float vtx_x[max_vtx], vtx_y[max_vtx], vtx_z[max_vtx];
  float met;
  float pu_weight;

  map<string,float> in_pot_accept;
  for ( const auto& pot : pots_accept ) {
    in_pot_accept[pot.first] = 0.; // first initialise the counter
  }
  cout << "at run start: " << endl;
  for ( const auto& pot : in_pot_accept ) {
    cout << "pot " << pot.first << ":: " << pot.second << " events" << endl;
  }

  for ( unsigned short i=0; i<num_samples; i++ ) {
    Sample s = dsh.sample( i );
    cout << ">>> Processing " << s.name() << endl;
    const float sample_weight = ( s.type()!=Sample::kData )
      ? s.cross_section()/s.num_events()*the_lumi
      : 1.;

    sample_types sample_type = the_data;
    if ( s.type()==Sample::kData ) sample_type = the_data;
    else if ( s.type()==Sample::kBackground ) { sample_type = mc_inclusive; }
    else if ( s.type()==Sample::kSignal ) { sample_type = mc_signal; }

    for ( unsigned short k=0; k<num_classes-2; k++ ) {
      for ( unsigned short j=0; j<num_regions; j++ ) {
        const char* region = kin_regions[j].c_str();
        h_ndiph[j][k][i] = new TH1D( Form( "h_%s_%d_ndiph_%d", region, k, i ), "Number of diphoton candidates\\Events", 7, 1., 8. );
        h_nvtx[j][k][i] = new TH1D( Form( "h_%s_%d_nvtx_%d", region, k, i ), "Number of primary vertices\\Events", 50, 0., 50. );
      }
      h_ptlead[nosel][k][i] = new TH1D( Form( "h_%s_%d_ptlead_%d", kin_regions[nosel].c_str(), k, i ), "Leading photon p_{T}\\Events\\GeV", 40, 0., 600. );
      h_ptlead[presel][k][i] = new TH1D( Form( "h_%s_%d_ptlead_%d", kin_regions[presel].c_str(), k, i ), "Leading photon p_{T}\\Events\\GeV", 40, 0., 600. );
      h_ptlead[elastic][k][i] = new TH1D( Form( "h_%s_%d_ptlead_%d", kin_regions[elastic].c_str(), k, i ), "Leading photon p_{T}\\Events\\GeV", 20, 0., 600. );
      h_ptlead[inelastic][k][i] = new TH1D( Form( "h_%s_%d_ptlead_%d", kin_regions[inelastic].c_str(), k, i ), "Leading photon p_{T}\\Events\\GeV", 40, 0., 600. );

      h_ptsublead[nosel][k][i] = new TH1D( Form( "h_%s_%d_ptsublead_%d", kin_regions[nosel].c_str(), k, i ), "Subleading photon p_{T}\\Events\\GeV", 40, 0., 600. );
      h_ptsublead[presel][k][i] = new TH1D( Form( "h_%s_%d_ptsublead_%d", kin_regions[presel].c_str(), k, i ), "Subleading photon p_{T}\\Events\\GeV", 40, 0., 600. );
      h_ptsublead[elastic][k][i] = new TH1D( Form( "h_%s_%d_ptsublead_%d", kin_regions[elastic].c_str(), k, i ), "Subleading photon p_{T}\\Events\\GeV", 20, 0., 600. );
      h_ptsublead[inelastic][k][i] = new TH1D( Form( "h_%s_%d_ptsublead_%d", kin_regions[inelastic].c_str(), k, i ), "Subleading photon p_{T}\\Events\\GeV", 40, 0., 600. );

      h_pt[nosel][k][i] = new TH1D( Form( "h_%s_%d_pt_%d", kin_regions[nosel].c_str(), k, i ), "Single photon p_{T}\\Events\\GeV", 40, 0., 600. );
      h_pt[presel][k][i] = new TH1D( Form( "h_%s_%d_pt_%d", kin_regions[presel].c_str(), k, i ), "Single photon p_{T}\\Events\\GeV", 40, 0., 600. );
      h_pt[elastic][k][i] = new TH1D( Form( "h_%s_%d_pt_%d", kin_regions[elastic].c_str(), k, i ), "Single photon p_{T}\\Events\\GeV", 20, 0., 600. );
      h_pt[inelastic][k][i] = new TH1D( Form( "h_%s_%d_pt_%d", kin_regions[inelastic].c_str(), k, i ), "Single photon p_{T}\\Events\\GeV", 40, 0., 600. );

      h_mass[nosel][k][i] = new TH1D( Form( "h_%s_%d_mass_%d", kin_regions[nosel].c_str(), k, i ), "Diphoton mass\\Events\\GeV", 33, 350., 2000. );
      h_mass[presel][k][i] = new TH1D( Form( "h_%s_%d_mass_%d", kin_regions[presel].c_str(), k, i ), "Diphoton mass\\Events\\GeV", 33, 350., 2000. );
      h_mass[elastic][k][i] = new TH1D( Form( "h_%s_%d_mass_%d", kin_regions[elastic].c_str(), k, i ), "Diphoton mass\\Events\\GeV", 16, 350., 1950. );
      h_mass[inelastic][k][i] = new TH1D( Form( "h_%s_%d_mass_%d", kin_regions[inelastic].c_str(), k, i ), "Diphoton mass\\Events\\GeV", 33, 350., 2000. );

      h_ptpair[nosel][k][i] = new TH1D( Form( "h_%s_%d_ptpair_%d", kin_regions[nosel].c_str(), k, i ), "Diphoton p_{T}\\Events\\GeV", 25, 0., 500. );
      h_ptpair[presel][k][i] = new TH1D( Form( "h_%s_%d_ptpair_%d", kin_regions[presel].c_str(), k, i ), "Diphoton p_{T}\\Events\\GeV", 25, 0., 500. );
      h_ptpair[elastic][k][i] = new TH1D( Form( "h_%s_%d_ptpair_%d", kin_regions[elastic].c_str(), k, i ), "Diphoton p_{T}\\Events\\GeV", 25, 0., 50. );
      h_ptpair[inelastic][k][i] = new TH1D( Form( "h_%s_%d_ptpair_%d", kin_regions[inelastic].c_str(), k, i ), "Diphoton p_{T}\\Events\\GeV", 25, 0., 500. );

      h_r9lead[nosel][k][i] = new TH1D( Form( "h_%s_%d_r9lead_%d", kin_regions[nosel].c_str(), k, i ), "Leading photon r_{9}\\Events\\?.2f", 50, 0.5, 1. );
      h_r9lead[presel][k][i] = new TH1D( Form( "h_%s_%d_r9lead_%d", kin_regions[presel].c_str(), k, i ), "Leading photon r_{9}\\Events\\?.2f", 30, 0.94, 1. );
      h_r9lead[elastic][k][i] = new TH1D( Form( "h_%s_%d_r9lead_%d", kin_regions[elastic].c_str(), k, i ), "Leading photon r_{9}\\Events\\?.2f", 15, 0.94, 1. );
      h_r9lead[inelastic][k][i] = new TH1D( Form( "h_%s_%d_r9lead_%d", kin_regions[inelastic].c_str(), k, i ), "Leading photon r_{9}\\Events\\?.2f", 30, 0.94, 1. );

      h_r9sublead[nosel][k][i] = new TH1D( Form( "h_%s_%d_r9sublead_%d", kin_regions[nosel].c_str(), k, i ), "Subleading photon r_{9}\\Events\\?.2f", 50, 0.5, 1. );
      h_r9sublead[presel][k][i] = new TH1D( Form( "h_%s_%d_r9sublead_%d", kin_regions[presel].c_str(), k, i ), "Subleading photon r_{9}\\Events\\?.2f", 30, 0.94, 1. );
      h_r9sublead[elastic][k][i] = new TH1D( Form( "h_%s_%d_r9sublead_%d", kin_regions[elastic].c_str(), k, i ), "Subleading photon r_{9}\\Events\\?.2f", 15, 0.94, 1. );
      h_r9sublead[inelastic][k][i] = new TH1D( Form( "h_%s_%d_r9sublead_%d", kin_regions[inelastic].c_str(), k, i ), "Subleading photon r_{9}\\Events\\?.2f", 30, 0.94, 1. );

      h_r9[nosel][k][i] = new TH1D( Form( "h_%s_%d_r9_%d", kin_regions[nosel].c_str(), k, i ), "Single photon r_{9}\\Events\\?.2f", 50, 0.5, 1. );
      h_r9[presel][k][i] = new TH1D( Form( "h_%s_%d_r9_%d", kin_regions[presel].c_str(), k, i ), "Single photon r_{9}\\Events\\?.2f", 30, 0.94, 1. );
      h_r9[elastic][k][i] = new TH1D( Form( "h_%s_%d_r9_%d", kin_regions[elastic].c_str(), k, i ), "Single photon r_{9}\\Events\\?.2f", 15, 0.94, 1. );
      h_r9[inelastic][k][i] = new TH1D( Form( "h_%s_%d_r9_%d", kin_regions[inelastic].c_str(), k, i ), "Single photon r_{9}\\Events\\?.2f", 30, 0.94, 1. );

      h_dphi[nosel][k][i] = new TH1D( Form( "h_%s_%d_dphi_%d", kin_regions[nosel].c_str(), k, i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 20, 0., 1. );
      h_dphi[presel][k][i] = new TH1D( Form( "h_%s_%d_dphi_%d", kin_regions[presel].c_str(), k, i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 40, 0., 1. );
      h_dphi[elastic][k][i] = new TH1D( Form( "h_%s_%d_dphi_%d", kin_regions[elastic].c_str(), k, i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 20, 0., 0.01 );
      h_dphi[inelastic][k][i] = new TH1D( Form( "h_%s_%d_dphi_%d", kin_regions[inelastic].c_str(), k, i ), "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.3f", 40, 0., 1. );

      h_met[nosel][k][i] = new TH1D( Form( "h_%s_%d_met_%d", kin_regions[nosel].c_str(), k, i ), "Event #slash{E}_{T}\\Events\\GeV", 50, 0., 150. );
      h_met[presel][k][i] = new TH1D( Form( "h_%s_%d_met_%d", kin_regions[presel].c_str(), k, i ), "Event #slash{E}_{T}\\Events\\GeV", 50, 0., 150. );
      h_met[elastic][k][i] = new TH1D( Form( "h_%s_%d_met_%d", kin_regions[elastic].c_str(), k, i ), "Event #slash{E}_{T}\\Events\\GeV", 25, 0., 150. );
      h_met[inelastic][k][i] = new TH1D( Form( "h_%s_%d_met_%d", kin_regions[inelastic].c_str(), k, i ), "Event #slash{E}_{T}\\Events\\GeV", 50, 0., 150. );

      h_diph_vtxz[nosel][k][i] = new TH1D( Form( "h_%s_%d_diph_vtxz_%d", kin_regions[nosel].c_str(), k, i ), "Diphoton vertex z\\Events\\cm", 40, -20., 20. );
      h_diph_vtxz[presel][k][i] = new TH1D( Form( "h_%s_%d_diph_vtxz_%d", kin_regions[presel].c_str(), k, i ), "Diphoton vertex z\\Events\\cm", 40, -20., 20. );
      h_diph_vtxz[elastic][k][i] = new TH1D( Form( "h_%s_%d_diph_vtxz_%d", kin_regions[elastic].c_str(), k, i ), "Diphoton vertex z\\Events\\cm", 20, -20., 20. );
      h_diph_vtxz[inelastic][k][i] = new TH1D( Form( "h_%s_%d_diph_vtxz_%d", kin_regions[inelastic].c_str(), k, i ), "Diphoton vertex z\\Events\\cm", 40, -20., 20. );

      h_xip[nosel][k][i] = new TH1D( Form( "h_%s_%d_xip_%d", kin_regions[nosel].c_str(), k, i ), "Diphoton #xi^{+}\\Events", 50, 0., 0.5 );
      h_xip[presel][k][i] = new TH1D( Form( "h_%s_%d_xip_%d", kin_regions[presel].c_str(), k, i ), "Diphoton #xi^{+}\\Events", 50, 0., 0.5 );
      h_xip[elastic][k][i] = new TH1D( Form( "h_%s_%d_xip_%d", kin_regions[elastic].c_str(), k, i ), "Diphoton #xi^{+}\\Events", 25, 0., 0.5 );
      h_xip[inelastic][k][i] = new TH1D( Form( "h_%s_%d_xip_%d", kin_regions[inelastic].c_str(), k, i ), "Diphoton #xi^{+}\\Events", 50, 0., 0.5 );

      h_xim[nosel][k][i] = new TH1D( Form( "h_%s_%d_xim_%d", kin_regions[nosel].c_str(), k, i ), "Diphoton #xi^{-}\\Events", 50, 0., 0.5 );
      h_xim[presel][k][i] = new TH1D( Form( "h_%s_%d_xim_%d", kin_regions[presel].c_str(), k, i ), "Diphoton #xi^{-}\\Events", 50, 0., 0.5 );
      h_xim[elastic][k][i] = new TH1D( Form( "h_%s_%d_xim_%d", kin_regions[elastic].c_str(), k, i ), "Diphoton #xi^{-}\\Events", 25, 0., 0.5 );
      h_xim[inelastic][k][i] = new TH1D( Form( "h_%s_%d_xim_%d", kin_regions[inelastic].c_str(), k, i ), "Diphoton #xi^{-}\\Events", 50, 0., 0.5 );
    }
    for ( unsigned short j=0; j<num_regions; j++ ) {
      const char* region = kin_regions[j].c_str();
    }

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
    tree->SetBranchAddress( "jet_dipho_match", jet_dipho_match );
    tree->SetBranchAddress( "num_vertex", &num_vtx );
    tree->SetBranchAddress( "vertex_x", vtx_x );
    tree->SetBranchAddress( "vertex_y", vtx_y );
    tree->SetBranchAddress( "vertex_z", vtx_z );
    tree->SetBranchAddress( "met", &met );
    if ( tree->SetBranchAddress( "pileup_weight", &pu_weight )!=0 ) pu_weight = 1.;
    const unsigned long long num_entries = tree->GetEntriesFast();

    const float pt_cut = 50.;
    const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;
    const float r9_cut = 0.94;
    const float mass_cut = 350.;

    // loop on events
    for ( unsigned long long j=0; j<num_entries; j++ ) {
      if ( fmod( j*1./num_entries, 0.1 )==0 ) { cout << "   event " << j << " / " << num_entries << endl; }

      tree->GetEntry( j );

      bool has_diphoton_cand = false,
           has_elastic_diphoton_cand = false,
           has_inelastic_diphoton_cand = false;

      if ( s.type()==Sample::kData && !ev_selector.isSelected( run_id, lumisection, event_number ) ) continue;


      //const float s_weight = ( s.type()==Sample::kData ) ? 1. : ( sample_weight * pu_weight );
      const double sp_weight = pu_weight;

      // loop on diphotons
      for ( unsigned short k=0; k<num_diph; k++ ) {

        // EB: 0 < |eta| < 1.4442
        // EE: |eta| > 1.566
        unsigned short ev_class = invalid;
        if ( diph_eta1[k]<min_etaveto && diph_eta2[k]>max_etaveto ) ev_class = ebee;
        if ( diph_eta1[k]>max_etaveto && diph_eta2[k]<min_etaveto ) ev_class = ebee;
        else if ( diph_eta1[k]<min_etaveto && diph_eta2[k]<min_etaveto ) ev_class = ebeb;
        else if ( diph_eta1[k]>max_etaveto && diph_eta2[k]>max_etaveto ) ev_class = eeee;

        //----- only keep EBEE and EBEB diphoton events

        if ( ev_class==invalid || ev_class==eeee ) continue;
        //cout << ">> evt class: " << classes[ev_class] << "  " << diph_eta1[k] << "\t" << diph_eta2[k] << endl;

        const float acop = 1.-fabs( diph_dphi[k]/M_PI );
        const float xip = ( diph_pt1[k]*exp( +diph_eta1[k] ) + diph_pt2[k]*exp( +diph_eta2[k] ) ) / sqrt_s,
                    xim = ( diph_pt1[k]*exp( -diph_eta1[k] ) + diph_pt2[k]*exp( -diph_eta2[k] ) ) / sqrt_s;

        const TVector3 diph_vtx( vtx_x[diph_vtx_id[k]], vtx_y[diph_vtx_id[k]], vtx_z[diph_vtx_id[k]] );
        const TVector3 sc_pho1( diph_sc_x_1[k], diph_sc_y_1[k], diph_sc_z_1[k] ),
                       sc_pho2( diph_sc_x_2[k], diph_sc_y_2[k], diph_sc_z_2[k] );

        float s_weight = 1.;
        if ( s.type()!=Sample::kData ) {
          const float eff_pho1 = pho_scales.efficiency( diph_pt1[k], sc_pho1.Eta() ),
                      eff_pho2 = pho_scales.efficiency( diph_pt2[k], sc_pho2.Eta() );

//std::cout << "--> " << eff_pho1 << "\t" << eff_pho2 << std::endl;

          s_weight = sp_weight * ( eff_pho1*eff_pho2 );
        }
        if ( sample_type==mc_signal ) s_weight *= scaling_signal;

        //--- no selection

        h_ptlead[nosel][ev_class][i]->Fill( diph_pt1[k], s_weight );
        h_ptsublead[nosel][ev_class][i]->Fill( diph_pt2[k], s_weight );
        h_pt[nosel][ev_class][i]->Fill( diph_pt1[k], s_weight );
        h_pt[nosel][ev_class][i]->Fill( diph_pt2[k], s_weight );

        h_r9lead[nosel][ev_class][i]->Fill( diph_r9_1[k], s_weight );
        h_r9sublead[nosel][ev_class][i]->Fill( diph_r9_2[k], s_weight );
        h_r9[nosel][ev_class][i]->Fill( diph_r9_1[k], s_weight );
        h_r9[nosel][ev_class][i]->Fill( diph_r9_2[k], s_weight );

        h_mass[nosel][ev_class][i]->Fill( diph_mass[k], s_weight );
        h_ptpair[nosel][ev_class][i]->Fill( diph_ptpair[k], s_weight );
        h_dphi[nosel][ev_class][i]->Fill( acop, s_weight );
        h_met[nosel][ev_class][i]->Fill( met, s_weight );

        h_xip[nosel][ev_class][i]->Fill( xip, s_weight );
        h_xim[nosel][ev_class][i]->Fill( xim, s_weight );

        h_diph_vtxz[nosel][ev_class][i]->Fill( diph_vtx.z(), s_weight );

        //--- preselection

        if ( diph_pt1[k]<pt_cut || diph_pt2[k]<pt_cut ) continue;
        if ( ( fabs( diph_eta1[k] )>min_etaveto && fabs( diph_eta1[k] )<max_etaveto ) ||
             ( fabs( diph_eta2[k] )>min_etaveto && fabs( diph_eta2[k] )<max_etaveto ) ||
             ( fabs( diph_eta1[k] )>eta_cut ) ||
             ( fabs( diph_eta2[k] )>eta_cut ) ) continue;
        if ( diph_r9_1[k]<r9_cut || diph_r9_2[k]<r9_cut ) continue;
        if ( diph_mass[k]<mass_cut ) continue;

        h_ptlead[presel][ev_class][i]->Fill( diph_pt1[k], s_weight );
        h_ptsublead[presel][ev_class][i]->Fill( diph_pt2[k], s_weight );
        h_pt[presel][ev_class][i]->Fill( diph_pt1[k], s_weight );
        h_pt[presel][ev_class][i]->Fill( diph_pt2[k], s_weight );

        h_r9lead[presel][ev_class][i]->Fill( diph_r9_1[k], s_weight );
        h_r9sublead[presel][ev_class][i]->Fill( diph_r9_2[k], s_weight );
        h_r9[presel][ev_class][i]->Fill( diph_r9_1[k], s_weight );
        h_r9[presel][ev_class][i]->Fill( diph_r9_2[k], s_weight );

        h_mass[presel][ev_class][i]->Fill( diph_mass[k], s_weight );
        h_ptpair[presel][ev_class][i]->Fill( diph_ptpair[k], s_weight );
        h_dphi[presel][ev_class][i]->Fill( acop, s_weight );
        h_met[presel][ev_class][i]->Fill( met, s_weight );

        h_xip[presel][ev_class][i]->Fill( xip, s_weight );
        h_xim[presel][ev_class][i]->Fill( xim, s_weight );

        h_diph_vtxz[presel][ev_class][i]->Fill( diph_vtx.z(), s_weight );

        //----- look at surrounding objects

        const float min_lep_vtx_dist = 2.0; // in cm
        float closest_lep_vtx_dist = 999.; // in cm

        unsigned short num_matched_ele = 0, num_matched_mu = 0, num_matched_jets = 0;
        for ( unsigned short l=0; l<num_electron; l++ ) {
          TLorentzVector ele; ele.SetPtEtaPhiE( ele_pt[l], ele_eta[l], ele_phi[l], ele_energy[l] );
          const TVector3 ele_vtx( ele_vtxx[l], ele_vtxy[l], ele_vtxz[l] );
          const float ele_dist = ( ele_vtx-diph_vtx ).Mag();
          if ( ele_dist<closest_lep_vtx_dist ) closest_lep_vtx_dist = ele_dist;
          if ( ele_dist<min_lep_vtx_dist ) num_matched_ele++;
        }
        for ( unsigned short l=0; l<num_muon; l++ ) {
          TLorentzVector mu; mu.SetPtEtaPhiE( muon_pt[l], muon_eta[l], muon_phi[l], muon_energy[l] );
          const TVector3 mu_vtx( muon_vtxx[l], muon_vtxy[l], muon_vtxz[l] );
          const float mu_dist = ( mu_vtx-diph_vtx ).Mag();
          if ( mu_dist<closest_lep_vtx_dist ) closest_lep_vtx_dist = mu_dist;
          if ( mu_dist<min_lep_vtx_dist ) num_matched_mu++;
        }
        for ( unsigned short l=0; l<num_jet; l++ ) {
          if ( jet_dipho_match[l]!=k ) continue;
          if ( jet_pt[l]<500. ) continue;
          TLorentzVector jet; jet.SetPtEtaPhiE( jet_pt[l], jet_eta[l], jet_phi[l], jet_energy[l] );
          num_matched_jets++;
        }

        has_diphoton_cand = true;

        //----- from that point on, we have a diphoton candidate

        bool is_elastic = false;
        //if ( diph_ptpair[k]<20. && acop<0.005 ) {
        //if ( diph_ptpair[k]<20. && acop<0.005 ) {
        if ( acop<0.005 ) {
          if ( num_matched_jets==0 && num_matched_ele==0 && num_matched_mu==0 ) {
            is_elastic = true;
          }
        }

        //--- elastic/non-elastic regions

        if ( is_elastic ) {
          h_ptlead[elastic][ev_class][i]->Fill( diph_pt1[k], s_weight );
          h_ptsublead[elastic][ev_class][i]->Fill( diph_pt2[k], s_weight );
          h_pt[elastic][ev_class][i]->Fill( diph_pt1[k], s_weight );
          h_pt[elastic][ev_class][i]->Fill( diph_pt2[k], s_weight );

          h_r9lead[elastic][ev_class][i]->Fill( diph_r9_1[k], s_weight );
          h_r9sublead[elastic][ev_class][i]->Fill( diph_r9_2[k], s_weight );
          h_r9[elastic][ev_class][i]->Fill( diph_r9_1[k], s_weight );
          h_r9[elastic][ev_class][i]->Fill( diph_r9_2[k], s_weight );

          h_mass[elastic][ev_class][i]->Fill( diph_mass[k], s_weight );
          h_ptpair[elastic][ev_class][i]->Fill( diph_ptpair[k], s_weight );
          h_dphi[elastic][ev_class][i]->Fill( acop, s_weight );
          h_met[elastic][ev_class][i]->Fill( met, s_weight );

          h_xip[elastic][ev_class][i]->Fill( xip, s_weight );
          h_xim[elastic][ev_class][i]->Fill( xim, s_weight );

          h_diph_vtxz[elastic][ev_class][i]->Fill( diph_vtx.z(), s_weight );

          has_elastic_diphoton_cand = true;
        }
        if ( !is_elastic ) {
          h_ptlead[inelastic][ev_class][i]->Fill( diph_pt1[k], s_weight );
          h_ptsublead[inelastic][ev_class][i]->Fill( diph_pt2[k], s_weight );
          h_pt[inelastic][ev_class][i]->Fill( diph_pt1[k], s_weight );
          h_pt[inelastic][ev_class][i]->Fill( diph_pt2[k], s_weight );

          h_r9lead[inelastic][ev_class][i]->Fill( diph_r9_1[k], s_weight );
          h_r9sublead[inelastic][ev_class][i]->Fill( diph_r9_2[k], s_weight );
          h_r9[inelastic][ev_class][i]->Fill( diph_r9_1[k], s_weight );
          h_r9[inelastic][ev_class][i]->Fill( diph_r9_2[k], s_weight );

          h_mass[inelastic][ev_class][i]->Fill( diph_mass[k], s_weight );
          h_ptpair[inelastic][ev_class][i]->Fill( diph_ptpair[k], s_weight );
          h_dphi[inelastic][ev_class][i]->Fill( acop, s_weight );
          h_met[inelastic][ev_class][i]->Fill( met, s_weight );

          h_xip[inelastic][ev_class][i]->Fill( xip, s_weight );
          h_xim[inelastic][ev_class][i]->Fill( xim, s_weight );

          h_diph_vtxz[inelastic][ev_class][i]->Fill( diph_vtx.z(), s_weight );

          has_inelastic_diphoton_cand = true;
        }

        if ( s.type()==Sample::kSignal ) {
          //cout << ">> MC signal candidate!" << endl;
        }

        //h2_excl_sel[i]->Fill( num_matched_ele+num_matched_mu+num_matched_jets, acop, 1./num_entries );
        h2_excl_sel[sample_type]->Fill( closest_lep_vtx_dist, acop, s_weight );
        h2_excl_acop_dpt[sample_type]->Fill( fabs( diph_pt1[k]-diph_pt2[k] ), acop, s_weight );
        h2_excl_jet[sample_type]->Fill( num_matched_jets, acop, s_weight );

        if ( sample_type==mc_inclusive && has_elastic_diphoton_cand ) {
          for ( const auto& pot : pots_accept ) {
            if ( pot.first.find( "45" )!=string::npos && xip>pot.second ) { in_pot_accept[pot.first] += s_weight; }
            if ( pot.first.find( "56" )!=string::npos && xim>pot.second ) { in_pot_accept[pot.first] += s_weight; }
//cout << "----> event with weight: " << s_weight << ":" << xim << "/" << xip << endl;
//cout << "new " << pot.first << " pot content: " << in_pot_accept[pot.first] << endl;
          }
        }

      } // loop on diphotons

      h_ndiph[nosel][0][i]->Fill( num_diph, sp_weight );
      h_nvtx[nosel][0][i]->Fill( num_vtx, sp_weight );
      if ( has_diphoton_cand ) {
        h_ndiph[presel][0][i]->Fill( num_diph, sp_weight );
        h_nvtx[presel][0][i]->Fill( num_vtx, sp_weight );
      }
      if ( has_elastic_diphoton_cand ) {
        h_ndiph[elastic][0][i]->Fill( num_diph, sp_weight );
        h_nvtx[elastic][0][i]->Fill( num_vtx, sp_weight );
      }
      if ( has_inelastic_diphoton_cand ) {
        h_ndiph[inelastic][0][i]->Fill( num_diph, sp_weight );
        h_nvtx[inelastic][0][i]->Fill( num_vtx, sp_weight );
      }
    } // loop on events

    for ( unsigned short j=0; j<num_regions; j++ ) {
      for ( unsigned short k=0; k<num_classes-2; k++ ) {
        h_nvtx[j][k][i]->Scale( sample_weight );
        h_ndiph[j][k][i]->Scale( sample_weight );
        h_mass[j][k][i]->Scale( sample_weight );
        h_ptpair[j][k][i]->Scale( sample_weight );
        h_met[j][k][i]->Scale( sample_weight );
        h_ptlead[j][k][i]->Scale( sample_weight );
        h_ptsublead[j][k][i]->Scale( sample_weight );
        h_pt[j][k][i]->Scale( sample_weight );
        h_r9lead[j][k][i]->Scale( sample_weight );
        h_r9sublead[j][k][i]->Scale( sample_weight );
        h_r9[j][k][i]->Scale( sample_weight );
        h_dphi[j][k][i]->Scale( sample_weight );
        h_xip[j][k][i]->Scale( sample_weight );
        h_xim[j][k][i]->Scale( sample_weight );
        h_diph_vtxz[j][k][i]->Scale( sample_weight );

        const string sample_name = ( sample_type==mc_signal && scaling_signal!=1. ) ? Form( "%s (#times%.0f)", s.name(), scaling_signal ) : s.name();

        hm_nvtx[j][k][sample_type].emplace_back( sample_name, h_nvtx[j][k][i] );
        hm_ndiph[j][k][sample_type].emplace_back( sample_name, h_ndiph[j][k][i] );
        hm_mass[j][k][sample_type].emplace_back( sample_name, h_mass[j][k][i] );
        hm_ptpair[j][k][sample_type].emplace_back( sample_name, h_ptpair[j][k][i] );
        hm_met[j][k][sample_type].emplace_back( sample_name, h_met[j][k][i] );
        hm_ptlead[j][k][sample_type].emplace_back( sample_name, h_ptlead[j][k][i] );
        hm_ptsublead[j][k][sample_type].emplace_back( sample_name, h_ptsublead[j][k][i] )  ;
        hm_pt[j][k][sample_type].emplace_back( sample_name, h_pt[j][k][i] )  ;
        hm_r9lead[j][k][sample_type].emplace_back( sample_name, h_r9lead[j][k][i] );
        hm_r9sublead[j][k][sample_type].emplace_back( sample_name, h_r9sublead[j][k][i] );
        hm_r9[j][k][sample_type].emplace_back( sample_name, h_r9[j][k][i] );
        hm_dphi[j][k][sample_type].emplace_back( sample_name, h_dphi[j][k][i] );
        hm_xip[j][k][sample_type].emplace_back( sample_name, h_xip[j][k][i] );
        hm_xim[j][k][sample_type].emplace_back( sample_name, h_xim[j][k][i] );
        hm_diph_vtxz[j][k][sample_type].emplace_back( sample_name, h_diph_vtxz[j][k][i] );

        if ( s.type()==Sample::kBackground ) {
          h_nvtx[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_ndiph[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_mass[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_ptpair[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_met[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_ptlead[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_ptsublead[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_pt[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_r9lead[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_r9sublead[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_r9[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_dphi[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_xip[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_xim[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
          h_diph_vtxz[j][k][i]->SetFillColorAlpha( s.colour(), 0.5 );
        }
      }
    }

  } // loop on samples

std::cout << "event finished" << std::endl;

  for ( const auto& pot : in_pot_accept ) {
    cout << " --> events in the acceptance of " << pot.first << ": " << pot.second << endl;
  }


  gStyle->SetOptStat( 0 );
  Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp", Form( "CMS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", the_lumi/1.e3 ) );
  for ( unsigned short i=0; i<num_regions; i++ ) {
    for ( unsigned short j=0; j<num_classes-2; j++ ) {
      const string class_name = Form( "#font[62]{%s}", classes[j].c_str() );
      plt.draw_multiplot( Form( "%s_diphoton_mass_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_mass[i][j][0], hm_mass[i][j][1], hm_mass[i][j][2], class_name.c_str(), false, false );
      plt.draw_multiplot( Form( "%s_diphoton_mass_logscale_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_mass[i][j][0], hm_mass[i][j][1], hm_mass[i][j][2], class_name.c_str(), false, true );
      plt.draw_multiplot( Form( "%s_diphoton_ptpair_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_ptpair[i][j][0], hm_ptpair[i][j][1], hm_ptpair[i][j][2], class_name.c_str(), false, false );
      plt.draw_multiplot( Form( "%s_diphoton_ptpair_logscale_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_ptpair[i][j][0], hm_ptpair[i][j][1], hm_ptpair[i][j][2], class_name.c_str(), false, true );
      plt.draw_multiplot( Form( "%s_diphoton_met_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_met[i][j][0], hm_met[i][j][1], hm_met[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_pt_singlepho_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_pt[i][j][0], hm_pt[i][j][1], hm_pt[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_pt_leadpho_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_ptlead[i][j][0], hm_ptlead[i][j][1], hm_ptlead[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_pt_subleadpho_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_ptsublead[i][j][0], hm_ptsublead[i][j][1], hm_ptsublead[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_r9_singlepho_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_r9[i][j][0], hm_r9[i][j][1], hm_r9[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_r9_leadpho_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_r9lead[i][j][0], hm_r9lead[i][j][1], hm_r9lead[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_r9_subleadpho_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_r9sublead[i][j][0], hm_r9sublead[i][j][1], hm_r9sublead[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_dphi_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_dphi[i][j][0], hm_dphi[i][j][1], hm_dphi[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_xip_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_xip[i][j][0], hm_xip[i][j][1], hm_xip[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_xim_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_xim[i][j][0], hm_xim[i][j][1], hm_xim[i][j][2], class_name.c_str(), false );
      plt.draw_multiplot( Form( "%s_diphoton_vtx_z_%s", kin_regions[i].c_str(), classes[j].c_str() ), hm_diph_vtxz[i][j][0], hm_diph_vtxz[i][j][1], hm_diph_vtxz[i][j][2], class_name.c_str(), false );
    }
    plt.draw_multiplot( Form( "%s_diphoton_mult", kin_regions[i].c_str() ), hm_ndiph[i][0][0], hm_ndiph[i][0][1], hm_ndiph[i][0][2], "", false, true );
    plt.draw_multiplot( Form( "%s_diphoton_nvtx", kin_regions[i].c_str() ), hm_nvtx[i][0][0], hm_nvtx[i][0][1], hm_nvtx[i][0][2], "", false, false );
  }

  /*FIXME FIXME FIXME FIXME FIXME FIXME FIXME
  plot_2ddiscrim( "presel_excl_selection", h2_excl_sel, true );
  plot_2ddiscrim( "presel_excl_dpt_vs_acop", h2_excl_acop_dpt, true );
  plot_2ddiscrim( "presel_excl_matched_selection", h2_excl_jet, false );
  */


  for ( unsigned short i=0; i<num_classes-2; i++ ) {
    Canvas c( "presel_diphoton_pt_sigovbckg", Form( "CMS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", the_lumi/1.e3 ) );
    TH1D* hist_bck = (TH1D*)hm_ptpair[presel][i][1][0].second->Clone( "bck" ),
         *hist_sig = (TH1D*)hm_ptpair[presel][i][2][0].second->Clone( "sig" );
    hist_bck->Scale( 0. ); hist_sig->Scale( 0. );
    for ( Plotter::HistsMap::const_iterator it=hm_ptpair[presel][i][1].begin(); it!=hm_ptpair[presel][i][1].end(); it++ ) { hist_bck->Add( it->second ); }
    for ( Plotter::HistsMap::const_iterator it=hm_ptpair[presel][i][2].begin(); it!=hm_ptpair[presel][i][2].end(); it++ ) { hist_sig->Add( it->second ); }

    TH1D* hist_signif = (TH1D*)hist_bck->Clone( "signif" ); hist_signif->Scale( 0. );
    for ( int j=1; j<hist_signif->GetNbinsX(); j++ ) {
      const double bck_value = hist_bck->GetBinContent( j ), sig_value = hist_sig->GetBinContent( j );
      const double signif = ( sig_value+bck_value!=0. ) ? sig_value / sqrt( sig_value+bck_value ) : 0.;
      cout << " bin " << j << ": background: " << bck_value << ", signal: " << sig_value << " --> significance: " << signif << endl;
      hist_signif->SetBinContent( j, signif );
    }
    hist_signif->Draw( "p" );
    hist_signif->SetMarkerStyle( 20 );
    c.Prettify( hist_signif );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/mc_comparison" );
  }
}


void logarithmicBins( TAxis* axis )
{
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = ( to-from )/bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for ( int i=0; i<=bins; i++ ) {
    new_bins[i] = TMath::Power( 10, from+i*width );
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}

void plot_2ddiscrim( const char* name, TH2D* h2[], bool logx )
{
  Canvas c( name, "CMS Simulation Preliminary, #sqrt{s} = 13 TeV" );
  THStack st;
  for ( unsigned short i=1; i<num_types; i++ ) {
    st.Add( h2[i] );
    st.SetTitle( h2[i]->GetTitle() );
    h2[i]->SetLineColor( 1+i );
    h2[i]->SetMarkerStyle( 19+i );
    h2[i]->SetMarkerColor( 1+i );
    h2[i]->SetMarkerSize( 0.5 );
  }
  st.Draw( "p,nostack" );
  //st.Draw( "box,nostack" );
  //c.AddLegendEntry( h2[the_data], "Data" );
  c.AddLegendEntry( h2[mc_inclusive], "Incl. backgrounds" );
  c.AddLegendEntry( h2[mc_signal], "Exclusive signal" );
  if ( logx ) c.SetLogx();
  c.SetLogy();
  TLine cut( st.GetHistogram()->GetXaxis()->GetXmin(), 0.005, st.GetHistogram()->GetXaxis()->GetXmax(), 0.005 );
  cut.SetLineColor( kBlack );
  cut.SetLineWidth( 3 );
  cut.SetLineStyle( 2 );
  cut.Draw();
  c.Prettify( st.GetHistogram() );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
