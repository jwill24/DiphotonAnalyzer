#include "DatasetHandler.h"
#include "Plotter.h"
#include "Canvas.h"

#include "THStack.h"
#include "TH1.h"
#include "TMath.h"

void
plotter()
{
  gSystem->Load( "libEventFilterUtilities.so" );
  DatasetHandler dsh( "utils/datasets_list.json" );

  const unsigned short num_samples = dsh.size();
  TH1D* h_presel_mass[num_samples], *h_presel_pt[num_samples], *h_presel_dphi[num_samples],
       *h_presel_xip[num_samples], *h_presel_xim[num_samples],
       *h_presel_nvtx[num_samples];
  Plotter::HistsMap hm_presel_mass_mc, hm_presel_mass_data,
                    hm_presel_pt_mc, hm_presel_pt_data,
                    hm_presel_dphi_mc, hm_presel_dphi_data,
                    hm_presel_xip_mc, hm_presel_xip_data,
                    hm_presel_xim_mc, hm_presel_xim_data,
                    hm_presel_nvtx_mc, hm_presel_nvtx_data;

  unsigned int num_diph;
  const unsigned short max_diph = 100;
  float diph_mass[max_diph], diph_pt[max_diph], diph_dphi[max_diph],
        diph_pt1[max_diph], diph_eta1[max_diph],
        diph_pt2[max_diph], diph_eta2[max_diph];
  unsigned int num_vtx;
  //const unsigned short max_vtx = 60;

  const float sqrt_s = 13.e3;
  const float lumi = ( 5.055026851041 + 1.470474265456 + 7.487318307770 ) * 1.e3;

  for ( unsigned short i=0; i<num_samples; i++ ) {
    Sample s = dsh.sample( i );
    cout << ">>> Processing " << s.name() << endl;
    const float s_weight = ( s.type()!=Sample::kData )
      ? s.cross_section()/s.num_events()*lumi
      : 1.;

    h_presel_mass[i] = new TH1D( Form( "h_presel_mass_%d", i ), "Diphoton mass\\Events\\GeV", 50, 300., 1300. );
    h_presel_pt[i] = new TH1D( Form( "h_presel_pt_%d", i ), "Diphoton p_{T}\\Events\\GeV", 50, 0., 500. );
    h_presel_dphi[i] = new TH1D( Form( "h_presel_dphi_%d", i ), "Diphoton 1-|#Delta#phi/#pi|\\Events", 50, 0., 1. );
    h_presel_xip[i] = new TH1D( Form( "h_presel_xip_%d", i ), "Diphoton #xi^{+}\\Events", 50, 0., 0.5 );
    h_presel_xim[i] = new TH1D( Form( "h_presel_xim_%d", i ), "Diphoton #xi^{-}\\Events", 50, 0., 0.5 );
    h_presel_nvtx[i] = new TH1D( Form( "h_presel_nvtx_%d", i ), "Number of primary vertices\\Events", 60, 0., 60. );

    TTree* tree = s.tree();
    tree->SetBranchAddress( "num_diphoton", &num_diph );
    tree->SetBranchAddress( "diphoton_mass", diph_mass );
    tree->SetBranchAddress( "diphoton_pt", diph_pt );
    tree->SetBranchAddress( "diphoton_dphi", diph_dphi );
    tree->SetBranchAddress( "diphoton_pt1", diph_pt1 );
    tree->SetBranchAddress( "diphoton_eta1", diph_eta1 );
    tree->SetBranchAddress( "diphoton_pt2", diph_pt2 );
    tree->SetBranchAddress( "diphoton_eta2", diph_eta2 );
    tree->SetBranchAddress( "num_vertex", &num_vtx );
    const unsigned long long num_entries = tree->GetEntries();
    for ( unsigned long long j=0; j<num_entries; j++ ) {
      if ( fmod( j*100./num_entries, 10. )==0. ) {
        cout << "   event " << j << " / " << num_entries << endl;
      }
      tree->GetEntry( j );
      for ( unsigned short k=0; k<num_diph; k++ ) {
        h_presel_mass[i]->Fill( diph_mass[k] );
        h_presel_pt[i]->Fill( diph_pt[k] );
        const float acop = 1.-fabs( diph_dphi[k]/TMath::Pi() );
        h_presel_dphi[i]->Fill( acop );
        const float xip = ( diph_pt1[k]*exp( +diph_eta1[k] ) + diph_pt2[k]*exp( +diph_eta2[k] ) ) / sqrt_s,
                    xim = ( diph_pt1[k]*exp( -diph_eta1[k] ) + diph_pt2[k]*exp( -diph_eta2[k] ) ) / sqrt_s;
        h_presel_xip[i]->Fill( xip );
        h_presel_xim[i]->Fill( xim );
      }
      h_presel_nvtx[i]->Fill( num_vtx );
    }

    h_presel_mass[i]->Scale( s_weight );
    h_presel_pt[i]->Scale( s_weight );
    h_presel_dphi[i]->Scale( s_weight );
    h_presel_xip[i]->Scale( s_weight );
    h_presel_xim[i]->Scale( s_weight );
    h_presel_nvtx[i]->Scale( s_weight );
    if ( s.type()==Sample::kData ) {
      hm_presel_mass_data.push_back( make_pair( s.name(), h_presel_mass[i] ) );
      hm_presel_pt_data.push_back( make_pair( s.name(), h_presel_pt[i] ) );
      hm_presel_dphi_data.push_back( make_pair( s.name(), h_presel_dphi[i] ) );
      hm_presel_xip_data.push_back( make_pair( s.name(), h_presel_xip[i] ) );
      hm_presel_xim_data.push_back( make_pair( s.name(), h_presel_xim[i] ) );
      hm_presel_nvtx_data.push_back( make_pair( s.name(), h_presel_nvtx[i] ) );
    }
    else {
      hm_presel_mass_mc.push_back( make_pair( s.name(), h_presel_mass[i] ) );
      hm_presel_pt_mc.push_back( make_pair( s.name(), h_presel_pt[i] ) );
      hm_presel_dphi_mc.push_back( make_pair( s.name(), h_presel_dphi[i] ) );
      hm_presel_xip_mc.push_back( make_pair( s.name(), h_presel_xip[i] ) );
      hm_presel_xim_mc.push_back( make_pair( s.name(), h_presel_xim[i] ) );
      hm_presel_nvtx_mc.push_back( make_pair( s.name(), h_presel_nvtx[i] ) );
    }
  }

  gStyle->SetOptStat( 0 );
  Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/test", "CMS Preliminary 2016, #sqrt{s} = 13 TeV" );
  {
    plt.draw_multiplot( "presel_diphoton_mass", hm_presel_mass_data, hm_presel_mass_mc, false );
    plt.draw_multiplot( "presel_diphoton_mass_logscale", hm_presel_mass_data, hm_presel_mass_mc, true );
    plt.draw_multiplot( "presel_diphoton_pt", hm_presel_pt_data, hm_presel_pt_mc, false );
    plt.draw_multiplot( "presel_diphoton_pt_logscale", hm_presel_pt_data, hm_presel_pt_mc, true );
    plt.draw_multiplot( "presel_diphoton_dphi", hm_presel_dphi_data, hm_presel_dphi_mc );
    plt.draw_multiplot( "presel_diphoton_xip", hm_presel_xip_data, hm_presel_xip_mc );
    plt.draw_multiplot( "presel_diphoton_xim", hm_presel_xim_data, hm_presel_xim_mc );
    plt.draw_multiplot( "presel_diphoton_nvtx", hm_presel_nvtx_data, hm_presel_nvtx_mc );
  }
}
