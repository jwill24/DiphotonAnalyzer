#include "Canvas.h"
#include "diproton_candidate.h"
#include "CTPPSAnalysisTools/Reconstruction/interface/LHCConditionsFactory.h"
#include "CTPPSAnalysisTools/Reconstruction/interface/XiReconstructor.h"
#include "CTPPSAnalysisTools/Alignment/interface/AlignmentsFactory.h"

#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"

#include <fstream>

#define MAX_PROTONS 30
#define MAX_DIPH 20
#define MAX_ELE 50
#define MAX_MUON 50
#define MAX_JET 200

using namespace std;

const float max_xi = 0.25;
struct track_t {
  track_t() : xi( 0. ), err_xi( 0. ), x( 0. ), y( 0. ) {}
  track_t( float xi, float err_xi, float x, float y ) : xi( xi ), err_xi( err_xi ), x( x ), y( y ) {}
  float xi, err_xi, x, y;
};

bool is_matched( int n_sigma, float x1, float x2, float err_x1, float err_x2 )
{
  const double combined_error = sqrt( err_x1*err_x1 + err_x2*err_x2 );
  const double delta = fabs( x1-x2 );

 return ( delta/combined_error<=n_sigma );
}
void plot_matching( const char* name, TGraphErrors&, TGraphErrors&, TGraphErrors&, TGraphErrors&, double min, double max, double xleg=0.5 );
vector<pair<float,float> > merge_nearfar( const vector<track_t>& near_tracks, const vector<track_t>& far_tracks, float xdiff_cut=0.01 );

map<string,float> pots_accept = { { "45N", 0.033 }, { "45F", 0.024 }, { "56N", 0.050 }, { "56F", 0.037 } };

void massrap_matcher()
{

  ofstream myfile;
  myfile.open ("events.txt");

  TFile f( "/afs/cern.ch/work/j/juwillia/CMSSW_9_4_5_cand1/src/DiphotonAnalyzer/ntp_2017BCDEF.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  ctpps::LHCConditionsFactory cond_fac;
  ostringstream cond_file1_path, cond_file2_path;
  cond_file1_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Reconstruction/data/2017/xangle_tillTS2.csv";
  cond_file2_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Reconstruction/data/2017/xangle_afterTS2.csv";
  cond_fac.feedConditions( cond_file1_path.str().c_str() );
  cond_fac.feedConditions( cond_file2_path.str().c_str() );
  
  ctpps::XiReconstructor reco;
  ostringstream disp_file_path;
  disp_file_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Reconstruction/data/2017/dispersions.txt";
  reco.feedDispersions( disp_file_path.str().c_str() );

  ctpps::AlignmentsFactory align_fac;
  ostringstream align_file_path;
  align_file_path << getenv( "CMSSW_BASE" ) << "/src/CTPPSAnalysisTools/Alignment/data/2017/alignments_30jan2017.txt";
  align_fac.feedAlignments( align_file_path.str().c_str() );
  
  const float num_sigma = 2.0;

  TreeEvent ev;
  ev.attach( tr, true );

  TGraphErrors gr_mass_massrapmatch, gr_mass_massmatch, gr_mass_rapmatch, gr_mass_nomatch;
  TGraphErrors gr_rap_massrapmatch, gr_rap_massmatch, gr_rap_rapmatch, gr_rap_nomatch;

  unsigned int num_massmatch = 0, num_rapmatch = 0, num_massrapmatch = 0, num_nomatch = 0;

  TCanvas* c1 = new TCanvas("c1", "", 500, 500 );
  TCanvas* c2 = new TCanvas("c2", "", 500, 500 );
  TCanvas* c3 = new TCanvas("c3", "", 500, 500 );
  TCanvas* c4 = new TCanvas("c4", "", 500, 500 );
  TCanvas* c5 = new TCanvas("c5", "", 500, 500 );
  TCanvas* c6 = new TCanvas("c6", "", 500, 500 );
  TCanvas* c7 = new TCanvas("c7", "", 500, 500 );

  TH1D* h_eta = new TH1D( "h_eta", "#eta of matching events", 100, -2.5, 2.5 );
  TH1D* h_dphi = new TH1D( "h_dphi", "#delta#phi of matching events", 100, -4., 4. );
  TH1D* h_acop = new TH1D( "h_acop", "acoplanarity", 100, 0., 1. );
  TH1D* h_pt = new TH1D( "h_pt", "p#_{T} of matching events", 100, 50., 500. );
  TH1D* h_numjets = new TH1D( "h_humjets", "num jets in matching events", 100, 0., 50. );
  TH1D* h_numleps = new TH1D( "h_numleps", "num leptons in matching events", 100, 0., 50. );
  TH1D* h_mass = new TH1D( "h_mass", "Diphoton Mass", 100., 0., 2500. );
  TH1D* h_mass_all = new TH1D( "mass_all", "Diproton missing mass@@Events@@GeV", 12, 200., 2000. );
  TH1D* h_rap_all = new TH1D( "rap_all", "Diproton rapidity@@Events", 20, -1., 1. );
  TH2D* h_pu_v_run = new TH2D( "pu_v_run", "Pile Up as a function of Run number", 100, 294645, 302663, 100, 0, 100 );

  const unsigned long long num_events = tr->GetEntriesFast();

  int count = 0;
  for ( unsigned long long i = 0; i < num_events; ++i ) {
    tr->GetEntry( i );

    // first loop to identify the tracks and their respective pot

    vector<track_t> xi_45n, xi_45f, xi_56n, xi_56f;

    double side_event[33] = {0};
    double xi_event[33] = {0};

    for ( unsigned short j = 0; j < ev.num_proton_track; ++j ) {
      const unsigned short pot_id = 100*ev.proton_track_side[j] + 10*ev.proton_track_station[j] + ev.proton_track_pot[j];
      const unsigned short raw_id = 100*ev.proton_track_side[j] + 10*ev.proton_track_station[j] + ev.proton_track_pot[j];
      
      if ( pot_id != 3 && pot_id != 23 && pot_id != 103 && pot_id != 123 ){
	count = count + 1;
	cout << "run: " << ev.run_id << " num tracks: " << ev.num_proton_track << endl; 
	cout << "side: " << ev.proton_track_side[j] << " station: " << ev.proton_track_station[j] << " pot: " << ev.proton_track_pot[j] << endl;
	cout << "" << endl;
	continue;
      }
      
      //----- reconstruct the kinematics
      double xi = 0.;
      double xi_err = 0.0;
      double trk_x_corr = 0.0;

      const auto cond = cond_fac.get( edm::EventID( ev.run_id, ev.lumisection, ev.event_number ) );
      const ctpps::alignment_t align = align_fac.get( edm::EventID( ev.run_id, ev.lumisection, ev.event_number ), pot_id );
      double align_quant = align.x_align;;
      double xangle = cond.crossing_angle;

      
      //----- associate each track to a RP
      if ( ev.proton_track_side[j] == 0 && ev.proton_track_station[j] == 0 ) {
	if (ev.fill_number < 303718 ) trk_x_corr = ev.proton_track_x[j] + align_quant;
	if (ev.fill_number > 303718 ) trk_x_corr = ev.proton_track_x[j]*cos( 16*M_PI/180. ) + ev.proton_track_y[j]*sin( 16*M_PI/180. ) + align_quant;
	reco.reconstruct( xangle, raw_id, trk_x_corr, xi, xi_err );
	xi_err = xi*0.10;
	xi_45n.emplace_back( xi, xi_err, ev.proton_track_x[j]+align_quant, ev.proton_track_y[j] );
      }
      else if ( ev.proton_track_side[j] == 0 && ev.proton_track_station[j] == 2 ) {
	trk_x_corr = ev.proton_track_x[j] + align_quant;
	reco.reconstruct( xangle, raw_id, trk_x_corr, xi, xi_err );
	xi_err = xi*0.10;
	xi_45f.emplace_back( xi, xi_err, ev.proton_track_x[j]+align_quant, ev.proton_track_y[j] );
      }
      else if ( ev.proton_track_side[j] == 1 && ev.proton_track_station[j] == 0 ) {
	if (ev.fill_number< 303718 ) trk_x_corr = ev.proton_track_x[j] + align_quant;
	if (ev.fill_number > 303718 ) trk_x_corr = ev.proton_track_x[j]*cos( 16*M_PI/180. ) + ev.proton_track_y[j]*sin( 16*M_PI/180. ) + align_quant;
	reco.reconstruct( xangle, raw_id, trk_x_corr, xi, xi_err );
	xi_err = xi*0.10;
	xi_56n.emplace_back( xi, xi_err, ev.proton_track_x[j]+align_quant, ev.proton_track_y[j] );
      }
      else if ( ev.proton_track_side[j] == 1 && ev.proton_track_station[j] == 2 ) {
	trk_x_corr = ev.proton_track_x[j] + align_quant;
	reco.reconstruct( xangle, raw_id, trk_x_corr, xi, xi_err );
	xi_err = xi*0.10;
	xi_56f.emplace_back( xi, xi_err, ev.proton_track_x[j]+align_quant, ev.proton_track_y[j] );
      }

      side_event[j] = ev.proton_track_side[j];
      xi_event[j] = xi;
    }

    //----- merge 2 tracks in one if N-F pot content is similar
    vector<pair<float,float> > xi_45, xi_56;

    const float xdiff_cut = 0.01;

    xi_45 = merge_nearfar( xi_45n, xi_45f, xdiff_cut );
    xi_56 = merge_nearfar( xi_56n, xi_56f, xdiff_cut );
    

    //---- identify the diproton candidates

    vector<diproton_candidate_t> candidates;
    //for ( const auto trk45 : xi_45 ) {
       //for ( const auto trk56 : xi_56 ) {	
	//candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
       
	//arr_45.emplace_back( trk45.first );
	//arr_56.emplace_back( trk56.first );
    
    float_t max_45 = 0.0;
    float_t max_56 = 0.0;
    for (unsigned int i = 0; i < xi_45.size(); i++) {
      for (unsigned int j = 0; j < xi_56.size(); j++) {
	if ( xi_45[i].first > max_45 ) max_45 = xi_45[i].first;
	if ( xi_56[j].first > max_56 ) max_56 = xi_56[j].first;
      }
    }
    
    float_t xi_45m = max_45;
    float_t xi_56m = max_56;
    
    
    candidates.emplace_back( xi_45m, xi_45m*0.1, xi_56m, xi_56m*0.1 );
    
    //cout << candidates.size() << " diproton candidate(s) in total!" << endl;

    //----- identify the diphoton candidates

    for ( unsigned short j = 0; j < ev.num_diphoton; ++j ) {

      //----- photon quality cuts


      if ( ev.diphoton_pt1[j] < 75. ) continue;
      if ( ev.diphoton_pt2[j] < 75. ) continue;
      if ( ev.diphoton_eta1[j] > 2.5 || ( ev.diphoton_eta1[j] > 1.4442 && ev.diphoton_eta1[j] < 1.566 ) ) continue;
      if ( ev.diphoton_eta2[j] > 2.5 || ( ev.diphoton_eta2[j] > 1.4442 && ev.diphoton_eta2[j] < 1.566 ) ) continue;
      if ( ev.diphoton_r91[j] < 0.94 ) continue;
      if ( ev.diphoton_r92[j] < 0.94 ) continue;
      if ( ev.diphoton_mass[j] < 350. ) continue;

      h_acop->Fill( 1.-fabs( ev.diphoton_dphi[j] )/M_PI );

      //----- back-to-back photons

      if ( 1.-fabs( ev.diphoton_dphi[j] )/M_PI > 0.005 ) continue;

      h_mass->Fill( ev.diphoton_mass[j] );
      
      const float xip = ( ev.diphoton_pt1[j]*exp( +ev.diphoton_eta1[j] ) + ev.diphoton_pt2[j]*exp( +ev.diphoton_eta2[j] ) ) / sqrt_s,
                  xim = ( ev.diphoton_pt1[j]*exp( -ev.diphoton_eta1[j] ) + ev.diphoton_pt2[j]*exp( -ev.diphoton_eta2[j] ) ) / sqrt_s;


      //----- search for associated jets

      TLorentzVector pho1, pho2;
      pho1.SetPtEtaPhiM( ev.diphoton_pt1[j], ev.diphoton_eta1[j], ev.diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( ev.diphoton_pt2[j], ev.diphoton_eta2[j], ev.diphoton_phi2[j], 0. );
      TLorentzVector cms = pho1+pho2;


      //----- another batch of "exclusivity" cuts

      /*if ( num_close_ele > 0 ) continue;
      if ( num_close_muon > 0 ) continue;*/
      //if ( num_associated_jet > 0 ) continue;

      //----- reconstruct the energy loss from central system

      float diphoton_mass_error = ev.diphoton_mass[j]*0.02;
      float diphoton_rapidity_error = fabs( ev.diphoton_rapidity[j] )*0.061;


      for ( const auto cand : candidates ) {
        h_mass_all->Fill( cand.mass() );
        h_rap_all->Fill( cand.rapidity() );
        bool mass_match = is_matched( num_sigma, cms.M(), cand.mass(), diphoton_mass_error, cand.mass_error() );
        bool rap_match = is_matched( num_sigma, cms.Rapidity(), cand.rapidity(), diphoton_rapidity_error, cand.rapidity_error() );
	
	// print statements for all events
	if ( cand.mass() > 100. && cand.mass() < 4000. ) {
	  myfile << "R:L:E ---> " << ev.run_id << ":" << ev.lumisection << ":" << ev.event_number << endl;                
	  myfile << "masses: central system: " << cms.M() << ", diphoton: " << ev.diphoton_mass[j] << " +/- " << diphoton_mass_error << ", diproton: " << cand.mass() << " +/- " << cand.mass_error() << endl; 
	  myfile << "rapidities: central system: " << cms.Rapidity() << ", diphoton: " << ev.diphoton_rapidity[j] << " +/- " << diphoton_rapidity_error << ", diproton: " << cand.rapidity() << " +/- " << cand.rapidity_error() << endl;        
	  myfile << "xip:" << xip << " xim: " << xim << endl;
	  for ( int t=0; t < 20; t++ ) {
	    if ( xi_event[t] > 0.0 ) {
	      myfile << "side: " << side_event[t] << " xi: " << xi_event[t] << endl;
	    }
	  }
	}
	
        if ( mass_match && rap_match ) {

	    cout << "@@@ DOUBLE TAGGING" << endl;
	    cout << "masses: central system: " << cms.M() << ", diphoton: " << ev.diphoton_mass[j] << " +/- " << diphoton_mass_error << ", diproton: " << cand.mass() << " +/- " << cand.mass_error() << endl;
	    cout << "rapidities: central system: " << cms.Rapidity() << ", diphoton: " << ev.diphoton_rapidity[j] << " +/- " << diphoton_rapidity_error << ", diproton: " << cand.rapidity() << " +/- " << cand.rapidity_error() << endl;
	    cout << "R:L:E ---> " << ev.run_id << ":" << ev.lumisection << ":" << ev.event_number << endl;
	    cout << "xip:" << xip << " xim: " << xim << " N tracks left: " << xi_45.size() << " N tracks right: " << xi_56.size() << endl;

          gr_mass_massrapmatch.SetPoint( num_massrapmatch, cand.mass(), cms.M() );
          gr_mass_massrapmatch.SetPointError( num_massrapmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_massrapmatch.SetPoint( num_massrapmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_massrapmatch.SetPointError( num_massrapmatch, cand.rapidity_error(), diphoton_rapidity_error );
	  
	  h_eta->Fill( ev.diphoton_eta1[j] );
	  h_eta->Fill( ev.diphoton_eta2[j] );
	  h_dphi->Fill( ev.diphoton_dphi[j] );
	  h_pt->Fill( ev.diphoton_pt1[j] );
	  h_pt->Fill( ev.diphoton_pt2[j] );
	  h_numjets->Fill( ev.num_jet );
	  h_numleps->Fill( ev.num_muon + ev.num_electron );
	  
          num_massrapmatch++;
        }
        else if ( mass_match ) { // only match in mass
          gr_mass_massmatch.SetPoint( num_massmatch, cand.mass(), cms.M() );
          gr_mass_massmatch.SetPointError( num_massmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_massmatch.SetPoint( num_massmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_massmatch.SetPointError( num_massmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_massmatch++;
        }
        else if ( rap_match ) { // only match in rapidity
          gr_mass_rapmatch.SetPoint( num_rapmatch, cand.mass(), cms.M() );
          gr_mass_rapmatch.SetPointError( num_rapmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_rapmatch.SetPoint( num_rapmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_rapmatch.SetPointError( num_rapmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_rapmatch++;
        }
        else { // no matching at all
          gr_mass_nomatch.SetPoint( num_nomatch, cand.mass(), cms.M() );
          gr_mass_nomatch.SetPointError( num_nomatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_nomatch.SetPoint( num_nomatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_nomatch.SetPointError( num_nomatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_nomatch++;
        }
      } // End candidates loop
    } // End diphoton loop
  } // End events loop

  cout << "Number of events without a real pot id: " << count << endl;
  myfile.close();
  
  cout << "in plot:\n\t" << "not matching: " << num_nomatch << "\n\tmass match: " << num_massmatch << "\n\trap match: " << num_rapmatch << "\n\tboth match: " << num_massrapmatch << endl;

  //----- plotting part

  c1->cd();
  h_mass->Draw();
  c1->SaveAs("/eos/user/j/juwillia/www/2017data/mass.pdf");

  c2->cd();
  h_eta->Draw();
  c2->SaveAs("/eos/user/j/juwillia/www/2017data/eta.pdf");

  c3->cd();
  h_pt->Draw();
  c3->SaveAs("/eos/user/j/juwillia/www/2017data/pt.pdf");

  c4->cd();
  h_dphi->Draw();
  c4->SaveAs("/eos/user/j/juwillia/www/2017data/dphi.pdf");

  c5->cd();
  h_numjets->Draw();
  c5->SaveAs("/eos/user/j/juwillia/www/2017data/numjets.pdf");

  c6->cd();
  h_numleps->Draw();
  c6->SaveAs("/eos/user/j/juwillia/www/2017data/numleps.pdf");
 
  c7->cd();
  h_acop->Draw();
  c7->SaveAs("/eos/user/j/juwillia/www/2017data/acoplanarity.pdf");


  gr_mass_massrapmatch.SetTitle( "Diproton system missing mass (GeV)@@Diphoton mass (GeV)" );
  gr_rap_massrapmatch.SetTitle( "Diproton system rapidity@@Diphoton rapidity" );

  plot_matching( "2d_massmatch", gr_mass_nomatch, gr_mass_rapmatch, gr_mass_massmatch, gr_mass_massrapmatch, 300., 1800. );
  plot_matching( "2d_rapmatch", gr_rap_nomatch, gr_rap_rapmatch, gr_rap_massmatch, gr_rap_massrapmatch, -3., 3., 0.15 );

  {
    Canvas c( "1d_massmatch", "CMS+TOTEM Preliminary 2017, #sqrt{s} = 13 TeV, L = 41.199 fb^{-1}" );
    h_mass_all->Sumw2();
    h_mass_all->Draw();
    c.Prettify( h_mass_all );
    //"/eos/user/j/juwillia/www/2017data"
    ///afs/cern.ch/work/j/juwillia/CMSSW_9_4_5_cand1/src/DiphotonAnalyzer
    c.Save( "pdf,png", "/eos/user/j/juwillia/www/2017data" );
  }
  {
    Canvas c( "1d_rapmatch", "CMS+TOTEM Preliminary 2017, #sqrt{s} = 13 TeV, L = 41.199 fb^{-1}" );
    h_rap_all->Sumw2();
    h_rap_all->Draw();
    c.Prettify( h_rap_all );
    c.Save( "pdf,png", "/eos/user/j/juwillia/www/2017data" );
  }

}

void plot_matching( const char* name, TGraphErrors& gr_nomatch, TGraphErrors& gr_rapmatch, TGraphErrors& gr_massmatch, TGraphErrors& gr_massrapmatch, double min, double max, double xleg )
{
  //Pre-TS2 = 17.541(94), Post-TS2 = 23.657(45), Total = 41.19939
  Canvas c( name, "CMS+TOTEM Preliminary 2017, #sqrt{s} = 13 TeV, L = 41.19939 fb^{-1}" );
  /*auto tmp = new TH2D( Form( "tmp_%s", name ), gr_massrapmatch.GetTitle(), 2, min, max, 2, min, max );
  tmp->Draw();*/
  auto diag = new TF1( "diag", "x", min, max );
  diag->SetLineColor( kGray );

  c.SetLegendX1( xleg );
  c.SetLegendY1( 0.7 );
  TMultiGraph mg;
  gr_nomatch.SetMarkerStyle( 24 );
  gr_rapmatch.SetMarkerStyle( 22 );
  gr_rapmatch.SetMarkerSize( 1.2 );
  gr_massmatch.SetMarkerStyle( 20 );
  gr_massmatch.SetMarkerSize( 1.0 );
  gr_massrapmatch.SetMarkerStyle( 21 );
  gr_massrapmatch.SetMarkerColor( kBlue+1 );
  gr_massrapmatch.SetMarkerSize( 1.3 );
  gr_massmatch.SetMarkerColor( kRed+1 );
  gr_rapmatch.SetMarkerColor( kGreen+2 );
  if ( gr_massrapmatch.GetN() > 0 ) c.AddLegendEntry( &gr_massrapmatch, "2D matching", "p" );
  c.AddLegendEntry( &gr_massmatch, "Mass matching", "lp" );
  c.AddLegendEntry( &gr_rapmatch, "Rapidity matching", "lp" );
  c.AddLegendEntry( &gr_nomatch, "No matching", "lp" );
  mg.Add( &gr_nomatch );
  mg.Add( &gr_massmatch );
  mg.Add( &gr_rapmatch );
  mg.Add( &gr_massrapmatch );

  diag->Draw( "l" );
  mg.Draw( "p" );

  //PaveText pt( 0.15, 0.8, 0.4, 0.95 );
  PaveText pt( xleg, 0.8, xleg+0.25, 0.95 );
  pt.SetTextAlign( kHAlignLeft+kVAlignBottom );
  pt.AddText( "1-|#Delta#phi(#gamma,#gamma)/#pi| < 0.005" );
  pt.Draw();

  diag->GetHistogram()->SetTitle( gr_massrapmatch.GetTitle() );
  c.Prettify( diag->GetHistogram() );
  //c.Prettify( tmp );
  //mg.GetXaxis()->SetLimits( min, max );
  diag->GetHistogram()->GetYaxis()->SetRangeUser( min, max );
  c.Save( "pdf,png", "/eos/user/j/juwillia/www/2017data" );
}

vector<pair<float,float> > merge_nearfar( const vector<track_t>& near_tracks, const vector<track_t>& far_tracks, float xdiff_cut ) {
  vector<pair<float,float> > out;
  set<unsigned short> matched_far_ids;

  //--- first loop to extract near tracks with matching
  for ( const auto& near : near_tracks ) {
    float min_xidiff = 999.;
    int matched_far_id = -1;
    pair<track_t,track_t> assoc;
    unsigned short id_far = 0;
    for ( const auto& far : far_tracks ) {
      if ( fabs( near.x-far.x ) < xdiff_cut ) {
        //--- near-far track association candidate
        float xidiff = fabs( near.xi-far.xi );
        if ( xidiff < min_xidiff ) {
          assoc = make_pair( near, far );
          min_xidiff = xidiff;
          matched_far_id = id_far;
        }
      }
      id_far++;
    }
    if ( min_xidiff < 999. ) {
      out.emplace_back( assoc.second.xi, assoc.second.err_xi ); // store the far tracks info
      matched_far_ids.insert( matched_far_id );
    }
    //--- store the near track if no matching found with far track
    else out.emplace_back( near.xi, near.err_xi );
  }
  //--- second loop to add the remaining, unmatched far tracks
  unsigned short id_far = 0;
  for ( const auto& far : far_tracks ) {
    if ( matched_far_ids.count( id_far ) == 0 ) {
      //--- discard far track if a mapping is already found
      out.emplace_back( far.xi, far.err_xi );
    }
  }
  return out;
}
