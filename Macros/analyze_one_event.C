#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "TFile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"

#include "Canvas.h"

void analyze_one_event()
{
  TFile file( "root:///eos/cms/store/group/phys_pps/diphoton/DoubleEG/crab_pickEvents/170320_152646/0000/pickevents_1.root" );
  fwlite::Event ev( &file );

  const char* event_num = "Event 279760:127:127079909";

  const TVector3 diphoton_vertex( 0.0517226, 0.100587, 3.57555 );

  TH1D* h_muon_dist_incl = new TH1D( "muon_dist_incl", "Distance to the track\\Events\\cm?.2f", 100, 0., 10. ),
       *h_electron_dist_incl = ( TH1D* )h_muon_dist_incl->Clone( "electron_dist_incl" );

  for ( ev.toBegin(); !ev.atEnd(); ++ev ) {
    fwlite::Handle< std::vector<reco::GsfElectron> > ele_handle;
    ele_handle.getByLabel( ev, "gedGsfElectrons" );
    const std::vector<reco::GsfElectron>* electrons = ele_handle.ptr();

    fwlite::Handle< std::vector<reco::Muon> > muo_handle;
    muo_handle.getByLabel( ev, "muons" );
    const std::vector<reco::Muon>* muons = muo_handle.ptr();

    fwlite::Handle< std::vector<reco::Vertex> > vtx_handle;
    vtx_handle.getByLabel( ev, "offlinePrimaryVertices" );
    const std::vector<reco::Vertex>* vertices = vtx_handle.ptr();

    for ( unsigned int i=0; i<vertices->size(); ++i ) {
      const reco::Vertex vtx = vertices->at( i );
      const TVector3 vtx_pos( vtx.x(), vtx.y(), vtx.z() );
      if ( (vtx_pos-diphoton_vertex).Perp()>0.001 ) continue;
      cout << " vtx at " << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << ", distance from diphoton vertex: " << (vtx_pos-diphoton_vertex).Perp() << endl;
      cout << "number of tracks associated to this vertex: " << vtx.tracksSize() << endl;
      for ( reco::Vertex::trackRef_iterator trk=vtx.tracks_begin(); trk!=vtx.tracks_end(); ++trk ) {
        const TVector3 track_kin( ( *trk )->px(), ( *trk )->py(), ( *trk )->pz() );
        for ( unsigned short j=0; j<electrons->size(); ++j ) {
          const reco::GsfElectron ele = electrons->at( j );
          const TVector3 ele_kin( ele.gsfTrack()->px(), ele.gsfTrack()->py(), ele.gsfTrack()->pz() );
          h_electron_dist_incl->Fill( ( ele_kin-track_kin ).Perp() );
          cout << "electron distance: " << (ele_kin-track_kin).Perp() << endl;
        }
        for ( unsigned short j=0; j<muons->size(); ++j ) {
          const reco::Muon mu = muons->at( j );
          const TVector3 mu_kin( mu.innerTrack()->px(), mu.innerTrack()->py(), mu.innerTrack()->pz() );
          h_muon_dist_incl->Fill( ( mu_kin-track_kin ).Perp() );
          cout << "muon distance: " << (mu_kin-track_kin).Perp() << endl;
        }
      }
    }
  }

  gStyle->SetOptStat( 0 );
  {
    Canvas c( "distance_track", Form( "%s - diphoton candidate vertex", event_num ) );
    h_electron_dist_incl->Draw( "p" );
    h_electron_dist_incl->SetMarkerStyle( 20 );
    h_muon_dist_incl->Draw( "p same" );
    h_muon_dist_incl->SetMarkerStyle( 24 );
    h_electron_dist_incl->GetYaxis()->SetRangeUser( 0.01, TMath::Max( h_electron_dist_incl->GetMaximum(), h_muon_dist_incl->GetMaximum() )*1.2 );
    c.AddLegendEntry( h_electron_dist_incl, "e^{#pm} GSF tracks" );
    c.AddLegendEntry( h_muon_dist_incl, "#mu^{#pm} tracker tracks" );
    c.Prettify( h_electron_dist_incl );
    c.Save( "png" );
  }
}
