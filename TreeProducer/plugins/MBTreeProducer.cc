// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/TreeProducer
// Class:      MBTreeProducer
// 
/**\class MBTreeProducer MBTreeProducer.cc DiphotonAnalyzer/TreeProducer/plugins/MBTreeProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Tue, 13 Sep 2016 03:57:43 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "DiphotonAnalyzer/TreeProducer/interface/SelectionUtils.h"
#include "DiphotonAnalyzer/TreeProducer/interface/FillNumberLUTHandler.h"

#include "TFile.h"
#include "TTree.h"

//
// class declaration
//
#define MAX_PROTON_TRK 20
#define MAX_VERTEX 100

class MBTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit MBTreeProducer( const edm::ParameterSet& );
    ~MBTreeProducer();

    static void fillDescriptions( edm::ConfigurationDescriptions& );


  private:
    virtual void beginJob() override;
    virtual void beginRun( const edm::Run&, const edm::EventSetup& );
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    void clearTree();
    
    edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > totemRPTracksToken_;
    edm::EDGetTokenT< edm::View<reco::Vertex> > vtxToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

    std::string filename_;

    std::unique_ptr<CTPPSAlCa::FillNumberLUTHandler> fillLUTHandler_;

    TFile* file_;
    TTree* tree_;

    // --- tree components ---

    unsigned int fBX, fFill, fRun, fLumiSection;
    unsigned long long fEventNum;

    unsigned int fProtonTrackNum;
    float fProtonTrackX[MAX_PROTON_TRK], fProtonTrackY[MAX_PROTON_TRK];
    float fProtonTrackTX[MAX_PROTON_TRK], fProtonTrackTY[MAX_PROTON_TRK];
    float fProtonTrackChi2[MAX_PROTON_TRK], fProtonTrackNormChi2[MAX_PROTON_TRK];
    unsigned int fProtonTrackSide[MAX_PROTON_TRK], fProtonTrackPot[MAX_PROTON_TRK];

    unsigned int fVertexNum;
    float fVertexX[MAX_VERTEX], fVertexY[MAX_VERTEX], fVertexZ[MAX_VERTEX];
    unsigned int fVertexTracks[MAX_VERTEX];
    unsigned int fVertexTracksWght0p75[MAX_VERTEX], fVertexTracksWght0p90[MAX_VERTEX], fVertexTracksWght0p95[MAX_VERTEX];
    float fBSX0, fBSY0, fBSZ0, fBSsigmaZ, fBSdxdz, fBSbeamWidthX, fBSbeamWidthY;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MBTreeProducer::MBTreeProducer( const edm::ParameterSet& iConfig ) :
  totemRPTracksToken_ ( consumes< edm::DetSetVector<TotemRPLocalTrack> >  ( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  vtxToken_           ( consumes< edm::View<reco::Vertex> >               ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
  beamSpotToken_      ( consumes<reco::BeamSpot>                          ( iConfig.getParameter<edm::InputTag>( "beamSpotLabel" ) ) ),
  filename_           ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  fillLUTHandler_ ( new CTPPSAlCa::FillNumberLUTHandler( iConfig.getParameter<edm::FileInPath>( "fillNumLUTFile" ).fullPath().c_str() ) ),
  file_( 0 ), tree_( 0 )
{

  file_ = new TFile( filename_.c_str(), "recreate" );
  file_->cd();

  tree_ = new TTree( "ntp", "diphoton ntuple" );
}


MBTreeProducer::~MBTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  if ( file_ ) {
    file_->Write();
    file_->Close();
    delete file_;
  }
  if ( tree_ ) delete tree_;
}


//
// member functions
//

void
MBTreeProducer::clearTree()
{
  fBX = fRun = fLumiSection = fEventNum = 0;
  fProtonTrackNum = 0;
  for ( unsigned int i=0; i<MAX_PROTON_TRK; i++ ) {
    fProtonTrackX[i] = fProtonTrackY[i] = -1.;
    fProtonTrackTX[i] = fProtonTrackTY[i] = -1.;
    fProtonTrackChi2[i] = fProtonTrackNormChi2[i] = -1.;
    fProtonTrackSide[i] = 2; //invalid
    fProtonTrackPot[i] = 0;
  }

  fVertexNum = 0;
  for ( unsigned int i=0; i<MAX_VERTEX; i++ ) {
    fVertexX[i] = fVertexY[i] = fVertexZ[i] = -999.;
    fVertexTracks[i] = fVertexTracksWght0p75[i] = fVertexTracksWght0p90[i] = fVertexTracksWght0p95[i] = 0;
  }
  fBSX0 = fBSY0 = fBSZ0 = fBSsigmaZ = fBSdxdz = fBSbeamWidthX = fBSbeamWidthY = -999.;

}

// ------------ method called for each event  ------------
void
MBTreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  clearTree();

  // Run and BX information
  fBX = iEvent.bunchCrossing();
  fRun = iEvent.id().run();
  fLumiSection = iEvent.luminosityBlock();
  fEventNum = iEvent.id().event();

  // get the fill number from the run id <-> fill number LUT
  fFill = ( fillLUTHandler_ ) ? fillLUTHandler_->getFillNumber( iEvent.id().run() ) : 1;

  //----- beam spot -----

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken( beamSpotToken_, beamSpot );
  if ( beamSpot.isValid() ) {
    fBSX0 = beamSpot->x0();
    fBSY0 = beamSpot->y0();
    fBSZ0 = beamSpot->z0();
    fBSsigmaZ = beamSpot->sigmaZ();
    fBSdxdz = beamSpot->dxdz();
    fBSbeamWidthX = beamSpot->BeamWidthX();
    fBSbeamWidthY = beamSpot->BeamWidthY();
  }

  //----- forward RP tracks -----

  edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rpLocalTracks;
  iEvent.getByToken( totemRPTracksToken_, rpLocalTracks );

  fProtonTrackNum = 0;
  for ( edm::DetSetVector<TotemRPLocalTrack>::const_iterator it=rpLocalTracks->begin(); it!=rpLocalTracks->end(); it++ ) {
    const TotemRPDetId detid( TotemRPDetId::decToRawId( it->detId()*10 ) );
    const unsigned short side = detid.arm(),
                         pot = detid.romanPot();

    for ( edm::DetSet<TotemRPLocalTrack>::const_iterator trk=it->begin(); trk!=it->end(); trk++ ) {
      if ( !trk->isValid() ) { continue; }

      fProtonTrackX[fProtonTrackNum] = trk->getX0() * 1.e-3; // store in m
      fProtonTrackY[fProtonTrackNum] = trk->getY0() * 1.e-3; // store in m
      fProtonTrackSide[fProtonTrackNum] = side; // 0 = left (45) ; 1 = right (56)
      fProtonTrackPot[fProtonTrackNum] = pot; // 2 = 210m ; 3 = 220m

      fProtonTrackChi2[fProtonTrackNum] = trk->getChiSquared();
      fProtonTrackNormChi2[fProtonTrackNum] = trk->getChiSquaredOverNDF();

      fProtonTrackNum++;
    }
  }

  //----- vertexing information -----

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );

  fVertexNum = 0;
  for ( unsigned int i=0; i<vertices->size() && i<MAX_VERTEX; i++ ) {
    const edm::Ptr<reco::Vertex> vtx = vertices->ptrAt( i );
    if ( fVertexNum>=MAX_VERTEX ) {
      edm::LogError("MBTreeProducer") << "More vertices than expected in this event (" << fVertexNum << ">MAX_VERTEX=" << MAX_VERTEX << "). Increase MAX_VERTEX for safety";
    }
    if ( !vtx->isValid() ) continue;

    fVertexX[fVertexNum] = vtx->x();
    fVertexY[fVertexNum] = vtx->y();
    fVertexZ[fVertexNum] = vtx->z();

    fVertexTracks[fVertexNum] = vtx->nTracks();
    fVertexTracksWght0p75[fVertexNum] = vtx->nTracks( 0.75 );
    fVertexTracksWght0p90[fVertexNum] = vtx->nTracks( 0.90 );
    fVertexTracksWght0p95[fVertexNum] = vtx->nTracks( 0.95 );
    fVertexNum++;
  }

  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
MBTreeProducer::beginJob()
{
  tree_->Branch( "run_id", &fRun, "run_id/i");
  tree_->Branch( "fill_number", &fFill, "fill_number/i");
  tree_->Branch( "lumisection", &fLumiSection, "lumisection/i");
  tree_->Branch( "bunch_crossing", &fBX, "bunch_crossing/i");
  tree_->Branch( "event_number", &fEventNum, "event_number/l");

  tree_->Branch( "num_proton_track", &fProtonTrackNum, "num_proton_track/i" );
  tree_->Branch( "proton_track_x", fProtonTrackX, "proton_track_x[num_proton_track]/F" );
  tree_->Branch( "proton_track_y", fProtonTrackY, "proton_track_y[num_proton_track]/F" );
  tree_->Branch( "proton_track_tx", fProtonTrackTX, "proton_track_tx[num_proton_track]/F" );
  tree_->Branch( "proton_track_ty", fProtonTrackTY, "proton_track_ty[num_proton_track]/F" );
  tree_->Branch( "proton_track_side", fProtonTrackSide, "proton_track_side[num_proton_track]/i" );
  tree_->Branch( "proton_track_chi2", fProtonTrackChi2, "proton_track_chi2[num_proton_track]/F" );
  tree_->Branch( "proton_track_normchi2", fProtonTrackNormChi2, "proton_track_normchi2[num_proton_track]/F" );
  tree_->Branch( "proton_track_pot", fProtonTrackPot, "proton_track_pot[num_proton_track]/i" );

  tree_->Branch( "num_vertex", &fVertexNum, "num_vertex/i" );
  tree_->Branch( "vertex_x", fVertexX, "vertex_x[num_vertex]/F" );
  tree_->Branch( "vertex_y", fVertexY, "vertex_y[num_vertex]/F" );
  tree_->Branch( "vertex_z", fVertexZ, "vertex_z[num_vertex]/F" );
  tree_->Branch( "vertex_tracks", fVertexTracks, "vertex_tracks[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_wgt0p75", fVertexTracksWght0p75, "vertex_tracks_wgt0p75[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_wgt0p90", fVertexTracksWght0p90, "vertex_tracks_wgt0p90[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_wgt0p95", fVertexTracksWght0p95, "vertex_tracks_wgt0p95[num_vertex]/i" );

  tree_->Branch( "bs_x0", &fBSX0, "bs_x0/F" );
  tree_->Branch( "bs_y0", &fBSY0, "bs_y0/F" );
  tree_->Branch( "bs_z0", &fBSZ0, "bs_z0/F" );
  tree_->Branch( "bs_sigma_z", &fBSsigmaZ, "bs_sigma_z/F" );
  tree_->Branch( "bs_dxdz", &fBSdxdz, "bs_dxdz/F" );
  tree_->Branch( "bs_beam_width_x", &fBSbeamWidthX, "bs_beam_width_x/F" );
  tree_->Branch( "bs_beam_width_y", &fBSbeamWidthY, "bs_beam_width_y/F" );
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MBTreeProducer::endJob() 
{}

void
MBTreeProducer::beginRun( const edm::Run&, const edm::EventSetup& )
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MBTreeProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( MBTreeProducer );
