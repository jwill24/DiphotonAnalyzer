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

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "DiphotonAnalyzer/TreeProducer/interface/MBTreeEvent.h"

#include "TFile.h"
#include "TTree.h"

class MBTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit MBTreeProducer( const edm::ParameterSet& );
    ~MBTreeProducer();

    static void fillDescriptions( edm::ConfigurationDescriptions& );


  private:
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack> > totemRPTracksToken_;

    std::string filename_;

    TFile* file_;
    TTree* tree_;
    MBTreeEvent ev_;
};

MBTreeProducer::MBTreeProducer( const edm::ParameterSet& iConfig ) :
  totemRPTracksToken_ ( consumes<edm::DetSetVector<TotemRPLocalTrack> >  ( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  filename_           ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  file_( 0 ), tree_( 0 )
{

  file_ = new TFile( filename_.c_str(), "recreate" );
  file_->cd();

  tree_ = new TTree( "ntp", "diphoton ntuple" );
  ev_.create( tree_ );
}


MBTreeProducer::~MBTreeProducer()
{
  if ( file_ ) {
    file_->Write();
    file_->Close();
    delete file_;
  }
  if ( tree_ ) delete tree_;
}

void
MBTreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  ev_.clear();

  // Run and BX information
  ev_.bunch_crossing = iEvent.bunchCrossing();
  ev_.run_id = iEvent.id().run();
  ev_.lumisection = iEvent.luminosityBlock();
  ev_.event_number = iEvent.id().event();

  //----- forward RP tracks -----

  edm::Handle<edm::DetSetVector<TotemRPLocalTrack> > rpLocalTracks;
  iEvent.getByToken( totemRPTracksToken_, rpLocalTracks );

  ev_.num_strips_track = 0;
  for ( const auto& dsv : *rpLocalTracks) {
    const TotemRPDetId detid( dsv.detId()*10 );
    const unsigned short arm = detid.arm(), pot = detid.rp();

    for ( const auto& trk : dsv ) {
      if ( !trk.isValid() ) continue;

      ev_.strips_track_x[ev_.num_strips_track] = trk.getX0() * 1.e-3; // store in m
      ev_.strips_track_y[ev_.num_strips_track] = trk.getY0() * 1.e-3; // store in m
      ev_.strips_track_arm[ev_.num_strips_track] = arm; // 0 = left (45) ; 1 = right (56)
      ev_.strips_track_pot[ev_.num_strips_track] = pot; // 2 = 210m ; 3 = 220m

      ev_.strips_track_chi2[ev_.num_strips_track] = trk.getChiSquared();
      ev_.strips_track_normchi2[ev_.num_strips_track] = trk.getChiSquaredOverNDF();

      ev_.num_strips_track++;
    }
  }
  if ( ev_.num_strips_track == 0 ) return; // do not store in the tree if no valid tracks are found

  tree_->Fill();
}

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
