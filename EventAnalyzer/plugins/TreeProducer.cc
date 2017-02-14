// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/EventAnalyzer
// Class:      TreeProducer
// 
/**\class TreeProducer TreeProducer.cc DiphotonAnalyzer/EventAnalyzer/plugins/TreeProducer.cc

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

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
//                               JW
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Jet.h"
//

#include "DiphotonAnalyzer/EventAnalyzer/interface/SelectionUtils.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/XiInterpolator.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/FillNumberLUTHandler.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/AlignmentLUTHandler.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//

#define MAX_PROTON_TRK 15
#define MAX_DIPHOTON 5
//                               JW
#define MAX_ELECTRON 10
#define MAX_MUON 10
#define MAX_JET 10
//

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit TreeProducer(const edm::ParameterSet&);
    ~TreeProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    void clearTree();
    
    edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > totemRPTracksToken_;
    edm::EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
    edm::EDGetTokenT< edm::View<reco::Vertex> > vtxToken_;
//                                 JW
    edm::EDGetTokenT< edm::View<flashgg::Electron> > electronToken_;
    edm::EDGetTokenT< edm::View<flashgg::Muon> > muonToken_; 
    edm::EDGetTokenT< edm::View<vector<flashgg::Jet> > >  jetToken_; 
//
    
    double sqrtS_;
    double singlePhotonMinPt_, singlePhotonMaxEta_, singlePhotonMinR9_;
    double photonPairMinMass_;
    std::string filename_;

    bool useXiInterp_;
    ProtonUtils::XiInterpolator* xiInterp_;

    CTPPSAlCa::FillNumberLUTHandler* fillLUTHandler_;
    CTPPSAlCa::AlignmentLUTHandler* alignmentLUTHandler_;

    TFile* file_;
    TTree* tree_;

    // --- tree components ---

    unsigned int fBX, fRun, fLumiSection;
    unsigned long long fEventNum;

    unsigned int fProtonTrackNum;
    float fProtonTrackXi[MAX_PROTON_TRK], fProtonTrackXiError[MAX_PROTON_TRK];
    unsigned int fProtonTrackSide[MAX_PROTON_TRK], fProtonTrackPot[MAX_PROTON_TRK];
 //                   JW
    unsigned int fElectronNum;
    unsigned int fMuonNum;
    unsigned int fJetNum;
    float fElectronPt[MAX_ELECTRON], fElectronEta[MAX_ELECTRON], fElectronPhi[MAX_ELECTRON], fElectronE[MAX_ELECTRON];
    float fMuonPt[MAX_MUON], fMuonEta[MAX_MUON], fMuonPhi[MAX_MUON], fMuonE[MAX_MUON];
    float fJetPt[MAX_JET], fJetEta[MAX_JET], fJetPhi[MAX_JET], fJetE[MAX_JET], fJetMass[MAX_JET];
 //

    unsigned int fDiphotonNum;
    float fDiphotonPt1[MAX_DIPHOTON], fDiphotonPt2[MAX_DIPHOTON];
    float fDiphotonEta1[MAX_DIPHOTON], fDiphotonEta2[MAX_DIPHOTON];
    float fDiphotonPhi1[MAX_DIPHOTON], fDiphotonPhi2[MAX_DIPHOTON];
    float fDiphotonR91[MAX_DIPHOTON], fDiphotonR92[MAX_DIPHOTON];
    float fDiphotonM[MAX_DIPHOTON], fDiphotonY[MAX_DIPHOTON];
    float fDiphotonPt[MAX_DIPHOTON], fDiphotonDphi[MAX_DIPHOTON];

    unsigned int fDiphotonVertexTracks[MAX_DIPHOTON];
    unsigned int fDiphotonVerticesAt1mmDist[MAX_DIPHOTON], fDiphotonVerticesAt2mmDist[MAX_DIPHOTON], fDiphotonVerticesAt5mmDist[MAX_DIPHOTON], fDiphotonVerticesAt1cmDist[MAX_DIPHOTON];
    float fDiphotonVertexX[MAX_DIPHOTON], fDiphotonVertexY[MAX_DIPHOTON], fDiphotonVertexZ[MAX_DIPHOTON];
    float fDiphotonNearestDist[MAX_DIPHOTON];

    float fMET, fMETPhi;

    unsigned int fVertexNum;

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
TreeProducer::TreeProducer(const edm::ParameterSet& iConfig) :
  totemRPTracksToken_( mayConsume< edm::DetSetVector<TotemRPLocalTrack> >( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  diphotonToken_     ( consumes< edm::View<flashgg::DiPhotonCandidate> > ( iConfig.getParameter<edm::InputTag>( "diphotonLabel" ) ) ),
  metToken_          ( mayConsume< edm::View<pat::MET> >                 ( iConfig.getParameter<edm::InputTag>( "metLabel") ) ),
  vtxToken_          ( mayConsume< edm::View<reco::Vertex> >             ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
//                               JW
  electronToken_     ( mayConsume< edm::View<flashgg::Electron> >        ( iConfig.getParameter<edm::InputTag>( "electronLabel") ) ),
  muonToken_         ( mayConsume< edm::View<flashgg::Muon> >            ( iConfig.getParameter<edm::InputTag>( "muonLabel") ) ),
  jetToken_          ( consumes< edm::View< std::vector<flashgg::Jet> > >( iConfig.getParameter<edm::InputTag>( "jetLabel") ) ),
//
  sqrtS_             ( iConfig.getParameter<double>( "sqrtS")),
  singlePhotonMinPt_ ( iConfig.getParameter<double>( "minPtSinglePhoton" ) ),
  singlePhotonMaxEta_( iConfig.getParameter<double>( "maxEtaSinglePhoton" ) ),
  singlePhotonMinR9_ ( iConfig.getParameter<double>( "minR9SinglePhoton" ) ),
  photonPairMinMass_ ( iConfig.getParameter<double>( "minMassDiPhoton" ) ),
  filename_          ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  useXiInterp_       ( iConfig.getParameter<bool>( "useXiInterpolation" ) ),
  xiInterp_( 0 ),
  file_( 0 ), tree_( 0 )

{
  xiInterp_ = new ProtonUtils::XiInterpolator;
  if ( useXiInterp_ ) { xiInterp_->loadInterpolationGraphs( iConfig.getParameter<edm::FileInPath>( "xiInterpolationFile" ).fullPath().c_str() ); }
  fillLUTHandler_ = new CTPPSAlCa::FillNumberLUTHandler( iConfig.getParameter<edm::FileInPath>( "fillNumLUTFile" ).fullPath().c_str() );
  alignmentLUTHandler_ = new CTPPSAlCa::AlignmentLUTHandler( iConfig.getParameter<edm::FileInPath>( "alignmentLUTFile" ).fullPath().c_str() );

  file_ = new TFile( filename_.c_str(), "recreate" );
  file_->cd();

  tree_ = new TTree( "ntp", "diphoton ntuple" );
}


TreeProducer::~TreeProducer()
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
TreeProducer::clearTree()
{
  fProtonTrackNum = 0;
  for ( unsigned int i=0; i<MAX_PROTON_TRK; i++ ) {
    fProtonTrackXi[i] = fProtonTrackXiError[i] = -1.;
    fProtonTrackSide[i] = 2; //invalid
    fProtonTrackPot[i] = 0;
  }

  fDiphotonNum = 0;
  for ( unsigned int i=0; i<MAX_DIPHOTON; i++ ) {
    fDiphotonPt1[i] = fDiphotonPt2[i] = -1.;
    fDiphotonEta1[i] = fDiphotonEta2[i] = -1.;
    fDiphotonPhi1[i] = fDiphotonPhi2[i] = -1.;
    fDiphotonR91[i] = fDiphotonR92[i] = -1.;
    fDiphotonM[i] = fDiphotonY[i] = fDiphotonPt[i] = fDiphotonDphi[i] = -1.;
    fDiphotonVertexTracks[i] = 0;
    fDiphotonVerticesAt1mmDist[i] = fDiphotonVerticesAt2mmDist[i] = fDiphotonVerticesAt5mmDist[i] = fDiphotonVerticesAt1cmDist[i] = 0;
    fDiphotonVertexX[i] = fDiphotonVertexY[i] = fDiphotonVertexZ[i] = -1.;
    fDiphotonNearestDist[i] = 999.;
  }

  fElectronNum = 0;
  for ( unsigned int i=0; i<MAX_ELECTRON; i++ ) {
    fElectronPt[i] = -1.;
    fElectronEta[i] = -1.;
    fElectronPhi[i] = -1.;
    fElectronE[i] = -1.;
  }

  fMuonNum = 0;
  for ( unsigned int i=0; i<MAX_MUON; i++ ) {
    fMuonPt[i] = -1.;
    fMuonEta[i] = -1.;
    fMuonPhi[i] = -1.;
    fMuonE[i] = -1.;
  }

  fJetNum = 0;
  for ( unsigned int i=0; i<MAX_JET; i++ ) {
    fJetPt[i] = -1.;
    fJetEta[i] = -1.;
    fJetPhi[i] = -1.;
    fJetE[i] = -1.;
    fJetMass[i] = -1.;
  }

  fMET = fMETPhi = -1.;

  fVertexNum = 0;

}

// ------------ method called for each event  ------------
void
TreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  #include <iostream> // for debugging purposes

  clearTree();

  // Run and BX information
  fBX = iEvent.bunchCrossing();
  fRun = iEvent.id().run();
  fLumiSection = iEvent.luminosityBlock();
  fEventNum = iEvent.id().event();

  // fetch the jet collection from EDM file
  edm::Handle< edm::View< std::vector<flashgg::Jet> > > jetsColls;
  iEvent.getByToken( jetToken_, jetsColls );

  // fetch the diphoton collection from EDM file
  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken(diphotonToken_, diphotons);

  fDiphotonNum = 0;
  edm::Ptr<reco::Vertex> diphoton_vtx[MAX_DIPHOTON];

  fJetNum = 0;

  for ( unsigned int i=0; i<diphotons->size() && fDiphotonNum<MAX_DIPHOTON; i++ ) {
    edm::Ptr<flashgg::DiPhotonCandidate> diphoton = diphotons->ptrAt( i );

    if ( diphoton->leadPhotonId()<-0.9 ) continue;
    if ( diphoton->subLeadPhotonId()<-0.9 ) continue;

    if ( !passSinglePhotonCuts( diphoton->leadingPhoton() ) ) continue;
    if ( !passSinglePhotonCuts( diphoton->subLeadingPhoton() ) ) continue;

    if ( fabs( diphoton->leadingPhoton()->eta() )>=singlePhotonMaxEta_ or fabs( diphoton->subLeadingPhoton()->eta() )>=singlePhotonMaxEta_ ) continue;
    if ( diphoton->leadingPhoton()->pt()<singlePhotonMinPt_ or diphoton->subLeadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( diphoton->leadingPhoton()->full5x5_r9()<singlePhotonMinR9_ or diphoton->subLeadingPhoton()->full5x5_r9()<singlePhotonMinR9_ ) continue;

    if ( diphoton->mass()<photonPairMinMass_ ) continue;

    fDiphotonVertexTracks[fDiphotonNum] = diphoton->vtx()->tracksSize();
    fDiphotonVertexX[fDiphotonNum] = diphoton->vtx()->x();
    fDiphotonVertexY[fDiphotonNum] = diphoton->vtx()->y();
    fDiphotonVertexZ[fDiphotonNum] = diphoton->vtx()->z();
    diphoton_vtx[fDiphotonNum] = diphoton->vtx();

    fDiphotonPt1[fDiphotonNum] = diphoton->leadingPhoton()->pt();
    fDiphotonEta1[fDiphotonNum] = diphoton->leadingPhoton()->eta();
    fDiphotonPhi1[fDiphotonNum] = diphoton->leadingPhoton()->phi();
    fDiphotonR91[fDiphotonNum] = diphoton->leadingPhoton()->full5x5_r9();

    fDiphotonPt2[fDiphotonNum] = diphoton->subLeadingPhoton()->pt();
    fDiphotonEta2[fDiphotonNum] = diphoton->subLeadingPhoton()->eta();
    fDiphotonPhi2[fDiphotonNum] = diphoton->subLeadingPhoton()->phi();
    fDiphotonR92[fDiphotonNum] = diphoton->subLeadingPhoton()->full5x5_r9();

    fDiphotonM[fDiphotonNum] = diphoton->mass();
    fDiphotonY[fDiphotonNum] = diphoton->rapidity();
    fDiphotonPt[fDiphotonNum] = diphoton->pt();

    float dphi = diphoton->leadingPhoton()->phi()-diphoton->subLeadingPhoton()->phi();
    while ( dphi<-TMath::Pi() ) dphi += 2.*TMath::Pi();
    while ( dphi> TMath::Pi() ) dphi -= 2.*TMath::Pi();
    fDiphotonDphi[fDiphotonNum] = dphi;

    //           JW

    edm::Ptr< std::vector<flashgg::Jet> > jets = jetsColls->ptrAt( i );
    for ( unsigned int j=0; j<jets->size() && fJetNum<MAX_JET; j++ ) {
      const flashgg::Jet jet = jets->at( j );
      fJetPt[fJetNum]   = jet.pt();
      fJetEta[fJetNum]  = jet.eta();
      fJetPhi[fJetNum]  = jet.phi();
      fJetE[fJetNum]    = jet.energy();
      fJetMass[fJetNum] = jet.mass();

      fJetNum++;
    }
    
    //

    //std::cout << fDiphotonPt1[fDiphotonNum] << " --- " << fDiphotonPt2[fDiphotonNum] << " --- " << fDiphotonM[fDiphotonNum] << std::endl;

    fDiphotonNum++;
  }

  if ( fDiphotonNum<1 ) return;

  // fetch the proton collection from EDM file
  edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rpLocalTracks;
  iEvent.getByToken( totemRPTracksToken_, rpLocalTracks );

  const unsigned short fill_num = fillLUTHandler_->getFillNumber( iEvent.id().run() );
  xiInterp_->setAlignmentConstants( alignmentLUTHandler_->getAlignmentConstants( fill_num ) ); // fill-based alignment corrections
  xiInterp_->setCalibrationConstants( iEvent.id().run() ); // run-based calibration parameters

  fProtonTrackNum = 0;
  for ( edm::DetSetVector<TotemRPLocalTrack>::const_iterator it=rpLocalTracks->begin(); it!=rpLocalTracks->end(); it++ ) {
    const TotemRPDetId detid( it->detId() );
    for ( edm::DetSet<TotemRPLocalTrack>::const_iterator trk=it->begin(); trk!=it->end(); trk++ ) {
      if ( !trk->isValid() ) { continue; }
      float xi = -1., err_xi = -1.;
      if ( useXiInterp_ ) { xiInterp_->computeXiSpline( detid, *trk, &xi, &err_xi ); }
      else                { xiInterp_->computeXiLinear( detid, *trk, &xi, &err_xi ); }
      fProtonTrackXi[fProtonTrackNum] = xi;
      fProtonTrackXiError[fProtonTrackNum] = err_xi;
      fProtonTrackSide[fProtonTrackNum] = detid.arm(); // 0 = left ; 1 = right
      fProtonTrackPot[fProtonTrackNum] = detid.romanPot();
      fProtonTrackNum++;
    }
  }
 
 //                               JW
 
 //Implementing electron 4 vector
 
 
 // fetch the electron collection from EDM file
  edm::Handle< edm::View<flashgg::Electron> > electrons;
  iEvent.getByToken( electronToken_, electrons );
 
 fElectronNum=0;
 //  edm::Ptr<reco::Vertex> electron_vtx[MAX_ELECTRON];

  for ( unsigned int i=0; i<electrons->size() && fElectronNum<MAX_ELECTRON; i++ ) {
  edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( i );
  fElectronPt[fElectronNum] = electron->pt();
  fElectronEta[fElectronNum] = electron->eta();
  fElectronPhi[fElectronNum] = electron->phi();
  fElectronE[fElectronNum] = electron->energy();

  //  fElectronVertexX[fElectronNum] = electron->vtx()->x();
  //  fElectronVertexY[fElectronNum] = electron->vtx()->y();
  //  fElectronVertexZ[fElectronNum] = electron->vtx()->z();
  //  electron_vtx[fElectronNum] = electron->vtx();

 fElectronNum++;
 }

 
 //Implementing muon 4 vector
 
 
 edm::Handle< edm::View<flashgg::Muon> > muons;
 iEvent.getByToken(muonToken_,muons);

 // Add muon vertex here
 
 fMuonNum=0;
 //edm::Ptr<reco::Vertex> muon_vtx[MAX_MUON];


 for ( unsigned int i=0; i<muons->size() && fMuonNum<MAX_MUON; i++ ) {
  edm::Ptr<flashgg::Muon> muon = muons->ptrAt( i );
  fMuonPt[fMuonNum] = muon->pt();
  fMuonEta[fMuonNum] = muon->eta();
  fMuonPhi[fMuonNum] = muon->phi();
  fMuonE[fMuonNum] = muon->energy();
 
 fMuonNum++;
 }
 
 //Implementing jet 4 vector
 
 /*
 edm::Handle< edm::View<vector<flashgg::Jet> > > jets;
 iEvent.getByToken(jetToken_,jets);
 
 // Add jet vertex here

 fJetNum=0;
 //edm::Ptr<reco::Vertex> jet_vtx{MAX_JET];

 for ( unsigned int i=0; i<jets->size() && fJetNum<MAX_JET; i++ ) {
   edm::Ptr<vector<flashgg::Jet> > jet = jets->ptrAt( i );
   fJetPt[fJetNum] = jet->ptrAt(i)->pt();
   fJetEta[fJetNum] = jet->eta();
   fJetPhi[fJetNum] = jet->phi();
   fJetE[fJetNum] = jet->energy();
   fJetMass[fJetNum] = jet->mass();
 
 fJetNum++;
 }
 */
//
 //
 
  std::cout << "# found " << fDiphotonNum << " diphoton candidate(s) with " << fProtonTrackNum << " proton track(s) (all pots)!" << std::endl;
  // retrieve the missing ET
  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken( metToken_, mets );

  const edm::View<pat::MET>* metColl = mets.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();
  fMET = met->sumEt();
  fMETPhi = met->phi();

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );
  fVertexNum = vertices->size();

  for ( unsigned int i=0; i<vertices->size(); i++ ) {
    edm::Ptr<reco::Vertex> vtx = vertices->ptrAt( i );
    if ( !vtx->isValid() ) continue;

    // loop over all the diphoton candidates to find the closest vertices
    for ( unsigned int j=0; j<fDiphotonNum; j++ ) {
      if ( diphoton_vtx[j]->position()==vtx->position() ) continue; // found the diphoton vertex
      const float vtx_dist = sqrt( pow( diphoton_vtx[j]->x()-vtx->x(), 2 )+pow( diphoton_vtx[j]->y()-vtx->y(), 2 )+pow( diphoton_vtx[j]->z()-vtx->z(), 2 ) );
      if ( vtx_dist<fDiphotonNearestDist[j] ) fDiphotonNearestDist[j] = vtx_dist;
      if ( vtx_dist<=0.1 ) fDiphotonVerticesAt1mmDist[j]++;
      if ( vtx_dist<=0.2 ) fDiphotonVerticesAt2mmDist[j]++;
      if ( vtx_dist<=0.5 ) fDiphotonVerticesAt5mmDist[j]++;
      if ( vtx_dist<=1.0 ) fDiphotonVerticesAt1cmDist[j]++;
    }
  }

  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{
  tree_->Branch( "run_id", &fRun, "run_id/i");
  tree_->Branch( "lumisection", &fLumiSection, "lumisection/i");
  tree_->Branch( "bunch_crossing", &fBX, "bunch_crossing/i");
  tree_->Branch( "event_number", &fEventNum, "event_number/l");

  tree_->Branch( "num_proton_track", &fProtonTrackNum, "num_proton_track/i" );
  tree_->Branch( "proton_track_xi", fProtonTrackXi, "proton_track_xi[num_proton_track]/F" );
  tree_->Branch( "proton_track_xi_error", fProtonTrackXiError, "proton_track_xi_error[num_proton_track]/F" );
  tree_->Branch( "proton_track_side", fProtonTrackSide, "proton_track_side[num_proton_track]/i" );
  tree_->Branch( "proton_track_pot", fProtonTrackPot, "proton_track_pot[num_proton_track]/i" );

  tree_->Branch( "num_diphoton", &fDiphotonNum, "num_diphoton/i" );
  tree_->Branch( "diphoton_pt1", fDiphotonPt1, "diphoton_pt1[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt2", fDiphotonPt2, "diphoton_pt2[num_diphoton]/F" );
  tree_->Branch( "diphoton_eta1", fDiphotonEta1, "diphoton_eta1[num_diphoton]/F" );
  tree_->Branch( "diphoton_eta2", fDiphotonEta2, "diphoton_eta2[num_diphoton]/F" );
  tree_->Branch( "diphoton_phi1", fDiphotonPhi1, "diphoton_phi1[num_diphoton]/F" );
  tree_->Branch( "diphoton_phi2", fDiphotonPhi2, "diphoton_phi2[num_diphoton]/F" );
  tree_->Branch( "diphoton_r91", fDiphotonR91, "diphoton_r91[num_diphoton]/F" );
  tree_->Branch( "diphoton_r92", fDiphotonR92, "diphoton_r92[num_diphoton]/F" );
  tree_->Branch( "diphoton_mass", fDiphotonM, "diphoton_mass[num_diphoton]/F" );
  tree_->Branch( "diphoton_rapidity", fDiphotonY, "diphoton_rapidity[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt", fDiphotonPt, "diphoton_pt[num_diphoton]/F" );
  tree_->Branch( "diphoton_dphi", fDiphotonDphi, "diphoton_dphi[num_diphoton]/F" );

  tree_->Branch( "diphoton_vertex_tracks", fDiphotonVertexTracks, "diphoton_vertex_tracks[num_diphoton]/i" );
  tree_->Branch( "diphoton_vertex_x", fDiphotonVertexX, "diphoton_vertex_x[num_diphoton]/F" );
  tree_->Branch( "diphoton_vertex_y", fDiphotonVertexY, "diphoton_vertex_y[num_diphoton]/F" );
  tree_->Branch( "diphoton_vertex_z", fDiphotonVertexZ, "diphoton_vertex_z[num_diphoton]/F" );
  tree_->Branch( "diphoton_vertex_nearestvtxdist", fDiphotonNearestDist, "diphoton_vertex_nearestvtxdist[num_diphoton]/F" );
  tree_->Branch( "diphoton_vertex_vtx1mmdist", fDiphotonVerticesAt1mmDist, "diphoton_vertex_vtx1mmdist[num_diphoton]/i" );
  tree_->Branch( "diphoton_vertex_vtx2mmdist", fDiphotonVerticesAt2mmDist, "diphoton_vertex_vtx2mmdist[num_diphoton]/i" );
  tree_->Branch( "diphoton_vertex_vtx5mmdist", fDiphotonVerticesAt5mmDist, "diphoton_vertex_vtx5mmdist[num_diphoton]/i" );
  tree_->Branch( "diphoton_vertex_vtx1cmdist", fDiphotonVerticesAt1cmDist, "diphoton_vertex_vtx1cmdist[num_diphoton]/i" );

  tree_->Branch( "num_electron", &fElectronNum, "num_electron/i" );
  tree_->Branch( "electron_pt", fElectronPt, "electron_pt[num_electron]/F" );
  tree_->Branch( "electron_eta", fElectronEta, "electron_eta[num_electron]/F" );
  tree_->Branch( "electron_phi", fElectronPhi, "electron_phi[num_electron]/F" );
  tree_->Branch( "electron_energy", fElectronE, "electron_energy[num_electron]/F" );

  tree_->Branch( "num_muon", &fMuonNum, "num_muon/i" );
  tree_->Branch( "muon_pt", fMuonPt, "muon_pt[num_muon]/F" );
  tree_->Branch( "muon_eta", fMuonEta, "muon_eta[num_muon]/F" );
  tree_->Branch( "muon_phi", fMuonPhi, "muon_phi[num_muon]/F" );
  tree_->Branch( "muon_energy", fMuonE, "muon_energy[num_muon]/F" );

  tree_->Branch( "num_jet", &fJetNum, "num_jet/i" );
  tree_->Branch( "jet_pt", fJetPt, "jet_pt[num_jet]/F" );
  tree_->Branch( "jet_eta", fJetEta, "jet_eta[num_jet]/F" );
  tree_->Branch( "jet_phi", fJetPhi, "jet_phi[num_jet]/F" );
  tree_->Branch( "jet_energy", fJetE, "jet_energy[num_jet]/F" );
  tree_->Branch( "jet_mass", fJetMass, "jet_mass[num_jet]/F" );

  tree_->Branch( "num_vertex", &fVertexNum, "num_vertex/i" );

  tree_->Branch( "met", &fMET );
  tree_->Branch( "met_phi", &fMETPhi );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer);
