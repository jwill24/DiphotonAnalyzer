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
#include "flashgg/DataFormats/interface/Met.h"

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

#define MAX_PROTON_TRK 20
#define MAX_DIPHOTON 10
//                               JW
#define MAX_ELECTRON 20
#define MAX_MUON 20
#define MAX_JET 50
//
#define MAX_VERTEX 100

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit TreeProducer( const edm::ParameterSet& );
    ~TreeProducer();

    static void fillDescriptions( edm::ConfigurationDescriptions& );


  private:
    virtual void beginJob() override;
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    void clearTree();
    
    edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > totemRPTracksToken_;
    edm::EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT< edm::View<flashgg::Met> > metToken_;
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

    unsigned int fBX, fFill, fRun, fLumiSection;
    unsigned long long fEventNum;

    unsigned int fProtonTrackNum;
    float fProtonTrackX[MAX_PROTON_TRK], fProtonTrackY[MAX_PROTON_TRK];
    float fProtonTrackXi[MAX_PROTON_TRK], fProtonTrackXiError[MAX_PROTON_TRK];
    unsigned int fProtonTrackSide[MAX_PROTON_TRK], fProtonTrackPot[MAX_PROTON_TRK], fProtonTrackLinkNF[MAX_PROTON_TRK];
    float fProtonTrackMinLinkDist[MAX_PROTON_TRK];

    unsigned int fElectronNum;
    float fElectronPt[MAX_ELECTRON], fElectronEta[MAX_ELECTRON], fElectronPhi[MAX_ELECTRON], fElectronE[MAX_ELECTRON];
    float fElectronVertexX[MAX_ELECTRON], fElectronVertexY[MAX_ELECTRON], fElectronVertexZ[MAX_ELECTRON];

    unsigned int fMuonNum;
    float fMuonPt[MAX_MUON], fMuonEta[MAX_MUON], fMuonPhi[MAX_MUON], fMuonE[MAX_MUON];
    float fMuonVertexX[MAX_MUON], fMuonVertexY[MAX_MUON], fMuonVertexZ[MAX_MUON];

    unsigned int fJetNum;
    float fJetPt[MAX_JET], fJetEta[MAX_JET], fJetPhi[MAX_JET], fJetE[MAX_JET], fJetMass[MAX_JET];
    float fJetVertexX[MAX_JET], fJetVertexY[MAX_JET], fJetVertexZ[MAX_JET];

    unsigned int fDiphotonNum;
    float fDiphotonPt1[MAX_DIPHOTON], fDiphotonPt2[MAX_DIPHOTON];
    float fDiphotonEta1[MAX_DIPHOTON], fDiphotonEta2[MAX_DIPHOTON];
    float fDiphotonPhi1[MAX_DIPHOTON], fDiphotonPhi2[MAX_DIPHOTON];
    float fDiphotonR91[MAX_DIPHOTON], fDiphotonR92[MAX_DIPHOTON];
    float fDiphotonM[MAX_DIPHOTON], fDiphotonY[MAX_DIPHOTON];
    float fDiphotonPt[MAX_DIPHOTON], fDiphotonDphi[MAX_DIPHOTON];

    int fDiphotonVertex[MAX_DIPHOTON];
    unsigned int fDiphotonVertexTracks[MAX_DIPHOTON];
    unsigned int fDiphotonVerticesAt1mmDist[MAX_DIPHOTON], fDiphotonVerticesAt2mmDist[MAX_DIPHOTON], fDiphotonVerticesAt5mmDist[MAX_DIPHOTON], fDiphotonVerticesAt1cmDist[MAX_DIPHOTON];
    float fDiphotonVertexX[MAX_DIPHOTON], fDiphotonVertexY[MAX_DIPHOTON], fDiphotonVertexZ[MAX_DIPHOTON];
    float fDiphotonNearestDist[MAX_DIPHOTON];

    float fMET, fMETPhi;

    unsigned int fVertexNum;
    float fVertexX[MAX_VERTEX], fVertexY[MAX_VERTEX], fVertexZ[MAX_VERTEX];
    //unsigned int fVertexTracks[MAX_VERTEX], fVertexTracksWght0p75[MAX_VERTEX], fVertexTracksWght0p9[MAX_VERTEX], fVertexTracksWght0p95[MAX_VERTEX];
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
TreeProducer::TreeProducer( const edm::ParameterSet& iConfig ) :
  totemRPTracksToken_( consumes< edm::DetSetVector<TotemRPLocalTrack> >  ( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  diphotonToken_     ( consumes< edm::View<flashgg::DiPhotonCandidate> > ( iConfig.getParameter<edm::InputTag>( "diphotonLabel" ) ) ),
  metToken_          ( consumes< edm::View<flashgg::Met> >               ( iConfig.getParameter<edm::InputTag>( "metLabel") ) ),
  vtxToken_          ( consumes< edm::View<reco::Vertex> >               ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
//                               JW
  electronToken_     ( consumes< edm::View<flashgg::Electron> >          ( iConfig.getParameter<edm::InputTag>( "electronLabel") ) ),
  muonToken_         ( consumes< edm::View<flashgg::Muon> >              ( iConfig.getParameter<edm::InputTag>( "muonLabel") ) ),
  jetToken_          ( consumes< edm::View< std::vector<flashgg::Jet> > >( iConfig.getParameter<edm::InputTag>( "jetLabel") ) ),
//
  sqrtS_             ( iConfig.getParameter<double>     ( "sqrtS")),
  singlePhotonMinPt_ ( iConfig.getParameter<double>     ( "minPtSinglePhoton" ) ),
  singlePhotonMaxEta_( iConfig.getParameter<double>     ( "maxEtaSinglePhoton" ) ),
  singlePhotonMinR9_ ( iConfig.getParameter<double>     ( "minR9SinglePhoton" ) ),
  photonPairMinMass_ ( iConfig.getParameter<double>     ( "minMassDiPhoton" ) ),
  filename_          ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  useXiInterp_       ( iConfig.getParameter<bool>       ( "useXiInterpolation" ) ),
  xiInterp_( 0 ),
  file_( 0 ), tree_( 0 )

{
  xiInterp_ = new ProtonUtils::XiInterpolator;
  if ( useXiInterp_ ) {
    xiInterp_->loadInterpolationGraphs( iConfig.getParameter<edm::FileInPath>( "xiInterpolationFile" ).fullPath().c_str() );
  }
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
  if ( xiInterp_ ) delete xiInterp_;
  if ( fillLUTHandler_ ) delete fillLUTHandler_;
  if ( alignmentLUTHandler_ ) delete alignmentLUTHandler_;
}


//
// member functions
//

void
TreeProducer::clearTree()
{
  fProtonTrackNum = 0;
  for ( unsigned int i=0; i<MAX_PROTON_TRK; i++ ) {
    fProtonTrackX[i] = fProtonTrackY[i] = -1.;
    fProtonTrackXi[i] = fProtonTrackXiError[i] = -1.;
    fProtonTrackSide[i] = 2; //invalid
    fProtonTrackPot[i] = 0;
    fProtonTrackLinkNF[i] = 999;
    fProtonTrackMinLinkDist[i] = -1.;
  }

  fDiphotonNum = 0;
  for ( unsigned int i=0; i<MAX_DIPHOTON; i++ ) {
    fDiphotonPt1[i] = fDiphotonPt2[i] = -1.;
    fDiphotonEta1[i] = fDiphotonEta2[i] = -1.;
    fDiphotonPhi1[i] = fDiphotonPhi2[i] = -1.;
    fDiphotonR91[i] = fDiphotonR92[i] = -1.;
    fDiphotonM[i] = fDiphotonY[i] = fDiphotonPt[i] = fDiphotonDphi[i] = -1.;
    fDiphotonVertexTracks[i] = 0;
    fDiphotonVertex[i] = -1;
    fDiphotonVerticesAt1mmDist[i] = fDiphotonVerticesAt2mmDist[i] = fDiphotonVerticesAt5mmDist[i] = fDiphotonVerticesAt1cmDist[i] = 0;
    fDiphotonVertexX[i] = fDiphotonVertexY[i] = fDiphotonVertexZ[i] = -1.;
    fDiphotonNearestDist[i] = 999.;
  }

  fElectronNum = 0;
  for ( unsigned int i=0; i<MAX_ELECTRON; i++ ) {
    fElectronPt[i] = fElectronEta[i] = fElectronPhi[i] = fElectronE[i] = -1.;
    fElectronVertexX[i] = fElectronVertexY[i] = fElectronVertexZ[i] = -1.;
  }

  fMuonNum = 0;
  for ( unsigned int i=0; i<MAX_MUON; i++ ) {
    fMuonPt[i] = fMuonEta[i] = fMuonPhi[i] = fMuonE[i] = -1.;
    fMuonVertexX[i] = fMuonVertexY[i] = fMuonVertexZ[i] = -1.;
  }

  fJetNum = 0;
  for ( unsigned int i=0; i<MAX_JET; i++ ) {
    fJetPt[i] = fJetEta[i] = fJetPhi[i] = fJetE[i] = fJetMass[i] = -1.;
    fJetVertexX[i] = fJetVertexY[i] = fJetVertexZ[i] = -1.;
  }

  fMET = fMETPhi = -1.;

  fVertexNum = 0;
  for ( unsigned int i=0; i<MAX_VERTEX; i++ ) {
    fVertexX[i] = fVertexY[i] = fVertexZ[i] = -999.;
    //fVertexTracks[i] = fVertexTracksWght0p75[i] = fVertexTracksWght0p9[i] = fVertexTracksWght0p95[i] = 0;
  }

}

// ------------ method called for each event  ------------
void
TreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& )
{
  #include <iostream> // for debugging purposes

  clearTree();

  // Run and BX information
  fBX = iEvent.bunchCrossing();
  fRun = iEvent.id().run();
  fLumiSection = iEvent.luminosityBlock();
  fEventNum = iEvent.id().event();

  // get the fill number from the run id <-> fill number LUT
  fFill = fillLUTHandler_->getFillNumber( iEvent.id().run() );

  //----- diphoton candidates -----

  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken(diphotonToken_, diphotons);

  fDiphotonNum = 0;
  edm::Ptr<reco::Vertex> diphoton_vtx[MAX_DIPHOTON];
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

    fDiphotonVertexTracks[fDiphotonNum] = diphoton->vtx()->nTracks();
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

    fDiphotonNum++;
  }

  //----- only store events with > 0 diphoton candidate(s) -----

  if ( fDiphotonNum<1 ) return;

  //----- jets collection -----

  edm::Handle< edm::View< std::vector<flashgg::Jet> > > jetsColls;
  iEvent.getByToken( jetToken_, jetsColls );

  fJetNum = 0;
  for ( unsigned int i=0; i<jetsColls->size(); i++ ) {
    //           JW
    edm::Ptr< std::vector<flashgg::Jet> > jets = jetsColls->ptrAt( i );
    if ( fJetNum>MAX_JET ) { std::cerr << ">> More jets than expected in this event (" << fJetNum << ">MAX_JET=" << MAX_JET << "). Increase MAX_JET for safety" << std::endl; }

    for ( unsigned int j=0; j<jets->size() && fJetNum<MAX_JET; j++ ) {
      const flashgg::Jet jet = jets->at( j );
      fJetPt[fJetNum]   = jet.pt();
      fJetEta[fJetNum]  = jet.eta();
      fJetPhi[fJetNum]  = jet.phi();
      fJetE[fJetNum]    = jet.energy();
      fJetMass[fJetNum] = jet.mass();

      fJetVertexX[fJetNum] = jet.vertex().x();
      fJetVertexY[fJetNum] = jet.vertex().y();
      fJetVertexZ[fJetNum] = jet.vertex().z();
      fJetNum++;
    }
  }

  //----- forward RP tracks -----

  edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rpLocalTracks;
  iEvent.getByToken( totemRPTracksToken_, rpLocalTracks );

  const CTPPSAlCa::RPAlignmentConstants align = alignmentLUTHandler_->getAlignmentConstants( fFill ); // fill-based alignment corrections
  xiInterp_->setAlignmentConstants( align );
  xiInterp_->setCalibrationConstants( fRun ); // run-based calibration parameters

  typedef std::pair<unsigned int, const TotemRPLocalTrack&> localtrack_t; // RP id -> local track object
  std::map<unsigned int,localtrack_t> map_near, map_far; // local track id in the tree -> localtrack_t object

  fProtonTrackNum = 0;
  for ( edm::DetSetVector<TotemRPLocalTrack>::const_iterator it=rpLocalTracks->begin(); it!=rpLocalTracks->end(); it++ ) {
    const TotemRPDetId detid( TotemRPDetId::decToRawId( it->detId()*10 ) );
    const unsigned short side = detid.arm(),
                         pot = detid.romanPot();
    const CTPPSAlCa::RPAlignmentConstants::Quantities align_quant = align.quantities( it->detId() );

    for ( edm::DetSet<TotemRPLocalTrack>::const_iterator trk=it->begin(); trk!=it->end(); trk++ ) {
      if ( !trk->isValid() ) { continue; }

      fProtonTrackX[fProtonTrackNum] = ( trk->getX0() + align_quant.x ) * 1.e-3; // store in m
      fProtonTrackY[fProtonTrackNum] = ( trk->getY0() - align_quant.y ) * 1.e-3; // store in m
      fProtonTrackSide[fProtonTrackNum] = side; // 0 = left (45) ; 1 = right (56)
      fProtonTrackPot[fProtonTrackNum] = pot; // 2 = 210m ; 3 = 220m

      // x-to-xi interpolation
      float xi = -1., err_xi = -1.;
      if ( useXiInterp_ ) { xiInterp_->computeXiSpline( detid, *trk, xi, err_xi ); }
      else                { xiInterp_->computeXiLinear( detid, *trk, xi, err_xi ); }
      fProtonTrackXi[fProtonTrackNum] = xi;
      fProtonTrackXiError[fProtonTrackNum] = err_xi;

      switch ( pot ) {
        case 2: { map_near.insert( std::make_pair( fProtonTrackNum, localtrack_t( it->detId(), *trk ) ) ); } break;
        case 3: {  map_far.insert( std::make_pair( fProtonTrackNum, localtrack_t( it->detId(), *trk ) ) ); } break;
      }

      fProtonTrackNum++;
    }
  }

  // second loop to associate near-far tracks
  for ( std::map<unsigned int,localtrack_t>::const_iterator it_n=map_near.begin(); it_n!=map_near.end(); it_n++ ) {
    float min_dist = 999.999;
    unsigned int cand = 999;
    for ( std::map<unsigned int,localtrack_t>::const_iterator it_f=map_far.begin(); it_f!=map_far.end(); it_f++ ) {
      const float dist = ProtonUtils::tracksDistance( align, it_n->second, it_f->second ); // in cm
      if ( dist<0. ) continue; // skip the comparison if opposite side tracks
      if ( dist<min_dist ) {
        min_dist = dist;
        cand = it_f->first;
      }
    }
    if ( cand!=999 ) { // near-far match found
      fProtonTrackLinkNF[it_n->first] = cand;
      fProtonTrackLinkNF[cand] = it_n->first;
      fProtonTrackMinLinkDist[it_n->first] = fProtonTrackMinLinkDist[cand] = min_dist;
    }
  }
 
  //                               JW
 
  //----- electrons collection -----
 
  edm::Handle< edm::View<flashgg::Electron> > electrons;
  iEvent.getByToken( electronToken_, electrons );
 
  fElectronNum=0;
  for ( unsigned int i=0; i<electrons->size() && fElectronNum<MAX_ELECTRON; i++ ) {
    const edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( i );
    fElectronPt[fElectronNum] = electron->pt();
    fElectronEta[fElectronNum] = electron->eta();
    fElectronPhi[fElectronNum] = electron->phi();
    fElectronE[fElectronNum] = electron->energy();

    fElectronVertexX[fElectronNum] = electron->vertex().x();
    fElectronVertexY[fElectronNum] = electron->vertex().y();
    fElectronVertexZ[fElectronNum] = electron->vertex().z();
    fElectronNum++;
  }
 
  //----- muons collection -----

  edm::Handle< edm::View<flashgg::Muon> > muons;
  iEvent.getByToken(muonToken_,muons);

  fMuonNum=0;
  for ( unsigned int i=0; i<muons->size() && fMuonNum<MAX_MUON; i++ ) {
    const edm::Ptr<flashgg::Muon> muon = muons->ptrAt( i );
    fMuonPt[fMuonNum] = muon->pt();
    fMuonEta[fMuonNum] = muon->eta();
    fMuonPhi[fMuonNum] = muon->phi();
    fMuonE[fMuonNum] = muon->energy();
 
    fMuonVertexX[fMuonNum] = muon->vertex().x();
    fMuonVertexY[fMuonNum] = muon->vertex().y();
    fMuonVertexZ[fMuonNum] = muon->vertex().z();
    fMuonNum++;
  }

  //----- event-wide information -----

  std::cout << "# found " << fDiphotonNum << " diphoton candidate(s) with " << fProtonTrackNum << " proton track(s) (all pots)!" << std::endl;
  // retrieve the missing ET
  edm::Handle< edm::View<flashgg::Met> > mets;
  iEvent.getByToken( metToken_, mets );

  const edm::View<flashgg::Met>* metColl = mets.product();
  edm::View<flashgg::Met>::const_iterator met = metColl->begin();
  fMET = met->pt();
  fMETPhi = met->phi();

  //----- vertexing information -----

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );

  fVertexNum = 0;
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
      fDiphotonVertex[j] = fVertexNum;
    }

    fVertexX[fVertexNum] = vtx->x();
    fVertexY[fVertexNum] = vtx->y();
    fVertexZ[fVertexNum] = vtx->z();

    // tracks content not stored in the miniAOD event format...
    /*fVertexTracks[fVertexNum] = vtx->nTracks();
    fVertexTracksWght0p75[fVertexNum] = vtx->nTracks( 0.75 );
    fVertexTracksWght0p9[fVertexNum] = vtx->nTracks( 0.9 );
    fVertexTracksWght0p95[fVertexNum] = vtx->nTracks( 0.95 );*/
    fVertexNum++;
  }

  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{
  tree_->Branch( "run_id", &fRun, "run_id/i");
  tree_->Branch( "fill_number", &fFill, "fill_number/i");
  tree_->Branch( "lumisection", &fLumiSection, "lumisection/i");
  tree_->Branch( "bunch_crossing", &fBX, "bunch_crossing/i");
  tree_->Branch( "event_number", &fEventNum, "event_number/l");

  tree_->Branch( "num_proton_track", &fProtonTrackNum, "num_proton_track/i" );
  tree_->Branch( "proton_track_x", fProtonTrackX, "proton_track_x[num_proton_track]/F" );
  tree_->Branch( "proton_track_y", fProtonTrackY, "proton_track_y[num_proton_track]/F" );
  tree_->Branch( "proton_track_xi", fProtonTrackXi, "proton_track_xi[num_proton_track]/F" );
  tree_->Branch( "proton_track_xi_error", fProtonTrackXiError, "proton_track_xi_error[num_proton_track]/F" );
  tree_->Branch( "proton_track_side", fProtonTrackSide, "proton_track_side[num_proton_track]/i" );
  tree_->Branch( "proton_track_pot", fProtonTrackPot, "proton_track_pot[num_proton_track]/i" );
  tree_->Branch( "proton_track_link_nearfar", fProtonTrackLinkNF, "proton_track_link_nearfar[num_proton_track]/i" );
  tree_->Branch( "proton_track_link_mindist", fProtonTrackMinLinkDist, "proton_track_link_mindist[num_proton_track]/F" );

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
  tree_->Branch( "diphoton_vertex_id", fDiphotonVertex, "diphoton_vertex_id[num_diphoton]/I" );
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
  tree_->Branch( "electron_vtx_x", fElectronVertexX, "electron_vtx_x[num_electron]/F" );
  tree_->Branch( "electron_vtx_y", fElectronVertexY, "electron_vtx_y[num_electron]/F" );
  tree_->Branch( "electron_vtx_z", fElectronVertexZ, "electron_vtx_z[num_electron]/F" );

  tree_->Branch( "num_muon", &fMuonNum, "num_muon/i" );
  tree_->Branch( "muon_pt", fMuonPt, "muon_pt[num_muon]/F" );
  tree_->Branch( "muon_eta", fMuonEta, "muon_eta[num_muon]/F" );
  tree_->Branch( "muon_phi", fMuonPhi, "muon_phi[num_muon]/F" );
  tree_->Branch( "muon_energy", fMuonE, "muon_energy[num_muon]/F" );
  tree_->Branch( "muon_vtx_x", fMuonVertexX, "muon_vtx_x[num_muon]/F" );
  tree_->Branch( "muon_vtx_y", fMuonVertexY, "muon_vtx_y[num_muon]/F" );
  tree_->Branch( "muon_vtx_z", fMuonVertexZ, "muon_vtx_z[num_muon]/F" );

  tree_->Branch( "num_jet", &fJetNum, "num_jet/i" );
  tree_->Branch( "jet_pt", fJetPt, "jet_pt[num_jet]/F" );
  tree_->Branch( "jet_eta", fJetEta, "jet_eta[num_jet]/F" );
  tree_->Branch( "jet_phi", fJetPhi, "jet_phi[num_jet]/F" );
  tree_->Branch( "jet_energy", fJetE, "jet_energy[num_jet]/F" );
  tree_->Branch( "jet_mass", fJetMass, "jet_mass[num_jet]/F" );
  tree_->Branch( "jet_vtx_x", fJetVertexX, "jet_vtx_x[num_jet]/F" );
  tree_->Branch( "jet_vtx_y", fJetVertexY, "jet_vtx_y[num_jet]/F" );
  tree_->Branch( "jet_vtx_z", fJetVertexZ, "jet_vtx_z[num_jet]/F" );

  tree_->Branch( "num_vertex", &fVertexNum, "num_vertex/i" );
  tree_->Branch( "vertex_x", fVertexX, "vertex_x[num_vertex]/F" );
  tree_->Branch( "vertex_y", fVertexY, "vertex_y[num_vertex]/F" );
  tree_->Branch( "vertex_z", fVertexZ, "vertex_z[num_vertex]/F" );
  /*tree_->Branch( "vertex_tracks", fVertexTracks, "vertex_tracks[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_weight0p75", fVertexTracksWght0p75, "vertex_tracks_weight0p75[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_weight0p9", fVertexTracksWght0p9, "vertex_tracks_weight0p9[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_weight0p95", fVertexTracksWght0p95, "vertex_tracks_weight0p95[num_vertex]/i" );*/

  tree_->Branch( "met", &fMET );
  tree_->Branch( "met_phi", &fMETPhi );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer::endJob() 
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( TreeProducer );
