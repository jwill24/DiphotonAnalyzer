// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/TreeProducer
// Class:      TreeProducer
// 
/**\class TreeProducer TreeProducer.cc DiphotonAnalyzer/TreeProducer/plugins/TreeProducer.cc

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

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
//                               JW
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Jet.h"
//
#include "flashgg/DataFormats/interface/Met.h"

#include "DiphotonAnalyzer/TreeProducer/interface/HLTMatcher.h"
#include "DiphotonAnalyzer/TreeProducer/interface/SelectionUtils.h"
#include "DiphotonAnalyzer/TreeProducer/interface/XiInterpolator.h"
#include "DiphotonAnalyzer/TreeProducer/interface/FillNumberLUTHandler.h"
#include "DiphotonAnalyzer/TreeProducer/interface/AlignmentLUTHandler.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//
#define MAX_HLT 10
#define MAX_PROTON_TRK 20
#define MAX_DIPHOTON 10
//                               JW
#define MAX_ELECTRON 20
#define MAX_MUON 20
#define MAX_JET 50
//
#define MAX_VERTEX 100
#define MAX_GEN_PHOTON 10
#define MAX_GEN_PART 20

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit TreeProducer( const edm::ParameterSet& );
    ~TreeProducer();

    static void fillDescriptions( edm::ConfigurationDescriptions& );


  private:
    virtual void beginJob() override;
    virtual void beginRun( const edm::Run&, const edm::EventSetup& );
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    void clearTree();
    void analyzeTriggers( const edm::Event&, const edm::EventSetup& );
    
    edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > totemRPTracksToken_;
    edm::EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT< edm::View<flashgg::Met> > metToken_;
    edm::EDGetTokenT< edm::View<reco::Vertex> > vtxToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
//                                 JW
    edm::EDGetTokenT< edm::View<flashgg::Electron> > electronToken_;
    edm::EDGetTokenT< edm::View<flashgg::Muon> > muonToken_; 
    edm::EDGetTokenT< edm::View<vector<flashgg::Jet> > > jetToken_;
//
    edm::EDGetTokenT< edm::View<reco::GenParticle> > genPartToken_;
    edm::EDGetTokenT< edm::View<pat::PackedGenParticle> > genPhoToken_;
    edm::EDGetTokenT< edm::View<PileupSummaryInfo> > pileupToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

    bool isData_;    
    double sqrtS_;
    double singlePhotonMinPt_, singlePhotonMaxEta_, singlePhotonMinR9_;
    double photonPairMinMass_;
    std::string filename_;
    double maxGenLevelDR_;
    edm::FileInPath puMCfile_, puDatafile_;
    std::string puMCpath_, puDatapath_;

    std::string hltMenuLabel_;
    std::vector<std::string> triggersList_;
    HLTConfigProvider hlt_config_;
    HLTPrescaleProvider hlt_prescale_;
    HLTMatcher hlt_matcher_;

    bool useXiInterp_;
    std::unique_ptr<ProtonUtils::XiInterpolator> xiInterp_;

    std::unique_ptr<CTPPSAlCa::FillNumberLUTHandler> fillLUTHandler_;
    std::unique_ptr<CTPPSAlCa::AlignmentLUTHandler> alignmentLUTHandler_;
    std::unique_ptr<edm::LumiReWeighting> lumiReweighter_;

    TFile* file_;
    TTree* tree_;

    // --- tree components ---

    unsigned int fBX, fFill, fRun, fLumiSection;
    unsigned long long fEventNum;

    int fHLTNum;
    int fHLTAccept[MAX_HLT], fHLTPrescl[MAX_HLT];

    unsigned int fProtonTrackNum;
    float fProtonTrackX[MAX_PROTON_TRK], fProtonTrackY[MAX_PROTON_TRK];
    float fProtonTrackXCorr[MAX_PROTON_TRK], fProtonTrackYCorr[MAX_PROTON_TRK];
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

    unsigned int fGenPhotonNum;
    float fGenPhotonPt[MAX_GEN_PHOTON], fGenPhotonEta[MAX_GEN_PHOTON], fGenPhotonPhi[MAX_GEN_PHOTON], fGenPhotonE[MAX_GEN_PHOTON];
    float fGenPhotonVertexX[MAX_GEN_PHOTON], fGenPhotonVertexY[MAX_GEN_PHOTON], fGenPhotonVertexZ[MAX_GEN_PHOTON];

    unsigned int fGenPartNum;
    int fGenPartPDGId[MAX_GEN_PART], fGenPartStatus[MAX_GEN_PART];
    float fGenPartPt[MAX_GEN_PART], fGenPartEta[MAX_GEN_PART], fGenPartPhi[MAX_GEN_PART], fGenPartE[MAX_GEN_PART];
    float fGenPartVertexX[MAX_GEN_PART], fGenPartVertexY[MAX_GEN_PART], fGenPartVertexZ[MAX_GEN_PART];

    unsigned int fDiphotonNum;
    float fDiphotonPt1[MAX_DIPHOTON], fDiphotonPt2[MAX_DIPHOTON];
    float fDiphotonEta1[MAX_DIPHOTON], fDiphotonEta2[MAX_DIPHOTON];
    float fDiphotonPhi1[MAX_DIPHOTON], fDiphotonPhi2[MAX_DIPHOTON];
    float fDiphotonE1[MAX_DIPHOTON], fDiphotonE2[MAX_DIPHOTON];
    float fDiphotonR91[MAX_DIPHOTON], fDiphotonR92[MAX_DIPHOTON];
    float fDiphotonId1[MAX_DIPHOTON], fDiphotonId2[MAX_DIPHOTON];
    float fDiphotonSigEOverE1[MAX_DIPHOTON], fDiphotonSigEOverE2[MAX_DIPHOTON];
    float fDiphotonM[MAX_DIPHOTON], fDiphotonY[MAX_DIPHOTON];
    float fDiphotonPt[MAX_DIPHOTON], fDiphotonDphi[MAX_DIPHOTON];

    float fDiphotonSCX1[MAX_DIPHOTON], fDiphotonSCY1[MAX_DIPHOTON], fDiphotonSCZ1[MAX_DIPHOTON];
    float fDiphotonSCX2[MAX_DIPHOTON], fDiphotonSCY2[MAX_DIPHOTON], fDiphotonSCZ2[MAX_DIPHOTON];

    float fDiphotonGenPt1[MAX_DIPHOTON], fDiphotonGenPt2[MAX_DIPHOTON];
    float fDiphotonGenEta1[MAX_DIPHOTON], fDiphotonGenEta2[MAX_DIPHOTON];
    float fDiphotonGenPhi1[MAX_DIPHOTON], fDiphotonGenPhi2[MAX_DIPHOTON];
    float fDiphotonGenE1[MAX_DIPHOTON], fDiphotonGenE2[MAX_DIPHOTON];

    int fDiphotonVertex[MAX_DIPHOTON];
    unsigned int fDiphotonVertexTracks[MAX_DIPHOTON];
    unsigned int fDiphotonVerticesAt1mmDist[MAX_DIPHOTON], fDiphotonVerticesAt2mmDist[MAX_DIPHOTON], fDiphotonVerticesAt5mmDist[MAX_DIPHOTON], fDiphotonVerticesAt1cmDist[MAX_DIPHOTON];
    float fDiphotonVertexX[MAX_DIPHOTON], fDiphotonVertexY[MAX_DIPHOTON], fDiphotonVertexZ[MAX_DIPHOTON];
    float fDiphotonNearestDist[MAX_DIPHOTON];

    float fDiphotonGenVertexX, fDiphotonGenVertexY, fDiphotonGenVertexZ;
    float fDiphotonGenVertexSmearX, fDiphotonGenVertexSmearY, fDiphotonGenVertexSmearZ;

    float fMET, fMETPhi, fMETsignif;

    unsigned int fVertexNum;
    float fVertexX[MAX_VERTEX], fVertexY[MAX_VERTEX], fVertexZ[MAX_VERTEX];
    int fJetDiphoMatch[MAX_VERTEX];
    //unsigned int fVertexTracks[MAX_VERTEX], fVertexTracksWght0p75[MAX_VERTEX], fVertexTracksWght0p9[MAX_VERTEX], fVertexTracksWght0p95[MAX_VERTEX];
    float fBSX0, fBSY0, fBSZ0, fBSsigmaZ, fBSdxdz, fBSbeamWidthX, fBSbeamWidthY;

    float fPileupWeight;
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
  totemRPTracksToken_ ( consumes< edm::DetSetVector<TotemRPLocalTrack> >  ( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  diphotonToken_      ( consumes< edm::View<flashgg::DiPhotonCandidate> > ( iConfig.getParameter<edm::InputTag>( "diphotonLabel" ) ) ),
  metToken_           ( consumes< edm::View<flashgg::Met> >               ( iConfig.getParameter<edm::InputTag>( "metLabel") ) ),
  vtxToken_           ( consumes< edm::View<reco::Vertex> >               ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
  beamSpotToken_      ( consumes<reco::BeamSpot>                          ( iConfig.getParameter<edm::InputTag>( "beamSpotLabel" ) ) ),
//                               JW
  electronToken_      ( consumes< edm::View<flashgg::Electron> >          ( iConfig.getParameter<edm::InputTag>( "electronLabel") ) ),
  muonToken_          ( consumes< edm::View<flashgg::Muon> >              ( iConfig.getParameter<edm::InputTag>( "muonLabel") ) ),
  jetToken_           ( consumes< edm::View< std::vector<flashgg::Jet> > >( iConfig.getParameter<edm::InputTag>( "jetLabel") ) ),
//
  genPartToken_       ( consumes< edm::View<reco::GenParticle> >          ( iConfig.getParameter<edm::InputTag>( "genPartLabel" ) ) ),
  genPhoToken_        ( consumes< edm::View<pat::PackedGenParticle> >     ( iConfig.getParameter<edm::InputTag>( "genPhoLabel" ) ) ),
  pileupToken_        ( consumes< edm::View<PileupSummaryInfo> >          ( iConfig.getUntrackedParameter<edm::InputTag>( "pileupLabel", edm::InputTag( "slimmedAddPileupInfo" ) ) ) ),
  triggerResultsToken_( consumes<edm::TriggerResults>                     ( iConfig.getParameter<edm::InputTag>( "triggerResults" ) ) ),
  isData_             ( iConfig.getParameter<bool>       ( "isData" ) ),
  sqrtS_              ( iConfig.getParameter<double>     ( "sqrtS" ) ),
  singlePhotonMinPt_  ( iConfig.getParameter<double>     ( "minPtSinglePhoton" ) ),
  singlePhotonMaxEta_ ( iConfig.getParameter<double>     ( "maxEtaSinglePhoton" ) ),
  singlePhotonMinR9_  ( iConfig.getParameter<double>     ( "minR9SinglePhoton" ) ),
  photonPairMinMass_  ( iConfig.getParameter<double>     ( "minMassDiPhoton" ) ),
  filename_           ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  maxGenLevelDR_      ( iConfig.getParameter<double>     ( "maxGenLevelDeltaR" ) ),
  puMCpath_           ( iConfig.getUntrackedParameter<std::string>( "pileupMCPath", "pileup" ) ),
  puDatapath_         ( iConfig.getUntrackedParameter<std::string>( "pileupDataPath", "pileup" ) ),
  hltMenuLabel_       ( iConfig.getParameter<std::string>( "hltMenuLabel" ) ),
  triggersList_       ( iConfig.getParameter<std::vector<std::string> >( "triggersList" ) ),
  hlt_prescale_       ( iConfig, consumesCollector(), *this ),
  hlt_matcher_        ( iConfig.getParameter<std::vector<std::string> >( "triggersList" ) ),
  useXiInterp_        ( iConfig.getParameter<bool>       ( "useXiInterpolation" ) ),
  file_( 0 ), tree_( 0 )

{
  if ( isData_ ) {
    xiInterp_ = std::make_unique<ProtonUtils::XiInterpolator>();
    if ( useXiInterp_ ) {
      xiInterp_->loadInterpolationGraphs( iConfig.getParameter<edm::FileInPath>( "xiInterpolationFile" ).fullPath().c_str() );
    }
    alignmentLUTHandler_ = std::make_unique<CTPPSAlCa::AlignmentLUTHandler>( iConfig.getParameter<edm::FileInPath>( "alignmentLUTFile" ).fullPath().c_str() );
  }
  if ( !isData_ ) {
    puMCfile_ =iConfig.getParameter<edm::FileInPath>( "pileupMCFile" );
    puDatafile_ = iConfig.getParameter<edm::FileInPath>( "pileupDataFile" );
    std::cout << ">> Pileup reweighting will be used with files:\n\t" << puMCfile_.fullPath() << " for MC ;\n\t" << puDatafile_.fullPath() << " for data" << std::endl;
    lumiReweighter_ = std::make_unique<edm::LumiReWeighting>( puMCfile_.fullPath(), puDatafile_.fullPath(), puMCpath_, puDatapath_ );
  }
  else {
    fillLUTHandler_ = std::make_unique<CTPPSAlCa::FillNumberLUTHandler>( iConfig.getParameter<edm::FileInPath>( "fillNumLUTFile" ).fullPath().c_str() );
  }

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
  fBX = fRun = fLumiSection = fEventNum = 0;
  fHLTNum = 0;
  for ( unsigned int i=0; i<MAX_HLT; i++ ) {
    fHLTAccept[i] = fHLTPrescl[i] = -1;
  }
  fProtonTrackNum = 0;
  for ( unsigned int i=0; i<MAX_PROTON_TRK; i++ ) {
    fProtonTrackX[i] = fProtonTrackY[i] = -1.;
    fProtonTrackXCorr[i] = fProtonTrackYCorr[i] = -1.;
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
    fDiphotonE1[i] = fDiphotonE2[i] = -1.;
    fDiphotonR91[i] = fDiphotonR92[i] = -1.;
    fDiphotonId1[i] = fDiphotonId2[i] = -1.;
    fDiphotonSigEOverE1[i] = fDiphotonSigEOverE2[i] = -1.;
    fDiphotonM[i] = fDiphotonY[i] = fDiphotonPt[i] = fDiphotonDphi[i] = -1.;
    fDiphotonVertexTracks[i] = 0;
    fDiphotonVertex[i] = -1;
    fDiphotonVerticesAt1mmDist[i] = fDiphotonVerticesAt2mmDist[i] = fDiphotonVerticesAt5mmDist[i] = fDiphotonVerticesAt1cmDist[i] = 0;
    fDiphotonVertexX[i] = fDiphotonVertexY[i] = fDiphotonVertexZ[i] = -1.;
    fDiphotonNearestDist[i] = 999.;

    fDiphotonGenPt1[i] = fDiphotonGenPt2[i] = -1.;
    fDiphotonGenEta1[i] = fDiphotonGenEta2[i] = -1.;
    fDiphotonGenPhi1[i] = fDiphotonGenPhi2[i] = -1.;
    fDiphotonGenE1[i] = fDiphotonGenE2[i] = -1.;

    fDiphotonSCX1[i] = fDiphotonSCY1[i] = fDiphotonSCZ1[i] = -999.;
    fDiphotonSCX2[i] = fDiphotonSCY2[i] = fDiphotonSCZ2[i] = -999.;
  }
  fDiphotonGenVertexX = fDiphotonGenVertexY = fDiphotonGenVertexZ = -999.;
  fDiphotonGenVertexSmearX = fDiphotonGenVertexSmearY = fDiphotonGenVertexSmearZ = -999.;

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

  fGenPhotonNum = 0;
  for ( unsigned int i=0; i<MAX_GEN_PHOTON; i++ ) {
    fGenPhotonPt[i] = fGenPhotonEta[i] = fGenPhotonPhi[i] = fGenPhotonE[i] = -1.;
    fGenPhotonVertexX[i] = fGenPhotonVertexY[i] = fGenPhotonVertexZ[i] = -999.;
  }

  fGenPartNum = 0;
  for ( unsigned int i=0; i<MAX_GEN_PART; i++ ) {
    fGenPartPDGId[i] = 0;
    fGenPartStatus[i] = -999;
    fGenPartPt[i] = fGenPartEta[i] = fGenPartPhi[i] = fGenPartE[i] = -1.;
    fGenPartVertexX[i] = fGenPartVertexY[i] = fGenPartVertexZ[i] = -999.;
  }

  fMET = fMETPhi = fMETsignif = -1.;

  fVertexNum = 0;
  for ( unsigned int i=0; i<MAX_VERTEX; i++ ) {
    fVertexX[i] = fVertexY[i] = fVertexZ[i] = -999.;
    //fVertexTracks[i] = fVertexTracksWght0p75[i] = fVertexTracksWght0p9[i] = fVertexTracksWght0p95[i] = 0;
  }
  fBSX0 = fBSY0 = fBSZ0 = fBSsigmaZ = fBSdxdz = fBSbeamWidthX = fBSbeamWidthY = -999.;

  fPileupWeight = 1.;
}

// ------------ method called for each event  ------------
void
TreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  #include <iostream> // for debugging purposes

  analyzeTriggers( iEvent, iSetup );

  //std::cout << "---> passing the trigger!" << std::endl;

  clearTree();
  const reco::Candidate::Point orig( -999., -999, -999. );

  // Run and BX information
  fBX = iEvent.bunchCrossing();
  fRun = iEvent.id().run();
  fLumiSection = iEvent.luminosityBlock();
  fEventNum = iEvent.id().event();

  // get the fill number from the run id <-> fill number LUT
  if ( fillLUTHandler_ ) {
    fFill = fillLUTHandler_->getFillNumber( iEvent.id().run() );
  }
  else fFill = 1;

  //----- gen-level information -----

  edm::Handle< edm::View<pat::PackedGenParticle> > genPhotons;
  edm::Handle< edm::View<reco::GenParticle> > genParts;

  if ( !isData_ ) {
    iEvent.getByToken( genPhoToken_, genPhotons );
    for ( unsigned int i=0; i<genPhotons->size() && fGenPhotonNum<MAX_GEN_PHOTON; i++ ) {
      const edm::Ptr<pat::PackedGenParticle> genPho = genPhotons->ptrAt( i );

      //----- gen-level information

      if ( genPho->pdgId()!=22 ) continue; // only keep gen-level photons

      fGenPhotonPt[fGenPhotonNum] = genPho->pt();
      fGenPhotonEta[fGenPhotonNum] = genPho->eta();
      fGenPhotonPhi[fGenPhotonNum] = genPho->phi();
      fGenPhotonE[fGenPhotonNum] = genPho->energy();

      fGenPhotonVertexX[fGenPhotonNum] = genPho->vx();
      fGenPhotonVertexY[fGenPhotonNum] = genPho->vy();
      fGenPhotonVertexZ[fGenPhotonNum] = genPho->vz();

      fGenPhotonNum++;
    }

    //----- smeared-level information

    iEvent.getByToken( genPartToken_, genParts );
    for ( unsigned int i=0; i<genParts->size() && fGenPartNum<MAX_GEN_PART; i++ ) {
      const edm::Ptr<reco::GenParticle> genPart = genParts->ptrAt( i );

      fGenPartPDGId[i] = genPart->pdgId();
      fGenPartStatus[i] = genPart->status();

      fGenPartPt[i] = genPart->pt();
      fGenPartEta[i] = genPart->eta();
      fGenPartPhi[i] = genPart->phi();
      fGenPartE[i] = genPart->energy();

      fGenPartVertexX[i] = genPart->vx();
      fGenPartVertexY[i] = genPart->vy();
      fGenPartVertexZ[i] = genPart->vz();

      fGenPartNum++;
    }

  }
  //std::cout << "---> passing the mc matching!" << std::endl;

  //----- diphoton candidates -----

  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken( diphotonToken_, diphotons );

  edm::Handle< edm::View< std::vector<flashgg::Jet> > > jetsColls;
  iEvent.getByToken( jetToken_, jetsColls );

  fDiphotonNum = 0;
  fJetNum = 0;
  edm::Ptr<reco::Vertex> diphoton_vtx[MAX_DIPHOTON];
  for ( unsigned int i=0; i<diphotons->size() && fDiphotonNum<MAX_DIPHOTON; i++ ) {
    const edm::Ptr<flashgg::DiPhotonCandidate> diphoton = diphotons->ptrAt( i );

    if ( diphoton->leadPhotonId()<-0.9 ) continue;
    if ( diphoton->subLeadPhotonId()<-0.9 ) continue;

    const flashgg::Photon* lead_pho = diphoton->leadingPhoton(),
                          *sublead_pho = diphoton->subLeadingPhoton();
    if ( !passSinglePhotonCuts( lead_pho ) ) continue;
    if ( !passSinglePhotonCuts( sublead_pho ) ) continue;

    if ( fabs( lead_pho->eta() )>=singlePhotonMaxEta_ or fabs( sublead_pho->eta() )>=singlePhotonMaxEta_ ) continue;
    if ( lead_pho->pt()<singlePhotonMinPt_ or sublead_pho->pt()<singlePhotonMinPt_ ) continue;
    if ( lead_pho->full5x5_r9()<singlePhotonMinR9_ or sublead_pho->full5x5_r9()<singlePhotonMinR9_ ) continue;

    if ( diphoton->mass()<photonPairMinMass_ ) continue;

    fDiphotonVertexTracks[fDiphotonNum] = diphoton->vtx()->nTracks();
    fDiphotonVertexX[fDiphotonNum] = diphoton->vtx()->x();
    fDiphotonVertexY[fDiphotonNum] = diphoton->vtx()->y();
    fDiphotonVertexZ[fDiphotonNum] = diphoton->vtx()->z();
    diphoton_vtx[fDiphotonNum] = diphoton->vtx();

    fDiphotonPt1[fDiphotonNum] = lead_pho->pt();
    fDiphotonEta1[fDiphotonNum] = lead_pho->eta();
    fDiphotonPhi1[fDiphotonNum] = lead_pho->phi();
    fDiphotonE1[fDiphotonNum] = lead_pho->energy();
    fDiphotonR91[fDiphotonNum] = lead_pho->full5x5_r9();
    fDiphotonId1[fDiphotonNum] = lead_pho->phoIdMvaDWrtVtx( diphoton->vtx() );
    fDiphotonSigEOverE1[fDiphotonNum] = lead_pho->sigEOverE();

    fDiphotonPt2[fDiphotonNum] = sublead_pho->pt();
    fDiphotonEta2[fDiphotonNum] = sublead_pho->eta();
    fDiphotonPhi2[fDiphotonNum] = sublead_pho->phi();
    fDiphotonE2[fDiphotonNum] = sublead_pho->energy();
    fDiphotonR92[fDiphotonNum] = sublead_pho->full5x5_r9();
    fDiphotonId2[fDiphotonNum] = sublead_pho->phoIdMvaDWrtVtx( diphoton->vtx() );
    fDiphotonSigEOverE2[fDiphotonNum] = sublead_pho->sigEOverE();

    fDiphotonM[fDiphotonNum] = diphoton->mass();
    fDiphotonY[fDiphotonNum] = diphoton->rapidity();
    fDiphotonPt[fDiphotonNum] = diphoton->pt();

    float dphi = lead_pho->phi()-sublead_pho->phi();
    while ( dphi<-TMath::Pi() ) dphi += 2.*TMath::Pi();
    while ( dphi> TMath::Pi() ) dphi -= 2.*TMath::Pi();
    fDiphotonDphi[fDiphotonNum] = dphi;

    //----- retrieve the SC information

    fDiphotonSCX1[fDiphotonNum] = lead_pho->superCluster()->x();
    fDiphotonSCY1[fDiphotonNum] = lead_pho->superCluster()->y();
    fDiphotonSCZ1[fDiphotonNum] = lead_pho->superCluster()->z();
    fDiphotonSCX2[fDiphotonNum] = sublead_pho->superCluster()->x();
    fDiphotonSCY2[fDiphotonNum] = sublead_pho->superCluster()->y();
    fDiphotonSCZ2[fDiphotonNum] = sublead_pho->superCluster()->z();

    //----- retrieve the gen-level information

    if ( !isData_ ) {
      reco::Candidate::Point vtx1( orig ), vtx2( orig ),
                             vtx1_smear( orig ), vtx2_smear( orig );
      if ( lead_pho->matchedGenPhoton() ) {
        fDiphotonGenPt1[fDiphotonNum] = lead_pho->matchedGenPhoton()->pt();
        fDiphotonGenEta1[fDiphotonNum] = lead_pho->matchedGenPhoton()->eta();
        fDiphotonGenPhi1[fDiphotonNum] = lead_pho->matchedGenPhoton()->phi();
        fDiphotonGenE1[fDiphotonNum] = lead_pho->matchedGenPhoton()->energy();
        vtx1 = lead_pho->matchedGenPhoton()->vertex();
      }
      if ( sublead_pho->matchedGenPhoton() ) {
        fDiphotonGenPt2[fDiphotonNum] = sublead_pho->matchedGenPhoton()->pt();
        fDiphotonGenEta2[fDiphotonNum] = sublead_pho->matchedGenPhoton()->eta();
        fDiphotonGenPhi2[fDiphotonNum] = sublead_pho->matchedGenPhoton()->phi();
        fDiphotonGenE2[fDiphotonNum] = sublead_pho->matchedGenPhoton()->energy();
        vtx2 = sublead_pho->matchedGenPhoton()->vertex();
      }
      if ( vtx2!=vtx1 ) {
        std::cerr << "-> 2 different gen-level vertices: " << std::endl
                  << "      leading photon vertex: (" << vtx1.x() << ", " << vtx1.y() << ", " << vtx1.z() << ")" << std::endl
                  << "   subleading photon vertex: (" << vtx2.x() << ", " << vtx2.y() << ", " << vtx2.z() << ")" << std::endl;
      }
      if ( vtx1!=orig ) {
        fDiphotonGenVertexX = vtx1.x();
        fDiphotonGenVertexY = vtx1.y();
        fDiphotonGenVertexZ = vtx1.z();
      }
      else {
        fDiphotonGenVertexX = vtx2.x();
        fDiphotonGenVertexY = vtx2.y();
        fDiphotonGenVertexZ = vtx2.z();
      }

      //----- retrieve the smeared-level generated vertex

      for ( unsigned int j=0; j<genParts->size(); j++ ) {
        const edm::Ptr<reco::GenParticle> genPart = genParts->ptrAt( j );
        if ( sqrt( pow( genPart->eta()-   lead_pho->eta(), 2 ) + pow( genPart->phi()-   lead_pho->phi(), 2 ) )<maxGenLevelDR_ ) { vtx1_smear = genPart->vertex(); continue; }
        if ( sqrt( pow( genPart->eta()-sublead_pho->eta(), 2 ) + pow( genPart->phi()-sublead_pho->phi(), 2 ) )<maxGenLevelDR_ ) { vtx2_smear = genPart->vertex(); continue; }
      }
      /*if ( vtx2_smear!=vtx1_smear ) {
        std::cerr << "-> 2 different smeared-level vertices: " << std::endl
                  << "      leading photon vertex: (" << vtx1_smear.x() << ", " << vtx1_smear.y() << ", " << vtx1_smear.z() << ")" << std::endl
                  << "   subleading photon vertex: (" << vtx2_smear.x() << ", " << vtx2_smear.y() << ", " << vtx2_smear.z() << ")" << std::endl;
      }*/
      if ( vtx1_smear!=orig ) {
        fDiphotonGenVertexSmearX = vtx1_smear.x();
        fDiphotonGenVertexSmearY = vtx1_smear.y();
        fDiphotonGenVertexSmearZ = vtx1_smear.z();
      }
      else {
        fDiphotonGenVertexSmearX = vtx2_smear.x();
        fDiphotonGenVertexSmearY = vtx2_smear.y();
        fDiphotonGenVertexSmearZ = vtx2_smear.z();
      }
    }

    //----- retrieve the associated jets

    if ( diphoton->jetCollectionIndex()<jetsColls->size() ) {
      for ( unsigned int j=0; j<jetsColls->at( diphoton->jetCollectionIndex() ).size(); j++ ) {
        if ( fJetNum>=MAX_JET ) { std::cerr << ">> More jets than expected in this event (" << fJetNum << ">MAX_JET=" << MAX_JET << "). Increase MAX_JET for safety" << std::endl; }

        const flashgg::Jet jet = jetsColls->at( diphoton->jetCollectionIndex() ).at( j );
        fJetPt[fJetNum]   = jet.pt();
        fJetEta[fJetNum]  = jet.eta();
        fJetPhi[fJetNum]  = jet.phi();
        fJetE[fJetNum]    = jet.energy();
        fJetMass[fJetNum] = jet.mass();

        fJetVertexX[fJetNum] = jet.vertex().x();
        fJetVertexY[fJetNum] = jet.vertex().y();
        fJetVertexZ[fJetNum] = jet.vertex().z();
        fJetDiphoMatch[fJetNum] = i;
        fJetNum++;
      }
    }

    fDiphotonNum++;
  }

  //----- only store events with > 0 diphoton candidate(s) -----

  if ( fDiphotonNum<1 ) return;

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

  if ( isData_ ) {
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

        fProtonTrackX[fProtonTrackNum] = trk->getX0() * 1.e-3; // store in m
        fProtonTrackY[fProtonTrackNum] = trk->getY0() * 1.e-3; // store in m
        fProtonTrackXCorr[fProtonTrackNum] = ( trk->getX0() + align_quant.x ) * 1.e-3; // store in m
        fProtonTrackYCorr[fProtonTrackNum] = ( trk->getY0() - align_quant.y ) * 1.e-3; // store in m
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
  } 
  //                               JW
 
  //----- electrons collection -----
 
  edm::Handle< edm::View<flashgg::Electron> > electrons;
  iEvent.getByToken( electronToken_, electrons );
 
  fElectronNum=0;
  for ( unsigned int i=0; i<electrons->size() && fElectronNum<MAX_ELECTRON; i++ ) {
    const edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( i );
    if ( fElectronNum>=MAX_ELECTRON ) { std::cerr << ">> More electrons than expected in this event (" << fElectronNum << ">MAX_ELECTRON=" << MAX_ELECTRON << "). Increase MAX_ELECTRON for safety" << std::endl; }

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
    if ( fMuonNum>=MAX_MUON ) { std::cerr << ">> More muons than expected in this event (" << fMuonNum << ">MAX_MUON=" << MAX_MUON << "). Increase MAX_MUON for safety" << std::endl; }

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
  fMET = met->et();
  fMETPhi = met->phi();
  fMETsignif = met->significance();

  //std::cout << "passed MET" << std::endl;
  //--- pileup information ---
  if ( !isData_ ) {
    edm::Handle< edm::View<PileupSummaryInfo> > pileupInfo;
    iEvent.getByToken( pileupToken_, pileupInfo );
    int npv0true = 0;
    for ( unsigned int i=0; i<pileupInfo->size(); i++ ) {
      const edm::Ptr<PileupSummaryInfo> PVI = pileupInfo->ptrAt( i );
      const int beamXing = PVI->getBunchCrossing(),
                npvtrue = PVI->getTrueNumInteractions();
      if ( beamXing==0 ) { npv0true += npvtrue; }
    }
    fPileupWeight = lumiReweighter_->weight( npv0true );
  }

  //std::cout << "passed PU" << std::endl;
  //----- vertexing information -----

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );

  fVertexNum = 0;
  for ( unsigned int i=0; i<vertices->size() && i<MAX_VERTEX; i++ ) {
    const edm::Ptr<reco::Vertex> vtx = vertices->ptrAt( i );
    if ( fVertexNum>=MAX_VERTEX ) { std::cerr << ">> More vertices than expected in this event (" << fVertexNum << ">MAX_VERTEX=" << MAX_VERTEX << "). Increase MAX_VERTEX for safety" << std::endl; }
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

  std::cout << ">> Filling the tree" << std::endl;
  tree_->Fill();
}

void
TreeProducer::analyzeTriggers( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // Get the trigger information from the event
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken( triggerResultsToken_, hltResults );

  const edm::TriggerNames& trigNames = iEvent.triggerNames( *hltResults );

  std::ostringstream os;
  os << "Trigger names: " << std::endl;
  for ( unsigned int i=0; i<trigNames.size(); i++ ) {

    const int trigNum = hlt_matcher_.TriggerNum( trigNames.triggerNames().at( i ) );
    if ( trigNum<0 ) continue; // Trigger didn't match the interesting ones
    fHLTAccept[trigNum] = hltResults->accept( i );

    // extract prescale value for this path
    fHLTPrescl[trigNum] = 1;
    if ( isData_ ) { //FIXME
      int prescale_set = hlt_prescale_.prescaleSet( iEvent, iSetup );
      fHLTPrescl[trigNum] = ( prescale_set<0 ) ? 0. : hlt_config_.prescaleValue( prescale_set, trigNames.triggerNames().at( i ) ); //FIXME
    }
    os << "--> (fired? " << fHLTAccept[trigNum] << ") " << trigNames.triggerNames().at( i ) << "\tprescale: " << fHLTPrescl[trigNum] << std::endl;
  }
  LogDebug( "TreeProducer" ) << os.str();
  //#include <iostream>
  //std::cout << os.str() << std::endl;
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

  tree_->Branch("num_hlt", &fHLTNum, "num_hlt/I");
  tree_->Branch("hlt_accept", fHLTAccept, "hlt_accept[num_hlt]/I");
  tree_->Branch("hlt_prescale", fHLTPrescl, "hlt_prescale[num_hlt]/I");

  /*std::vector<std::string>* HLT_Name;
  tree_->Branch("hlt_name", &HLT_Name);
  *HLT_Name = triggersList_;*/

  if ( isData_ ) {
    tree_->Branch( "num_proton_track", &fProtonTrackNum, "num_proton_track/i" );
    tree_->Branch( "proton_track_x", fProtonTrackX, "proton_track_x[num_proton_track]/F" );
    tree_->Branch( "proton_track_y", fProtonTrackY, "proton_track_y[num_proton_track]/F" );
    tree_->Branch( "proton_track_x_corr", fProtonTrackXCorr, "proton_track_x_corr[num_proton_track]/F" );
    tree_->Branch( "proton_track_y_corr", fProtonTrackYCorr, "proton_track_y_corr[num_proton_track]/F" );
    tree_->Branch( "proton_track_xi", fProtonTrackXi, "proton_track_xi[num_proton_track]/F" );
    tree_->Branch( "proton_track_xi_error", fProtonTrackXiError, "proton_track_xi_error[num_proton_track]/F" );
    tree_->Branch( "proton_track_side", fProtonTrackSide, "proton_track_side[num_proton_track]/i" );
    tree_->Branch( "proton_track_pot", fProtonTrackPot, "proton_track_pot[num_proton_track]/i" );
    tree_->Branch( "proton_track_link_nearfar", fProtonTrackLinkNF, "proton_track_link_nearfar[num_proton_track]/i" );
    tree_->Branch( "proton_track_link_mindist", fProtonTrackMinLinkDist, "proton_track_link_mindist[num_proton_track]/F" );
  }
  if ( !isData_ ) {
    tree_->Branch( "num_gen_photon", &fGenPhotonNum, "num_gen_photon/i" );
    tree_->Branch( "gen_photon_pt", fGenPhotonPt, "gen_photon_pt[num_gen_photon]/F" );
    tree_->Branch( "gen_photon_eta", fGenPhotonEta, "gen_photon_eta[num_gen_photon]/F" );
    tree_->Branch( "gen_photon_phi", fGenPhotonPhi, "gen_photon_phi[num_gen_photon]/F" );
    tree_->Branch( "gen_photon_energy", fGenPhotonE, "gen_photon_energy[num_gen_photon]/F" );
    tree_->Branch( "gen_photon_vertex_x", fGenPhotonVertexX, "gen_photon_vertex_x[num_gen_photon]/F" );
    tree_->Branch( "gen_photon_vertex_y", fGenPhotonVertexY, "gen_photon_vertex_y[num_gen_photon]/F" );
    tree_->Branch( "gen_photon_vertex_z", fGenPhotonVertexZ, "gen_photon_vertex_z[num_gen_photon]/F" );

    tree_->Branch( "num_gen_part", &fGenPartNum, "num_gen_part/i" );
    tree_->Branch( "gen_part_pdgid", fGenPartPDGId, "gen_part_pdgid[num_gen_part]/I" );
    tree_->Branch( "gen_part_status", fGenPartStatus, "gen_part_status[num_gen_part]/I" );
    tree_->Branch( "gen_part_pt", fGenPartPt, "gen_part_pt[num_gen_part]/F" );
    tree_->Branch( "gen_part_eta", fGenPartEta, "gen_part_eta[num_gen_part]/F" );
    tree_->Branch( "gen_part_phi", fGenPartPhi, "gen_part_phi[num_gen_part]/F" );
    tree_->Branch( "gen_part_energy", fGenPartE, "gen_part_energy[num_gen_part]/F" );
    tree_->Branch( "gen_part_vertex_x", fGenPartVertexX, "gen_part_vertex_x[num_gen_part]/F" );
    tree_->Branch( "gen_part_vertex_y", fGenPartVertexY, "gen_part_vertex_y[num_gen_part]/F" );
    tree_->Branch( "gen_part_vertex_z", fGenPartVertexZ, "gen_part_vertex_z[num_gen_part]/F" );
  }

  tree_->Branch( "num_diphoton", &fDiphotonNum, "num_diphoton/i" );
  tree_->Branch( "diphoton_pt1", fDiphotonPt1, "diphoton_pt1[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt2", fDiphotonPt2, "diphoton_pt2[num_diphoton]/F" );
  tree_->Branch( "diphoton_eta1", fDiphotonEta1, "diphoton_eta1[num_diphoton]/F" );
  tree_->Branch( "diphoton_eta2", fDiphotonEta2, "diphoton_eta2[num_diphoton]/F" );
  tree_->Branch( "diphoton_phi1", fDiphotonPhi1, "diphoton_phi1[num_diphoton]/F" );
  tree_->Branch( "diphoton_phi2", fDiphotonPhi2, "diphoton_phi2[num_diphoton]/F" );
  tree_->Branch( "diphoton_energy1", fDiphotonE1, "diphoton_energy1[num_diphoton]/F" );
  tree_->Branch( "diphoton_energy2", fDiphotonE2, "diphoton_energy2[num_diphoton]/F" );
  tree_->Branch( "diphoton_r91", fDiphotonR91, "diphoton_r91[num_diphoton]/F" );
  tree_->Branch( "diphoton_r92", fDiphotonR92, "diphoton_r92[num_diphoton]/F" );
  tree_->Branch( "diphoton_mass", fDiphotonM, "diphoton_mass[num_diphoton]/F" );
  tree_->Branch( "diphoton_rapidity", fDiphotonY, "diphoton_rapidity[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt", fDiphotonPt, "diphoton_pt[num_diphoton]/F" );
  tree_->Branch( "diphoton_dphi", fDiphotonDphi, "diphoton_dphi[num_diphoton]/F" );
  tree_->Branch( "diphoton_id1", fDiphotonId1, "diphoton_id1[num_diphoton]/F" );
  tree_->Branch( "diphoton_id2", fDiphotonId2, "diphoton_id2[num_diphoton]/F" );
  tree_->Branch( "diphoton_sigeove1", fDiphotonSigEOverE1, "diphoton_sigeove1[num_diphoton]/F" );
  tree_->Branch( "diphoton_sigeove2", fDiphotonSigEOverE2, "diphoton_sigeove2[num_diphoton]/F" );

  tree_->Branch( "diphoton_supercluster_x1", fDiphotonSCX1, "diphoton_supercluster_x1[num_diphoton]/F" );
  tree_->Branch( "diphoton_supercluster_y1", fDiphotonSCY1, "diphoton_supercluster_y1[num_diphoton]/F" );
  tree_->Branch( "diphoton_supercluster_z1", fDiphotonSCZ1, "diphoton_supercluster_z1[num_diphoton]/F" );
  tree_->Branch( "diphoton_supercluster_x2", fDiphotonSCX2, "diphoton_supercluster_x2[num_diphoton]/F" );
  tree_->Branch( "diphoton_supercluster_y2", fDiphotonSCY2, "diphoton_supercluster_y2[num_diphoton]/F" );
  tree_->Branch( "diphoton_supercluster_z2", fDiphotonSCZ2, "diphoton_supercluster_z2[num_diphoton]/F" );

  if ( !isData_ ) {
    tree_->Branch( "diphoton_genpt1", fDiphotonGenPt1, "diphoton_genpt1[num_diphoton]/F" );
    tree_->Branch( "diphoton_genpt2", fDiphotonGenPt2, "diphoton_genpt2[num_diphoton]/F" );
    tree_->Branch( "diphoton_geneta1", fDiphotonGenEta1, "diphoton_geneta1[num_diphoton]/F" );
    tree_->Branch( "diphoton_geneta2", fDiphotonGenEta2, "diphoton_geneta2[num_diphoton]/F" );
    tree_->Branch( "diphoton_genphi1", fDiphotonGenPhi1, "diphoton_genphi1[num_diphoton]/F" );
    tree_->Branch( "diphoton_genphi2", fDiphotonGenPhi2, "diphoton_genphi2[num_diphoton]/F" );
    tree_->Branch( "diphoton_genenergy1", fDiphotonGenE1, "diphoton_genenergy1[num_diphoton]/F" );
    tree_->Branch( "diphoton_genenergy2", fDiphotonGenE2, "diphoton_genenergy2[num_diphoton]/F" );
  }

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

  if ( !isData_ ) {
    tree_->Branch( "diphoton_genvertex_x", &fDiphotonGenVertexX, "diphoton_genvertex_x/F" );
    tree_->Branch( "diphoton_genvertex_y", &fDiphotonGenVertexY, "diphoton_genvertex_y/F" );
    tree_->Branch( "diphoton_genvertex_z", &fDiphotonGenVertexZ, "diphoton_genvertex_z/F" );
    tree_->Branch( "diphoton_genvertex_smeared_x", &fDiphotonGenVertexSmearX, "diphoton_genvertex_smeared_x/F" );
    tree_->Branch( "diphoton_genvertex_smeared_y", &fDiphotonGenVertexSmearY, "diphoton_genvertex_smeared_y/F" );
    tree_->Branch( "diphoton_genvertex_smeared_z", &fDiphotonGenVertexSmearZ, "diphoton_genvertex_smeared_z/F" );
  }

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
  tree_->Branch( "jet_dipho_match", fJetDiphoMatch, "jet_dipho_match[num_jet]/I" );

  tree_->Branch( "num_vertex", &fVertexNum, "num_vertex/i" );
  tree_->Branch( "vertex_x", fVertexX, "vertex_x[num_vertex]/F" );
  tree_->Branch( "vertex_y", fVertexY, "vertex_y[num_vertex]/F" );
  tree_->Branch( "vertex_z", fVertexZ, "vertex_z[num_vertex]/F" );
  /*tree_->Branch( "vertex_tracks", fVertexTracks, "vertex_tracks[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_weight0p75", fVertexTracksWght0p75, "vertex_tracks_weight0p75[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_weight0p9", fVertexTracksWght0p9, "vertex_tracks_weight0p9[num_vertex]/i" );
  tree_->Branch( "vertex_tracks_weight0p95", fVertexTracksWght0p95, "vertex_tracks_weight0p95[num_vertex]/i" );*/

  tree_->Branch( "met", &fMET, "met/F" );
  tree_->Branch( "met_phi", &fMETPhi, "met_phi/F" );
  tree_->Branch( "met_significance", &fMETsignif, "met_significance/F" );

  tree_->Branch( "bs_x0", &fBSX0, "bs_x0/F" );
  tree_->Branch( "bs_y0", &fBSY0, "bs_y0/F" );
  tree_->Branch( "bs_z0", &fBSZ0, "bs_z0/F" );
  tree_->Branch( "bs_sigma_z", &fBSsigmaZ, "bs_sigma_z/F" );
  tree_->Branch( "bs_dxdz", &fBSdxdz, "bs_dxdz/F" );
  tree_->Branch( "bs_beam_width_x", &fBSbeamWidthX, "bs_beam_width_x/F" );
  tree_->Branch( "bs_beam_width_y", &fBSbeamWidthY, "bs_beam_width_y/F" );

  tree_->Branch( "pileup_weight", &fPileupWeight, "pileup_weight/F" );
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer::endJob() 
{}

void
TreeProducer::beginRun( const edm::Run& iRun, const edm::EventSetup& iSetup )
{
  std::string processName = "";
  bool changed = true;
  if ( !hlt_prescale_.init( iRun, iSetup, hltMenuLabel_, changed ) ) {
    throw cms::Exception("TreeProducer") << " prescales extraction failure with process name " << hltMenuLabel_;
  }
  hlt_config_ = hlt_prescale_.hltConfigProvider();
  if ( !hlt_config_.init( iRun, iSetup, hltMenuLabel_, changed ) ) {
    throw cms::Exception("TreeProducer") << " HLT config extraction failure with process name " << processName;
  }
  if ( hlt_config_.size()<=0 ) {
    edm::LogError("TreeProducer") << "HLT config size error";
  }
  if (changed) {
    // The HLT config has actually changed wrt the previous Run, hence rebook your
    // histograms or do anything else dependent on the revised HLT config
  }
}

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
