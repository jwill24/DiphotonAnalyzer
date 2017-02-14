import FWCore.ParameterSet.Config as cms

treeProducer = cms.EDAnalyzer('TreeProducer',
    sqrtS = cms.double(13.e3),
    metLabel = cms.InputTag('slimmedMETs'),
    diphotonLabel = cms.InputTag('flashggDiPhotons'),
    vertexLabel = cms.InputTag('offlineSlimmedPrimaryVertices'),
    electronLabel = cms.InputTag('flashggSelectedElectrons'),
    muonLabel = cms.InputTag('flashggSelectedMuons'),
    jetLabel = cms.InputTag('flashggFinalJets'),
    minPtSinglePhoton = cms.double(50.),
    minR9SinglePhoton = cms.double(0.94),
    maxEtaSinglePhoton = cms.double(2.5),
    minMassDiPhoton = cms.double(500.),
    outputFilename = cms.string('output.root'),
    # totem RP information extraction
    totemRPTracksLabel = cms.InputTag('totemRPLocalTrackFitter'),
    useXiInterpolation = cms.bool(True),
    xiInterpolationFile = cms.FileInPath('DiphotonAnalyzer/EventAnalyzer/data/ctpps_optics_v1.root'),
    fillNumLUTFile = cms.FileInPath('DiphotonAnalyzer/EventAnalyzer/data/fill_run_lut.dat'),
    alignmentLUTFile = cms.FileInPath('DiphotonAnalyzer/EventAnalyzer/data/alignment_collection_v2.out'),
)
