import FWCore.ParameterSet.Config as cms

mbTreeProducer = cms.EDAnalyzer('MBTreeProducer',

    # general parameters
    outputFilename = cms.string('output.root'),

    # input collections
    vertexLabel = cms.InputTag('offlinePrimaryVertices'),
    beamSpotLabel = cms.InputTag('offlineBeamSpot'),
    totemRPTracksLabel = cms.InputTag('totemRPLocalTrackFitter'),

    # totem RP information extraction
    fillNumLUTFile = cms.FileInPath('DiphotonAnalyzer/TreeProducer/data/fill_run_lut_v2.dat'),
)
