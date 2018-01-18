import FWCore.ParameterSet.Config as cms

mbTreeProducer = cms.EDAnalyzer('MBTreeProducer',

    # general parameters
    outputFilename = cms.string('output.root'),

    # input collections
    totemRPTracksLabel = cms.InputTag('totemRPLocalTrackFitter'),
)
