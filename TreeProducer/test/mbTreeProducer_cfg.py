import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'',
    )
)

process.load('DiphotonAnalyzer.TreeProducer.mbTreeProducer_cfi')

# set some parameters to the run

process.p = cms.Path(
    process.mbTreeProducer
)
