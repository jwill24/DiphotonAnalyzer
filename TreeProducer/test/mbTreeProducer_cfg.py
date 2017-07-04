import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/data/Run2016B/MinimumBias/AOD/23Sep2016-v1/50000/040C2700-A388-E611-AF56-00259073E390.root',
    )
)
# Trigger

#from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltHighLevel.HLTPaths = ['HLT_L1MinimumBias*']
process.hltHighLevel.throw = cms.bool(False)

process.load('DiphotonAnalyzer.TreeProducer.mbTreeProducer_cfi')

# set some parameters to the run

process.p = cms.Path(
    process.hltHighLevel*
    process.mbTreeProducer
)
