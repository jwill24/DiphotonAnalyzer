import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'/store/group/phys_pps/diphoton/DoubleEG/lforthom-microAOD-ctpps_Run2016C-23Sep2016_v3/170303_022624/0000/myMicroAODOutputFile_64.root',
'/store/group/phys_pps/diphoton/DoubleEG/lforthom-microAOD-ctpps_Run2016H_postTS2-07Aug2017_v1/180117_094639/0000/myMicroAODOutputFile_120.root'
    )
)

# Trigger
#from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltHighLevel.HLTPaths = ['HLT_DoublePhoton60*', 'HLT_DoublePhoton85*']
process.hltHighLevel.throw = cms.bool(False)

process.load('DiphotonAnalyzer.TreeProducer.TreeProducer_cfi')

# set some parameters to the run
process.treeProducer.minPtSinglePhoton = cms.double(50.)
process.treeProducer.minMassDiPhoton = cms.double(50.)
#process.treeProducer.maxMassDiPhoton = cms.untracked.double(200.)
process.treeProducer.minR9SinglePhoton = cms.double(0.)
process.treeProducer.triggersList = process.hltHighLevel.HLTPaths

process.p = cms.Path(
    process.hltHighLevel*
    process.treeProducer
)
