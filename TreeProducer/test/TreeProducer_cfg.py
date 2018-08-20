import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
#'/store/group/phys_higgs/cmshgg/lforthom/flashgg/pps_run2016/Moriond16WSFinal-106-g90923ae/DoubleEG/pps_run2016-Moriond16WSFinal-106-g90923ae-v0-Run2016B-PromptReco-v2/160624_003754/0000/myMicroAODOutputFile_1.root',
#'/store/group/phys_higgs/cmshgg/lforthom/flashgg/pps_run2016/Moriond16WSFinal-106-g90923ae/DoubleEG/pps_run2016-Moriond16WSFinal-106-g90923ae-v0-Run2016B-PromptReco-v2/160624_003754/0000/myMicroAODOutputFile_2.root',
#'/store/group/phys_higgs/cmshgg/lforthom/flashgg/DoubleEG/pps_lforthom-miniAOD_run2016B_v2/160831_083550/0000/myMicroAODOutputFile_881.root',
#'/store/group/phys_higgs/cmshgg/lforthom/flashgg/DoubleEG/pps_lforthom-miniAOD_run2016B_v5/161208_215044/0000/myMicroAODOutputFile_881.root'
#'/store/group/phys_pps/diphoton/DoubleEG/lforthom-microAOD-ctpps_Run2016C-23Sep2016_v3/170303_022624/0000/myMicroAODOutputFile_64.root',
#'/store/group/phys_pps/diphoton/juwillia/GGToGG_bSM-AA1e-14_Pt-50_M-300_13TeV-fpmc-herwig6/GGtoGG_e-14_realistic_microAOD/180426_225340/0000/myMicroAODOutputFile_24.root'
'/store/group/phys_pps/diphoton/DoubleEG/DoubleEG_microAOD_run2017B/180517_175342/0000/myMicroAODOutputFile_2.root'
#'file:/afs/cern.ch/work/j/juwillia/CMSSW_9_4_5_cand1/src/DiphotonAnalyzer/myMicroAODOutputFile_1.root',
#'file:/afs/cern.ch/work/j/juwillia/CMSSW_9_4_5_cand1/src/DiphotonAnalyzer/myMicroAODOutputFile_2.root',
#'file:/afs/cern.ch/work/j/juwillia/CMSSW_9_4_5_cand1/src/DiphotonAnalyzer/myMicroAODOutputFile_3.root',
#'file:/afs/cern.ch/work/j/juwillia/CMSSW_9_4_5_cand1/src/DiphotonAnalyzer/myMicroAODOutputFile_4.root'
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
process.treeProducer.minMassDiPhoton = cms.double(350.)
process.treeProducer.minR9SinglePhoton = cms.double(0.)
process.treeProducer.triggersList = process.hltHighLevel.HLTPaths

process.p = cms.Path(
    process.hltHighLevel*
    process.treeProducer
)
