#cmsRun runTauEfficiency.py inputFiles=/store/relval/CMSSW_8_1_0_pre12/RelValZTT_14TeV/MINIAODSIM/81X_mcRun2_asymptotic_v8_2023D1-v1/00000/4CC5E70D-8B8F-E611-8079-0CC47A4D76D0.root

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')

options.outputFile = "MiniAOD_effi_80x_DYtoLL.root"
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '81X_mcRun2_asymptotic_v8', '')

#how many events to run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)


##################################################
# Main
process.cutBased = cms.EDAnalyzer("phase2TausRECO",
    tracks             = cms.InputTag("generalTracks"),#checkme
    times              = cms.InputTag(""),
    timesResos         = cms.InputTag(""),
    vertices           = cms.InputTag("offlineSlimmedPrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discrimininator    = cms.InputTag("byCombinedIsolationDeltaBetaCorrRaw3Hits")
#    genParticles       = cms.InputTag("prunedGenParticles")
)

#process.MVA = cms.EDAnalyzer("phase2Taus",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    taus = cms.InputTag("slimmedTaus"),
#    tauID = cms.string("byLooseIsolationMVArun2v1DBoldDMwLT"),
#    pruned = cms.InputTag("prunedGenParticles")
#)

###################################################
#Global sequence

process.p = cms.Path(
         process.cutBased
 #        process.MVA
                     )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)

#print out all processes used when running- useful check to see if module ran
#UNCOMMENT BELOW
#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
