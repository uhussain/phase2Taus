#cmsRun runTauTiming.py inputFiles=/store/relval/CMSSW_8_1_0_pre12/RelValZTT_14TeV/MINIAODSIM/81X_mcRun2_asymptotic_v8_2023D1-v1/00000/4CC5E70D-8B8F-E611-8079-0CC47A4D76D0.root

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')

options.outputFile = "timing-ttbar.root"
#options.register('outputFile','timing.root',VarParsing.multiplicity.singleton, VarParsing.varType.string,'outputfile name is timing.root if one is not given')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")
process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '81X_mcRun2_asymptotic_v8', '')

#how many events to run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

if len(options.inputFileList) > 0 :
    with open(options.inputFileList) as f :
        inputFiles = list((line.strip() for line in f))
else :
    inputFiles = cms.untracked.vstring(options.inputFiles)
    
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles),
)


##################################################
# Main
process.VLooseMVA = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.LooseMVA = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.MediumMVA = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.TightMVA = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.LooseCutBased = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.MediumCutBased = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.TightCutBased = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
    cutByDiscriminator = cms.untracked.bool(True),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

process.PFChargedBased = cms.EDAnalyzer("phase2TausRECO",
    vertices           = cms.InputTag("offlinePrimaryVertices4D"),#check me
    taus               = cms.InputTag("hpsPFTauProducer"),
    jets               = cms.InputTag("ak4PFJetsCHS"),
    discriminator      = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    cutByDiscriminator = cms.untracked.bool(False),
    hpsPFTauChargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
    dmf                = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs"),
    genParticles       = cms.InputTag("genParticles")
)

###################################################
#Global sequence

process.p = cms.Path(
         process.VLooseMVA*
         process.LooseMVA*
         process.MediumMVA*
         process.TightMVA*
         process.LooseCutBased*
         process.MediumCutBased*
         process.TightCutBased*
         process.PFChargedBased
                     )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)

#print out all processes used when running- useful check to see if module ran
#UNCOMMENT BELOW
#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
