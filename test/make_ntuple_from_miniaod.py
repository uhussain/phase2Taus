import FWCore.ParameterSet.Config as cms

'''

Configuration to produce PAT ntuples from miniaod for tau upgrade validation purposes.

Author: Martin Flechl, Vienna

Created: Sep 11, 2016


'''

process = cms.Process("nt")

local=2 #0 for crab

pu='pu200'
#geo='flat'
geo='2023D1'
rel='810_pre12'
ver='v6'
proc='ztt'

nevents=-1
fstr='%s_%s_%s_%s_%s.root' % ( proc,pu,geo,rel,ver )

infile='file:ZTT-RelVal-810-pre12.root'

outfile='outFile.root'
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring( infile )
    )


print 'Output file: ',outfile


#infile='file:../make_pat_new/ptuple_ztt_pu140_tilted_v6.root'

process.out = cms.OutputModule(
    "PoolOutputModule", 
    fileName = cms.untracked.string( outfile ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ), # save PAT Layer 1 output; you need a '*' to # unpack the list of commands 'patEventContent' 
    outputCommands = cms.untracked.vstring('drop *', 'keep *_*_*_*nt')
    )
process.out.outputCommands.append('drop *_prunedGenParticles2_*_*')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents) )

process.MessageLogger = cms.Service("MessageLogger")
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.truemu = cms.EDProducer(
    "PileupSummaryInfoSlimmer",
    src=cms.InputTag('addPileupInfo'),
#    keepDetailedInfoFor=cms.vint32( 1 , 2 )
    keepDetailedInfoFor=cms.vint32( )
    )


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.allprimvtx = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(0.0), maxZ = cms.double(1000000.0) ),
    src=cms.InputTag('offlineSlimmedPrimaryVertices')
    )

process.jet = cms.EDProducer( 
    "CandViewNtpProducer", 
#    src = cms.InputTag("selectedPatJets"), 
    src = cms.InputTag("slimmedJets:"), 
    lazyParser = cms.untracked.bool(True), 
    prefix = cms.untracked.string(""), 
    eventInfo = cms.untracked.bool(False), 
    variables = cms.VPSet( cms.PSet( tag = cms.untracked.string("pt"), 
                                     quantity = cms.untracked.string("pt") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("eta"), 
                                     quantity = cms.untracked.string("eta") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("phi"), 
                                     quantity = cms.untracked.string("phi") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("mass"), 
                                     quantity = cms.untracked.string("mass") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("partonFlavour"), 
                                     quantity = cms.untracked.string("partonFlavour") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("chargedHadronMultiplicity"), 
                                     quantity = cms.untracked.string("chargedHadronMultiplicity") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("neutralHadronMultiplicity"), 
                                     quantity = cms.untracked.string("neutralHadronMultiplicity") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isPFJet"), 
                                     quantity = cms.untracked.string("isPFJet") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("hadronFlavour"), 
                                     quantity = cms.untracked.string("hadronFlavour") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("photonMultiplicity"), 
                                     quantity = cms.untracked.string("photonMultiplicity") 
                                     ), 
                           ) 
    )

process.mu = cms.EDProducer( 
    "CandViewNtpProducer", 
#    src = cms.InputTag("selectedPatMuons"), 
    src = cms.InputTag("slimmedMuons"), 
    lazyParser = cms.untracked.bool(True), 
    prefix = cms.untracked.string(""), 
    eventInfo = cms.untracked.bool(False), 
    variables = cms.VPSet( cms.PSet( tag = cms.untracked.string("pt"), 
                                     quantity = cms.untracked.string("pt") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("eta"), 
                                     quantity = cms.untracked.string("eta") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("phi"), 
                                     quantity = cms.untracked.string("phi") 
                                     ), 
                           ) 
    )

process.prunedGenParticles2 = cms.EDProducer(
    "GenParticlePruner",
     src = cms.InputTag("prunedGenParticles"),
     select = cms.vstring(
         "drop  *",
         #keeps all particles from the hard matrix element
         #         "keep status = 3"
         #keeps all stable muons and electrons and their (direct) mothers. #keep+ ?
         "keep (abs(pdgId) = 11 | abs(pdgId) = 13) & status = 1",
         #keeps all taus after PS.
         "keep+ abs(pdgId) = 15",
         "keep abs(pdgId) = 23",
         "keep abs(pdgId) = 24",
         "keep abs(pdgId) = 25",
         "keep abs(pdgId) = 5",
         "drop abs(pdgId) = 15 & status = 3",
         #keep neutrinos.
         "keep (abs(pdgId) = 12 | abs(pdgId) = 14 | abs(pdgId) = 16 ) & status = 1"
         )
    )


process.gen = cms.EDProducer( 
    "CandViewNtpProducer", 
    src = cms.InputTag("prunedGenParticles2"), 
    lazyParser = cms.untracked.bool(True), 
    prefix = cms.untracked.string(""), 
    eventInfo = cms.untracked.bool(False), 
    variables = cms.VPSet( cms.PSet( tag = cms.untracked.string("pt"), 
                                     quantity = cms.untracked.string("pt") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("eta"), 
                                     quantity = cms.untracked.string("eta") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("phi"), 
                                     quantity = cms.untracked.string("phi") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("m"), 
                                     quantity = cms.untracked.string("mass") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("pdgId"), 
                                     quantity = cms.untracked.string("pdgId") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("status"), 
                                     quantity = cms.untracked.string("status") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("nDau"), 
                                     quantity = cms.untracked.string("numberOfDaughters") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("nMom"), 
                                     quantity = cms.untracked.string("numberOfMothers") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isPromptDecayed"), 
                                     quantity = cms.untracked.string("isPromptDecayed") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isPromptFinalState"), 
                                     quantity = cms.untracked.string("isPromptFinalState") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isLastCopy"), 
                                     quantity = cms.untracked.string("isLastCopy") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isFromHardProcessDecayed"), 
                                     quantity = cms.untracked.string("fromHardProcessDecayed") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isFromHardProcessFinalState"), 
                                     quantity = cms.untracked.string("fromHardProcessFinalState") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isDirectHardProcessTauDecayProductFinalState"), 
                                     quantity = cms.untracked.string("isDirectHardProcessTauDecayProductFinalState") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("isDirectPromptTauDecayProductFinalState"), 
                                     quantity = cms.untracked.string("isDirectPromptTauDecayProductFinalState") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("charge"), 
                                     quantity = cms.untracked.string("charge")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("vx"),
                                     quantity = cms.untracked.string("vx")
                                     ),
                           cms.PSet( tag = cms.untracked.string("vy"),
                                     quantity = cms.untracked.string("vy")
                                     ),
                           cms.PSet( tag = cms.untracked.string("vz"),
                                     quantity = cms.untracked.string("vz")
                                     ),

                           cms.PSet( tag = cms.untracked.string("mother0pdgId"), 
                                     quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).pdgId : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("mother0status"), 
                                     quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).status : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("mother0pt"), 
                                     quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("mother0eta"), 
                                     quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).eta : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("mother0phi"), 
                                     quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).phi : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("mother0m"), 
                                     quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).mass : -900")
                                     ), 

                          cms.PSet( tag = cms.untracked.string("daughter0pdgId"), 
                                     quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).pdgId : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("daughter0pt"), 
                                     quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter0eta"), 
                                     quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).eta : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("daughter0phi"), 
                                     quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).phi : -900")
                                     ), 
                           cms.PSet( tag = cms.untracked.string("daughter0m"), 
                                     quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).mass : -900")
                                     ), 

                           cms.PSet( tag = cms.untracked.string("daughter1pdgId"),
                                     quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).pdgId : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter1pt"),
                                     quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter1eta"),
                                     quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).eta : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter1phi"),
                                     quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).phi : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter1m"),
                                     quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).mass : -900")
                                     ),

                           cms.PSet( tag = cms.untracked.string("daughter2pdgId"),
                                     quantity = cms.untracked.string("?numberOfDaughters>2 ? daughter(2).pdgId : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter2pt"),
                                     quantity = cms.untracked.string("?numberOfDaughters>2 ? daughter(2).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter2eta"),
                                     quantity = cms.untracked.string("?numberOfDaughters>2 ? daughter(2).eta : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter2phi"),
                                     quantity = cms.untracked.string("?numberOfDaughters>2 ? daughter(2).phi : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter2m"),
                                     quantity = cms.untracked.string("?numberOfDaughters>2 ? daughter(2).mass : -900")
                                     ),

                           cms.PSet( tag = cms.untracked.string("daughter3pdgId"),
                                     quantity = cms.untracked.string("?numberOfDaughters>3 ? daughter(3).pdgId : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter3pt"),
                                     quantity = cms.untracked.string("?numberOfDaughters>3 ? daughter(3).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter3eta"),
                                     quantity = cms.untracked.string("?numberOfDaughters>3 ? daughter(3).eta : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter3phi"),
                                     quantity = cms.untracked.string("?numberOfDaughters>3 ? daughter(3).phi : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter3m"),
                                     quantity = cms.untracked.string("?numberOfDaughters>3 ? daughter(3).mass : -900")
                                     ),


                           cms.PSet( tag = cms.untracked.string("daughter4pdgId"),
                                     quantity = cms.untracked.string("?numberOfDaughters>4 ? daughter(4).pdgId : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter4pt"),
                                     quantity = cms.untracked.string("?numberOfDaughters>4 ? daughter(4).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter4eta"),
                                     quantity = cms.untracked.string("?numberOfDaughters>4 ? daughter(4).eta : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter4phi"),
                                     quantity = cms.untracked.string("?numberOfDaughters>4 ? daughter(4).phi : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter4m"),
                                     quantity = cms.untracked.string("?numberOfDaughters>4 ? daughter(4).mass : -900")
                                     ),


                           cms.PSet( tag = cms.untracked.string("daughter5pdgId"),
                                     quantity = cms.untracked.string("?numberOfDaughters>5 ? daughter(5).pdgId : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter5pt"),
                                     quantity = cms.untracked.string("?numberOfDaughters>5 ? daughter(5).pt : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter5eta"),
                                     quantity = cms.untracked.string("?numberOfDaughters>5 ? daughter(5).eta : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter5phi"),
                                     quantity = cms.untracked.string("?numberOfDaughters>5 ? daughter(5).phi : -900")
                                     ),
                           cms.PSet( tag = cms.untracked.string("daughter5m"),
                                     quantity = cms.untracked.string("?numberOfDaughters>5 ? daughter(5).mass : -900")
                                     ),

                           ) 
    )


process.el = cms.EDProducer( 
    "CandViewNtpProducer", 
    src = cms.InputTag("slimmedElectrons"), 
    lazyParser = cms.untracked.bool(True), 
    prefix = cms.untracked.string(""), 
    eventInfo = cms.untracked.bool(False), 
    variables = cms.VPSet( cms.PSet( tag = cms.untracked.string("pt"), 
                                     quantity = cms.untracked.string("pt") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("eta"), 
                                     quantity = cms.untracked.string("eta") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("phi"), 
                                     quantity = cms.untracked.string("phi") 
                                     ), 
                           ) 
    )

process.tau = cms.EDProducer( 
    "CandViewNtpProducer", 
    src = cms.InputTag("slimmedTaus"), 
    lazyParser = cms.untracked.bool(True), 
    prefix = cms.untracked.string(""), 
    eventInfo = cms.untracked.bool(False), 
    variables = cms.VPSet( cms.PSet( tag = cms.untracked.string("pt"), 
                                     quantity = cms.untracked.string("pt") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("eta"), 
                                     quantity = cms.untracked.string("eta") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("phi"), 
                                     quantity = cms.untracked.string("phi") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("genSize"), 
                                     quantity = cms.untracked.string("genParticlesSize") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("pdgId"), 
                                     quantity = cms.untracked.string("pdgId") 
                                     ), 
                           cms.PSet( tag = cms.untracked.string("dxy"),
                                     quantity = cms.untracked.string("dxy")
                                     ),
                           cms.PSet( tag = cms.untracked.string("vx"),
                                     quantity = cms.untracked.string("vx")
                                     ),
                           cms.PSet( tag = cms.untracked.string("vy"),
                                     quantity = cms.untracked.string("vy")
                                     ),
                           cms.PSet( tag = cms.untracked.string("vz"),
                                     quantity = cms.untracked.string("vz")
                                     ),
                           cms.PSet( tag = cms.untracked.string("decayMode"), 
                                     quantity = cms.untracked.string("decayMode")
                                     ),
cms.PSet( tag = cms.untracked.string("againstElectronLooseMVA6"),quantity = cms.untracked.string("tauID('againstElectronLooseMVA6')") ),
cms.PSet( tag = cms.untracked.string("againstElectronMVA6Raw"),quantity = cms.untracked.string("tauID('againstElectronMVA6Raw')") ),
cms.PSet( tag = cms.untracked.string("againstElectronMVA6category"),quantity = cms.untracked.string("tauID('againstElectronMVA6category')") ),
cms.PSet( tag = cms.untracked.string("againstElectronMediumMVA6"),quantity = cms.untracked.string("tauID('againstElectronMediumMVA6')") ),
cms.PSet( tag = cms.untracked.string("againstElectronTightMVA6"),quantity = cms.untracked.string("tauID('againstElectronTightMVA6')") ),
cms.PSet( tag = cms.untracked.string("againstElectronVLooseMVA6"),quantity = cms.untracked.string("tauID('againstElectronVLooseMVA6')") ),
cms.PSet( tag = cms.untracked.string("againstElectronVTightMVA6"),quantity = cms.untracked.string("tauID('againstElectronVTightMVA6')") ),
cms.PSet( tag = cms.untracked.string("againstMuonLoose3"),quantity = cms.untracked.string("tauID('againstMuonLoose3')") ),
cms.PSet( tag = cms.untracked.string("againstMuonTight3"),quantity = cms.untracked.string("tauID('againstMuonTight3')") ),
cms.PSet( tag = cms.untracked.string("byCombinedIsolationDeltaBetaCorrRaw3Hits"),quantity = cms.untracked.string("tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits')") ),
cms.PSet( tag = cms.untracked.string("byIsolationMVArun2v1DBdR03oldDMwLTraw"),quantity = cms.untracked.string("tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw')") ),
cms.PSet( tag = cms.untracked.string("byIsolationMVArun2v1DBnewDMwLTraw"),quantity = cms.untracked.string("tauID('byIsolationMVArun2v1DBnewDMwLTraw')") ),
cms.PSet( tag = cms.untracked.string("byIsolationMVArun2v1DBoldDMwLTraw"),quantity = cms.untracked.string("tauID('byIsolationMVArun2v1DBoldDMwLTraw')") ),
cms.PSet( tag = cms.untracked.string("byIsolationMVArun2v1PWdR03oldDMwLTraw"),quantity = cms.untracked.string("tauID('byIsolationMVArun2v1PWdR03oldDMwLTraw')") ),
cms.PSet( tag = cms.untracked.string("byIsolationMVArun2v1PWnewDMwLTraw"),quantity = cms.untracked.string("tauID('byIsolationMVArun2v1PWnewDMwLTraw')") ),
cms.PSet( tag = cms.untracked.string("byIsolationMVArun2v1PWoldDMwLTraw"),quantity = cms.untracked.string("tauID('byIsolationMVArun2v1PWoldDMwLTraw')") ),
cms.PSet( tag = cms.untracked.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"),quantity = cms.untracked.string("tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits')") ),
cms.PSet( tag = cms.untracked.string("byLooseIsolationMVArun2v1DBdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byLooseIsolationMVArun2v1DBdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byLooseIsolationMVArun2v1DBnewDMwLT"),quantity = cms.untracked.string("tauID('byLooseIsolationMVArun2v1DBnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byLooseIsolationMVArun2v1DBoldDMwLT"),quantity = cms.untracked.string("tauID('byLooseIsolationMVArun2v1DBoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byLooseIsolationMVArun2v1PWdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byLooseIsolationMVArun2v1PWdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byLooseIsolationMVArun2v1PWnewDMwLT"),quantity = cms.untracked.string("tauID('byLooseIsolationMVArun2v1PWnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byLooseIsolationMVArun2v1PWoldDMwLT"),quantity = cms.untracked.string("tauID('byLooseIsolationMVArun2v1PWoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byMediumCombinedIsolationDeltaBetaCorr3Hits"),quantity = cms.untracked.string("tauID('byMediumCombinedIsolationDeltaBetaCorr3Hits')") ),
cms.PSet( tag = cms.untracked.string("byMediumIsolationMVArun2v1DBdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byMediumIsolationMVArun2v1DBdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byMediumIsolationMVArun2v1DBnewDMwLT"),quantity = cms.untracked.string("tauID('byMediumIsolationMVArun2v1DBnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byMediumIsolationMVArun2v1DBoldDMwLT"),quantity = cms.untracked.string("tauID('byMediumIsolationMVArun2v1DBoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byMediumIsolationMVArun2v1PWdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byMediumIsolationMVArun2v1PWdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byMediumIsolationMVArun2v1PWnewDMwLT"),quantity = cms.untracked.string("tauID('byMediumIsolationMVArun2v1PWnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byMediumIsolationMVArun2v1PWoldDMwLT"),quantity = cms.untracked.string("tauID('byMediumIsolationMVArun2v1PWoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byPhotonPtSumOutsideSignalCone"),quantity = cms.untracked.string("tauID('byPhotonPtSumOutsideSignalCone')") ),
cms.PSet( tag = cms.untracked.string("byTightCombinedIsolationDeltaBetaCorr3Hits"),quantity = cms.untracked.string("tauID('byTightCombinedIsolationDeltaBetaCorr3Hits')") ),
cms.PSet( tag = cms.untracked.string("byTightIsolationMVArun2v1DBdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byTightIsolationMVArun2v1DBdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byTightIsolationMVArun2v1DBnewDMwLT"),quantity = cms.untracked.string("tauID('byTightIsolationMVArun2v1DBnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byTightIsolationMVArun2v1DBoldDMwLT"),quantity = cms.untracked.string("tauID('byTightIsolationMVArun2v1DBoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byTightIsolationMVArun2v1PWdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byTightIsolationMVArun2v1PWdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byTightIsolationMVArun2v1PWnewDMwLT"),quantity = cms.untracked.string("tauID('byTightIsolationMVArun2v1PWnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byTightIsolationMVArun2v1PWoldDMwLT"),quantity = cms.untracked.string("tauID('byTightIsolationMVArun2v1PWoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byVLooseIsolationMVArun2v1DBdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVLooseIsolationMVArun2v1DBnewDMwLT"),quantity = cms.untracked.string("tauID('byVLooseIsolationMVArun2v1DBnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVLooseIsolationMVArun2v1DBoldDMwLT"),quantity = cms.untracked.string("tauID('byVLooseIsolationMVArun2v1DBoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVLooseIsolationMVArun2v1PWdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byVLooseIsolationMVArun2v1PWdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVLooseIsolationMVArun2v1PWnewDMwLT"),quantity = cms.untracked.string("tauID('byVLooseIsolationMVArun2v1PWnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVLooseIsolationMVArun2v1PWoldDMwLT"),quantity = cms.untracked.string("tauID('byVLooseIsolationMVArun2v1PWoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVTightIsolationMVArun2v1DBdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byVTightIsolationMVArun2v1DBdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVTightIsolationMVArun2v1DBnewDMwLT"),quantity = cms.untracked.string("tauID('byVTightIsolationMVArun2v1DBnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVTightIsolationMVArun2v1DBoldDMwLT"),quantity = cms.untracked.string("tauID('byVTightIsolationMVArun2v1DBoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVTightIsolationMVArun2v1PWdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byVTightIsolationMVArun2v1PWdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVTightIsolationMVArun2v1PWnewDMwLT"),quantity = cms.untracked.string("tauID('byVTightIsolationMVArun2v1PWnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVTightIsolationMVArun2v1PWoldDMwLT"),quantity = cms.untracked.string("tauID('byVTightIsolationMVArun2v1PWoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVVTightIsolationMVArun2v1DBdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byVVTightIsolationMVArun2v1DBdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVVTightIsolationMVArun2v1DBnewDMwLT"),quantity = cms.untracked.string("tauID('byVVTightIsolationMVArun2v1DBnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVVTightIsolationMVArun2v1DBoldDMwLT"),quantity = cms.untracked.string("tauID('byVVTightIsolationMVArun2v1DBoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVVTightIsolationMVArun2v1PWdR03oldDMwLT"),quantity = cms.untracked.string("tauID('byVVTightIsolationMVArun2v1PWdR03oldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVVTightIsolationMVArun2v1PWnewDMwLT"),quantity = cms.untracked.string("tauID('byVVTightIsolationMVArun2v1PWnewDMwLT')") ),
cms.PSet( tag = cms.untracked.string("byVVTightIsolationMVArun2v1PWoldDMwLT"),quantity = cms.untracked.string("tauID('byVVTightIsolationMVArun2v1PWoldDMwLT')") ),
cms.PSet( tag = cms.untracked.string("chargedIsoPtSum"),quantity = cms.untracked.string("tauID('chargedIsoPtSum')") ),
cms.PSet( tag = cms.untracked.string("chargedIsoPtSumdR03"),quantity = cms.untracked.string("tauID('chargedIsoPtSumdR03')") ),
cms.PSet( tag = cms.untracked.string("decayModeFinding"),quantity = cms.untracked.string("tauID('decayModeFinding')") ),
cms.PSet( tag = cms.untracked.string("decayModeFindingNewDMs"),quantity = cms.untracked.string("tauID('decayModeFindingNewDMs')") ),
cms.PSet( tag = cms.untracked.string("footprintCorrection"),quantity = cms.untracked.string("tauID('footprintCorrection')") ),
cms.PSet( tag = cms.untracked.string("footprintCorrectiondR03"),quantity = cms.untracked.string("tauID('footprintCorrectiondR03')") ),
cms.PSet( tag = cms.untracked.string("neutralIsoPtSum"),quantity = cms.untracked.string("tauID('neutralIsoPtSum')") ),
cms.PSet( tag = cms.untracked.string("neutralIsoPtSumWeight"),quantity = cms.untracked.string("tauID('neutralIsoPtSumWeight')") ),
cms.PSet( tag = cms.untracked.string("neutralIsoPtSumWeightdR03"),quantity = cms.untracked.string("tauID('neutralIsoPtSumWeightdR03')") ),
cms.PSet( tag = cms.untracked.string("neutralIsoPtSumdR03"),quantity = cms.untracked.string("tauID('neutralIsoPtSumdR03')") ),
cms.PSet( tag = cms.untracked.string("photonPtSumOutsideSignalCone"),quantity = cms.untracked.string("tauID('photonPtSumOutsideSignalCone')") ),
cms.PSet( tag = cms.untracked.string("photonPtSumOutsideSignalConedR03"),quantity = cms.untracked.string("tauID('photonPtSumOutsideSignalConedR03')") ),
cms.PSet( tag = cms.untracked.string("puCorrPtSum"),quantity = cms.untracked.string("tauID('puCorrPtSum')") ),
                           ) 
    )

process.p = cms.Path(process.allprimvtx*process.jet*process.mu*process.el*process.tau*process.prunedGenParticles2*process.gen)

process.outpath = cms.EndPath(process.out)
