import FWCore.ParameterSet.Config as cms

QGTagger = cms.EDProducer('QGTagger'
    , QGFrames = cms.InputTag('QGFrames','QGFrames')
    , reducedHBHERecHitCollection = cms.InputTag('reducedHcalRecHits:hbhereco')
    , ak4PFJetCollection = cms.InputTag('ak4PFJetsCHS')
    , ak4GenJetCollection = cms.InputTag('ak4GenJets')
    , ak4RecoJetsForBTagging = cms.InputTag("ak4PFJetsCHS")
    , QGModelName = cms.string('ResNet.pb')
    , doECALstitched = cms.bool(True)
    , doTracksAtECALstitchedPt = cms.bool(True)
    , doTracksAtECALadjPt = cms.bool(True)
    , doHBHEenergy = cms.bool(True)
    )
