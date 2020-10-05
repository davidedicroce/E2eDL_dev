import FWCore.ParameterSet.Config as cms

EGTagger = cms.EDProducer('EGTagger'
    , photonCollection = cms.InputTag('gedPhotons')
    , EGFrames = cms.InputTag('EGFrameProducer','EGFrames')
    )
