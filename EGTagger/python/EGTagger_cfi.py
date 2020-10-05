import FWCore.ParameterSet.Config as cms

EGTagger = cms.EDProducer('EGTagger'
    , photonCollection = cms.InputTag('gedPhotons')
    , EGframes = cms.InputTag('EGFrameProducer','EGframes')
    )
