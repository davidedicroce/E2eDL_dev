import FWCore.ParameterSet.Config as cms

EGTagger = cms.EDProducer('EGTagger'
    , photonCollection = cms.InputTag('gedPhotons')
    , EGFrames = cms.InputTag('EGFrames','EGFrames')
    , EGModelName = cms.string('e_vs_ph_model.pb')
    )
