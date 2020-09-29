import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')
options.register('skipEvents', 
    default=0, 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
# TODO: put this option in cmsRun scripts
options.register('processMode', 
    default='JetLevel', 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process mode: JetLevel or EventLevel")
# Set doEBenergy to 1 to produce EGSeeds and EGFrames.
options.register('doEBenergy',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doEBenergy")
# Name of the EGInference model to be used for inference.
options.register('EGModelName',
    default='e_vs_ph_model.pb',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "EGInference Model name")
# Name of the QGInference model to be used for inference.
options.register('QGModelName',
    default='ResNet.pb',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "QGInference Model name")
# Name of the TopInference model to be used for inference.
options.register('TopQuarksModelName',
    default='ResNet.pb',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "TopQuarks Inference Model name")
options.parseArguments()

process = cms.Process("Classifier")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.GlobalTag.globaltag = cms.string('80X_dataRun2_HLT_v12')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents) 
    #input = cms.untracked.int32(1000) 
    #input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(1000000) 
    )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.inputFiles
      )#"file:/afs/cern.ch/user/s/schaudha/public/CMSSW_10_6_8/src/demo/ZprimeToTT_M-2000_W-20_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_AODSIM_PUMoriond17.root"
    , skipEvents = cms.untracked.uint32(0)#options.skipEvents
    )
print (" >> Loaded",len(options.inputFiles),"input files from list.")

#process.load("ProdTutorial.ProducerTest.DetImg_cfi")
process.load("E2eDL.FrameProducers.DetFrameProducer_cfi")
process.load("E2eDL.FrameProducers.EGFrameProducer_cfi")

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('myOutputFile.root')
#    ,outputCommands = cms.untracked.vstring('drop *',
#      "keep *_generalTracks_*_*",
#      "keep *_globalMuons_*_*",
#       "keep *_MuonTrackPoints_*_*",
#      "keep *_TrackTrackPoints_*_*")
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'))
#print " >> Processing as:",(process.fevt_tf.mode)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("myoutput.root")#options.outputFile
   )

process.p = cms.Path(process.DetFrames + process.EGFrames)
process.ep=cms.EndPath(process.out)
process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(False),
  useJobReport = cms.untracked.bool(True)
)
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    ignoreTotal = cms.untracked.int32(1)
)
