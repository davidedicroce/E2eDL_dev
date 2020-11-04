import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# TODO: put this option in cmsRun scripts
options.register('processMode', 
    default='JetLevel', 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process mode: JetLevel or EventLevel")
# Skip Events.
options.register('skipEvents',
    default=0,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
# Set doECALstitched to 1 to produce JetSeeds and JetFrames.
options.register('doECALstitched',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doECALstitched")
# Set doTracksAtECALstitchedPt to 1 to produce JetSeeds and JetFrames.
options.register('doTracksAtECALstitchedPt',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doTracksAtECALstitchedPt")
# Set doTracksAtECALadjPt to 1 to produce JetSeeds and JetFrames.
options.register('doTracksAtECALadjPt',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doTracksAtECALadjPt")
# Set doHBHEenergy to 1 to produce JetSeeds and JetFrames.
options.register('doHBHEenergy',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doHBHEenergy")
# Name of the TopInference model to be used for inference.
options.register('TopModelName',
    default='ResNet.pb',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "TopInference Model name")
options.parseArguments()

process = cms.Process("TopClassifier")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.connect = cms.string('frontier://(backupproxyurl=http://cmsbpfrontier.cern.ch:3128)(backupproxyurl=http://cmsbpfrontier1.cern.ch:3128)(backupproxyurl=http://cmsbpfrontier2.cern.ch:3128)(backupproxyurl=http://cmsbproxy.fnal.gov:3128)(serverurl=http://cmsfrontier.cern.ch:8000/FrontierProd)(serverurl=http://cmsfrontier1.cern.ch:8000/FrontierProd)(serverurl=http://cmsfrontier2.cern.ch:8000/FrontierProd)(serverurl=http://cmsfrontier3.cern.ch:8000/FrontierProd)(serverurl=http://cmsfrontier4.cern.ch:8000/FrontierProd)/CMS_CONDITIONS')
#process.GlobalTag.connect   = cms.string('frontier://(proxyurl=http://localhost:8888)(serverurl=http://localhost:10800/FrontierOnProd)/CMS_CONDITIONS')
#process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://(proxyurl=http://localhost:8888)(serverurl=http://localhost:10800/FrontierOnProd)')
process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v15')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    #input = cms.untracked.int32(10)
    )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.inputFiles
      #"file:myOutputFile.root"#SinglePhotonPt50_noPU_AODSIM.root
      )
    , skipEvents = cms.untracked.uint32(0)#options.skipEvents
    )
print (" >> Loaded",len(options.inputFiles),"input files from list.")

process.load("E2eDL.FrameProducers.DetFrameProducer_cfi")
process.load("E2eDL.FrameProducers.JetFrameProducer_cfi")
process.load("E2eDL.TopTagger.TopTagger_cfi")
#process.EGTagger.EGModelName = options.EGModelName
process.JetFrames.jetCollection = cms.string("ak8")
process.JetFrames.minJetPt = cms.double(35.)
process.JetFrames.maxJetEta = cms.double(2.4)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('TopPt+TopFrames.root') 
    )
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root")#options.outputFile
    )

process.p = cms.Path(process.DetFrames + process.JetFrames+process.TopTagger)
process.ep=cms.EndPath(process.out)

#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
