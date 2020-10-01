#include "E2eDL/FrameProducers/plugins/JetFrameProducer.h"

JetFrameProducer::JetFrameProducer(const edm::ParameterSet& iConfig)
{
  photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
  EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
  ECALstitched_energy_token=consumes<e2e::Frame1D>(iConfig.getParameter<edm::InputTag>("ECALstitchedenergy"));
  TracksAtECALstitchedPt_token=consumes<e2e::Frame1D>(iConfig.getParameter<edm::InputTag>("TracksAtECALstitchedPt"));
  HBHEenergy_token = consumes<e2e::Frame1D>>(iConfig.getParameter<edm::InputTag>("HBHEenergy"));
  vertexCollectionT_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  
  TRKRecHitCollectionT_   = consumes<TrackingRecHitCollection>(iConfig.getParameter<edm::InputTag>("trackRecHitCollection"));
  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  
  pfCollectionT_          = consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCollection"));
  
  jetTagCollectionT_      = consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("jetTagCollection"));
  ipTagInfoCollectionT_   = consumes<std::vector<reco::CandIPTagInfo> > (iConfig.getParameter<edm::InputTag>("ipTagInfoCollection"));
  
  // Jet Collection switches
  jetCollection_sel = iConfig.getParameter<std::string>("jetCollection");//
  if (jetCollection_sel == "ak4"){
    jetCollectionT_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak4PFJetCollection"));
    genJetCollectionT_      = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak4GenJetCollection"));
    recoJetsT_              = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("ak4RecoJetsForBTagging"));
  }
  else if (jetCollection_sel == "ak8"){
    jetCollectionT_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak8PFJetCollection"));
    genJetCollectionT_      = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak8GenJetCollection"));
    recoJetsT_              = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("ak8RecoJetsForBTagging"));
    TracksAtECALadjPt_token=consumes<e2e::Frame1D>(iConfig.getParameter<edm::InputTag>("TracksAtECALadjPt"));
  }
  
  mode_      = iConfig.getParameter<std::string>("mode");
  minJetPt_  = iConfig.getParameter<double>("minJetPt");
  maxJetEta_ = iConfig.getParameter<double>("maxJetEta");
  z0PVCut_   = iConfig.getParameter<double>("z0PVCut");
  modelName = iConfig.getParameter<std::string>("QGModelName");
  std::cout << " >> Mode set to " << mode_ << std::endl;	
  if ( mode_ == "JetLevel" ) {
     doJets_ = true;
     nJets_ = iConfig.getParameter<int>("nJets");
     std::cout << "\t>> nJets set to " << nJets_ << std::endl;
  } else if ( mode_ == "EventLevel" ) {
     doJets_ = false;
  } else {
     std::cout << " >> Assuming EventLevel Config. " << std::endl;
     doJets_ = false;
  }
  
  // Output collections to be produced
  produces<e2e::PhoSeedCollection>   ("QGTracksAtECALstitchedJetCollectionPtSeeds");
  produces<e2e::PhoFrame3DCollection>("QGtracksAtECALstitchedJetCollectionPtFrames");
  produces<e2e::PhoSeedCollection>   ("QGECALStitchedJetCollectionSeeds");
  produces<e2e::PhoFrame3DCollection>("QGECALStitchedJetCollectionFrames");
  produces<e2e::PhoSeedCollection>   ("QGHBHEJetCollectionSeeds");
  produces<e2e::PhoFrame3DCollection>("QGHBHEJetCollectionFrames");
}

JetFrameProducer::~JetFrameProducer()
{
}

void
JetFrameProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::LogInfo("JetFrameProducer") << " >> Running JetFrameProducer...";
  
  nTotal++;
  std::vector<e2e::pred>    vPhProbs ( nPhos, defaultVal );
  std::vector<e2e::seed>    vJetSeeds ( nPhos, e2e::seed(nSeedCoords, int(defaultVal)) );
  std::vector<e2e::Frame3D> vPhoFrames( nPhos,
                                        e2e::Frame3D(nFrameD,
                                        e2e::Frame2D(nFrameH,
                                        e2e::Frame1D(nFrameW, 0.))) );
  
  bool passedSelection = false;
  // Selecting Jet Seeds (ak8 / ak4) and storing them in edm root file.
  if ( doJets_ ) {
    edm::LogInfo("JetFrameProducer") << " >> doJets set";
    passedSelection = runEventSel_jet( iEvent, iSetup );
    std::cout<<" >> Number of Jets: "<<vJetSeeds.size()<<std::endl;
    std::cout<<" >> The jet seeds are (ieta,iphi): ";
    for (int idx=0;idx<int(vJetSeeds.size());idx++){
    	std::cout<<"("<<vJetSeeds[idx][0]<<","<<vJetSeeds[idx][1]<<") ";
     }
     std::cout<<std::endl;
     for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     	if(vJetSeeds[idx][0]>=0){vJetSeeds[idx][0]=int(vJetSeeds[idx][0]*5+2);}  //5 EB xtals per HB tower
	if(vJetSeeds[idx][1]>=0){vJetSeeds[idx][1]=int(vJetSeeds[idx][1]*5+2);}  //5 EB xtals per HB tower
     }
     std::unique_ptr<e2e::seed> JetSeeds_edm (new e2e::seed(vJetSeeds));
     if (jetCollection_sel == "ak4"){
     	iEvent.put(std::move(JetSeeds_edm),"ak8JetSeeds");
     }
     else if (jetCollection_sel == "ak8"){
	iEvent.put(std::move(JetSeeds_edm, "ak4JetSeeds");
     }
  } 
  else {
     edm::LogInfo("JetFrameProducer") << " >> doJets not set";
     passedSelection = runEvtSel( iEvent, iSetup );
     std::cout<<" >> Number of Jets: "<<vJetSeeds.size()<<std::endl;
     std::cout<<" The jet seeds are (ieta,iphi): ";
     for (int idx=0;idx<int(vJetSeeds.size());idx++){
     	std::cout<<"("<<vJetSeeds[idx][0]<<","<<vJetSeeds[idx][1]<<") ";
     }
     std::cout<<std::endl;
     for (int idx=0;idx<int(vJetSeeds.size());idx++){
     	if(vJetSeeds[idx][0]>=0){vJetSeeds[idx][0]=int(vJetSeeds[idx][0]*5+2);}  //5 EB xtals per HB tower
	if(vJetSeeds[idx][1]>=0){vJetSeeds[idx][1]=int(vJetSeeds[idx][1]*5+2);}  //5 EB xtals per HB tower
     }
     std::unique_ptr<e2e::seed> JetSeeds_edm (new e2e::seed(vJetSeed_ieta_));
     if (jetCollection_sel == "ak4"){
     	iEvent.put(std::move(JetSeeds_edm),"ak4JetSeeds");
     }
     else if (jetCollection_sel == "ak8"){
	iEvent.put(std::move(JetSeeds_edm, "ak8JetSeeds");
     }
   }

   edm::Handle<e2e::Frame1D> ECALstitched_energy_handle;
   iEvent.getByToken(ECALstitched_energy_token, ECALstitched_energy_handle);
   edm::Handle<e2e::Frame1D> TracksAtECALstitchedPt_handle;
   iEvent.getByToken(TracksAtECALstitchedPt_token, TracksAtECALstitchedPt_handle);
   if (jetCollection_sel == "ak8"){
   edm::Handle<e2e::Frame1D> TracksAtECALadjPt_handle;
   	iEvent.getByToken(TracksAtECALadjPt_token, TracksAtECALadjPt_handle);
   }
   edm::Handle<std::vector<float>> HBHEenergy_handle;
   iEvent.getByToken(HBHEenergy_token, HBHEenergy_handle);
  // Load required tokens into input collection handles
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetFrameProducer);
