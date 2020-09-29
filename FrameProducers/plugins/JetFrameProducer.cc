#include "E2eDL/FrameProducers/plugins/JetFrameProducer.h"

JetFrameProducer::JetFrameProducer(const edm::ParameterSet& iConfig)
{
  photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
  EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
  ECALstitched_energy_token=consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("ECALstitchedenergy"));
  TracksAtECALstitchedPt_token=consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("TracksAtECALstitchedPt"));
  HBHEenergy_token = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("HBHEenergy"));
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
  vJetSeed_ieta_.clear(); vJetSeed_iphi_.clear();
  bool passedSelection = false;
  // Selecting Jet Seeds (ak8 / ak4) and storing them in edm root file.
  if ( doJets_ ) {
    edm::LogInfo("JetFrameProducer") << " >> doJets set";
    passedSelection = runEventSel_jet( iEvent, iSetup );
    std::cout<<" >> Size of JetSeed vector (JetSeed_eta_size, JetSeed_phi_size) is: ("<<vJetSeed_ieta_.size()<<", "<<vJetSeed_iphi_.size()<<")"<<std::endl;
    std::cout<<" >> The jet seeds are (ieta,iphi): ";
    if (vJetSeed_ieta_.size()==0){vJetSeed_ieta_.push_back(-1); vJetSeed_iphi_.push_back(-1); std::cout<<"(-1, -1)"<<std::endl;}
    else{
     	for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     	  std::cout<<"("<<vJetSeed_ieta_[idx]<<","<<vJetSeed_iphi_[idx]<<") ";
     	}
     	std::cout<<std::endl;
     }
     if (vJetSeed_ieta_.size()==vJetSeed_iphi_.size()){
     	for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		if(vJetSeed_ieta_[idx]>=0){vJetSeed_ieta_[idx]=int(vJetSeed_ieta_[idx]*5+2);}  //5 EB xtals per HB tower
		if(vJetSeed_iphi_[idx]>=0){vJetSeed_iphi_[idx]=int(vJetSeed_iphi_[idx]*5+2);}  //5 EB xtals per HB tower
		//std::cout<<vJetSeed_ieta_[idx]<<" "<<vJetSeed_iphi_[idx];
     	}
     }
     std::unique_ptr<std::vector<int>> JetSeedieta_edm (new std::vector<int>(vJetSeed_ieta_));
     std::unique_ptr<std::vector<int>> JetSeediphi_edm (new std::vector<int>(vJetSeed_iphi_));
     iEvent.put(std::move(JetSeedieta_edm),"ak8JetSeedieta");
     iEvent.put(std::move(JetSeediphi_edm),"ak8JetSeediphi");
     //vJetSeed_ieta_.clear(); vJetSeed_iphi_.clear();
   } else {
     std::cout<<" >> doJets not set"<<std::endl;
     passedSelection = runEvtSel( iEvent, iSetup );
     std::cout<<" >> Size of JetSeed vector (ak8JetSeed_eta_size, ak8JetSeed_phi_size) is: ("<<vJetSeed_ieta_.size()<<", "<<vJetSeed_iphi_.size()<<")"<<std::endl;
     std::cout<<" The jet seeds are (ieta,iphi): ";
     if (vJetSeed_ieta_.size()==0){vJetSeed_ieta_.push_back(-1); vJetSeed_iphi_.push_back(-1); std::cout<<"(-1, -1)"<<std::endl;}
     else{
	   for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		std::cout<<"("<<vJetSeed_ieta_[idx]<<","<<vJetSeed_iphi_[idx]<<") ";
     	}
     	std::cout<<std::endl;
     }
     if (vJetSeed_ieta_.size()==vJetSeed_iphi_.size()){
     	for (int idx=0;idx<int(vJetSeed_ieta_.size());idx++){
     		if(vJetSeed_ieta_[idx]>=0){vJetSeed_ieta_[idx]=int(vJetSeed_ieta_[idx]*5+2);}  //5 EB xtals per HB tower
		if(vJetSeed_iphi_[idx]>=0){vJetSeed_iphi_[idx]=int(vJetSeed_iphi_[idx]*5+2);}  //5 EB xtals per HB tower
		//std::cout<<vJetSeed_ieta_[idx]<<" "<<vJetSeed_iphi_[idx];
     	}
     }
     std::unique_ptr<std::vector<int>> JetSeedieta_edm (new std::vector<int>(vJetSeed_ieta_));
     std::unique_ptr<std::vector<int>> JetSeediphi_edm (new std::vector<int>(vJetSeed_iphi_));
     iEvent.put(std::move(JetSeedieta_edm),"ak8JetSeedieta");
     iEvent.put(std::move(JetSeediphi_edm),"ak8JetSeediphi");
   }

   if ( !passedSelection ) {
     h_sel->Fill( 0. );;
     return;
   }  
  // Load required tokens into input collection handles
  
}
