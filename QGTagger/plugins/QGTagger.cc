#include "E2eDL/QGTagger/interface/QGTagger.h"

QGTagger::QGTagger(const edm::ParameterSet& iConfig)
{
  // Input tokens
  HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
  tQGframeCollection = consumes<std::vector<e2e::Frame3D> >(iConfig.getParameter<edm::InputTag>("QGFrames"));
  jetCollectionT_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak4PFJetCollection"));
  genJetCollectionT_      = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak4GenJetCollection"));
  recoJetsT_              = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("ak4RecoJetsForBTagging"));
  TracksAtECALadjPt_token = consumes<e2e::Frame1D>(iConfig.getParameter<edm::InputTag>("TracksAtECALadjPt"));
  JetFramesT_ = consumes<e2e::Frame4D>(iConfig.getParameter<edm::InputTag>("JetFrames"));
  //tEGframeCollection = consumes<e2e::PhoFrame3DCollection>(iConfig.getParameter<edm::InputTag>("EGFrames"));

  mode_      = iConfig.getParameter<std::string>("mode");
  minJetPt_  = iConfig.getParameter<double>("minJetPt");
  maxJetEta_ = iConfig.getParameter<double>("maxJetEta");
  z0PVCut_   = iConfig.getParameter<double>("z0PVCut");
  
  // DL inference model
  modelName = iConfig.getParameter<std::string>("QGModelName");

  // Output collections to be produced
  produces<e2e::PhoPredCollection>("QGProbs");
}

QGTagger::~QGTagger()
{
}

void
QGTagger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::LogInfo("QGTagger") << " >> Running QGTagger...";

  // Load required tokens into input collection handles
  iEvent.getByToken( tQGframeCollection, hQGframe );
  iEvent.getByToken( jetCollectionT_, jets );
  iEvent.getByToken( JetFramesT_, hJetFrames );
  assert( hJetFrames->size() == jets->size() );
  
  nJets = jets->size();
  std::vector<e2e::pred>    vJetProbs ( nJets, defaultVal );
  if (hJetFrames->size()>0) {
    // Get pointer to input EG frames
    const std::vector<e2e::Frame3D>* pJetFrame = hJetFrames.product();
    nFrameD = pJetFrame->front().size(); // get size of depth dimension

    // Initialize product values to be stored with default values at start of every event
    // Each object is a vector over the no. of photons in the event
  
    std::vector<e2e::Frame3D> vJetFrames( nJets,
                                        e2e::Frame3D(nFrameD,
                                        e2e::Frame2D(nFrameH,
                                        e2e::Frame1D(nFrameW, 0.))) );

    //_____ Load EG frame collection into `vPhoFrames` for each photon _____//

    for ( unsigned int iJ = 0; iJ < jets->size(); iJ++ ) {
      // Get EG frame for this photon
      vJetFrames[iP] = pJetFrame->at(iJ);
    } // photons

    //_____ Run DL inference _____//

    // Run inference on `vPhoFrames` batch of size nPhos*nFrameD*nFrameH*nFrameW: store output in `vPhoProbs`
    // Running on entire batch at once maximizes computing parellization
    // runInference( vPhoProbs, vPhoFrames, modelName );
  
    e2e::Frame2D tmp_out = e2e::predict_tf(vJetFrames, "ResNet.pb", "inputs","outputs");
  
    //_____ Store products associated with each photon _____//

    // Initialize pointers to edm::AssociationVector (key,val) collections
    // These collections create explicit associations between the photon object (key) and the stored product (val)
  }
  cJetProbs  = std::make_unique<e2e::Frame2D>   ( tmp_out );
  // Set association between photon ref (key) and products to be stored (val)
  /*for ( unsigned int iP = 0; iP < hPhoton->size(); iP++ ) {
    PhotonRef iRecoPho( hPhoton, iP );
    cPhoProbs->setValue( iP, vPhoProbs[iP] );
  } */ // photons
    
    
  // Put collections into output EDM file
  iEvent.put( std::move(cJetProbs), "JetProbs" );

  return;
} // EGTagger::produce()

void
QGTagger::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
QGTagger::endStream()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QGTagger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(QGTagger);
