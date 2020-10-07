#ifndef RecoE2E_QGTagger_h
#define RecoE2E_QGTagger_h

#include <memory>
#include <iostream>
#include <vector>
#include <cassert>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here
//#include "DataFormats/PatCandidates/interface/Photon.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "E2eDL/DataFormats/interface/FrameCollections.h"
#include "E2eDL/FrameProducers/interface/JetFrameProducer.h"
#include "E2eDL/FrameProducers/interface/predict_tf.h"

using namespace std;

//using pat::PhotonCollection;
//using pat::PhotonRef;

class QGTagger : public edm::stream::EDProducer<> {

   public:

      explicit QGTagger(const edm::ParameterSet&);
      ~QGTagger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:

      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      // Input tokens
      edm::EDGetTokenT<reco::PFJetCollection> jetCollectionT_;
      edm::EDGetTokenT<e2e::Frame4D> JetFramesT_;
      //edm::EDGetTokenT<e2e::PhoFrame3DCollection> tEGframeCollection;
      // Handles
      edm::Handle<reco::PFJetCollection> jets;
      edm::Handle<e2e::Frame4D> hJetFrames;
      //edm::Handle<e2e::PhoFrame3DCollection> hEGframe;
   
      // DL inference model
      std::string modelName;
      void runInference( std::vector<e2e::pred>&, const std::vector<e2e::Frame3D>&, const std::string );

      // Vector to hold input EG frames for inference
      std::vector<e2e::Frame3D> vJetFrames;

      // Frame dimensions determined at runtime
      int nJets;   // frame batch size in no. of photons
      int nFrameD; // frame depth in no. of detector layers

      // Output collections to be produced and values stored in them
      std::unique_ptr<e2e::Frame2D> cPhoProbs;
      std::vector<e2e::pred> vJetProbs;

}; // EGTagger

#endif
