#include "E2eDL/FrameProducers/interface/JetFrameProducer.h"

// Run event selection ////////////////////////////////

extern unsigned long long eventId_;
extern unsigned int runId_;
extern unsigned int lumiId_;
extern float m0_;
//float nJet_;
extern float diPhoE_;
extern float diPhoPt_;
extern std::vector<float> vFC_inputs_;

extern float m0cut; //= 90.;
//float m0cut = 80.;

// Run event selection _______________________________________________________________//
bool TopProducer::runEvtSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::PhotonCollection> photons;
  //edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken( photonCollectionT_, photons );

  int nPhoTrg = 0;

  // Perform photon pre-selection
  float m0; //dR;
  float leadPhoPt = 0.;
  math::PtEtaPhiELorentzVectorD vDiPho;
  std::vector<int> vPhoIdxs;
  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    reco::PhotonRef iPho( photons, iP );
    //pat::PhotonRef iPho( photons, iP );
    if ( std::abs(iPho->pt()) < 18. ) continue;
    //std::cout << iPho->full5x5_sigmaIetaIeta() << std::endl;
    if ( std::abs(iPho->eta()) > 1.44 ) continue;
    if ( iPho->r9() < 0.5 ) continue;
    if ( iPho->hadTowOverEm() > 0.07 ) continue;
    if ( iPho->full5x5_sigmaIetaIeta() > 0.0105 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;
    //if ( std::abs(iPho->eta()) > 2.1 ) continue;
    //if ( std::abs(iPho->eta()) > 1.44 && std::abs(iPho->eta()) < 1.57 ) continue;
    if (debug) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
    nPhoTrg++;
    if ( std::abs(iPho->pt()) > leadPhoPt ) leadPhoPt = std::abs(iPho->pt()); 
    vDiPho += iPho->p4();
    vPhoIdxs.push_back( iP );

  } // reco photons 
  m0 = vDiPho.mass();
  if ( m0 < m0cut ) return false;
  if ( nPhoTrg != 2 ) return false;
  if ( leadPhoPt < 30. ) return false;

  // Apply selection
  int nPho = 0;
  int leadPho = -1;
  leadPhoPt = 0.;
  for ( int iP = 0; iP < nPhoTrg; iP++ ) {
    reco::PhotonRef iPho( photons, vPhoIdxs[iP] );
    //pat::PhotonRef iPho( photons, vPhoIdxs[iP] );
    // Get leading photon pt
    if ( std::abs(iPho->pt()) > leadPhoPt ) {
      leadPhoPt = std::abs(iPho->pt()); 
      leadPho = iP;
    }
    // Minimum pt/m0 cut
    if ( std::abs(iPho->pt()) < m0/4. ) continue;
    nPho++;
  }
  if ( nPho != 2 ) return false;
  if ( leadPhoPt < m0/3 ) return false;
  nPho = nPhoTrg;
  
  // Get photon order
  int ptOrder[2] = {0, 1};
  if ( leadPho == 1 ) {
    ptOrder[0] = 1;
    ptOrder[1] = 0;
  }
  
  // Fill kinematic variables
  //h_nJet->Fill( nJet );
  //h_m0->Fill( m0 );
  diPhoE_  = 0.;
  diPhoPt_ = 0.;
  float dphi[2] = {0., 0.};
  vFC_inputs_.clear();
  for ( int iP = 0; iP < nPho; iP++ ) {
    reco::PhotonRef iPho( photons, vPhoIdxs[ptOrder[iP]] );
    diPhoE_  += std::abs( iPho->energy() );
    diPhoPt_ += std::abs( iPho->pt() );
    vFC_inputs_.push_back( iPho->pt()/m0 );
    vFC_inputs_.push_back( iPho->eta() );
    dphi[iP] = iPho->phi();
  }
  vFC_inputs_.push_back( TMath::Cos(reco::deltaPhi(dphi[0], dphi[1])) );
  // Write out event
  m0_ = m0;
  //nJet_ = nJet;
  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();

  return true;

} // runEvtSel()
