// -*- C++ -*-
//
// Package:    RecoTauTag/phase2Taus
// Class:      phase2Taus
// 
/**\class phase2Taus phase2Taus.cc RecoTauTag/phase2Taus/plugins/phase2Taus.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isabel Ojalvo
//         Created:  Tue, 15 Nov 2016 16:00:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "iostream"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class phase2Taus : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2Taus(const edm::ParameterSet&);
      ~phase2Taus();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      std::string tauID_;
  edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;

  TTree* tree;
  double tauPt_;
  double tauEta_;
  double tauMass_;
  double genTauPt_;
  double genTauEta_;

  int nvtx_;
  int dmf_;
  int goodReco_;
  int genTauMatch_;

  reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
  void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
  bool isNeutrino(const reco::Candidate* daughter);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
phase2Taus::phase2Taus(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  prunedGenToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   tauID_    = iConfig.getParameter<std::string>("tauID");
   edm::Service<TFileService> fs;

   tree = fs->make<TTree>("Ntuple", "Ntuple");
   tree->Branch("tauPt", &tauPt_,"tauPt_/D");
   tree->Branch("tauEta", &tauEta_,"tauEta_/D");
   tree->Branch("genTauPt", &genTauPt_,"genTauPt_/D");
   tree->Branch("genTauEta", &genTauEta_,"genTauEta_/D");
   tree->Branch("genTauMatch", &genTauMatch_,"genTauMatch_/I");
   tree->Branch("nvtx",&nvtx_,"nvtx_/I");
   tree->Branch("dmf",&dmf_,"dmf_/I");
   tree->Branch("goodReco",&goodReco_,"goodReco_/I");
   tree->Branch("tauMass",&tauMass_,"tauMass_/D");

}


phase2Taus::~phase2Taus()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2Taus::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   nvtx_=vertices->size();
   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(prunedGenToken_, genParticles);
   std::vector<const reco::GenParticle*> GenTaus;

   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
     if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
     //if(TMath::Abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
     //if(TMath::Abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
   }

   for(auto genTau : GenTaus){
      tauPt_=-10;
      tauEta_=-10;
      tauMass_=-10;
      genTauPt_=-10;
      genTauEta_=-10;
     
     nvtx_=-10;
     dmf_=-10;
     goodReco_=0;
     genTauMatch_=0;

     std::vector<const reco::GenParticle*> genTauDaughters;
     findDaughters(genTau, genTauDaughters);
     reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
     genTauPt_  = (float) genTauVis.pt();
     genTauEta_ = (float) genTauVis.eta();

     //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
     for(const pat::Tau &tau : *taus){
       if (reco::deltaR(tau.eta(),tau.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5 && tau.tauID("decayModeFinding")>0.5){
	 genTauMatch_ = 1;
	 tauPt_  = tau.pt();
	 tauEta_ = tau.eta();
	 dmf_ = tau.decayMode();

	 goodReco_ = (bool) tau.tauID(tauID_) >0.5;

	 break;
       }
     }
     tree->Fill(); 
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2Taus::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2Taus::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2Taus::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}


//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector phase2Taus::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
  reco::Candidate::LorentzVector p4_vis(0,0,0,0);
  for(size_t i = 0; i < daughters.size(); ++i){
    if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
      p4_vis += daughters[i]->p4();
    }
  }
  return p4_vis;
}


//Creates a vector of all (including intermediate) daughters for a given mother particle
void phase2Taus::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
{
  unsigned numDaughters = mother->numberOfDaughters();
  if (numDaughters == 0) std::cout << " none ";
  for (unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
    if (daughter->status() == 1){  //status = 1 is a final state daughter
      daughters.push_back(daughter); 
    }
    if (daughter->status() == 2){  //status = 2 is an intermediate daughter; will decay further
      daughters.push_back(daughter); 
      findDaughters(daughter, daughters);
    }
  }
}

bool phase2Taus::isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}


//define this as a plug-in
DEFINE_FWK_MODULE(phase2Taus);
