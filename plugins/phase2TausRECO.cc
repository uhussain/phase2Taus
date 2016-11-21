// -*- C++ -*-
//
// Package:    RecoTauTag/phase2TausRECO
// Class:      phase2TausRECO
// 
/**\class phase2TausRECO phase2TausRECO.cc RecoTauTag/phase2TausRECO/plugins/phase2TausRECO.cc

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
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// user include files
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

/*
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
*/

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include <unordered_map>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class phase2TausRECO : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2TausRECO(const edm::ParameterSet&);
      ~phase2TausRECO();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  edm::InputTag  tracksTag_;
  edm::InputTag  timesTag_;
  edm::InputTag  timeResosTag_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<reco::PFTauCollection>      tauToken_;
  edm::EDGetTokenT<reco::PFJetCollection>    jetSrc_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorSrc_;
  //edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;

  std::string tauID_;

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
  bool initializeTrackTiming( edm::Handle<edm::View<reco::Track> > *tracks, reco::VertexCollection *vtxs, edm::ValueMap<float> *times, edm::ValueMap<float> *timeResos);

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
phase2TausRECO::phase2TausRECO(const edm::ParameterSet& iConfig):
  tracksTag_       ("tracks"),
  timesTag_        ("times"),
  timeResosTag_    ("timesResos"),
  vtxToken_        (consumes<reco::VertexCollection>          (iConfig.getParameter<edm::InputTag>("vertices")      )),
  tauToken_        (consumes<reco::PFTauCollection>             (iConfig.getParameter<edm::InputTag>("taus")          )),
  jetSrc_          (consumes<reco::PFJetCollection>           (iConfig.getParameter<edm::InputTag>("jets")          )),
  discriminatorSrc_(consumes<reco::PFTauDiscriminator>        (iConfig.getParameter<edm::InputTag>("discriminator") ))
  //prunedGenToken_  (consumes<std::vector<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genParticles")  ))
{
  consumes<edm::View<reco::Track> >(tracksTag_);
  consumes<edm::ValueMap<float> >(timesTag_);
  consumes<edm::ValueMap<float> >(timeResosTag_);
  //now do what ever initialization is needed
  usesResource("TFileService");

   //tauID_    = iConfig.getParameter<std::string>("tauID");
   edm::Service<TFileService> fs;

   tree = fs->make<TTree>(      "Ntuple",      "Ntuple"          );
   tree->Branch("tauPt",        &tauPt_,       "tauPt_/D"        );
   tree->Branch("tauEta",       &tauEta_,      "tauEta_/D"       );
   tree->Branch("genTauPt",     &genTauPt_,    "genTauPt_/D"     );
   tree->Branch("genTauEta",    &genTauEta_,   "genTauEta_/D"    );
   tree->Branch("genTauMatch",  &genTauMatch_, "genTauMatch_/I"  );
   tree->Branch("nvtx",         &nvtx_,        "nvtx_/I"         );
   tree->Branch("dmf",          &dmf_,         "dmf_/I"          );
   tree->Branch("goodReco",     &goodReco_,    "goodReco_/I"     );
   tree->Branch("tauMass",      &tauMass_,     "tauMass_/D"      );

}


phase2TausRECO::~phase2TausRECO()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2TausRECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //get track colletions
   edm::Handle<edm::View<reco::Track> > tracks;
   iEvent.getByLabel(tracksTag_,tracks);

   edm::Handle<edm::ValueMap<float> > times;
   iEvent.getByLabel(timesTag_,times);

   edm::Handle<edm::ValueMap<float> > timeResos;
   iEvent.getByLabel(timeResosTag_,timeResos);

   Handle<reco::VertexCollection> vtxs;
   iEvent.getByToken(vtxToken_, vtxs);
   nvtx_ = vtxs->size();

   //get objects
   Handle<reco::PFJetCollection> jetObjects;
   iEvent.getByToken(jetSrc_, jetObjects);

   edm::Handle<pat::PFTauCollection> taus;
   iEvent.getByToken(tauToken_, taus);

   Handle<reco::PFTauDiscriminator> discriminator;
   iEvent.getByToken(discriminatorSrc_, discriminator);

   Handle<reco::PFTauDiscriminator> DMF; 
   iEvent.getByToken("hpsPFTauDiscriminationByDecayModeFindingOldDMs",DMF);

   std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection(iEvent);
   if (GenObjects.size()!=0) return;

   //edm::Handle<std::vector<reco::GenParticle> > genParticles;
   //iEvent.getByToken(prunedGenToken_, genParticles);

   //initialize timing objects
   initializeTrackTiming( tracks, vtxs, times, timeResos);

   std::vector<const reco::GenParticle*> GenTaus;
   for(std::vector<reco::GenParticle>::const_iterator genParticle = GenObjects->begin(); genParticle != GenObjects->end(); genParticle++){
     if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
     //if(TMath::Abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
     //if(TMath::Abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
   }
   
   std::vector<reco::PFTau> goodTaus;
   for (unsigned int iTau = 0; iTau<taus->size() ; ++iTau){
     reco::PFTauRef tauCandidate(taus, iTau);
     if((*DMF)[tauCandidate] > -1)
       if((*discriminator)[tauCandidate] > 0.5)
	 goodTaus.push_back(*tauCandidate);
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
      //for(const pat::Tau &tau : *taus){
      for(const reco::PFTau tau : goodTaus){
	if (reco::deltaR(tau.eta(),tau.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5){
	  genTauMatch_ = 1;
	  tauPt_  = tau.pt();
	  tauEta_ = tau.eta();
	  dmf_ = tau.decayMode();
	  //goodReco_ = (bool) tau.tauID(tauID_) >0.5;

	  //get the matched vertex
	  int vtx_index = -1;

	  // find the 4D vertex this muon is best associated to..
	  float max_weight = 0.f;
	  for( unsigned i = 0; i < vtxs->size(); ++i ) {
	    const auto& vtx = (*vtxs)[i];      
	    const float weight = vtx.trackWeight(tau.leadPFChargedHadrCand());// check me -> Get track ref for charged hadron candidate
	    if( weight > max_weight ) {
	      max_weight = weight;
	      vtx_index = i;
	    }
	  }
    


	  break;
	}
      }
      tree->Fill(); 
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2TausRECO::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2TausRECO::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2TausRECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}


//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector phase2TausRECO::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
  reco::Candidate::LorentzVector p4_vis(0,0,0,0);
  for(size_t i = 0; i < daughters.size(); ++i){
    if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
      p4_vis += daughters[i]->p4();
    }
  }
  return p4_vis;
}


//Creates a vector of all (including intermediate) daughters for a given mother particle
void phase2TausRECO::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
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

bool phase2TausRECO::isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}

bool phase2TausRECO::initializeTrackTiming( edm::Handle<edm::View<reco::Track> >*tracks, 
						reco::VertexCollection *vtxs, 
						edm::ValueMap<float>   *times, 
						edm::ValueMap<float>   *timeResos)
{
  //make a map of vertices to track refs within cuts
  std::unordered_multimap<unsigned,reco::TrackBaseRef> vertices_to_tracks_z, vertices_to_tracks_zt3, vertices_to_tracks_zt4, vertices_to_tracks_zt5, vertices_to_tracks_zt6 ;
  for( unsigned i = 0; i < tracks->size(); ++i ) {
    auto ref = tracks->refAt(i);
    const float time = (*times)[ref];
    const float timeReso = (*timeResos)[ref] != 0.f ? (*timeResos)[ref] : 0.170f;
    for( unsigned ivtx = 0; ivtx < vtxs->size(); ++ivtx ) {
      const auto& thevtx = (*vtxs)[ivtx];
      const float dz = std::abs(ref->dz(thevtx.position()));
      const float dt = std::abs(time - thevtx.t());
      const bool useTime = (thevtx.t() != 0.);

      const float base_cut = std::sqrt(thevtx.tError()*thevtx.tError()
				       + timeReso*timeReso);

      const float time_cut3 = 3.f*base_cut;
      const float time_cut4 = 4.f*base_cut;
      const float time_cut5 = 5.f*base_cut;
      const float time_cut6 = 6.f*base_cut;

      const bool keepz = ( dz < 0.2f );
      const bool keept3 = (!useTime || dt < time_cut3);
      const bool keept4 = (!useTime || dt < time_cut4);
      const bool keept5 = (!useTime || dt < time_cut5);
      const bool keept6 = (!useTime || dt < time_cut6);
      
      if( keepz ) {
	vertices_to_tracks_z.emplace(ivtx,ref);
	if( keept6 ) {
	  vertices_to_tracks_zt6.emplace(ivtx,ref);
	}
	if( keept5 ) {
	  vertices_to_tracks_zt5.emplace(ivtx,ref);
	}
	if( keept4 ) {
	  vertices_to_tracks_zt4.emplace(ivtx,ref);
	}
	if( keept3 ) {
	  vertices_to_tracks_zt3.emplace(ivtx,ref);
	}
      }      
    }
  }
  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(phase2TausRECO);
