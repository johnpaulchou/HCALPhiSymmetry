// *************F********* for reco data ********************** system
// include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"


#include "CondFormats/HcalObjects/interface/HcalQIEShape.h"
#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CondFormats/HcalObjects/interface/HcalPedestals.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TROOT.h>
#include <TSystem.h>
#include "TFile.h"
#include <TCanvas.h>
#include <cmath>
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <utility>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <set>
#include <memory>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

// class declaration
namespace cms {

class phiSym : public edm::one::EDAnalyzer<> {
 public:
  explicit phiSym (const edm::ParameterSet&);
  ~phiSym();

 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::Service<TFileService> fs;
//  std::set<std::tuple<int, int, int>> vetoEvents;
//  void loadVetoList(const std::string& jsonFile);


  std::string histfile;
  std::string textfile;

  TFile* mFile;
  FILE* tFile;

  Int_t runNumb, EventN, NvtxEv, NumNvtxEv1, NumNvtxEv2, NumNvtxEv3;
  Int_t i_HLT_HcalNZS_v9;
  Int_t i_HLT_HcalPhiSym_v10;
  Int_t i_HLT_L1SingleEG5_v1;
  Int_t i_HLT_L1SingleEG12_v6;
  Int_t i_HLT_ZeroBias_v7;
  Int_t i_HLT_Mu5_v21;
  Int_t i_HLT_Mu8_v19;
  Int_t i_HLT_Mu12_eta2p1_L1Mu10erJetC12WdEtaPhi1DiJetsC_v8 ;

  static const int N =20;
  TH1F *hcounter,*herun, *hlumi, *heventn, *hBX, *hvertex, *htrigger1, *htrigger2, *htrigger3, *htrigger4, *htrigger5, *htrigger6, *htrigger7;
  TH1F *hen[26][36][2][N];
  TH1F *henhbp[16][72][4][N], *henhbm[16][72][4][N];
  TH1F *henhep[14][72][7][N], *henhem[14][72][7][N];

  InputTag HBHENoiseFilterResultLabel_;
  HcalCalibrations calibs_;

  edm::EDGetTokenT<HFRecHitCollection> mhfreco;
  edm::EDGetTokenT<HBHERecHitCollection> mhbhereco;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
//  std::vector<std::string> triggerNamesSingleMu_;
//  std::vector<std::string> triggerNamesDoubleMu_;
//  std::string triggerType_; // "single" or "double"

};


// constructors and destructor

phiSym::phiSym(const edm::ParameterSet& iConfig)
{
  triggerResultsToken_ = (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults")));
//  triggerNamesSingleMu_ = iConfig.getUntrackedParameter<std::vector<std::string>>("triggerNamesSingleMu");
//  triggerNamesDoubleMu_ = iConfig.getUntrackedParameter<std::vector<std::string>>("triggerNamesDoubleMu");
//  triggerType_ = iConfig.getUntrackedParameter<std::string>("triggerType", "single");
  //metFiltersToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilters"));
  mhfreco  = consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("hfreco"));//for RECO
  mhbhereco = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbhereco"));//for RECO
  textfile = iConfig.getUntrackedParameter<string>("textFile");
  //std::string vetoJson = iConfig.getUntrackedParameter<std::string>("vetoJson");
//  loadVetoList(vetoJson);
}

phiSym::~phiSym()
{

}

//void phiSym::loadVetoList(const std::string& jsonFile) {
  //  boost::property_tree::ptree pt;
//    boost::property_tree::read_json(jsonFile, pt);

    // Iterate over all top-level keys
  //  for (const auto& key : pt) {
  //      const std::string& filter_name = key.first;
    //    std::cout << "Processing filter: " << filter_name << std::endl;

        // Iterate over each array item (each event/run pair)
    //    for (const auto& item : key.second) {
    //        try {
    //            int run = item.second.get<int>("run");
  //              int event = item.second.get<int>("event");
	//	int lumi = item.second.get<int>("lumi");
	//	vetoEvents.emplace(run, event, lumi);

                //std::cout << "  Run: " << run << ", Event: " << event << ", Lumi: "<< lumi << std::endl;
  //          }
  //          catch (boost::property_tree::ptree_error &e) {
//                std::cerr << "  Parsing error: " << e.what() << std::endl;
  //          }
  //      }
//    }
//}


// member functions

// ------------ method called to for each event  ------------
 void phiSym::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   printf("Starting :\n");
   edm::EventID eventId = iEvent.id();
   int runNumber = eventId.run ();
   int eventNumber = eventId.event ();
   int lumi = eventId.luminosityBlock();


   float eventWeight = 1;

   cout<<"weight = "<<eventWeight<<endl;
   runNumb=runNumber;
   EventN++;
   cout <<"event number----------------------------------------------------------------------------------------------- "<<EventN<<endl;

   NvtxEv = 0;
   i_HLT_HcalNZS_v9 = 0;
   i_HLT_HcalPhiSym_v10 = 0;
   i_HLT_L1SingleEG5_v1=0;
   i_HLT_L1SingleEG12_v6=0;
   i_HLT_ZeroBias_v7=0;
   i_HLT_Mu5_v21=0;
   i_HLT_Mu8_v19=0;
   i_HLT_Mu12_eta2p1_L1Mu10erJetC12WdEtaPhi1DiJetsC_v8 = 0;


   NvtxEv++;

   edm::Handle<edm::TriggerResults> triggerResults;
   edm::InputTag trigResultsTag("TriggerResults","","HLT");

   iEvent.getByToken(triggerResultsToken_, triggerResults);

//   if (!triggerResults.isValid()) {
//     edm::LogWarning("phiSym") << "TriggerResults not valid!";
//     return;
//   }

//   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);

//   const std::vector<std::string>& allowList = (triggerType_ == "single") ? triggerNamesSingleMu_ : triggerNamesDoubleMu_;
//   const std::vector<std::string>& vetoList  = (triggerType_ == "single") ? triggerNamesDoubleMu_ : triggerNamesSingleMu_;

//   bool passAllowed = false;
//   bool vetoed = false;

   // Check allowed triggers
//   for (const std::string& trig : allowList) {
//     unsigned int index = trigNames.triggerIndex(trig);
//     if (index < triggerResults->size() && triggerResults->accept(index)) {
//       passAllowed = true;
//       break;
//     }
//   }

   // Check veto triggers
  // for (const std::string& trig : vetoList) {
  //   unsigned int index = trigNames.triggerIndex(trig);
  //   if (index < triggerResults->size() && triggerResults->accept(index)) {
  //     vetoed = true;
  //     break;
  //   }
  // }

  // if (!passAllowed || vetoed) return;  // Reject if not in allow list or in veto list



   // MET filter handling
   //edm::Handle<edm::TriggerResults> metFiltersHandle;
   //iEvent.getByToken(metFiltersToken_, metFiltersHandle);

   //const edm::TriggerNames &metFilterNames = iEvent.triggerNames(*metFiltersHandle);


   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
   //std::vector<std::string> metFiltersToCheck = {
  //     "Flag_goodVertices",
    //   "Flag_globalSuperTightHalo2016Filter",
  //     "Flag_HBHENoiseFilter",//wiki says not needed in run 3
  //     "Flag_HBHENoiseIsoFilter",// wiki says not needed in run 3
  //     "Flag_EcalDeadCellTriggerPrimitiveFilter",
  //     "Flag_BadPFMuonFilter",
  //     "Flag_ecalBadCalibFilter",
  //     "Flag_eeBadScFilter",
  //     "Flag_BadChargedCandidateFilter",
//       "Flag_BadPFMuonDzFilter", // Mini v2
//       "Flag_hfNoisyHitsFilter", // Mini v2
//   };

//   if (vetoEvents.count(std::make_tuple(runNumber, eventNumber, lumi))) {
//        return; // Veto this event
//    }


   //bool passMETFilters = true;
   //for (const auto &filter : metFiltersToCheck)
   // {
   //  unsigned int index = metFilterNames.triggerIndex(filter);
       //cout<<metFiltersHandle->size()<<endl;
   //  if (index < metFiltersHandle->size() && !metFiltersHandle->accept(index))
   //  {
   //      passMETFilters = false;
   //      metFilterFailCounts_[filter]++; // Increment the counter for this filter
           // Do not break; we want to check all filters to see which ones failed
           // break;
   //     }

   // }

    // Skip event if it fails MET filters
    //if (!passMETFilters)
   //{
  //   return; // Skip this event
   //}




   hcounter->Fill(0); // total

   Int_t nBX=0;// iBX=1, nORBIT=0;
   nBX = iEvent.bunchCrossing();

   hcounter->Fill(1);

   //hf rechits
   Handle<HFRecHitCollection> hf_hits_h;
   iEvent.getByToken(mhfreco, hf_hits_h);
   const HFRecHitCollection* hf_hits = hf_hits_h.failedToGet () ? 0 : &*hf_hits_h;

   //hcal rechits
   Handle<HBHERecHitCollection> hbhe_hits_h;
   iEvent.getByToken(mhbhereco, hbhe_hits_h);
   const HBHERecHitCollection* hbhe_hits = hbhe_hits_h.failedToGet () ? 0 : &*hbhe_hits_h;



   double Etot=0;
   int jeta;//,jphi,jdepth;

   // FILE *input_file_1;
   Float_t FACTOR;
   int eventNum = abs(eventNumber);
   // ------------ HF -----------
   if (hf_hits_h.isValid()) {

     for (HFRecHitCollection::const_iterator hfhit=hf_hits->begin(); hfhit!=hf_hits->end(); hfhit++) {

       int ieta = hfhit->id().ieta();
       int iphi = hfhit->id().iphi();
       int depth = hfhit->id().depth();
       //cout<<"depth: "<<endl;
       //cout<<depth<<endl;
       FACTOR = 1.;

       double energy = hfhit->energy()*FACTOR;
       if (ieta>0) jeta = ieta-29;
       else jeta = 13-ieta-29;
       //     cout<<eventNumber<<endl;
       //cout<<eventNumber%N<<endl;
       hen[jeta][(iphi-1)/2][depth-1][eventNum%N]->Fill(energy, eventWeight );//****
       if (energy>10) Etot += energy;
     }
   }
   else printf("No HF RecHits: run= %d  ev= %d :\n",runNumber,eventNumber); // ------------

   // ------------ HBHE -----------

    if (hbhe_hits_h.isValid()) {
      hcounter->Fill(3);

     for (HBHERecHitCollection::const_iterator hbhehit=hbhe_hits->begin(); hbhehit!=hbhe_hits->end(); hbhehit++) {
       int ieta = hbhehit->id().ieta();
       int iphi = hbhehit->id().iphi();
       int depth = hbhehit->id().depth();

       FACTOR = 1.;

       double energy = hbhehit->energy()*FACTOR;

        if (abs(ieta)>16 || (abs(ieta)==16 && depth==4)) { // HE
        if (ieta>0) henhep[ieta-16][iphi-1][depth-1][eventNum%N]->Fill(energy, eventWeight);
        else  henhem[-ieta-16][iphi-1][depth-1][eventNum%N]->Fill(energy, eventWeight);
         if (energy>4) Etot += energy;
       }


		if (abs(ieta)<16){
	 if (ieta>0) henhbp[ieta-1][iphi-1][depth-1][eventNum%N]->Fill(energy, eventWeight);
	 else henhbm[-ieta-1][iphi-1][depth-1][eventNum%N]->Fill(energy, eventWeight);
	 if (energy>4) Etot += energy;
		}
		if ((abs(ieta)==16)&&(depth<4)){
                  if (ieta>0) henhbp[ieta-1][iphi-1][depth-1][eventNum%N]->Fill(energy, eventWeight);
		  else henhbm[-ieta-1][iphi-1][depth-1][eventNum%N]->Fill(energy, eventWeight);
		  if (energy>4) Etot += energy;
		}
	    }



   }


   else printf("No RecHits: run= %d  ev= %d :\n",runNumber,eventNumber); // ------------


   herun->Fill(runNumber);
   heventn->Fill(EventN);
   hlumi->Fill(lumi);
   hBX->Fill(nBX,Etot);
   hvertex->Fill(NvtxEv);


  return;
}


void phiSym::beginJob() {

  EventN=0;
  char htit[64];
  if ((tFile = fopen(textfile.c_str(),"w"))==NULL) printf("\nNo textfile open\n\n");

  herun = fs->make< TH1F>("herun","E(HF) vs Nrun;Nrun;GeV",5000,246000,251000);
  heventn = fs->make< TH1F>("heventn","E(HF) vs EventN;EventN;GeV",100000,0,100000);
  hlumi = fs->make< TH1F>("hlumi","E(HF) vs Lumi;EventN;GeV",100000,0,100000);
  hBX = fs->make< TH1F>("hBX","E(HF) vs nBX;BX;GeV",4096,0,4096);
  hcounter = fs->make< TH1F>("hcounter","hcounter",101,-0.5,100.5);
  hvertex = fs->make< TH1F>("hvertex","hvertex",100, 0, 99.5);
  htrigger1 = fs->make< TH1F>("hHLT_HcalPhiSym","hHLT_HcalPhiSym",100, 0, 100);
  htrigger2 = fs->make< TH1F>("hHLT_L1SingleEG5","hHLT_L1SingleEG5",100, 0, 100);
  htrigger3 = fs->make< TH1F>("hHLT_L1SingleEG12","hHLT_L1SingleEG12",100, 0, 100);
  htrigger4 = fs->make< TH1F>("hHLT_ZeroBias","hHLT_ZeroBias",100, 0, 100);
  htrigger5 = fs->make< TH1F>("hHLT_Mu5","hHLT_Mu5",100, 0, 100);
  htrigger6 = fs->make< TH1F>("hHLT_Mu8","hHLT_Mu8",100, 0, 100);
  htrigger7 = fs->make< TH1F>("hHLT_Mu12_eta2p1","hHLT_Mu12_eta2p1",100, 0, 100);
  // Initialize MET filter fail counts
  //for (const auto &filter : {
  //       "Flag_goodVertices",
  //       "Flag_globalSuperTightHalo2016Filter",
  //       "Flag_HBHENoiseFilter",
  //       "Flag_HBHENoiseIsoFilter",
  //       "Flag_EcalDeadCellTriggerPrimitiveFilter",
  //       "Flag_BadPFMuonFilter",
  //       "Flag_ecalBadCalibFilter",
  //       "Flag_eeBadScFilter",
  //       "Flag_BadChargedCandidateFilter",
  //       "Flag_BadPFMuonDzFilter", // Mini v2
  //       "Flag_hfNoisyHitsFilter", // Mini v2
  //   })
  //{
    //    metFilterFailCounts_[filter] = 0;



  TFileDirectory ESpec = fs->mkdir( "espec" );
  for (int i=0;i<13;i++) for (int j=0;j<36;j++) for (int k=0;k<2;k++) for (int l=0; l<N;l++){
    if (i>10 && j%2==0) continue;
    sprintf(htit,"E_+%d_%d_%d_%d",i+29,j*2+1,k+1,l);
    hen[i][j][k][l] = ESpec.make< TH1F>(htit,htit,1000,0,250); // E Rec
    sprintf(htit,"E_-%d_%d_%d_%d",i+29,j*2+1,k+1,l);
    hen[i+13][j][k][l] = ESpec.make< TH1F>(htit,htit,1000,0,250);
  }

  TFileDirectory EBSpec = fs->mkdir( "eHBspec" );
  for (int i=0;i<16;i++) for (int j=0;j<72;j++) for (int k=0;k<4;k++) for (int l=0; l<N;l++) {
    sprintf(htit,"E_+%d_%d_%d_%d",i+1,j+1,k+1,l);
    henhbp[i][j][k][l] = EBSpec.make< TH1F>(htit,htit,1000,0,250); // E Rec
    sprintf(htit,"E_-%d_%d_%d_%d",i+1,j+1,k+1,l);
    henhbm[i][j][k][l] = EBSpec.make< TH1F>(htit,htit,1000,0,250);
  }

  TFileDirectory EESpec = fs->mkdir( "eHEspec" );
  for (int i=0;i<14;i++) for (int j=0;j<72;j++) for (int k=0;k<7;k++) for (int l=0;l<N;l++) {
	if (i+16==16 && k!=3) continue;
	if (i+16==17 && k<1) continue;
	if (i+16==17 && k>2) continue;
	if (i+16==18 && k<1) continue;
	if (i+16==18 && k>4) continue;
	if (i+16==19 && k>5) continue;
	if (i+16==20 && k>5) continue;


	if ((i+16>20)&&((j+1)%2==0))continue;
	if ((i+16<26)&&(i+16>20)&&(k>5))continue;
	if ((i+16==29)&&(k>2))continue;


    sprintf(htit,"E_+%d_%d_%d_%d",i+16,j+1,k+1,l);
    henhep[i][j][k][l] = EESpec.make< TH1F>(htit,htit,1000,0,250); // E Rec
    sprintf(htit,"E_-%d_%d_%d_%d",i+16,j+1,k+1,l);
    henhem[i][j][k][l] = EESpec.make< TH1F>(htit,htit,1000,0,250);

  }

  std::cout<<std::endl<<"beginJob: histfile="<<histfile.c_str()<<"  textfile="<<textfile.c_str()<<std::endl;
}

void phiSym::endJob() {

    fprintf(tFile,"#RunN %d   Events processed %d\n",runNumb,EventN);
  std::cout<<"endJob: histos processing..."<<std::endl;
  std::cout<<"RunN= "<<runNumb<<"  Events processed= "<<EventN<<std::endl;
  // Output MET filter fail counts
  //std::cout << "---------------------->" << std::endl;
  //std::cout << "MET Filter Fail Counts:" << std::endl;
  //int totalFilteredEvents = 0;
  // for (const auto &filterCount : metFilterFailCounts_)
  //{
  //   std::cout << "  " << filterCount.first << ": " << filterCount.second << std::endl;
  //   totalFilteredEvents += filterCount.second;
  //}
  //std::cout << "Total events filtered by MET filters: " << totalFilteredEvents << std::endl;
  //std::cout << "----------------------/" << std::endl;

  fclose(tFile);

  //std::cout<<std::endl<<" --endJob-- done"<<std::endl;
  //return;
}
}
//define this as a plug-in
DEFINE_FWK_MODULE(phiSym);
