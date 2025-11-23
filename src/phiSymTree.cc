#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"
#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"


#include <fnmatch.h>

#include "TTree.h"

class phiSymTree : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:
  explicit phiSymTree (const edm::ParameterSet&);
  ~phiSymTree();

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void endRun(const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // tree and its contents
  TTree* tree_;
  unsigned long int run_, lumi_, event_;
  const static int MAXNHITS=20000;
  const static int MAXNTRIGS=100;
  int nHits_;
  float hitEnergy_[MAXNHITS];
  int ieta_[MAXNHITS], iphi_[MAXNHITS], depth_[MAXNHITS];
  bool isHF_[MAXNHITS];
  int nTrigs_;
  float trigPt_[MAXNTRIGS], trigEta_[MAXNTRIGS], trigPhi_[MAXNTRIGS];
  std::vector<std::string> trigNames_;
  
  // parameters
  bool undoRespCorr_;
  double HBHEthreshold_, HFthreshold_;
  std::vector<std::string> hltPaths_;
  
  // collection tokens
  edm::EDGetTokenT<HBHERecHitCollection> hbheRecHitToken_;
  edm::EDGetTokenT<HFRecHitCollection> hfRecHitToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryToken_;
  edm::ESGetToken<HcalRespCorrs, HcalRespCorrsRcd> hcalRespCorrToken_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> hcalTopologyToken_;

  // other stuff
  std::map<std::string, int> triggerMap_;
  HLTConfigProvider hltConfig_;
  edm::Service<TFileService> fs;
  
};


////////////////////////////////////////////////////////
// constructor and destructor
////////////////////////////////////////////////////////

phiSymTree::phiSymTree(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");
  
  hbheRecHitToken_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheRecHits"));
  hfRecHitToken_ = consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("hfRecHits"));
  triggerResultsToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
  triggerSummaryToken_ = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerSummary"));
  hcalRespCorrToken_ = esConsumes<HcalRespCorrs, HcalRespCorrsRcd>();
  hcalTopologyToken_ = esConsumes<HcalTopology, HcalRecNumberingRecord>();
						       

  undoRespCorr_ = iConfig.getParameter<bool>("undoRespCorr");
  HBHEthreshold_ = iConfig.getParameter<double>("HBHEthreshold");
  HFthreshold_ = iConfig.getParameter<double>("HFthreshold");
  hltPaths_ = iConfig.getParameter<std::vector<std::string>>("hltPaths");
}

phiSymTree::~phiSymTree() {}



////////////////////////////////////////////////////////
// beginRun:
// creates the tree
////////////////////////////////////////////////////////


void phiSymTree::beginRun(const edm::Run& run, const edm::EventSetup& setup) {
  bool changed = true;
  if (!hltConfig_.init(run, setup, "HLT", changed)) {
    edm::LogError("phiSymTree") << "HLTConfigProvider::init failed";
    return;
  }
}
void phiSymTree::endRun(const edm::Run& run, const edm::EventSetup& setup) {
}



////////////////////////////////////////////////////////
// beginJob:
// creates the tree
////////////////////////////////////////////////////////

void phiSymTree::beginJob() {

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("phisymtree", "A tree for information for the HCAL phi symmetry calibration");

  // Create branches. Names = how they'll appear in the ROOT file
  tree_->Branch("run",  &run_,  "run/l");
  tree_->Branch("lumi",  &lumi_,  "lumi/l");
  tree_->Branch("event",  &event_,  "event/l");
  tree_->Branch("nHits",  &nHits_,  "nHits/I");
  tree_->Branch("hitEnergy", hitEnergy_, "hitEnergy[nHits]/F");
  tree_->Branch("ieta", ieta_, "ieta[nHits]/I");
  tree_->Branch("iphi", iphi_, "iphi[nHits]/I");
  tree_->Branch("depth", depth_, "depth[nHits]/I");
  tree_->Branch("isHF", isHF_, "isHF[nHits]/O");
  tree_->Branch("nTrigs", &nTrigs_, "nTrigs/I");
  tree_->Branch("trigPt", trigPt_, "trigPt[nTrigs]/F");
  tree_->Branch("trigEta", trigEta_, "trigEta[nTrigs]/F");
  tree_->Branch("trigPhi", trigPhi_, "trigPhi[nTrigs]/F");
  tree_->Branch("trigNames", &trigNames_);
}

////////////////////////////////////////////////////////
// analyzer:
// fills the tree
////////////////////////////////////////////////////////

void phiSymTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<HBHERecHitCollection> hbheRecHits;
  edm::Handle<HFRecHitCollection> hfRecHits;
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<trigger::TriggerEvent> triggerSummary;


  iEvent.getByToken(hbheRecHitToken_, hbheRecHits);
  iEvent.getByToken(hfRecHitToken_, hfRecHits);
  iEvent.getByToken(triggerResultsToken_, triggerResults);
  iEvent.getByToken(triggerSummaryToken_, triggerSummary);
 
  if (!hbheRecHits.isValid()) {
    edm::LogWarning("phiSymTree") << "HBHERecHitCollection not found!";
    return;
  }
  if (!hfRecHits.isValid()) {
    edm::LogWarning("phiSymTree") << "HFRecHitCollection not found!";
    return;
  }

  if (!triggerResults.isValid()) {
    edm::LogWarning("MyAnalyzer") << "Trigger results not found!";
    return;
  }

  if (!triggerSummary.isValid()) {
    edm::LogWarning("MyAnalyzer") << "Trigger summary not found!";
    return;
  }

  
  const HcalTopology* theTopology = &iSetup.getData(hcalTopologyToken_);
  const HcalRespCorrs* constresp = &iSetup.getData(hcalRespCorrToken_);
  std::unique_ptr<HcalRespCorrs>  respCorrs = std::make_unique<HcalRespCorrs>(*constresp);
  respCorrs->setTopo(theTopology);

  // get the trigger information
  trigNames_.clear();
  nTrigs_=0;

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {

    // require that the trigger is fired
    if(!triggerResults->accept(i)) continue;

    // get the path name
    std::string pathName=triggerNames.triggerName(i);

    // store the trigger names and the map count for later
    triggerMap_[pathName]++;

    // select those that match the given hltPaths
    bool foundPath=false;
    for(const auto& path : hltPaths_)
      if(fnmatch(path.c_str(), pathName.c_str(), 0)==0) {
	foundPath=true;
	break;
      }
    if(!foundPath) continue;
    
    
    // find the last module label that is saved
    unsigned int pathIndex = hltConfig_.triggerIndex(pathName);
    if (pathIndex >= hltConfig_.size()) continue;
    const std::vector<std::string>& moduleLabels = hltConfig_.moduleLabels(pathIndex);
    std::vector<std::string> savedLabels;
    std::string lastLabel;
    for (const auto& mod : moduleLabels) if(hltConfig_.saveTags(mod)) savedLabels.push_back(mod);
    if(savedLabels.size()>0) lastLabel=savedLabels.back();
    else continue;

    // find the trigger object based on the label
    auto filterIndex = triggerSummary->filterIndex(edm::InputTag(lastLabel, "", "HLT"));
    if(filterIndex >= triggerSummary->sizeFilters()) continue;
    //    const trigger::Vids& ids  = triggerSummary->filterIds(filterIndex);
    const trigger::Keys& keys = triggerSummary->filterKeys(filterIndex);
    const trigger::TriggerObjectCollection& objects = triggerSummary->getObjects();

    for (size_t i = 0; i < keys.size(); ++i) {
      const trigger::TriggerObject& obj = objects[keys[i]];

      trigPt_[nTrigs_]=obj.pt();
      trigEta_[nTrigs_]=obj.eta();
      trigPhi_[nTrigs_]=obj.phi();
      trigNames_.push_back(pathName);
      nTrigs_++;
    }
  }
  
  
  // start filling the other tree variables
  nHits_=0;
  run_   = iEvent.id().run();
  lumi_  = iEvent.luminosityBlock();
  event_ = iEvent.id().event();
  
  
  for (const auto& hit : *hbheRecHits) {
    float energy = hit.energy();
    if(undoRespCorr_) energy = energy / respCorrs->getValues(hit.detid())->getValue();
    if(energy<=HBHEthreshold_) continue;
    hitEnergy_[nHits_]= energy;
    ieta_[nHits_]=hit.id().ieta();
    iphi_[nHits_]=hit.id().iphi();
    depth_[nHits_]=hit.id().depth();
    isHF_[nHits_]=false; // denote this is _not_ in the HF
    ++nHits_;
  }

  for (const auto& hit : *hfRecHits) {
    float energy = hit.energy();
    if(undoRespCorr_) energy = energy / respCorrs->getValues(hit.detid())->getValue();
    if(energy<=HFthreshold_) continue;
    hitEnergy_[nHits_]= energy;
    ieta_[nHits_]=hit.id().ieta();
    iphi_[nHits_]=hit.id().iphi();
    depth_[nHits_]=hit.id().depth();
    isHF_[nHits_]=true; // denote this _is_ in the HF
    ++nHits_;
  }

  // fill the tree
  tree_->Fill();
  
}

////////////////////////////////////////////////////////
// endjob
////////////////////////////////////////////////////////

void phiSymTree::endJob() {
  // Copy map entries into a vector to sort
  std::vector<std::pair<std::string, int>> vec(triggerMap_.begin(), triggerMap_.end());
  std::sort(vec.begin(), vec.end(),
	    [](const std::pair<std::string,int>& a, const std::pair<std::string,int>& b) {
	      return a.second > b.second; // descending
	    });
  for (const auto& entry : vec) {
    std::cout << "Trigger name: " << entry.first << ", Count: " << entry.second << std::endl;
  }
  return;
}






//define this as a plug-in
DEFINE_FWK_MODULE(phiSymTree);
