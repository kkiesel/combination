#include <math.h>
#include <iostream>
#include <map>
#include <utility>
#include <string>
#include <regex>
#include <vector>
#include <time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TTree.h"

#include "TreeParticles.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
#include "cxxopts/include/cxxopts.hpp"

using namespace std;

enum class Selection: char {
  original, di_cleaned, lep_cleaned, st_cleaned, bjet_cleaned, dilepst_cleaned, all_cleaned,
  original_ee
};

const map<Selection, string> selectionNames {
  {Selection::original, "original"},
  {Selection::di_cleaned, "di_cleaned"},
  {Selection::lep_cleaned, "lep_cleaned"},
  {Selection::st_cleaned, "st_cleaned"},
  {Selection::bjet_cleaned, "bjet_cleaned"},
  {Selection::dilepst_cleaned, "dilepst_cleaned"},
  {Selection::all_cleaned, "all_cleaned"}
};

enum class Histogram {
  gen,
  isr, isrUp, isrDn,
  nopu, nopuUp, nopuDn, puUp, puDn,
  jesUp, jesDn, jerUp, jerDn,
};

const map<Histogram, string> histogramNames {
  {Histogram::gen, "gen"},
  {Histogram::isr, "isr"},
  {Histogram::isrUp, "isrUp"},
  {Histogram::isrDn, "isrDn"},
  {Histogram::nopu, "nopu"},
  {Histogram::nopuUp, "nopuUp"},
  {Histogram::nopuDn, "nopuDn"},
  {Histogram::puUp, "puUp"},
  {Histogram::puDn, "puDn"},
  {Histogram::jesUp, "jesUp"},
  {Histogram::jesDn, "jesDn"},
  {Histogram::jerUp, "jerUp"},
  {Histogram::jerDn, "jerDn"},
};

enum class Region {
  sR, eCR, jCR, genE
};

typedef pair<UShort_t,UShort_t> SignalPoint;

class CombinationHistogramProducer : public TSelector {
 public:

  CombinationHistogramProducer();
  virtual ~CombinationHistogramProducer() { }

  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree){};
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Bool_t Notify();
  virtual Int_t Version() const { return 2; }

  void initHistograms();
  void fillHistograms(Selection, Region, bool);


  void resetSelection();
  void defaultSelection();
  void fillSelection(string const& s, float w, bool);
  void fillSignalSelection(string const& s, float w);
  void fillSelectedPhotons(const tree::Particle& p);
  tree::Jet* matchedJet(const tree::Particle& p);
  float getPhotonWeight(const tree::Photon& p);
  void printEventId();

  void initTriggerStudies();
  void fillTriggerStudies();

  void initUncut();
  void fillUncut();

  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Muon>> muons;
  TTreeReaderValue<std::vector<tree::Particle>> genJets;
  TTreeReaderValue<std::vector<tree::GenParticle>> genParticles;
  TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<tree::MET> metGen;
  TTreeReaderValue<tree::MET> metRaw;
  TTreeReaderValue<tree::MET> met_JESu;
  TTreeReaderValue<tree::MET> met_JESd;
  TTreeReaderValue<tree::MET> met_JERu;
  TTreeReaderValue<tree::MET> met_JERd;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<std::vector<Float_t>> pdf_weights;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Int_t> nTracksPV;
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<Float_t> rho;
  TTreeReaderValue<Int_t> nTruePV;
  TTreeReaderValue<Int_t> nISR;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<ULong64_t> evtNo;
  TTreeReaderValue<Bool_t> hlt_photon90_ht600;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_ht600;
  TTreeReaderValue<Bool_t> hlt_ht800;
  TTreeReaderValue<Bool_t> hlt_el27;
  TTreeReaderValue<Int_t> hlt_ht600_pre;

  // signal scan
  TTreeReaderValue<UShort_t> signal_nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

  bool electroweak;
  bool isData;
  bool isSignal;
  bool genHt600;
  bool noPromptPhotons;
  string inputName;
  map<string,Weighter> weighters;

  double startTime;
  TRandom2 rand;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;


  SignalPoint sp_; // signal point

  float weight_ = 1.;
  float htg_ = 0;
  float ptmiss_ = 0;

  map<SignalPoint, map<Selection, TH2F>> nominalHists_; // data, mc, signal
  map<SignalPoint, map<Selection, TH2F>> eControlHists_; // data, mc, signal
  map<SignalPoint, map<Selection, TH2F>> genEHists_;
  map<SignalPoint, map<Selection, map<float, TH2F>>> jControlHists_; // data, mc, signal
  map<SignalPoint, map<Selection, map<Histogram, TH2F>>> mcHists_; // signal, mc
  map<SignalPoint, map<Selection, map<unsigned, TH2F>>> weightedHists_; // signal, mc
  map<SignalPoint, int> nGens_;

//  ClassDef(CombinationHistogramProducer, 1) // not for root compilation
};

