#include "CombinationHistogramProducer.h"
using namespace std;

void CombinationHistogramProducer::printEventId() {
  cout << *runNo << ":" << *lumNo << ":" << *evtNo << " in " << inputName << endl;
}

void CombinationHistogramProducer::initHistograms() {
  TH2F hist("", ";#it{p}_{T}^{miss} (GeV);#it{H}_{T}^{#gamma} (GeV)", 200, 0, 2000, 1, 700, 2000);
  nGens_[sp_] = getNgenFromFile(fReader, *signal_m1, *signal_m2);
  nominalHists_[sp_] = map<Selection, TH2F>();
  eControlHists_[sp_] = map<Selection, TH2F>();
  jControlHists_[sp_] = map<Selection, map<float, TH2F>>();
  mcHists_[sp_] = map<Selection, map<Histogram, TH2F>>();
  weightedHists_[sp_] = map<Selection, map<unsigned, TH2F>>();
  for (const auto& sel : selectionNames) {
    nominalHists_.at(sp_)[sel.first] = hist;
    eControlHists_.at(sp_)[sel.first] = hist;

    mcHists_.at(sp_)[sel.first] = map<Histogram, TH2F>();
    for (const auto& histo : histogramNames) {
      mcHists_.at(sp_).at(sel.first)[histo.first] = hist;
    }
    weightedHists_.at(sp_)[sel.first] = map<unsigned, TH2F>();
    unsigned nPdfs = isSignal ? 9 : pdf_weights->size();
    for (unsigned i=0; i<nPdfs; i++) {
      weightedHists_.at(sp_).at(sel.first)[i] = hist;
    }
    jControlHists_.at(sp_)[sel.first] = map<float, TH2F>();
    for (float i=0.85; i<1.2; i+= 0.01) {
      jControlHists_.at(sp_).at(sel.first)[i] = hist;
    }
  }
}

void CombinationHistogramProducer::fillHistograms(Selection sel, Region region, bool isSel) {
  ptmiss_ = met->p.Pt();
  if (region == Region::eCR) {
    eControlHists_.at(sp_).at(sel).Fill(isSel?ptmiss_:-1, htg_, weight_);
  } else if (region == Region::jCR) {
    for (auto &it : jControlHists_.at(sp_).at(sel)) {
      it.second.Fill(isSel?it.first*ptmiss_:-1, htg_, weight_);
    }
  } else if (region == Region::sR) {
    nominalHists_.at(sp_).at(sel).Fill(isSel?ptmiss_:-1, htg_, weight_);
    if (!isData) {
      for (auto &it : weightedHists_.at(sp_).at(sel)) {
        it.second.Fill(isSel?ptmiss_:-1, htg_, weight_*pdf_weights->at(it.first));
      }

      // other mc uncertainties
      float isrPt = 0;
      if (electroweak) {
        isrPt = getPropagatorPt(*genParticles);
      }
      auto isrWeight = electroweak ? isrReweightingEWK(isrPt) : isrReweighting(*nISR);
      auto isrWeightE = electroweak ? isrReweightingEWK(isrPt, true) : isrReweighting(*nISR, true);

      auto hs = &mcHists_.at(sp_).at(sel);
      hs->at(Histogram::gen).Fill(isSel?metGen->p.Pt():-1, htg_, weight_);
      hs->at(Histogram::isr).Fill(isSel?ptmiss_:-1, htg_, weight_*isrWeight);
      hs->at(Histogram::isrUp).Fill(isSel?ptmiss_:-1, htg_, weight_*(isrWeight+isrWeightE));
      hs->at(Histogram::isrDn).Fill(isSel?ptmiss_:-1, htg_, weight_*(isrWeight-isrWeightE));
      hs->at(Histogram::nopu).Fill(isSel?ptmiss_:-1, htg_, weight_/ *pu_weight);
      if (*nTruePV>=23) hs->at(Histogram::nopuUp).Fill(isSel?ptmiss_:-1, htg_, weight_/ *pu_weight);
      else              hs->at(Histogram::nopuDn).Fill(isSel?ptmiss_:-1, htg_, weight_/ *pu_weight);
      hs->at(Histogram::puUp).Fill(isSel?ptmiss_:-1, htg_, weight_*weighters.at("puWeightUp").getWeight(*nTruePV)/ *pu_weight);
      hs->at(Histogram::puDn).Fill(isSel?ptmiss_:-1, htg_, weight_*weighters.at("puWeightDn").getWeight(*nTruePV)/ *pu_weight);
      hs->at(Histogram::jesUp).Fill(isSel?met_JESu->p.Pt():-1, htg_, weight_);
      hs->at(Histogram::jesDn).Fill(isSel?met_JESd->p.Pt():-1, htg_, weight_);
      hs->at(Histogram::jerUp).Fill(isSel?met_JERu->p.Pt():-1, htg_, weight_);
      hs->at(Histogram::jerDn).Fill(isSel?met_JERd->p.Pt():-1, htg_, weight_);
    }
  }
}

CombinationHistogramProducer::CombinationHistogramProducer():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  muons(fReader, "muons"),
  genJets(fReader, "genJets"),
  genParticles(fReader, "genParticles"),
  met(fReader, "met"),
  metGen(fReader, "met_gen"),
  met_JESu(fReader, "met_JESu"),
  met_JESd(fReader, "met_JESd"),
  met_JERu(fReader, "met_JERu"),
  met_JERd(fReader, "met_JERd"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
  pdf_weights(fReader, "pdf_weights"),
  nGoodVertices(fReader, "nGoodVertices"),
  genHt(fReader, "genHt"),
  nTruePV(fReader, "true_nPV"),
  nISR(fReader, "nISR"),
  runNo(fReader, "runNo"),
  lumNo(fReader, "lumNo"),
  evtNo(fReader, "evtNo"),
  hlt_photon90_ht600(fReader, "HLT_Photon90_CaloIdL_PFHT600_v"),
  hlt_ht600(fReader, "HLT_PFHT600_v"),
  hlt_ht600_pre(fReader, "HLT_PFHT600_v_pre"),
  signal_nBinos(fReader, "signal_nBinos"),
  signal_m1(fReader, "signal_m1"),
  signal_m2(fReader, "signal_m2"),
  startTime(time(NULL)),
  rand()
{
}

void CombinationHistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  tree->GetEntry(0); // to get current file
  inputName = tree->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;
  isSignal = inputName.find("SMS") != string::npos ||
             inputName.find("GMSB") != string::npos;
  electroweak = inputName.find("TChi") != string::npos ||
                inputName.find("GMSB") != string::npos;

  bool puSummer16 = inputName.find("T5bbbbZg") != string::npos ||
                    inputName.find("T5ttttZg") != string::npos ||
                    inputName.find("T5qqqqHg") != string::npos ||
                    inputName.find("GMSB") != string::npos;

  weighters["sf_photon_id_loose"] = Weighter("../../phd/plotter/data/dataMcScaleFactors_80X.root", "EGamma_SF2D");
  weighters["sf_photon_pixel"] = Weighter("../../phd/plotter/data/ScalingFactors_80X_Summer16.root", "Scaling_Factors_HasPix_R9 Inclusive");
  string puUp = "pileupWeightUp_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
  string puDn = "pileupWeightDown_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
  if (puSummer16) {
    puUp = "pileupWeightUp_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
    puDn = "pileupWeightDown_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
  }
  weighters["puWeightUp"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puUp);
  weighters["puWeightDn"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puDn);
  weighters.at("sf_photon_id_loose").fillOverflow2d();
  weighters.at("sf_photon_pixel").fillOverflow2d();

  genHt600 = inputName.find("HT-0to600") != string::npos;
  noPromptPhotons = inputName.find("QCD_HT") != string::npos
    || inputName.find("TTJets") != string::npos
    || inputName.find("WJetsToLNu") != string::npos
    || inputName.find("ZJetsToNuNu") != string::npos;
}

void CombinationHistogramProducer::defaultSelection()
{
  for (auto& mu : *muons) {
    if (mu.p.Pt() < 15) continue;
    if (indexOfMatchedParticle<tree::Photon*>(mu, selPhotons, .3) >= 0) continue;
    selMuons.push_back(&mu);
  }
  for (auto& el : *electrons) {
    if (!el.isLoose || el.p.Pt() < 15) continue;
    if (indexOfMatchedParticle<tree::Photon*>(el, selPhotons, .3) >= 0) continue;
    selElectrons.push_back(&el);
  }
  for (auto& jet : *jets) {
    if (!jet.isLoose
      || indexOfMatchedParticle<tree::Photon*>(jet, selPhotons, .4) >= 0
      || jet.p.Pt() < 30 || fabs(jet.p.Eta()) > 3) continue;
    selJets.push_back(&jet);
    if (jet.bDiscriminator > bTaggingWorkingPoints.at(CSVv2M) && fabs(jet.p.Eta()) < 2.5)
      selBJets.push_back(&jet);
  }
}

float CombinationHistogramProducer::getPhotonWeight(const tree::Photon& p) {
  float weight = 1.;
  if (!isData) {
    auto pt = p.p.Pt();
    auto eta = p.p.Eta();
    // sf for id and electron veto
    weight = weighters.at("sf_photon_id_loose").getWeight(eta, pt) * weighters.at("sf_photon_pixel").getWeight(fabs(eta), pt);
    // trigger efficiency
    weight *= fabs(eta)<photonsEtaMaxBarrel ? 0.964 : 0.94;
  }
  return weight;
}

Bool_t CombinationHistogramProducer::Process(Long64_t entry)
{
  fReader.SetLocalEntry(entry);

  if (genHt600 && *genHt>600) {
    return kTRUE;
  }

  sp_.first = *signal_m1;
  sp_.second = *signal_m2;
  if (!nominalHists_.count(sp_)) {
    initHistograms();
  }

  if (isSignal && unMatchedSuspiciousJet(*jets, *genJets)) {
    nGens_.at(sp_) --; // this is not considered in nGen for tree weights
    return kTRUE;
  }

  bool cutPrompt = noPromptPhotons && count_if(genParticles->begin(), genParticles->end(), [] (const tree::GenParticle& p) { return p.pdgId==22 && p.promptStatus == DIRECTPROMPT;});

  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  weight_ = *mc_weight * *pu_weight;
  if (selPhotons.size()) weight_ *= getPhotonWeight(*selPhotons.at(0));

  htg_ = 0;
  for (auto& p : selPhotons) htg_ += p->p.Pt();
  for (auto& p : selJets) htg_ += p->p.Pt();

  bool orthogonal = true;
  bool isGenEclean = true;
  if (selPhotons.size()) {
    auto g = selPhotons.at(0);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    orthogonal = .3<dPhi && dPhi<2.84;
    bool genE = genMatchToId(*g, *genParticles, 11);
    isGenEclean = isData || isSignal || !genE;
  }
  bool signalSel = !cutPrompt && selPhotons.size() && htg_ > 700 && (*hlt_photon90_ht600 || !isData) && orthogonal && isGenEclean;
  fillHistograms(Selection::original, Region::sR, signalSel);
  // todo : genE for electron closure
/*
  float mt = selPhotons.size() ? TMath::Sqrt( 2*selPhotons.at(0)->p.Pt()*met->p.Pt()*(1-TMath::Cos(selPhotons.at(0)->p.DeltaPhi(met->p)))) : 0;
  float st = met->p.Pt();
  for (auto& photon: *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.isEB() && photon.p.Pt()> 15) st += photon.p.Pt();
  }
  bool st_selection = selPhotons.size() && selPhotons.at(0)->p.Pt()>180 && mt > 300 && st > 600;
  bool lep_selection = mt > 100
    && count_if(photons->begin(), photons->end(), [] (const tree::Photon& p) { return p.isLoose && p.p.Pt()>35 && p.r9 >.5 && !p.hasPixelSeed && p.isEB();})
    && (count_if(electrons->begin(), electrons->end(), [] (const tree::Electron& p) { return p.isMedium && p.p.Pt()>25 && p.rIso < .1;}) // TOOD: r9 cut
       || count_if(muons->begin(), muons->end(), [] (const tree::Muon& p) { return p.p.Pt()>25 && p.rIso < .2;})); // TOOD: medium id
  bool hasBtags = count_if(jets->begin(), jets->end(), [] (const tree::Jet& p) { return p.bDiscriminator>0.8 && p.p.Pt()>30 && fabs(p.p.Eta())<2.5;}); // TODO: correct b-discriminator
  bool diPhoton = 1 < count_if(photons->begin(), photons->end(), [] (const tree::Photon& p) { return p.isMedium && p.p.Pt()>40 && p.r9 >.5 && !p.hasPixelSeed && p.isEB();});

  fillHistograms(Selection::di_cleaned, Region::sR, signalSel && !diPhoton);
  fillHistograms(Selection::lep_cleaned, Region::sR, signalSel && !lep_selection);
  fillHistograms(Selection::st_cleaned, Region::sR, signalSel && !st_selection);
  fillHistograms(Selection::bjet_cleaned, Region::sR, signalSel && !hasBtags);
  fillHistograms(Selection::dilepst_cleaned, Region::sR, signalSel && !lep_selection && !diPhoton && !st_selection);
  fillHistograms(Selection::all_cleaned, Region::sR, signalSel && !lep_selection && !diPhoton && !st_selection && !hasBtags);
*/
  if (!selPhotons.size() && htg_ > 700 && (*hlt_ht600 || !isData)) {
    fillHistograms(Selection::original, Region::jCR, false);
  }
  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // electron sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();
  weight_ = *mc_weight * *pu_weight;
  if (selPhotons.size()) weight_ *= getPhotonWeight(*selPhotons.at(0));

  htg_ = 0;
  for (auto& p : selPhotons) htg_ += p->p.Pt();
  for (auto& p : selJets) htg_ += p->p.Pt();
  if (selPhotons.size() && htg_ > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto g = selPhotons.at(0);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    orthogonal = .3<dPhi && dPhi<2.84;
    if (!cutPrompt && orthogonal) fillHistograms(Selection::original, Region::eCR, false);
  }
  resetSelection();
  return kTRUE;
}

template<typename T>
void save2File(const map<string,map<string,T>>& hMaps, TFile& file)
{
  for (auto& hMapIt : hMaps) {
    if (!file.Get(hMapIt.first.c_str())) {
      file.mkdir(hMapIt.first.c_str());
    }
    file.cd(hMapIt.first.c_str());
    for (auto& h : hMapIt.second) {
      h.second.Write(h.first.c_str(), TObject::kWriteDelete);
    }
    file.cd();
  }
}

void CombinationHistogramProducer::Terminate()
{
  auto outputName = getOutputFilename(inputName, "combiHists");
  TFile file(outputName.c_str(), "RECREATE");

  for (auto& spIt : nominalHists_) {
    auto dir = (to_string(spIt.first.first) + "_" + to_string(spIt.first.second));
    float scale = isData ? 1 : 1./nGens_.at(spIt.first);
    file.mkdir(dir.c_str());
    for (auto& selIt : spIt.second) {
      auto subDir = selectionNames.at(selIt.first).c_str();
      file.cd(dir.c_str());
      gDirectory->mkdir(subDir);
      gDirectory->cd(subDir);
      nominalHists_.at(spIt.first).at(selIt.first).Scale(scale);
      nominalHists_.at(spIt.first).at(selIt.first).Write("nominal", TObject::kWriteDelete);
      eControlHists_.at(spIt.first).at(selIt.first).Scale(scale);
      eControlHists_.at(spIt.first).at(selIt.first).Write("eControl", TObject::kWriteDelete);
      gDirectory->mkdir("jCR");
      gDirectory->mkdir("weights");
      gDirectory->mkdir("systematics");
      gDirectory->cd("jCR");
      for (auto scaleIt : jControlHists_.at(spIt.first).at(selIt.first)) {
        scaleIt.second.Scale(scale);
        scaleIt.second.Write(to_string(scaleIt.first).c_str(), TObject::kWriteDelete);
      }
      gDirectory->cd("../weights");
      for (auto scaleIt : weightedHists_.at(spIt.first).at(selIt.first)) {
        scaleIt.second.Scale(scale);
        scaleIt.second.Write(to_string(scaleIt.first).c_str(), TObject::kWriteDelete);
      }
      gDirectory->cd("../systematics");
      for (auto scaleIt : mcHists_.at(spIt.first).at(selIt.first)) {
        scaleIt.second.Scale(scale);
        scaleIt.second.Write(histogramNames.at(scaleIt.first).c_str(), TObject::kWriteDelete);
      }
    }
  }

  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

void CombinationHistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
}

Bool_t CombinationHistogramProducer::Notify() {
  for (auto& it: nGens_) {
    it.second += getNgenFromFile(fReader, it.first.first, it.first.second);
  }
  return kTRUE;
}


int main(int argc, char** argv) {

  /*
  cxxopts::Options options("CombinationHistogramProducer", "Writes histograms from ROOT TTrees");
  options.add_options()
    ("f,file", "File name", cxxopts::value<std::string>());
  auto result = options.parse(argc, argv);
  if (not result.count("file")) {
    cout << "Please provide an input file" << endl;
    return -1;
  }
  auto file = result["file"].as<std::string>();
  */
  CombinationHistogramProducer chp;
  TChain ch("TreeWriter/eventTree");
  for (int i=1;i<argc;i++) {
    cout << "+ Adding file " << argv[i] << endl;
    ch.AddFile(argv[i]);
  }
  chp.Init(&ch);
  chp.SlaveBegin(&ch);
  for(unsigned i=0; i<ch.GetEntries(); i++) {
    chp.Process(i);
  }
  chp.SlaveTerminate();
  chp.Terminate();
}



