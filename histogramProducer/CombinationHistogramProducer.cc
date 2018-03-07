#include "CombinationHistogramProducer.h"
using namespace std;

void CombinationHistogramProducer::printEventId() {
  string fname = fReader.GetTree()->GetCurrentFile()->GetName();
  smatch sm;
  regex e(".*2016(.)-.*");
  regex_match(fname, sm, e);
  string r = sm.size() > 1 ? sm[1] : fname;
  cout << *runNo << ":" << *lumNo << ":" << *evtNo << " in Run " << r << endl;
}

string CombinationHistogramProducer::id() {
  return to_string(*runNo) + ":" + to_string(*lumNo) + ":" + to_string(*evtNo);
}

void CombinationHistogramProducer::initHistograms() {
  TH2F hist("", ";#it{p}_{T}^{miss} (GeV);#it{H}_{T}^{#gamma} (GeV)", 200, 0, 2000, 1, 700, 2000);
  nGens_[sp_] = 0;
  nominalHists_[sp_] = map<Selection, TH2F>();
  nominalHistsGG_[sp_] = map<Selection, TH2F>();
  nominalHistsWG_[sp_] = map<Selection, TH2F>();
  eControlHists_[sp_] = map<Selection, TH2F>();
  genEHists_[sp_] = map<Selection, TH2F>();
  jControlHists_[sp_] = map<Selection, map<float, TH2F>>();
  mcHists_[sp_] = map<Selection, map<Histogram, TH2F>>();
  weightedHists_[sp_] = map<Selection, map<unsigned, TH2F>>();
  for (const auto& sel : selectionNames) {
    nominalHists_.at(sp_)[sel.first] = hist;
    nominalHistsGG_.at(sp_)[sel.first] = hist;
    nominalHistsWG_.at(sp_)[sel.first] = hist;
    eControlHists_.at(sp_)[sel.first] = hist;
    genEHists_.at(sp_)[sel.first] = hist;

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
    for (float i=0.60; i<1.2; i+= 0.01) {
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
    if (sel==Selection::original and not sp_.first and not sp_.second) {
      isLowEmht_ = htg_<2000;
      jCRTree_original->Fill();
    }
  } else if (region == Region::genE) {
    genEHists_.at(sp_).at(sel).Fill(isSel?ptmiss_:-1, htg_, weight_);
  } else if (region == Region::sR) {
    nominalHists_.at(sp_).at(sel).Fill(isSel?ptmiss_:-1, htg_, weight_);
    if (*signal_nBinos == 2) nominalHistsGG_.at(sp_).at(sel).Fill(isSel?ptmiss_:-1, htg_, weight_);
    if (*signal_nBinos == 1) nominalHistsWG_.at(sp_).at(sel).Fill(isSel?ptmiss_:-1, htg_, weight_);
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
  nGoodVertices(fReader, "nGoodVertices"), // TODO: cut on nGoodVertices?
  genHt(fReader, "genHt"),
  nTruePV(fReader, "true_nPV"),
  nISR(fReader, "nISR"),
  runNo(fReader, "runNo"),
  lumNo(fReader, "lumNo"),
  evtNo(fReader, "evtNo"),
  hlt_photon90_ht600(fReader, "HLT_Photon90_CaloIdL_PFHT600_v"),
  hlt_ht600(fReader, "HLT_PFHT600_v"),
  hlt_diphoton(fReader, "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v"),
  hlt_ht600_pre(fReader, "HLT_PFHT600_v_pre"),
  signal_nBinos(fReader, "signal_nBinos"),
  signal_m1(fReader, "signal_m1"),
  signal_m2(fReader, "signal_m2"),
  eventIdsGamGam("diPhoton.txt"),
  eventIdsGamLep("lepPhoton.txt"),
  eventIdsGamStg("../../other_photon_limits/SUS-16-046/signalRegion_st.txt"),
  eventIdsGamHtg("../SUS-16-047_ids.txt"),
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
  weighters["puWeightUp"] = Weighter("../../CMSSW/treewriter/v22/CMSSW_8_0_26_patch1/src/TreeWriter/PUreweighting/data/puWeights.root", puUp);
  weighters["puWeightDn"] = Weighter("../../CMSSW/treewriter/v22/CMSSW_8_0_26_patch1/src/TreeWriter/PUreweighting/data/puWeights.root", puDn);
  weighters.at("sf_photon_id_loose").fillOverflow2d();
  weighters.at("sf_photon_pixel").fillOverflow2d();

  genHt600 = inputName.find("HT-0to600") != string::npos;
  noPromptPhotons = inputName.find("QCD_HT") != string::npos
    || inputName.find("TTJets") != string::npos
    || inputName.find("WJetsToLNu") != string::npos
    || inputName.find("ZJetsToNuNu") != string::npos;

  jCRTree_original = new TTree("jCRTree", "");
  jCRTree_original->Branch("met", &ptmiss_);
  jCRTree_original->Branch("weight", &weight_);
  jCRTree_original->Branch("lowEMHT", &isLowEmht_, "lowEMHT/O");
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

bool CombinationHistogramProducer::isDiPhotonSel(bool pho=true, bool eb=true) {
  if (isData and !*hlt_diphoton) return false;
  float m = met->p.Pt();
  if (m<100) return false;
  vector<tree::Photon*> phos;
  if (pho and eb ) {
    for (auto &p : *photons) {
      if (p.isMedium && p.p.Pt()>40 && p.r9>.5 && p.r9<1 && p.sigmaIetaIeta>0.005 && !p.hasPixelSeed && p.isEB()) {
        phos.push_back(&p);
      }
    }
  }
  if (pho and !eb) {
    int nPho = 0;
    for (auto &p : *photons) {
      if (p.isMedium && p.p.Pt()>40 && p.r9>.5 && p.r9<1 && p.sigmaIetaIeta>0.005 && !p.hasPixelSeed) {
        if (p.isEB()) nPho ++;
        phos.push_back(&p);
      }
    }
    if (nPho>1) return false;
  }
  if (!pho and eb) {
    int nPho = 0;
    for (auto &p : *photons) {
      if (p.isMedium && p.p.Pt()>40 && p.r9>.5 && p.r9<1 && p.sigmaIetaIeta>0.005 && p.isEB()) {
        if (!p.hasPixelSeed) nPho ++;
        phos.push_back(&p);
      }
    }
    if (nPho>1) return false;
  }
  if (!pho and !eb) {
    int nPho = 0;
    for (auto &p : *photons) {
      if (p.isMedium && p.p.Pt()>40 && p.r9>.5 && p.r9<1 && p.sigmaIetaIeta>0.005) {
        if (!p.hasPixelSeed && p.isEB()) nPho ++;
        phos.push_back(&p);
      }
    }
    if (nPho>1) return false;
  }
  if (phos.size()<2) return false;
  if (phos.at(0)->p.DeltaR(phos.at(1)->p) < .3 or minv(phos.at(0)->p,phos.at(1)->p) < 105) return false;
  // lepton veto
  for (auto& lep : *electrons) {
    if (!lep.isLoose and lep.p.Pt() < 25) continue;
    bool overlap = false;
    for (auto& p : phos) { if (p->p.DeltaR(lep.p)>.3) overlap = true;}
    if (overlap) continue;
    return false;
  }
  for (auto& lep : *muons) {
    if (lep.p.Pt() < 25 and lep.PFminiIso < 0.2) continue;
    bool overlap = false;
    for (auto& p : phos) { if (p->p.DeltaR(lep.p)>.3) overlap = true;}
    if (overlap) continue;
    return false;
  }
  return true;
}


bool CombinationHistogramProducer::isLepSel(bool pho=true, bool eb=true) {
  // code copied from danilo

  std::vector<tree::Electron const *> el_comb;
  std::vector<tree::Muon const *> mu_comb;
  tree::Photon* lead_pho;


  for (auto g: *photons) {
    if (g.p.Pt()>40 and ((eb&&g.isEB()) || (!eb&&g.isEE())) and g.isLoose && (pho ^ g.hasPixelSeed)) {
      lead_pho = &g;
      break;
    }
  }
  if (!lead_pho) return false; // no photon found

  bool leptoVeto = false;
  for (tree::Electron const &ele: *electrons) {
     if (ele.p.Pt() < 25 || ele.isMedium == 0) { continue; }
     if (fabs(ele.p.Eta()) < 1.4442) {
        if (ele.r9 > 0.5 && ele.SigmaIEtaIEtaFull5x5 < 0.00998 &&
              fabs(ele.dEtaAtVtx) < 0.00311 && fabs(ele.dPhiAtVtx) < 0.103 &&
              ele.HoverE < 0.253 && ele.MissHits <= 1 && fabs(ele.EoverPInv) < 0.134 &&
              ele.ConvVeto == 1 && ele.PFminiIso < 0.1) {
                 leptoVeto = true;
                 el_comb.push_back(&ele);
        }
     }
     if (fabs(ele.p.Eta()) > 1.56 && fabs(ele.p.Eta()) < 2.5) {
        if (ele.r9 > 0.8 && ele.SigmaIEtaIEtaFull5x5 < 0.0298 &&
              fabs(ele.dEtaAtVtx) < 0.00609 && fabs(ele.dPhiAtVtx) < 0.045 &&
              ele.HoverE < 0.0878 && ele.MissHits <= 1 && fabs(ele.EoverPInv) < 0.13 &&
              ele.ConvVeto == 1 && ele.PFminiIso < 0.1) {
                 leptoVeto = true;
                 el_comb.push_back(&ele);
        }
     }
  }

  for (tree::Muon const &mu: *muons) {
     if (mu.p.Pt() < 25 || mu.isMedium == 0) { continue; }
     if (fabs(mu.p.Eta()) < 2.4 && mu.PFminiIso < 0.2 &&
           fabs(mu.d0) < 0.05 && fabs(mu.dZ) < 0.1) {
              leptoVeto = true;
              mu_comb.push_back(&mu);
     }
  }
  if (mu_comb.empty() and el_comb.empty()) {
    return false; // no lepton found
  }

  //leading lepton
  tree::Particle const *leadLep;
  bool leadLep_ele = false;
  bool MuAndEle = false;
  if (leptoVeto) {
    if (el_comb.size() > 0) {
      if (mu_comb.size() > 0) {
        if (el_comb[0]->p.Pt() > mu_comb[0]->p.Pt()) {
          leadLep = el_comb[0];
          leadLep_ele = true;
        } else {
          leadLep = mu_comb[0];
        }
        MuAndEle = true;
      } else {
        leadLep = el_comb[0];
        leadLep_ele = true;
      }
    } else {
      leadLep = mu_comb[0];
    }
  }
  //additional cuts from lepton analysis
  float invmassPhoLep = 0;
  float MT_LepMet = 0;
  if (leptoVeto) {
    if (leadLep->p.DeltaR(lead_pho->p) < 0.8) {
      leptoVeto = false;
    } else {
      for (tree::Electron const &ele: *electrons) {
        if (ele.p.DeltaR(lead_pho->p) < 0.3 && ele.p.Pt() > 2.0) leptoVeto = false;
      }
      for (tree::Muon const &mu: *muons) {
        if (mu.p.DeltaR(lead_pho->p) < 0.3 && mu.p.Pt() > 2.0) leptoVeto = false;
      }
    }
    if (leptoVeto) {
      invmassPhoLep = minv(lead_pho->p, leadLep->p);
      if (leadLep_ele && fabs(invmassPhoLep-91.1876) < 10) leptoVeto = false;
    }
    if (leptoVeto) {
      if (MuAndEle) {
        if (mt(el_comb[0]->p, met->p) < 100 && mt(mu_comb[0]->p, met->p) < 100) leptoVeto = false;
      } else {
        MT_LepMet = mt(leadLep->p, met->p);
        if (*evtNo == 570123420) {MT_LepMet = 99.7227;} // due to different phi definition of electron
        if (MT_LepMet < 100) leptoVeto = false;
      }
    }
  }
  return leptoVeto;
}

bool CombinationHistogramProducer::isStSel(bool pho=true, bool eb=true) {
  float st = met->p.Pt();
  if (st<300) return false; // ptmiss here

  // selection photons
  vector<tree::Photon*> joPhotons;
  for (auto& photon: *photons) {
    if (photon.isLoose && (pho ^ photon.hasPixelSeed) && ((eb && photon.isEB()) || (!eb && photon.isEE())) && photon.p.Pt()> 15
      && photon.seedCrystalE/photon.p.Pt()>.3 && photon.sigmaIetaIeta>0.001
      && photon.sigmaIphiIphi>0.001) {
        joPhotons.push_back(&photon);
        st += photon.p.Pt();
      }
  }
  if (joPhotons.empty()) return false;

  // jet selection
  bool cleanMet = true;
  bool drPhoJet = true;
  for (auto& j : *jets) {
    if (j.isLoose && j.p.Pt()>100 && fabs(j.p.Eta())<3 && !j.hasPhotonMatch
      && !j.hasElectronMatch && !j.hasMuonMatch && fabs(met->p.DeltaPhi(j.p))<0.3) cleanMet = false;
    auto dr = joPhotons.at(0)->p.DeltaR(j.p);
    auto dPt = fabs(joPhotons.at(0)->p.Pt()-j.p.Pt())/joPhotons.at(0)->p.Pt();
    if (j.isLoose && j.p.Pt()>30 && dr<.5 && (dr>.1 || dPt>.5)) drPhoJet = false;
  }
  if (!cleanMet || !drPhoJet) return false;

  if (joPhotons.at(0)->p.Pt()<180 || st<600) return false;

  float mtg = joPhotons.size() ? mt(met->p, joPhotons.at(0)->p) : 0;
  if (mtg<300) return false;
  return true;
}

Bool_t CombinationHistogramProducer::Process(Long64_t entry)
{
  fReader.Next();
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

  resetSelection();
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
  bool genE = false;
  if (selPhotons.size()) {
    auto g = selPhotons.at(0);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    orthogonal = .3<dPhi && dPhi<2.84;
    genE = genMatchToId(*g, *genParticles, 11);
    isGenEclean = isData || isSignal || !genE;
  }
  bool signalSelwoGen = !cutPrompt && selPhotons.size() && htg_ > 700 && (*hlt_photon90_ht600 || !isData) && orthogonal;
  bool signalSel = signalSelwoGen && isGenEclean;
  bool signalGenE = signalSelwoGen && genE;
  fillHistograms(Selection::original, Region::sR, signalSel);
  if (signalGenE) fillHistograms(Selection::original, Region::genE, true);

  /* overlap check
  // combination stuff
  if (signalSel && met->p.Pt()>350) {
    printEventId();
    auto pt = selPhotons.size() ? selPhotons.at(0)->p.Pt() : 0;
    cout << "ptmiss = " << met->p.Pt() << "\thtg = " << htg_ << "\tphoton pt = " << pt << endl;
    cout << "selected by gam+gam: " << eventIdsGamGam.check(*runNo, *lumNo, *evtNo) << " " << isDiPhotonSel() << endl;
    cout << "selected by gam+lep: " << eventIdsGamLep.check(*runNo, *lumNo, *evtNo) << " " << isLepSel() << endl;
    cout << "selected by gam+stg: " << eventIdsGamStg.check(*runNo, *lumNo, *evtNo) << " " << isStSel() << endl;
    cout << "selected by gam+htg: " << eventIdsGamHtg.check(*runNo, *lumNo, *evtNo) << " " << "1" << endl;
  }
  */
  bool isDiPhoton = isDiPhotonSel();
  bool isLepPhoton = isLepSel();
  bool isStPhoton = isStSel();

  fillHistograms(Selection::di_cleaned, Region::sR, signalSel && !isDiPhoton);
  fillHistograms(Selection::lep_cleaned, Region::sR, signalSel && !isLepPhoton);
  fillHistograms(Selection::st_cleaned, Region::sR, signalSel && !isStPhoton);
  fillHistograms(Selection::dilep_cleaned, Region::sR, signalSel && !isDiPhoton && !isLepPhoton);
  fillHistograms(Selection::all_cleaned, Region::sR, signalSel && !isDiPhoton && !isLepPhoton && !isStPhoton);
  if (signalGenE && !isDiPhoton) fillHistograms(Selection::di_cleaned, Region::genE, true);
  if (signalGenE && !isLepPhoton) fillHistograms(Selection::lep_cleaned, Region::genE, true);
  if (signalGenE && !isStPhoton) fillHistograms(Selection::st_cleaned, Region::genE, true);
  if (signalGenE && !isDiPhoton && !isLepPhoton) fillHistograms(Selection::dilep_cleaned, Region::genE, true);
  if (signalGenE && !isDiPhoton && !isLepPhoton && !isStPhoton) fillHistograms(Selection::all_cleaned, Region::genE, true);

  weight_ = *mc_weight * *pu_weight * *hlt_ht600_pre;
  if (!selPhotons.size() && htg_ > 700 && (*hlt_ht600 || !isData)) {
    fillHistograms(Selection::original, Region::jCR, true);
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
    if (!cutPrompt && orthogonal) {
      bool isDiPhotonEl = isDiPhotonSel(false);
      bool isLepPhotonEl = isLepSel(false);
      bool isStPhotonEl = isStSel(false);
      fillHistograms(Selection::original, Region::eCR, true);
      if (!isDiPhotonEl) fillHistograms(Selection::di_cleaned, Region::eCR, true);
      if (!isLepPhotonEl) fillHistograms(Selection::lep_cleaned, Region::eCR, true);
      if (!isStPhotonEl) fillHistograms(Selection::st_cleaned, Region::eCR, true);
      if (!isDiPhotonEl && !isLepPhotonEl) fillHistograms(Selection::dilep_cleaned, Region::eCR, true);
      if (!isDiPhotonEl && !isLepPhotonEl && !isStPhotonEl) fillHistograms(Selection::all_cleaned, Region::eCR, true);
    }
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
      nominalHistsGG_.at(spIt.first).at(selIt.first).Scale(scale);
      nominalHistsGG_.at(spIt.first).at(selIt.first).Write("nominalGG", TObject::kWriteDelete);
      nominalHistsWG_.at(spIt.first).at(selIt.first).Scale(scale);
      nominalHistsWG_.at(spIt.first).at(selIt.first).Write("nominalWG", TObject::kWriteDelete);
      eControlHists_.at(spIt.first).at(selIt.first).Scale(scale);
      eControlHists_.at(spIt.first).at(selIt.first).Write("eControl", TObject::kWriteDelete);
      genEHists_.at(spIt.first).at(selIt.first).Scale(scale);
      genEHists_.at(spIt.first).at(selIt.first).Write("genE", TObject::kWriteDelete);
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
  file.cd();
  jCRTree_original->Write("simpleTree", TObject::kWriteDelete);
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

void CombinationHistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
  weight_ = 0;
}

void CombinationHistogramProducer::FillNgen(const string& f) {
  for (auto& it: nGens_) {
    it.second += getNgenFromFile(f, it.first.first, it.first.second);
  }
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
    ch.AddFile(argv[i]);
  }
  chp.Init(&ch);
  chp.SlaveBegin(&ch);
  for(unsigned i=0; i<ch.GetEntries(); i++) {
    chp.Process(i);
  }
  for (int i=1;i<argc;i++) {
    chp.FillNgen(argv[i]);
  }
  chp.SlaveTerminate();
  chp.Terminate();
}



