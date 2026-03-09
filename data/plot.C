#include <iostream>
#include <string>

#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TSystem.h"

std::string make_path(const std::string& dir, const std::string& name) {
  if (dir.empty()) {
    return name;
  }
  if (dir.back() == '/') {
    return dir + name;
  }
  return dir + "/" + name;
}

void draw_single(TH1D& hist, const std::string& output_file, bool logy) {
  TCanvas canvas("canvas_single", "canvas_single", 900, 700);
  canvas.SetMargin(0.12, 0.04, 0.12, 0.06);
  canvas.SetLogy(logy);
  hist.SetLineWidth(2);
  hist.Draw("HIST");
  canvas.SaveAs(output_file.c_str());
}

void draw_pair(TH1D& hist1,
               TH1D& hist2,
               const char* label1,
               const char* label2,
               const std::string& output_file,
               bool logy) {
  TCanvas canvas("canvas_pair", "canvas_pair", 900, 700);
  canvas.SetMargin(0.12, 0.04, 0.12, 0.06);
  canvas.SetLogy(logy);

  hist1.SetLineWidth(2);
  hist2.SetLineWidth(2);
  hist2.SetLineColor(kRed + 1);

  hist1.Draw("HIST");
  hist2.Draw("HIST SAME");

  TLegend legend(0.72, 0.76, 0.92, 0.90);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&hist1, label1, "l");
  legend.AddEntry(&hist2, label2, "l");
  legend.Draw();

  canvas.SaveAs(output_file.c_str());
}

void plot() {
  const std::string kInputFile = "/eos/home-y/yiyangz/public/CMSDAS/data/2025G_Parking0.root";
  const std::string kTreePath = "rootuple/mm_tree";
  const std::string kOutputDir = "results";

  const bool kRequireTrigger = false;
  const bool kRequireOppositeCharge = false;
  const bool kRequireVertexProb = false;
  const bool kRequireSingleCandidate = false;
  const double kMinVertexProb = 0.01;

  const int kMassBins = 70;
  const double kMassMin = 8.0;
  const double kMassMax = 12.0;

  const int kDimuonPtBins = 130;
  const double kDimuonPtMin = 0.0;
  const double kDimuonPtMax = 130.0;

  const int kMuonPtBins = 100;
  const double kMuonPtMin = 0.0;
  const double kMuonPtMax = 100.0;

  const int kEtaBins = 48;
  const double kEtaMin = -2.4;
  const double kEtaMax = 2.4;

  const int kPhiBins = 64;
  const double kPhiMin = -3.2;
  const double kPhiMax = 3.2;

  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutputDir.c_str(), true);

  TChain chain(kTreePath.c_str());
  const int files_added = chain.Add(kInputFile.c_str());
  if (files_added == 0) {
    std::cerr << "No input file matched: " << kInputFile << std::endl;
    return;
  }

  TLorentzVector* dimuon_p4 = nullptr;
  TLorentzVector* muonP_p4 = nullptr;
  TLorentzVector* muonM_p4 = nullptr;
  Float_t vProb = 0.0;
  UInt_t trigger = 0;
  Int_t charge = 0;
  UInt_t nonia = 0;

  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("dimuon_p4", 1);
  chain.SetBranchStatus("muonP_p4", 1);
  chain.SetBranchStatus("muonM_p4", 1);
  chain.SetBranchAddress("dimuon_p4", &dimuon_p4);
  chain.SetBranchAddress("muonP_p4", &muonP_p4);
  chain.SetBranchAddress("muonM_p4", &muonM_p4);

  const bool has_vProb = (chain.GetBranch("vProb") != nullptr);
  const bool has_trigger = (chain.GetBranch("trigger") != nullptr);
  const bool has_charge = (chain.GetBranch("charge") != nullptr);
  const bool has_nonia = (chain.GetBranch("nonia") != nullptr);

  if (has_vProb) {
    chain.SetBranchStatus("vProb", 1);
    chain.SetBranchAddress("vProb", &vProb);
  }
  if (has_trigger) {
    chain.SetBranchStatus("trigger", 1);
    chain.SetBranchAddress("trigger", &trigger);
  }
  if (has_charge) {
    chain.SetBranchStatus("charge", 1);
    chain.SetBranchAddress("charge", &charge);
  }
  if (has_nonia) {
    chain.SetBranchStatus("nonia", 1);
    chain.SetBranchAddress("nonia", &nonia);
  }

  if (kRequireVertexProb && !has_vProb) {
    std::cerr << "Requested vProb selection, but branch vProb is missing." << std::endl;
    return;
  }
  if (kRequireTrigger && !has_trigger) {
    std::cerr << "Requested trigger selection, but branch trigger is missing." << std::endl;
    return;
  }
  if (kRequireOppositeCharge && !has_charge) {
    std::cerr << "Requested charge selection, but branch charge is missing." << std::endl;
    return;
  }
  if (kRequireSingleCandidate && !has_nonia) {
    std::cerr << "Requested nonia selection, but branch nonia is missing." << std::endl;
    return;
  }

  TH1D h_dimuon_mass("h_dimuon_mass", "Dimuon mass;M_{#mu#mu} (GeV);Events", kMassBins, kMassMin, kMassMax);
  TH1D h_dimuon_pt("h_dimuon_pt", "Dimuon p_{T};p_{T} (GeV);Events", kDimuonPtBins, kDimuonPtMin, kDimuonPtMax);
  TH1D h_dimuon_eta("h_dimuon_eta", "Dimuon #eta;#eta;Events", kEtaBins, kEtaMin, kEtaMax);
  TH1D h_dimuon_phi("h_dimuon_phi", "Dimuon #phi;#phi (rad);Events", kPhiBins, kPhiMin, kPhiMax);

  TH1D h_muonP_pt("h_muonP_pt", "Muon p_{T};p_{T} (GeV);Events", kMuonPtBins, kMuonPtMin, kMuonPtMax);
  TH1D h_muonM_pt("h_muonM_pt", "Muon p_{T};p_{T} (GeV);Events", kMuonPtBins, kMuonPtMin, kMuonPtMax);
  TH1D h_muonP_eta("h_muonP_eta", "Muon #eta;#eta;Events", kEtaBins, kEtaMin, kEtaMax);
  TH1D h_muonM_eta("h_muonM_eta", "Muon #eta;#eta;Events", kEtaBins, kEtaMin, kEtaMax);
  TH1D h_muonP_phi("h_muonP_phi", "Muon #phi;#phi (rad);Events", kPhiBins, kPhiMin, kPhiMax);
  TH1D h_muonM_phi("h_muonM_phi", "Muon #phi;#phi (rad);Events", kPhiBins, kPhiMin, kPhiMax);

  const Long64_t n_entries = chain.GetEntries();
  Long64_t n_filled = 0;

  for (Long64_t i = 0; i < n_entries; ++i) {
    chain.GetEntry(i);

    if (dimuon_p4 == nullptr || muonP_p4 == nullptr || muonM_p4 == nullptr) {
      continue;
    }
    if (kRequireTrigger && trigger == 0) {
      continue;
    }
    if (kRequireOppositeCharge && charge != 0) {
      continue;
    }
    if (kRequireVertexProb && vProb <= kMinVertexProb) {
      continue;
    }
    if (kRequireSingleCandidate && nonia != 1) {
      continue;
    }

    h_dimuon_mass.Fill(dimuon_p4->M());
    h_dimuon_pt.Fill(dimuon_p4->Pt());
    h_dimuon_eta.Fill(dimuon_p4->Eta());
    h_dimuon_phi.Fill(dimuon_p4->Phi());

    h_muonP_pt.Fill(muonP_p4->Pt());
    h_muonM_pt.Fill(muonM_p4->Pt());
    h_muonP_eta.Fill(muonP_p4->Eta());
    h_muonM_eta.Fill(muonM_p4->Eta());
    h_muonP_phi.Fill(muonP_p4->Phi());
    h_muonM_phi.Fill(muonM_p4->Phi());

    ++n_filled;
  }

  std::cout << "Entries in tree: " << n_entries << std::endl;
  std::cout << "Entries used: " << n_filled << std::endl;

  draw_single(h_dimuon_mass, make_path(kOutputDir, "dimuon_mass.pdf"), false);
  draw_single(h_dimuon_pt, make_path(kOutputDir, "dimuon_pt.pdf"), true);
  draw_single(h_dimuon_eta, make_path(kOutputDir, "dimuon_eta.pdf"), false);
  draw_single(h_dimuon_phi, make_path(kOutputDir, "dimuon_phi.pdf"), false);

  draw_pair(h_muonP_pt, h_muonM_pt, "#mu^{+}", "#mu^{-}", make_path(kOutputDir, "muon_pt.pdf"), true);
  draw_pair(h_muonP_eta, h_muonM_eta, "#mu^{+}", "#mu^{-}", make_path(kOutputDir, "muon_eta.pdf"), false);
  draw_pair(h_muonP_phi, h_muonM_phi, "#mu^{+}", "#mu^{-}", make_path(kOutputDir, "muon_phi.pdf"), false);
}
