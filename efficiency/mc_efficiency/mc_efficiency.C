#include <TChain.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <TEfficiency.h>

#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void mc_efficiency()
{
  const char* kTreeName   = "rootuple/mm_tree";
  const char* kInputRoot  = "/afs/ihep.ac.cn/users/y/yiyangzhao/Research/CMS_THU_Space/CMSDAS/mc/merged/Upsilon2MM*";
  const char* kOutDir     = "results";
  const char* kOutCSV     = "results/efficiency.csv";
  const char* kOutPDF     = "results/efficiency.pdf";

  const double kAbsYMax   = 2.4;
  const double kPtMin     = 0.0;
  const double kPtMax     = 130.0;

  const double kMuonAbsEtaMax = 2.0;
  const double kMuonPtMin     = 3.1;

  const int nPtBins = 38;
  const double PtBins[nPtBins + 1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,43,46,50,55,60,70,100,130};
  const int nYBins = 4;
  const double YBins[nYBins + 1] = {0.0, 0.6, 1.2, 1.8, 2.4};

  // I/O setup
  gSystem->mkdir(kOutDir, true);

  TChain* tree = new TChain(kTreeName);
  tree->Add(kInputRoot);

  TLorentzVector *gen_dimuon_p4 = nullptr, *gen_muonP_p4 = nullptr, *gen_muonM_p4 = nullptr, *dimuon_p4 = nullptr, *muonP_p4 = nullptr, *muonM_p4 = nullptr;
  UInt_t nonia, trigger;
  Int_t charge;
  Float_t vProb;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("gen_dimuon_p4", 1);
  tree->SetBranchStatus("gen_muonP_p4", 1);
  tree->SetBranchStatus("gen_muonM_p4", 1);
  tree->SetBranchStatus("dimuon_p4", 1);
  tree->SetBranchStatus("muonP_p4", 1);
  tree->SetBranchStatus("muonM_p4", 1);
  tree->SetBranchStatus("nonia", 1);
  tree->SetBranchStatus("trigger", 1);
  tree->SetBranchStatus("charge", 1);
  tree->SetBranchStatus("vProb", 1);

  tree->SetBranchAddress("gen_dimuon_p4", &gen_dimuon_p4);
  tree->SetBranchAddress("gen_muonP_p4", &gen_muonP_p4);
  tree->SetBranchAddress("gen_muonM_p4", &gen_muonM_p4);
  tree->SetBranchAddress("dimuon_p4", &dimuon_p4);
  tree->SetBranchAddress("muonP_p4", &muonP_p4);
  tree->SetBranchAddress("muonM_p4", &muonM_p4);
  tree->SetBranchAddress("nonia", &nonia);
  tree->SetBranchAddress("trigger", &trigger);
  tree->SetBranchAddress("charge", &charge);
  tree->SetBranchAddress("vProb", &vProb);

  double All[nYBins][nPtBins] = {{0}};
  double Passed[nYBins][nPtBins] = {{0}};

  TH2D* hEff = new TH2D("hEff", "Efficiency; #it{p}_{T} (GeV); |y|", nPtBins, PtBins, nYBins, YBins);
  TH2D* hEffErr = new TH2D("hEffErr", "Efficiency error; p_{T} (GeV); |y|", nPtBins, PtBins, nYBins, YBins);

  const Long64_t nEvents = tree->GetEntries();

  // Event loop
  for (Long64_t j = 0; j < nEvents; ++j)
  {
    tree->GetEntry(j);

    if (j % 1000000LL == 0) cout << "\rProcessing event " << j << " / " << nEvents << " (" << (100.0 * j / nEvents) << "%)" << std::flush;

    if (!gen_dimuon_p4 || !gen_muonP_p4 || !gen_muonM_p4 || !dimuon_p4 || !muonP_p4 || !muonM_p4) continue;

    const double pt = gen_dimuon_p4->Pt();
    const double absy = std::abs(gen_dimuon_p4->Rapidity());

    if (absy >= kAbsYMax) continue;
    if (pt <= kPtMin) continue;
    if (pt >= kPtMax) continue;

    int iPt = -1;
    for (int i = 0; i < nPtBins; ++i)
    {
      if (pt >= PtBins[i] && pt < PtBins[i + 1]) { iPt = i; break; }
    }

    int iY = -1;
    for (int i = 0; i < nYBins; ++i)
    {
      if (absy >= YBins[i] && absy < YBins[i + 1]) { iY = i; break; }
    }

    if (iPt < 0 || iY < 0) continue;

    // Acceptance cuts
    if (std::abs(gen_muonP_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (std::abs(gen_muonM_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (gen_muonP_p4->Pt() <= kMuonPtMin) continue;
    if (gen_muonM_p4->Pt() <= kMuonPtMin) continue;

    All[iY][iPt] += 1.0;

    if (std::abs(muonP_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (std::abs(muonM_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (muonP_p4->Pt() <= kMuonPtMin) continue;
    if (muonM_p4->Pt() <= kMuonPtMin) continue;
    if (nonia != 1 || !trigger || charge != 0 || vProb <= 0.01) continue;

    Passed[iY][iPt] += 1.0;
  }

  // Fill histograms + write CSV
  ofstream fout(kOutCSV);
  fout << "pt_min,pt_max,y_abs_min,y_abs_max,efficiency,efficiency_error\n";

  for (int iY = 0; iY < nYBins; ++iY)
  {
    for (int iPt = 0; iPt < nPtBins; ++iPt)
    {
      const double all = All[iY][iPt];
      const double passed = Passed[iY][iPt];

      double eff = 0.0;
      double eff_err = 0.0;

      if (all > 0.0)
      {
        eff = passed / all;

        double low = TEfficiency::ClopperPearson((int)all, (int)passed, 0.682689492, false);
        double up  = TEfficiency::ClopperPearson((int)all, (int)passed, 0.682689492, true);
        eff_err = std::max(eff - low, up - eff);
      }

      hEff->SetBinContent(iPt + 1, iY + 1, eff);
      hEffErr->SetBinContent(iPt + 1, iY + 1, eff_err);

      fout << PtBins[iPt] << "," << PtBins[iPt + 1] << ","
           << YBins[iY] << "," << YBins[iY + 1] << ","
           << eff << "," << eff_err << "\n";
    }
  }
  fout.close();

  // Plot
  gStyle->SetOptStat(0);

  const int NCont = 255;
  const int NStops = 3;
  double stops[NStops] = {0.00, 0.50, 1.00};
  double red[NStops]   = {0.90, 0.35, 0.05};
  double green[NStops] = {0.95, 0.60, 0.20};
  double blue[NStops]  = {1.00, 0.85, 0.55};
  TColor::CreateGradientColorTable(NStops, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  TCanvas* c = new TCanvas("c", "Efficiency", 1000, 500); // 2:1
  c->SetRightMargin(0.14);
  c->SetLeftMargin(0.09);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.05);

  hEff->SetMinimum(0.0);
  hEff->SetMaximum(1.0);
  hEff->SetTitle("");
  hEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
  hEff->GetYaxis()->SetTitle("|y|");
  hEff->GetZaxis()->SetTitle("#epsilon");
  hEff->GetXaxis()->SetTitleSize(0.04);
  hEff->GetYaxis()->SetTitleSize(0.04);
  hEff->GetZaxis()->SetTitleSize(0.04);
  hEff->GetYaxis()->SetNdivisions(nYBins, false);

  hEff->Draw("COLZ");

  c->SaveAs(kOutPDF);

  cout << "Wrote CSV: " << kOutCSV << "\n";
  cout << "Wrote PDF: " << kOutPDF << "\n";
}