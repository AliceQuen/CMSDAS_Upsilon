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

void acceptance()
{
  // Config
  const char* kTreeName   = "rootuple/oniaTree";
  const char* kInputRoot  = "/afs/ihep.ac.cn/users/y/yiyangzhao/Research/CMS_THU_Space/CMSDAS/mc/merged/ParticleGun.root";
  const char* kOutDir     = "results";
  const char* kOutCSV     = "results/acceptance.csv";
  const char* kOutPDF     = "results/acceptance.pdf";

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

  TLorentzVector *gen_dimuon_p4 = nullptr, *gen_muonP_p4 = nullptr, *gen_muonN_p4 = nullptr;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("gen_dimuon_p4", 1);
  tree->SetBranchStatus("gen_muonP_p4", 1);
  tree->SetBranchStatus("gen_muonN_p4", 1);

  tree->SetBranchAddress("gen_dimuon_p4", &gen_dimuon_p4);
  tree->SetBranchAddress("gen_muonP_p4", &gen_muonP_p4);
  tree->SetBranchAddress("gen_muonN_p4", &gen_muonN_p4);

  double All[nYBins][nPtBins] = {{0}};
  double Passed[nYBins][nPtBins] = {{0}};

  TH2D* hAcc = new TH2D("hAcc", "Acceptance; #it{p}_{T} (GeV); |y|",
                       nPtBins, PtBins, nYBins, YBins);
  TH2D* hAccErr = new TH2D("hAccErr", "Acceptance uncertainty; p_{T} (GeV); |y|",
                          nPtBins, PtBins, nYBins, YBins);

  const Long64_t nEvents = tree->GetEntries();

  // Event loop
  for (Long64_t j = 0; j < nEvents; ++j)
  {
    tree->GetEntry(j);

    if (j % 1000000LL == 0) cout << "\rProcessing event " << j << " / " << nEvents << " (" << (100.0 * j / nEvents) << "%)" << std::flush;

    if (!gen_dimuon_p4 || !gen_muonP_p4 || !gen_muonN_p4) continue;

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

    All[iY][iPt] += 1.0;

    // Acceptance cuts
    if (std::abs(gen_muonP_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (std::abs(gen_muonN_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (gen_muonP_p4->Pt() <= kMuonPtMin) continue;
    if (gen_muonN_p4->Pt() <= kMuonPtMin) continue;

    Passed[iY][iPt] += 1.0;
  }

  // Fill histograms + write CSV
  ofstream fout(kOutCSV);
  fout << "pt_min,pt_max,y_abs_min,y_abs_max,acceptance,acceptance_error\n";

  for (int iY = 0; iY < nYBins; ++iY)
  {
    for (int iPt = 0; iPt < nPtBins; ++iPt)
    {
      const double all = All[iY][iPt];
      const double passed = Passed[iY][iPt];

      double acc = 0.0;
      double acc_err = 0.0;

      if (all > 0.0)
      {
        acc = passed / all;

        double low = TEfficiency::ClopperPearson((int)all, (int)passed, 0.682689492, false);
        double up  = TEfficiency::ClopperPearson((int)all, (int)passed, 0.682689492, true);
        acc_err = std::max(acc - low, up - acc);
      }

      hAcc->SetBinContent(iPt + 1, iY + 1, acc);
      hAccErr->SetBinContent(iPt + 1, iY + 1, acc_err);

      fout << PtBins[iPt] << "," << PtBins[iPt + 1] << ","
           << YBins[iY] << "," << YBins[iY + 1] << ","
           << acc << "," << acc_err << "\n";
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

  TCanvas* c = new TCanvas("c", "Acceptance", 1000, 500); // 2:1
  c->SetRightMargin(0.14);
  c->SetLeftMargin(0.09);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.05);

  hAcc->SetMinimum(0.0);
  hAcc->SetMaximum(1.0);
  hAcc->SetTitle("");
  hAcc->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
  hAcc->GetYaxis()->SetTitle("|y|");
  hAcc->GetZaxis()->SetTitle("A");
  hAcc->GetXaxis()->SetTitleSize(0.04);
  hAcc->GetYaxis()->SetTitleSize(0.04);
  hAcc->GetZaxis()->SetTitleSize(0.04);
  hAcc->GetYaxis()->SetNdivisions(nYBins, false);

  hAcc->Draw("COLZ");

  c->SaveAs(kOutPDF);

  cout << "Wrote CSV: " << kOutCSV << "\n";
  cout << "Wrote PDF: " << kOutPDF << "\n";
}