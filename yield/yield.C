#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TColor.h>
#include <TGaxis.h>

#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooFormulaVar.h>
#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooGenericPdf.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

static const std::vector<std::string> INPUT_FILES = {
  "/eos/home-y/yiyangz/public/CMSDAS/data/2025G_Parking*"
};

static const char* TREE_PATH = "rootuple/mm_tree";

static const char* OUT_DIR   = "./results/2025G";
static const char* CSV_BASE  = "./results/2025G/yields.csv";
static const char* CSV_EXT   = "./results/2025G/results_ext.csv";

static const double M_MIN = 8.7;
static const double M_MAX = 11.5;
static const double M_MIN_HIPT = 8.1;
static const double M_MAX_HIPT = 11.9;
static const double M_MIN_GLOBAL = 8.1;
static const double M_MAX_GLOBAL = 11.9;
static const double PT_SKIP_GE = 130.0;
static const double YABS_MAX   = 2.4;
static const int MIN_COUNTS_TO_FIT = 10;
static const int NEDGE_BINS_SLOPE = 6;
static const double kMuonAbsEtaMax = 2.0;
static const double kMuonPtMin     = 3.1;

static const double Y_EDGES[] = {0.0, 0.6, 1.2, 1.8, 2.4};
static const int NY = (int)(sizeof(Y_EDGES)/sizeof(Y_EDGES[0])) - 1;

static std::vector<double> makePtEdges() {
  std::vector<double> e;
  for (int i = 0; i <= 20; ++i) e.push_back((double)i);
  const double tail[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 50, 55, 60, 70, 100, 130};
  for (double v : tail) if (v > e.back() + 1e-12) e.push_back(v);
  return e;
}

static int findBin(double x, const std::vector<double>& edges) {
  int n = (int)edges.size() - 1;
  if (n <= 0) return -1;
  if (x < edges.front()) return -1;
  if (x >= edges.back()) return -1;
  for (int i = 0; i < n; ++i) if (x >= edges[i] && x < edges[i+1]) return i;
  return -1;
}

static int findBinY(double yabs) {
  if (yabs < Y_EDGES[0]) return -1;
  if (yabs >= Y_EDGES[NY]) return -1;
  for (int i = 0; i < NY; ++i) if (yabs >= Y_EDGES[i] && yabs < Y_EDGES[i+1]) return i;
  return -1;
}

static double massBinWidthForPtBin(double pt_hi) {
  if (pt_hi <= 20.0 + 1e-12)  return 0.040;
  if (pt_hi <= 40.0 + 1e-12)  return 0.050;
  if (pt_hi <= 100.0 + 1e-12) return 0.080;
  return 0.100;
}

static void ensureDir(const char* path) {
  gSystem->mkdir(path, true);
}

static std::string fmtNum(double x, int ndp) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.*f", ndp, x);
  std::string s(buf);
  for (char& ch : s) {
    if (ch == '.') ch = 'p';
    if (ch == '-') ch = 'm';
  }
  return s;
}

static double meanEdgeCounts(const TH1D* h, bool left, int nEdge) {
  const int nb = h->GetNbinsX();
  const int nUse = std::max(1, std::min(nEdge, nb));
  double sum = 0.0;
  if (left) {
    for (int i = 1; i <= nUse; ++i) sum += h->GetBinContent(i);
  } else {
    for (int i = nb - nUse + 1; i <= nb; ++i) sum += h->GetBinContent(i);
  }
  return sum / (double)nUse;
}

static void massWindowForPtLo(double pt_lo, double& mmin, double& mmax) {
  if (pt_lo > 5.0) {  // requirement: pTlow > 40
    mmin = M_MIN_HIPT;
    mmax = M_MAX_HIPT;
  } else {
    mmin = M_MIN;
    mmax = M_MAX;
  }
}

static int nbinsFromBW(double bw, double mmin, double mmax) {
  return (int)std::lround((mmax - mmin) / bw);
}

int yield() {
  ensureDir(OUT_DIR);
  
  TGaxis::SetMaxDigits(4);

  std::ofstream fcsv(CSV_BASE);
  std::ofstream fcsv_ext(CSV_EXT);
  if (!fcsv.good() || !fcsv_ext.good()) {
    std::cerr << "ERROR: cannot open CSV outputs.\n";
    return 1;
  }

  fcsv
    << "pt_min,pt_max,y_abs_min,y_abs_max,"
    << "N_1S,N_1S_err,N_2S,N_2S_err,N_3S,N_3S_err,N_bkg,N_bkg_err\n";

  fcsv_ext
    << "pt_min,pt_max,y_abs_min,y_abs_max,"
    << "N_1S,N_1S_err,N_2S,N_2S_err,N_3S,N_3S_err,N_bkg,N_bkg_err,"
    << "fit_status,fit_covQual,fit_edm,"
    << "alpha1,alpha1_err,alpha2,alpha2_err,n1,n1_err,n2,n2_err,f,f_err,"
    << "m1,m1_err,sigma1,sigma1_err,a,a_err,bfactor,bfactor_err\n";

  TChain chain(TREE_PATH);
  for (const auto& f : INPUT_FILES) chain.Add(f.c_str());
  const Long64_t nEntries = chain.GetEntries();
  std::cout << "Total entries: " << nEntries << "\n";

  TLorentzVector *dimuon_p4 = nullptr, *muonP_p4 = nullptr, *muonM_p4 = nullptr;
  UInt_t nonia, trigger;
  Int_t charge;
  Float_t vProb;

  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("dimuon_p4", 1);
  chain.SetBranchStatus("muonP_p4", 1);
  chain.SetBranchStatus("muonM_p4", 1);
  chain.SetBranchStatus("nonia", 1);
  chain.SetBranchStatus("trigger", 1);
  chain.SetBranchStatus("charge", 1);
  chain.SetBranchStatus("vProb", 1);

  chain.SetBranchAddress("dimuon_p4", &dimuon_p4);
  chain.SetBranchAddress("muonP_p4", &muonP_p4);
  chain.SetBranchAddress("muonM_p4", &muonM_p4);
  chain.SetBranchAddress("nonia", &nonia);
  chain.SetBranchAddress("trigger", &trigger);
  chain.SetBranchAddress("charge", &charge);
  chain.SetBranchAddress("vProb", &vProb);

  std::vector<double> ptEdges = makePtEdges();
  const int NPT = (int)ptEdges.size() - 1;

  std::vector<std::vector<TH1D*>> hists(NY, std::vector<TH1D*>(NPT, nullptr));
  std::vector<std::vector<double>> bwMass(NY, std::vector<double>(NPT, 0.0));

  for (int iy = 0; iy < NY; ++iy) {
    for (int ipt = 0; ipt < NPT; ++ipt) {
      const double pt_lo = ptEdges[ipt];
      const double pt_hi = ptEdges[ipt+1];

      const double bw = massBinWidthForPtBin(pt_hi);
      bwMass[iy][ipt] = bw;

      double h_mmin, h_mmax;
      massWindowForPtLo(pt_lo, h_mmin, h_mmax);
      const int nb = nbinsFromBW(bw, h_mmin, h_mmax);

      char name[128];
      std::snprintf(name, sizeof(name), "h_m_iy%d_ipt%d", iy, ipt);
      hists[iy][ipt] = new TH1D(name, name, nb, h_mmin, h_mmax);
      hists[iy][ipt]->Sumw2(false);
    }
  }

  Long64_t filled = 0;
  for (Long64_t i = 0; i < nEntries; ++i) {
    
    chain.GetEntry(i);

    if (!dimuon_p4 || !muonP_p4 || !muonM_p4) continue;
    if (std::abs(muonP_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (std::abs(muonM_p4->Eta()) > kMuonAbsEtaMax) continue;
    if (muonP_p4->Pt() <= kMuonPtMin) continue;
    if (muonM_p4->Pt() <= kMuonPtMin) continue;
    if (nonia != 1 || !trigger || charge != 0 || vProb <= 0.01) continue;

    if (i % 1000000 == 0) {
      const double pct = (nEntries > 0) ? (100.0 * (double)i / (double)nEntries) : 0.0;
      std::cout << "Processing entry " << i << " / " << nEntries
                << " (" << pct << "%)\r" << std::flush;
    }

    const double pt = dimuon_p4->Pt();
    if (pt >= PT_SKIP_GE) continue;

    const double ya = std::fabs(dimuon_p4->Rapidity());
    if (ya >= YABS_MAX) continue;

    const double m = dimuon_p4->M();
    if (m < M_MIN_GLOBAL || m >= M_MAX_GLOBAL) continue;

    const int iy = findBinY(ya);
    const int ipt = findBin(pt, ptEdges);
    if (iy < 0 || ipt < 0) continue;

    hists[iy][ipt]->Fill(m);
    ++filled;
  }
  std::cout << "\nFilled entries: " << filled << "\n";

  const int col_bkg   = TColor::GetColor("#9c9ca1");
  const int col_3S    = TColor::GetColor("#5790fc");
  const int col_2S    = TColor::GetColor("#3a981d");
  const int col_1S    = TColor::GetColor("#e42536");
  const int col_total = TColor::GetColor("#717581");

  for (int iy = 0; iy < NY; ++iy) {
    for (int ipt = 0; ipt < NPT; ++ipt) {
      TH1D* h = hists[iy][ipt];
      if (!h) continue;

      const double Ntot = h->Integral();
      if (Ntot < MIN_COUNTS_TO_FIT) continue;

      const double y_lo  = Y_EDGES[iy];
      const double y_hi  = Y_EDGES[iy+1];
      const double pt_lo = ptEdges[ipt];
      const double pt_hi = ptEdges[ipt+1];
      const double bw    = bwMass[iy][ipt];

      double mmin_bin, mmax_bin;
      massWindowForPtLo(pt_lo, mmin_bin, mmax_bin);

      char suf[64];
      std::snprintf(suf, sizeof(suf), "_iy%d_ipt%d", iy, ipt);

      RooRealVar xvar((std::string("xvar")+suf).c_str(), "m(#mu#mu)", mmin_bin, mmax_bin);
      RooDataHist datahist((std::string("datahist")+suf).c_str(),
                           (std::string("datahist")+suf).c_str(),
                           RooArgList(xvar), h);

      RooRealVar alpha1((std::string("alpha1")+suf).c_str(), "alpha1", 2);
      RooRealVar alpha2((std::string("alpha2")+suf).c_str(), "alpha2", 2);
      RooRealVar n1    ((std::string("n1")+suf).c_str(),     "n1",  1.0);
      RooRealVar n2    ((std::string("n2")+suf).c_str(),     "n2",  2.0);
      RooRealVar fracDCB((std::string("f")+suf).c_str(),     "f",   0.4);

      RooRealVar m1    ((std::string("m1")+suf).c_str(),     "m1",  9.4604, 9.420, 9.485);
      RooRealVar sigma1((std::string("sigma1")+suf).c_str(), "sigma1", 0.052, 0.045, 0.15);

      RooFormulaVar sigma12((std::string("sigma12")+suf).c_str(),
                            "@0*1.55", RooArgList(sigma1));

      RooCBShape CB1 ((std::string("CB1")+suf).c_str(),  "CB1",  xvar, m1, sigma1,  alpha1, n1);
      RooCBShape CB12((std::string("CB12")+suf).c_str(), "CB12", xvar, m1, sigma12, alpha2, n2);
      RooAddPdf douCB1((std::string("douCB1")+suf).c_str(), "douCB1",
                       RooArgList(CB1, CB12), RooArgList(fracDCB));

      RooFormulaVar m2((std::string("m2")+suf).c_str(),
                       "@0 - 9.4604 + 10.0234", RooArgList(m1));
      RooFormulaVar sigma2((std::string("sigma2")+suf).c_str(),
                           "@0/@1*@2", RooArgList(sigma1, m1, m2));
      RooFormulaVar sigma22((std::string("sigma22")+suf).c_str(),
                            "@0/@1*@2", RooArgList(sigma12, m1, m2));

      RooCBShape CB2 ((std::string("CB2")+suf).c_str(),  "CB2",  xvar, m2, sigma2,  alpha1, n1);
      RooCBShape CB22((std::string("CB22")+suf).c_str(), "CB22", xvar, m2, sigma22, alpha2, n2);
      RooAddPdf douCB2((std::string("douCB2")+suf).c_str(), "douCB2",
                       RooArgList(CB2, CB22), RooArgList(fracDCB));

      RooFormulaVar m3((std::string("m3")+suf).c_str(),
                       "@0 - 9.4604 + 10.3501", RooArgList(m1));
      RooFormulaVar sigma3((std::string("sigma3")+suf).c_str(),
                           "@0/@1*@2", RooArgList(sigma1, m1, m3));
      RooFormulaVar sigma32((std::string("sigma32")+suf).c_str(),
                            "@0/@1*@2", RooArgList(sigma12, m1, m3));

      RooCBShape CB3 ((std::string("CB3")+suf).c_str(),  "CB3",  xvar, m3, sigma3,  alpha1, n1);
      RooCBShape CB32((std::string("CB32")+suf).c_str(), "CB32", xvar, m3, sigma32, alpha2, n2);
      RooAddPdf douCB3((std::string("douCB3")+suf).c_str(), "douCB3",
                       RooArgList(CB3, CB32), RooArgList(fracDCB));

      const double x0_val = 0.5*(mmin_bin + mmax_bin);
      RooConstVar x0((std::string("x0")+suf).c_str(), "x0", x0_val);

      const double yLseed = meanEdgeCounts(h, true,  NEDGE_BINS_SLOPE);
      const double yRseed = meanEdgeCounts(h, false, NEDGE_BINS_SLOPE);
      double r = 0.0;
      if ((yLseed + yRseed) > 0.0) r = (yRseed - yLseed) / (yRseed + yLseed);

      double beta_lo = -1.5;
      double beta_hi =  1.5;

      double c_lo    = -0.3;
      double c_hi    = 1e-4;

      if (pt_lo >= 40.0) {
        beta_lo = -0.2;
        beta_hi =  0.2;
        c_lo    = -0.05;
      }

      double beta_init = 0.8 * r;
      if (beta_init < beta_lo) beta_init = beta_lo;
      if (beta_init > beta_hi) beta_init = beta_hi;

      double c_init = 0.15 * std::fabs(r);
      if (c_init < c_lo) c_init = c_lo;
      if (c_init > c_hi) c_init = c_hi;

      RooRealVar a((std::string("a")+suf).c_str(), "a", beta_init, beta_lo, beta_hi);
      RooRealVar bfactor((std::string("bfactor")+suf).c_str(), "bfactor", c_init, c_lo, c_hi);

      RooGenericPdf bkg(
        (std::string("bkg")+suf).c_str(), "bkg",
        "exp(@0*(@1-@2))*(1+@3*(@1-@2)*(@1-@2))",
        RooArgList(a, xvar, x0, bfactor)
      );

      RooRealVar N1((std::string("N1")+suf).c_str(), "N1", 0.3*Ntot, 0.0, 0.8*Ntot);
      RooRealVar N2((std::string("N2")+suf).c_str(), "N2", 0.20*Ntot, 0.0, 0.5*Ntot);
      RooRealVar N3((std::string("N3")+suf).c_str(), "N3", 0.15*Ntot, 0.0, 0.3*Ntot);
      RooRealVar Nbkg((std::string("Nbkg")+suf).c_str(), "Nbkg", 0.5*Ntot, 0.0, 1.0*Ntot);

      RooAddPdf model((std::string("model")+suf).c_str(), "model",
                      RooArgList(douCB1, douCB2, douCB3, bkg),
                      RooArgList(N1, N2, N3, Nbkg));

      RooFitResult* fr = model.fitTo(
        datahist,
        RooFit::Extended(true),
        RooFit::Save(true),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(false),
        RooFit::Verbose(false)
      );

      const double N1v = N1.getVal(),   N1e = N1.getError();
      const double N2v = N2.getVal(),   N2e = N2.getError();
      const double N3v = N3.getVal(),   N3e = N3.getError();
      const double Nbv = Nbkg.getVal(), Nbe = Nbkg.getError();

      fcsv
        << pt_lo << "," << pt_hi << "," << y_lo << "," << y_hi << ","
        << N1v << "," << N1e << ","
        << N2v << "," << N2e << ","
        << N3v << "," << N3e << ","
        << Nbv << "," << Nbe << "\n";

      int fit_status = fr ? fr->status() : -999;
      int fit_covQual = fr ? fr->covQual() : -999;
      double fit_edm = fr ? fr->edm() : -1.0;

      fcsv_ext
        << pt_lo << "," << pt_hi << "," << y_lo << "," << y_hi << ","
        << N1v << "," << N1e << ","
        << N2v << "," << N2e << ","
        << N3v << "," << N3e << ","
        << Nbv << "," << Nbe << ","
        << fit_status << "," << fit_covQual << "," << fit_edm << ","
        << alpha1.getVal() << "," << alpha1.getError() << ","
        << alpha2.getVal() << "," << alpha2.getError() << ","
        << n1.getVal() << "," << n1.getError() << ","
        << n2.getVal() << "," << n2.getError() << ","
        << fracDCB.getVal() << "," << fracDCB.getError() << ","
        << m1.getVal() << "," << m1.getError() << ","
        << sigma1.getVal() << "," << sigma1.getError() << ","
        << a.getVal() << "," << a.getError() << ","
        << bfactor.getVal() << "," << bfactor.getError()
        << "\n";

      RooPlot* frame = xvar.frame();
      frame->SetTitle("");

      datahist.plotOn(
        frame,
        RooFit::Name("data"),
        RooFit::DrawOption("PE"),
        RooFit::MarkerStyle(20),
        RooFit::MarkerSize(1.1),
        RooFit::LineColor(kBlack),
        RooFit::MarkerColor(kBlack)
      );

      model.plotOn(frame, RooFit::Name("total"),
                   RooFit::LineColor(col_total),
                   RooFit::LineStyle(kSolid),
                   RooFit::LineWidth(3));

      model.plotOn(frame, RooFit::Name("bkg"),
                   RooFit::Components(bkg),
                   RooFit::LineColor(col_bkg),
                   RooFit::LineStyle(kDashed),
                   RooFit::LineWidth(3));

      model.plotOn(frame, RooFit::Name("y3"),
                   RooFit::Components(douCB3),
                   RooFit::LineColor(col_3S),
                   RooFit::LineStyle(kSolid),
                   RooFit::LineWidth(3));

      model.plotOn(frame, RooFit::Name("y2"),
                   RooFit::Components(douCB2),
                   RooFit::LineColor(col_2S),
                   RooFit::LineStyle(kSolid),
                   RooFit::LineWidth(3));

      model.plotOn(frame, RooFit::Name("y1"),
                   RooFit::Components(douCB1),
                   RooFit::LineColor(col_1S),
                   RooFit::LineStyle(kSolid),
                   RooFit::LineWidth(3));

      TCanvas c("c", "c", 900, 900);
      c.SetLeftMargin(0.14);
      c.SetRightMargin(0.04);
      c.SetBottomMargin(0.08);
      c.SetTopMargin(0.04);

      frame->GetXaxis()->SetTitle("m(#it{#mu}#it{#mu}) (GeV)");
      frame->GetYaxis()->SetTitle(Form("Events / %.0f MeV", bw*1000.0));
      frame->GetXaxis()->SetNoExponent(true);   // x axis never shows ×10^n
      frame->GetYaxis()->SetNoExponent(false);  // y axis may show ×10^n

      frame->Draw();

      TLegend leg(0.62, 0.56, 0.91, 0.85);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);

      std::string pt_line;
      if (pt_lo <= 0.0) pt_line = Form("#it{p}_{T} < %.0f GeV", pt_hi);
      else              pt_line = Form("%.0f < #it{p}_{T} < %.0f GeV", pt_lo, pt_hi);

      std::string y_line;
      if (y_lo <= 0.0) y_line = Form("|y| < %.1f", y_hi);
      else             y_line = Form("%.1f < |y| < %.1f", y_lo, y_hi);

      leg.SetHeader(Form("#splitline{#splitline{%s}{%s}}{}", pt_line.c_str(), y_line.c_str()), "C");

      leg.AddEntry(frame->findObject("data"),  "Data", "pe");
      leg.AddEntry(frame->findObject("y1"),    "#Upsilon(1S)", "l");
      leg.AddEntry(frame->findObject("y2"),    "#Upsilon(2S)", "l");
      leg.AddEntry(frame->findObject("y3"),    "#Upsilon(3S)", "l");
      leg.AddEntry(frame->findObject("bkg"),   "Background", "l");
      leg.AddEntry(frame->findObject("total"), "Total", "l");
      leg.Draw();

      std::string outPdf =
            std::string(OUT_DIR) +
            "/fit_pt_" + fmtNum(pt_lo, 0) + "-" + fmtNum(pt_hi, 0) +
            "_y_" + fmtNum(y_lo, 1) + "-" + fmtNum(y_hi, 1) +
            ".pdf";
      c.SaveAs(outPdf.c_str());

      delete frame;
      delete fr;
    }
  }

  for (int iy = 0; iy < NY; ++iy) {
    for (int ipt = 0; ipt < NPT; ++ipt) {
      delete hists[iy][ipt];
      hists[iy][ipt] = nullptr;
    }
  }

  fcsv.close();
  fcsv_ext.close();

  std::cout << "Saved: " << CSV_BASE << "\n";
  std::cout << "Saved: " << CSV_EXT  << "\n";
  std::cout << "PDFs in: " << OUT_DIR << "\n";
  return 0;
}