# Part 3 - Yield

In this part you fit the dimuon mass spectrum in each $(p_\mathrm{T},|y|)$ bin and extract $N_{1S}$, $N_{2S}$, and $N_{3S}$.
These are the numerator terms in the cross-section formula.

Bin definitions are in `yield.C`:
```cpp
static const double Y_EDGES[] = {0.0, 0.6, 1.2, 1.8, 2.4};
...
for (int i = 0; i <= 20; ++i) e.push_back((double)i);
const double tail[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 50, 55, 60, 70, 100, 130};
```

- binning in $(p_\mathrm{T},|y|)$ matches the final differential measurement,
- a per-bin fit allows direct propagation of statistical uncertainties.

### Fit Model
RooFit is a ROOT-based toolkit designed for fitting. It is widely used in high-energy physics to build probability density functions (PDFs), perform maximum-likelihood fits, and extract signal yields from data. In this part, we will use RooFit for modeling invariant-mass distributions and separating signal from background.

The fit model in this analysis is:

$F(m_{\mu\mu})=N_{1S}S_{1S}(m_{\mu\mu})+N_{2S}S_{2S}(m_{\mu\mu})+N_{3S}S_{3S}(m_{\mu\mu})+N_{\mathrm{bkg}}B(m_{\mu\mu})$.

Core implementation excerpt:
```cpp
RooCBShape CB1 ((std::string("CB1")+suf).c_str(),  "CB1",  xvar, m1, sigma1,  alpha1, n1);
RooCBShape CB12((std::string("CB12")+suf).c_str(), "CB12", xvar, m1, sigma12, alpha2, n2);
RooAddPdf douCB1((std::string("douCB1")+suf).c_str(), "douCB1",
                 RooArgList(CB1, CB12), RooArgList(fracDCB));
...
RooAddPdf douCB2((std::string("douCB2")+suf).c_str(), "douCB2",
                 RooArgList(CB2, CB22), RooArgList(fracDCB));
...
RooAddPdf douCB3((std::string("douCB3")+suf).c_str(), "douCB3",
                 RooArgList(CB3, CB32), RooArgList(fracDCB));
...
RooGenericPdf bkg(
  (std::string("bkg")+suf).c_str(), "bkg",
  "exp(@0*(@1-@2))*(1+@3*(@1-@2)*(@1-@2))",
  RooArgList(a, xvar, x0, bfactor)
);
...
RooAddPdf model((std::string("model")+suf).c_str(), "model",
  RooArgList(douCB1, douCB2, douCB3, bkg),
  RooArgList(N1, N2, N3, Nbkg));
```

Physical meaning:
- signal PDFs (double crystal-ball functions) shapes represent the three $\Upsilon$ resonances,
- background PDF (2nd order polynomial) absorbs non-resonant dimuon contributions under the peaks.

### Extended Likelihood
In the yield fit, we use
```cpp
RooFitResult* fr = model.fitTo(
  datahist,
  RooFit::Extended(true),
  RooFit::Save(true),
...
);
```

- `Extended(true)` means normalization parameters (`N1,N2,N3,Nbkg`) enter the Poisson counting term and can be interpreted directly as event yields.
- `Extended(false)` uses shape information only; normalization/yield interpretation is weaker.

## Yield Fit
**Note: DO NOT MODIFY the `INPUT_FILES` path in `yield.C`. The existing path can be accessed.**
```bash
cd /path/to/CMSDAS_Upsilon/yield
root -l yield.C
```

Outputs:
- `/path/to/CMSDAS_Upsilon/yield/results/2025G/yields.csv`
- `/path/to/CMSDAS_Upsilon/yield/results/2025G/results_ext.csv`
- per-bin `.pdf` files, for example `fit_pt_20-22_y_0p0-0p6.pdf`

What you obtain from this step:
- fitted values for $N_{1S}$, $N_{2S}$, $N_{3S}$,
- fit diagnostics used later for quality control and uncertainty handling.

### Fit Quality
Use both fit `.pdf` and `results_ext.csv` to check:
- stability of peak positions,
- stability of widths,
- sideband/background behavior.

> #### **Checkpoint**
> Before moving on, mark bins with very low statistics or unstable covariance quality.

> #### **Question**
> 1. Which bins are most likely statistics-limited, and why?
> 2. How is this reflected in both `.pdf` shapes and `results_ext.csv`, and does this match your expectation?
> 3. Why did we choose Crystall-Ball as the signal PDF? Why do we use double Crystall-Ball functions for each $\Upsilon$ state?
> #### **Task**
> Try modifying `yield.C` by yourself. Change the signal PDF used for fitting each $\Upsilon$ state to:
> 1. single Crystall-Ball function,
> 2. double-Gaussian function.
> Re-examine the fitting results (`.pdf` and `.csv`) to determine whether the fitting function fully describes the shape of the signal. Use these results to further elaborate on your answer to previous Question 3.

> #### **Question**
> 4. Why can we only observe the three signal peaks of $\Upsilon(\mathrm{nS})$ ($n=1,2,3$)? Is it because the mass of $\Upsilon(\mathrm{4S})$ is too large? 
> 5. Using PDG Live (https://pdglive.lbl.gov), search the mass and decay modes of $\Upsilon(\mathrm{4S})$. What is the branching ratio for $\Upsilon(\mathrm{4S})\to\mu^+\mu^-$? Why is it so small compared to $\Upsilon(\mathrm{nS})$ ($n=1,2,3$)?
