# Part 3 - Yield

In this part you fit the dimuon mass spectrum in each $(p_\mathrm{T},|y|)$ bin and extract $N_{1S}$, $N_{2S}$, and $N_{3S}$.
These are the numerator terms in the cross-section formula.

Tree and bin definitions are already in `yield.C`:
```cpp
static const char* TREE_PATH = "rootuple/mm_tree";
...
static const double Y_EDGES[] = {0.0, 0.6, 1.2, 1.8, 2.4};
...
for (int i = 0; i <= 20; ++i) e.push_back((double)i);
const double tail[] = {20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 50, 55, 60, 70, 100, 130};
```

- binning in $(p_\mathrm{T},|y|)$ matches the final differential measurement,
- a per-bin fit allows direct propagation of statistical uncertainties.

### Fit Model
The fit model is
$$
F(m_{\mu\mu})=N_{1S}S_{1S}(m_{\mu\mu})+N_{2S}S_{2S}(m_{\mu\mu})+N_{3S}S_{3S}(m_{\mu\mu})+N_{\mathrm{bkg}}B(m_{\mu\mu}).
$$

Core implementation excerpt:
```cpp
RooCBShape CB1 ((std::string("CB1")+suf).c_str(),  "CB1",  xvar, m1, sigma1,  alpha1, n1);
RooCBShape CB12((std::string("CB12")+suf).c_str(), "CB12", xvar, m1, sigma12, alpha2, n2);
RooAddPdf douCB1((std::string("douCB1")+suf).c_str(), "douCB1",
                 RooArgList(CB1, CB12), RooArgList(fracDCB));
...
RooFormulaVar m2((std::string("m2")+suf).c_str(),
                 "@0 - 9.4604 + 10.0234", RooArgList(m1));
RooFormulaVar m3((std::string("m3")+suf).c_str(),
                 "@0 - 9.4604 + 10.3501", RooArgList(m1));
...
RooGenericPdf bkg(
  (std::string("bkg")+suf).c_str(), "bkg",
  "exp(@0*(@1-@2))*(1+@3*(@1-@2)*(@1-@2))",
  RooArgList(a, xvar, x0, bfactor)
);
```

Physical meaning:
- signal shapes represent the three $\Upsilon$ resonances,
- parameter links between $1S/2S/3S$ encode known mass-spacing behavior,
- background absorbs non-resonant dimuon contributions under the peaks.

### Extended Likelihood
In yield fit, we used
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
```bash
cd /path/to/CMSDAS_Upsilon/yield
root -l yield.C
```

Outputs:
- `/path/to/CMSDAS_Upsilon/yield/results/2025G/yields.csv`
- `/path/to/CMSDAS_Upsilon/yield/results/2025G/results_ext.csv`
- per-bin PDF files, for example `fit_pt_20-22_y_0p0-0p6.pdf`

What you obtain from this step:
- fitted values for $N_{1S}$, $N_{2S}$, $N_{3S}$,
- fit diagnostics used later for quality control and uncertainty handling.

### Fit Quality
Use both fit PDFs and `results_ext.csv` to check:
- stability of peak positions,
- stability of widths,
- sideband/background behavior.

> #### **Checkpoint**
> Before moving on, mark bins with very low statistics or unstable covariance quality.

> #### **Question**
> 1. Which bins are most likely statistics-limited, and why?
> 2. How is this reflected in both PDF shapes and `results_ext.csv`, and does this match your expectation?
