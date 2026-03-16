# Part 7 - Results

This part combines $N$, $A$, $\epsilon$, and $\mathcal{L}$ to produce the differential cross sections and final figures.

## Python Setup
**Note: Before continuing, please open a new terminal or reconnect to `lxplus8`.**

```bash
ssh <cern_username>@lxplus8.cern.ch

cd /path/to/CMSDAS_Upsilon/data/CMSSW_15_0_18/src
cmsenv

python3 -c "import sys; print(sys.version)"
python3 -c "import pandas,numpy,matplotlib,mplhep; print('python packages ok')"
```

> #### **Checkpoint**
> Continue only after both `python3 -c` tests complete successfully.

## Luminosity
Luminosity, $\mathcal{L}$, is the absolute normalization of the measurement. Any luminosity bias directly rescales the extracted cross section. We will calculate the integrated luminosity from certified luminosity mask selection, `Cert_Collisions2025_391658_398903_Muon.json`, which was previously employed to generate the data sample.

`brilcalc` is a standard CMS tool for calculating integrated luminosity from certified luminosity mask selection (`Cert_Collisions2025_391658_398903_Muon.json`). In CMS workflows, `brilcalc` is the official and widely recommended tool for luminosity evaluation, and the calculation is typically performed with the appropriate normtag and stable-beam selection. 
- [Luminosity short exercise in CMS DAS 2024](https://cmsdas-cern-2024-short-ex-lumi.docs.cern.ch)

In this exercise, we will use `luminosity.py` basing on `brilcalc` to calculate the integrated luminosity.

```bash
cd /path/to/CMSDAS_Upsilon/luminosity

python3 luminosity.py
```

Outputs:
- `/path/to/CMSDAS_Upsilon/luminosity/results/2025G.csv`

## Cross Section

$\mathcal{B}(\Upsilon(\mathrm{nS})\to\mu^+\mu^-)\,\frac{d^2\sigma_n}{dp_\mathrm{T}\,d|y|}
=\frac{N_n}{\mathcal{L}\,A_n\,\epsilon_n\,\Delta p_\mathrm{T}\,2\Delta|y|}$.

**Note**: rapidity width is defined as $2\Delta|y|=\Delta y$.

```bash
cd /path/to/CMSDAS_Upsilon/cross_section
python3 cross_section.py
```

Output:
- `/path/to/CMSDAS_Upsilon/cross_section/results/cross_section.csv`

> #### **Question**
> 1. The same program also reports the integrated cross sections of the three $\Upsilon$ states in the region $|y|<2.4$ and $p_\mathrm{T}<130$ GeV; if your final target were only this integrated value, could you measure it with one very wide bin covering the full range?
> 2. What problems can appear when such a wide bin is used (for example: acceptance/efficiency variation inside the bin, model dependence, and loss of kinematic information)? Beyond resolving trends versus $p_\mathrm{T}$ and $y$, this is another advantage of differential measurements.

## Plots
```bash
python3 plot.py
python3 plot_comprison_13p6TeV.py
```

Outputs:
- `/path/to/CMSDAS_Upsilon/cross_section/results/cross_section.pdf`
- `/path/to/CMSDAS_Upsilon/cross_section/results/2025vs2022.pdf`

Plot meaning:
- `plot.py` provides absolute differential cross-section distributions,
- `plot_comprison_13p6TeV.py` provides 2025/2022 ratios for consistency and evolution checks.

> #### **Question**
> 3. After inspecting these `.pdf`, what are the observed differential cross-section trends versus $p_\mathrm{T}$ and rapidity?
> 4. Why do these trends appear, and are they consistent with your physical expectation?
