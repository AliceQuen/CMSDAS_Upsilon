# Part 7 - Results

This part combines $N$, $A$, $\epsilon$, and $\mathcal{L}$ to produce differential cross sections results and final figures.

## Python Setup

```bash
cd /path/to/CMSDAS_Upsilon

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el8-gcc11-opt/setup.sh

python3 -c "import sys; print(sys.version)"
python3 -c "import pandas,numpy,matplotlib,mplhep; print('python packages ok')"
```

> #### **Checkpoint**
> Continue only after both `python3 -c` tests complete successfully.

## Luminosity
Luminosity, $\mathcal{L}$, is the absolute normalization of the measurement. Any luminosity bias directly rescales the extracted cross section. We will use `Cert_Collisions2025_391658_398903_Muon.json`, which was previously employed to generate data samples, to compute the integral luminosity of the used dataset after mask selection.

```bash
cd /path/to/CMSDAS_Upsilon/luminosity
python3 luminosity.py
```

Output:
- `/path/to/CMSDAS_Upsilon/luminosity/results/2025G.csv`


## Cross Section

$\mathcal{B}(\Upsilon(nS)\to\mu^+\mu^-)\,\frac{d^2\sigma_n}{dp_\mathrm{T}\,d|y|}
=\frac{N_n}{\mathcal{L}\,A_n\,\epsilon_n\,\Delta p_\mathrm{T}\,2\Delta|y|}.$

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
> 1. After inspecting these PDFs, what are the observed differential cross-section trends versus $p_\mathrm{T}$ and rapidity?
> 2. Why do these trends appear, and are they consistent with your physical expectation?
