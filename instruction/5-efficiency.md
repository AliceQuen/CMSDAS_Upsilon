# Part 5 - Efficiency

This part computes efficiency $\epsilon(p_\mathrm{T},|y|)$, the probability that events inside the fiducial region are reconstructed and pass analysis selections.

Why efficiency is needed:
- similiar to acceptance, efficiency quantifies the probability that fiducial events survive the full analysis chain.
- the final selection also improves signal-to-background ratio, which stabilizes the mass fits and helps reduce statistical/fit uncertainty.

## Efficiency Definition
$$
\epsilon(p_\mathrm{T},|y|)=\frac{N^{\mathrm{sel}}_{\mathrm{reco}}(p_\mathrm{T},|y|)}{N^{\mathrm{fid}}_{\mathrm{gen}}(p_\mathrm{T},|y|)}.
$$
- denominator: generated fiducial events,
- numerator: reconstructed events that pass final analysis requirements.

Physical interpretation of the numerator cuts:
- kinematic cuts ensure consistency with fiducial definitions,
- `trigger` models online selection survival,
- `vProb` and candidate/charge requirements model analysis-level quality and topology constraints.

## Efficiency Calculation
This part uses a full-chain MC (generation, detector simulation, particle reconstruction, etc.), which hopes to match the process of collecting data in reality. In `mc_efficiency.C`, core logics are:
```cpp
if (std::abs(gen_muonP_p4->Eta()) > kMuonAbsEtaMax) continue;
if (std::abs(gen_muonM_p4->Eta()) > kMuonAbsEtaMax) continue;
if (gen_muonP_p4->Pt() <= kMuonPtMin) continue;
if (gen_muonM_p4->Pt() <= kMuonPtMin) continue;
All[iY][iPt] += 1.0;
...
if (std::abs(muonP_p4->Eta()) > kMuonAbsEtaMax) continue;
if (std::abs(muonM_p4->Eta()) > kMuonAbsEtaMax) continue;
if (muonP_p4->Pt() <= kMuonPtMin) continue;
if (muonM_p4->Pt() <= kMuonPtMin) continue;
if (nonia != 1 || !trigger || charge != 0 || vProb <= 0.01) continue;
Passed[iY][iPt] += 1.0;
```

```bash
cd /path/to/CMSDAS_Upsilon/efficiency/mc_efficiency
root -l mc_efficiency.C
```

Outputs are efficiency maps in $(p_\mathrm{T},|y|)$:
- `/path/to/CMSDAS_Upsilon/efficiency/mc_efficiency/results/efficiency.csv`
- `/path/to/CMSDAS_Upsilon/efficiency/mc_efficiency/results/efficiency.pdf`


> #### **Question**
> In low-$p_\mathrm{T}$ (high-$p_\mathrm{T}$) bins, which numerator requirement is most likely to dominate the efficiency loss (`trigger`, `vProb`, or single-candidate condition), and why?
> #### **Task**
> Verify your assumption by modifying the program to calculate the efficiency of each selection separately.
