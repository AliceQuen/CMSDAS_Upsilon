# Part 5 - Efficiency

This part computes efficiency $\epsilon(p_T,|y|)$, the probability that events inside the fiducial region are reconstructed and pass analysis selections.

Why efficiency is needed:
- similiar to acceptance, efficiency quantifies the probability that fiducial events survive the full analysis chain.
- the final selection also improves signal-to-background ratio, which stabilizes the mass fits and helps reduce statistical/fit uncertainty.

## Efficiency Definition
$$
\epsilon(p_T,|y|)=\frac{N^{\mathrm{sel}}_{\mathrm{reco}}(p_T,|y|)}{N^{\mathrm{fid}}_{\mathrm{gen}}(p_T,|y|)}.
$$
- denominator: generated fiducial events,
- numerator: reconstructed events that pass final analysis requirements.

Physical interpretation of the numerator cuts:
- kinematic cuts ensure consistency with fiducial definitions,
- `trigger` models online selection survival,
- `vProb` and candidate/charge requirements model analysis-level quality and topology constraints.

## Efficiency Calculation
This part uses a full-chain MC (generation, detector simulation, particle reconstruction, etc.), which hopes to match the process of collecting data in reality.

```bash
cd /path/to/CMSDAS/efficiency/mc_efficiency
root -l -b -q mc_efficiency.C
```

Core logic:
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

Outputs are efficiency maps in $(p_T,|y|)$:
- `efficiency/mc_efficiency/results/efficiency.csv`
- `efficiency/mc_efficiency/results/efficiency.pdf`

> #### **Checkpoint**
> Compare efficiency and acceptance maps directly.
> They should not have identical structures.

> #### **Question**
> In low-$p_T$ bins, which numerator requirement is most likely to dominate the efficiency loss (`trigger`, `vProb`, or single-candidate condition), and why?
