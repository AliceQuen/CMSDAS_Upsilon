# Part 1 - Introduction

## Overview
- This is the first analysis step that uses the new 2025 trigger path `HLT_Dimuon0_Upsilon`.
- Enabled by this trigger, the measurement covers both low-$p_\mathrm{T}$ and high-$p_\mathrm{T}$ regions in one workflow.
- The target phase space spans $p_\mathrm{T}=0-130$ GeV and $|y|<2.4$.
- This is the first $\Upsilon$ cross-section measurement based on 2025 data and also serves as a data-quality validation.

### References
- [Pre-Approval talk](https://indico.cern.ch/event/1505578/)
- [Approval talk](https://indico.cern.ch/event/1602931/#2-approval-of-bph-24-004-measu)
- CADI: [BPH-24-004](https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BPH-24-004)
- Analysis Note: [AN-23-142](https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2023/142)

## Cross Section
The (production) cross section of $\Upsilon$ characterizes its production probability in the collision. Because the $\Upsilon$ has a very short lifetime, it decays before it can be detected by the detector. Therefore, we can only infer whether an $\Upsilon$ meson was produced by analyzing the more stable particles from $\Upsilon$'s decay, thereby measuring the cross section.

More specifically, in this analysis, we will use the decay of $\Upsilon$ into two muons to do the measurement. 

> #### **Question**
> 1. Why we choose this decay channel to measure the cross section?

For each state $\Upsilon(\mathrm{nS})$ ($n=1,2,3$), the measured cross section is

$\mathcal{B}(\Upsilon(\mathrm{nS})\to\mu^+\mu^-)\,\frac{d^2\sigma_n}{dp_\mathrm{T}\,d|y|}
=\frac{N_n}{\mathcal{L}\,A_n\,\epsilon_n\,\Delta p_\mathrm{T}\,2\Delta|y|}$.

Here:
- $\mathcal{B}(\Upsilon(\mathrm{nS})\to\mu^+\mu^-)$ is the branching ratio (or probability) of $\Upsilon(\mathrm{nS})$ decaying into two muons.
- $N_n$ is the fitted signal yield.
- $A_n$ is the acceptance.
- $\epsilon_n$ is the efficiency.
- $\mathcal{L}$ is the integrated luminosity.
- The rapidity bins are defined in $|y|$, so the full rapidity width is $\Delta y = 2\Delta|y|$.

### Binning
$|y|$: $[0.0,0.6]$, $[0.6,1.2]$, $[1.2,1.8]$, $[1.8,2.4]$

$p_\mathrm{T}$: $0-20$ in $1$ GeV width, $20-40$ in $2$ GeV width, $[40,43]$, $[43,46]$, $[46,50]$, $[50,55]$, $[55,60]$, $[60,70]$, $[70,100]$, $[100,130]$

> #### **Question**
> 1. Is it always better to use finer bins?
> 2. What problems can appear if bins are too narrow? What problems can appear if bins are too wide?

## Trigger
The CMS trigger system has two levels:
- L1 (hardware level)
- HLT (software level)

`HLT_Dimuon0_Upsilon` is an HLT path, applied after L1 preselection.

To inspect the full trigger content, use [cmshltinfo](https://cmshltinfo-dev.app.cern.ch/summary):
1. Search `HLT_Dimuon0_Upsilon`. If not found, set year to `2025` and ensure `Parking` is selected in `Stream Select`.
2. Open the path and go to `filters`.
3. Choose a run range.
4. In the second filter (`hltL1s12ForUpsilonDimuon0Mass8to12`), read `L1SeedsLogicalExpression` to identify L1 seeds.
5. Use the `event` decision to understand how events propagate through the HLT path.

> #### **Task**
> Find the trigger path used in the 2022 $\Upsilon$ cross-section measurement (`HLT_Dimuon10_Upsilon_y1p4`), inspect its filters, and summarize the key differences with the 2025 setup.

> #### **Question**
> 1. Compare `HLT_Dimuon10_Upsilon_y1p4` and `HLT_Dimuon0_Upsilon`: what filter-level differences do you see?
> 2. What is the benefit of using `HLT_Dimuon0_Upsilon`?

## Tutorial Roadmap
- `Part 2` (`Data`): produce/inspect ntuples and validate baseline data distributions.
- `Part 3` (`Yield`): extract $N_{1S}$, $N_{2S}$, $N_{3S}$ with mass fits in each analysis bin.
- `Part 4` (`Acceptance`): evaluate geometric/kinematic acceptance from generated-level information.
- `Part 5` (`Efficiency`): quantify trigger/reconstruction/selection efficiency using full-chain MC.
- `Part 6` (`Systematics`): estimate uncertainty components and build a controlled acceptance example.
- `Part 7` (`Results`): combine all ingredients into differential cross sections and validate final outputs.
- `Part 8` (`Event Display`, optional): visualize selected events with `cmsShow`.
- `Part 9` (`Polarization`, optional): study polarization effects on acceptance variations with event-level reweighting.
