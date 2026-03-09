# Part 1 - Introduction

## Overview
- This is the first analysis step that uses the new 2025 trigger path `HLT_Dimuon0_Upsilon`.
- Enabled by this new trigger, the measurement covers low-$p_T$ and high-$p_T$ regions in one workflow.
- The target phase space is broad: $0$-$130$ GeV and $|y|<2.4$.
- This is the first $\Upsilon$ cross-section measurement based on 2025 data and also serves as a data-quality validation.

### References
- [Pre-Approval talk](https://indico.cern.ch/event/1505578/)
- [Approval talk](https://indico.cern.ch/event/1602931/#2-approval-of-bph-24-004-measu)
- CADI: [BPH-24-004](https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BPH-24-004)
- [AN-23-142](https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2023/142)

## Cross-section
For each state $\Upsilon(nS)$ ($n=1,2,3$), the measured quantity is
$$
\mathcal{B}(\Upsilon(nS)\to\mu^+\mu^-)\,\frac{d^2\sigma_n}{dp_T\,dy}
=\frac{N_n}{\mathcal{L}\,A_n\,\epsilon_n\,\Delta p_T\,\Delta y}.
$$

In this analysis you bin in $|y|$ and combine positive/negative rapidity, so $\Delta y=2\Delta|y|$ and
$$
\mathcal{B}(\Upsilon(nS)\to\mu^+\mu^-)\,\frac{d^2\sigma_n}{dp_T\,d|y|}
=\frac{N_n}{\mathcal{L}\,A_n\,\epsilon_n\,\Delta p_T\,2\Delta|y|}.
$$

Here:
- $N_n$ is the fitted signal yield.
- $A_n$ is the acceptance.
- $\epsilon_n$ is the efficiency.
- $\mathcal{L}$ is the integrated luminosity.

Each factor has a clear physical role:
- $N_n$ captures how many signal candidates are observed.
- $A_n$ corrects for geometric/kinematic phase-space loss.
- $\epsilon_n$ corrects for detector, trigger, and reconstruction effects.
- $\mathcal{L}$ normalizes event counts to cross section.

### Binning
$|y|$: $[0.0,0.6]$, $[0.6,1.2]$, $[1.2,1.8]$, $[1.8,2.4]$

$p_T$: $0$-$20$ in $1$ GeV width, $20$-$40$ in $2$ GeV width, $[40,43]$, $[43,46]$, $[46,50]$, $[50,55]$, $[55,60]$, $[60,70]$, $[70,100]$, $[100,130]$

> #### **Question**
> 1. Is it always better to use finer bins?
> 2. What problems can appear if bins are too narrow?
> 3. What problems can appear if bins are too wide?

## Trigger
The CMS trigger system has two levels:
- L1 (hardware level)
- HLT (software level)

`HLT_Dimuon0_Upsilon` is an HLT path, applied after L1 preselection.

To inspect the full trigger content, use [cmshltinfo](https://cmshltinfo-dev.app.cern.ch/summary):
1. Search `HLT_Dimuon0_Upsilon`.
2. If not found, set year to `2025`.
3. Ensure `Parking` is selected in `Stream Select`.
4. Open the path and go to `filters`.
5. Choose a run range.
6. In the second filter (`hltL1s12ForUpsilonDimuon0Mass8to12`), read `L1SeedsLogicalExpression` to identify L1 seeds.
7. Use the `event` decision to understand how events propagate through the HLT path.

- Trigger choices define which phase-space region is visible.
- Lower trigger thresholds are crucial for low-$p_T$ cross-section precision.

> #### **Task**
> Find the trigger path used in the 2022 $\Upsilon$ cross-section measurement, inspect its filters, and summarize the key differences with the 2025 setup.

> #### **Question**
> 1. Compare `HLT_Dimuon10_Upsilon` and `HLT_Dimuon0_Upsilon`: what filter-level differences do you see?
> 2. What is the benefit of using `HLT_Dimuon0_Upsilon`?

## Overview
- `Part 2` (`Data`): produce/inspect ntuples and validate baseline data distributions.
- `Part 3` (`Yield`): extract $N_{1S}$, $N_{2S}$, $N_{3S}$ with mass fits in each analysis bin.
- `Part 4` (`Acceptance`): evaluate geometric/kinematic acceptance from generated-level information.
- `Part 5` (`Efficiency`): quantify trigger/reconstruction/selection efficiency using full-chain MC.
- `Part 6` (`Systematics`): estimate uncertainty components and build a controlled acceptance example.
- `Part 7` (`Results`): combine all ingredients into differential cross sections and validate final outputs.
