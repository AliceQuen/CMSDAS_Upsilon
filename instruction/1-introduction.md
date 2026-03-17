# Part 1 - Introduction

## Overview
- This is the first analysis step that uses the new 2025 trigger path `HLT_Dimuon0_Upsilon`.
- Enabled by this trigger, the measurement covers both low-$p_\mathrm{T}$ and high-$p_\mathrm{T}$ regions in one workflow.
- The target phase space spans $p_\mathrm{T}=0-130$ GeV and $|y|<2.4$.
- This is the first $\Upsilon$ cross-section measurement based on 2025 data and also serves as a data-quality validation.


Note: in collider physics, the symbols for common physical varibles are defined as follows:
- The transverse momentum, $p_\mathrm{T}$, is the component of a particle momentum perpendicular to the beam axis;
- The rapidity, $y$ is defined as $\frac{1}{2}\ln\frac{E+p_z}{E-p_z}$;
- The pseudorapidity, $\eta$, is defined as $-\ln\tan(\theta/2)$, where $\theta$ is the polar angle with respect to the beam direction;
- For highly relativistic particles, $y$ and $\eta$ are close, but for massive particles they are not exactly the same.

> #### **Question**
> 1. In the analysis that follows, we use $\eta$ in selections of $\mu$ and $y$ in selections of $\Upsilon$. Given the previous definitions, consider why.


### References
#### Analysis References
- [Pre-Approval talk](https://indico.cern.ch/event/1505578/)
- [Approval talk](https://indico.cern.ch/event/1602931/#2-approval-of-bph-24-004-measu)
- CADI: [BPH-24-004](https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BPH-24-004)
- Analysis Note: [AN-23-142](https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2023/142)

#### Technical References
- [Linux command line for beginners](https://ubuntu.com/tutorials/command-line-for-beginners)
- [ROOT documentation](https://root.cern.ch/root/htmldoc/guides/users-guide/ROOTUsersGuide.html)
- [ROOT tutorials and examples](https://root.cern/doc/master/group__Tutorials.html)
- [RooFit documentation](https://root.cern.ch/doc/master/group__Roofitmain.html)
- [CMS Offline WorkBook (official CMS software introduction and workflow reference)](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBook)
- [CMSSW setup and first steps ](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter1)
- [CMS DAS pre-exercises](https://fnallpc.github.io/cms-das-pre-exercises/)
- [PDG live (the reference to the particle properties)](https://pdglive.lbl.gov/)


## Cross Section
The (production) cross section of $\Upsilon$ characterizes its production probability in the collision. Because the $\Upsilon$ has a very short lifetime, it decays before it can be detected by the detector. Therefore, we can only infer whether an $\Upsilon$ meson was produced by analyzing the more stable particles from $\Upsilon$'s decay, thereby measuring the cross section.

More specifically, in this analysis, we will use the decay of $\Upsilon\to\mu^+\mu^-$ to do the measurement. 
> #### **Question**
> 2. Using PDG Live (https://pdglive.lbl.gov), search common decay modes of $\Upsilon(\mathrm{nS})$ ($n=1,2,3$). Why we choose this decay channel ($\Upsilon\to\mu^+\mu^-$) to measure the cross section?

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
Trigger is the online event-selection system used by CMS to reduce the collision data size to a manageable level for storage and reconstruction. It selects potentially interesting events in real time based on detector signatures and physics objects. In CMS, the trigger system is usually divided into:
- L1 (hardware level);
- HLT (software level).

`HLT_Dimuon0_Upsilon` is an HLT trigger oath, applied after L1 preselection.

To inspect the full trigger content, use [cmshltinfo](https://cmshltinfo-dev.app.cern.ch/summary):
1. Search `HLT_Dimuon0_Upsilon`. If not found, set year to `2025` and ensure `Parking` is selected in `Stream Select`.
2. Open the path and go to `filters`.
3. Choose a run range.
4. Use the `event` decision to understand how events propagate through the HLT path.

In trigger names, the prefix usually tells the trigger level: `HLT_` denotes a High-Level Trigger path. The object part of the name describes the reconstructed candidate and its main requirements. 

For example, `HLT_Dimuon0_Upsilon` means an HLT path targeting dimuon candidates, with a low threshold encoded by `0`, and specifically designed for the $\Upsilon \to \mu^+\mu^-$ region. Likewise, the older `HLT_Dimuon10_Upsilon_y1p4` indicates a dimuon $\Upsilon$ path with $p_\mathrm{T}>10$ GeV and $|y|<1.4$ requirements.


> #### **Task**
> Find the trigger path used in the 2022 $\Upsilon$ cross-section measurement (`HLT_Dimuon10_Upsilon_y1p4`), inspect its filters that contain $p_\mathrm{T}>10$ GeV and $|y|<1.4$ requirements, and summarize the key differences with the 2025 setup.

> #### **Question**
> 3. Compare `HLT_Dimuon10_Upsilon_y1p4` and `HLT_Dimuon0_Upsilon`: what filter-level differences do you see?
> 4. What is the benefit of using `HLT_Dimuon0_Upsilon`?

## Analysis Overview
- `Part 2` (`Data`): produce/inspect ntuples and validate baseline data distributions.
- `Part 3` (`Yield`): extract $N_{1S}$, $N_{2S}$, $N_{3S}$ with mass fits in each analysis bin.
- `Part 4` (`Acceptance`): evaluate geometric/kinematic acceptance from generated-level information.
- `Part 5` (`Efficiency`): quantify trigger/reconstruction/selection efficiency using full-chain MC.
- `Part 6` (`Systematics`): estimate uncertainty components and build a controlled acceptance example.
- `Part 7` (`Results`): combine all ingredients into differential cross sections and validate final outputs.
- `Part 8` (`Event Display`, optional): visualize selected events with `cmsShow`.
- `Part 9` (`Polarization`, optional): study polarization effects on acceptance variations with event-level reweighting.
