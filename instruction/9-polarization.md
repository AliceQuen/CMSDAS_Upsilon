# Part 9 - Polarization Effects on Acceptance (Optional)

Quarkonium polarization changes the decay-angle distribution of the muons, which can modify kinematic acceptance.
A common angular variable is $\theta^*$ in the helicity frame, where the $z$ axis is defined by the $\Upsilon$ momentum direction in the laboratory frame.

<p align="center">
  <img src="figures/HX_frame.png" alt="Helicity frame definition" width="50%">
</p>

In the helicity frame, a standard parameterization is
$$
\frac{dN}{d\cos\theta^*} \propto 1 + \lambda_\vartheta \cos^2\theta^*.
$$

The default acceptance maps in this exercise assume unpolarized production, consistent with previous measurements of $\Upsilon$ polarization. If the true production is polarized, the muon angular and kinematic distributions change.
Therefore, acceptance changes, and the extracted cross section changes as well.

A practical strategy is to evaluate two extreme scenarios for $\theta^*$ polarization:
- transverse polarization: $\lambda_\vartheta = +1$,
- longitudinal polarization: $\lambda_\vartheta = -1$.

For any polarization scenario, you can decompose it into transverse/longitudinal components and propagate acceptance changes to cross-section corrections:
$$
R_A(\lambda_\vartheta) = \frac{A(\lambda_\vartheta)}{A(0)},
\qquad
R_\sigma(\lambda_\vartheta) = \frac{\sigma(\lambda_\vartheta)}{\sigma(0)} = \frac{1}{R_A(\lambda_\vartheta)}.
$$

## Polarized Particle-Gun Generation
To evaluate $A(\lambda_\vartheta=\pm 1)$, we need to generate polarized particle-gun samples.

Available files are:
- `polarization/ParticleGun-Upsilon2MM.py`
- `polarization/run_ntuple_gen.sh`
- `polarization/submit_ntuple_gen.jdl`

> #### **Task**
> Enable one polarization scenario (choose either transverse or longitudinal) in `ParticleGun-Upsilon2MM.py`.
> 
> Hints:
> 1. Focus on the `process.generator = cms.EDFilter("Pythia8PtGun", ...)` block.
> 2. The relevant place is in `processParameters`, where decay settings are defined.
> 3. Add polarization-related steering in the same block.
> 4. Due to time, you can generate one polarization case only.

### Condor Production
We will use HT-Condor to submit production jobs of polarized Particle-Gun MC.
- In `ParticleGun-Upsilon2MM.py`: `process.maxEvents.input` sets events per `cmsRun` job.
- In `submit_ntuple_gen.jdl`: `Queue JOBNUM from seq 1 1000 |` sets the number of submitted jobs.
- In `run_ntuple_gen.sh`: `seed` is derived from `JOBNUM`, and output transfer is controlled by `OUTPUT_DIR`.
- Note: on lxplus8, run Condor from AFS (not CERNBox).

> #### **Task**
> Submit polarized generation on Condor from AFS.
>
> 1. Go to your user AFS area and create a clean working directory:
>    ```bash
>    cd /afs/cern.ch/user/<u>/<username>
>    mkdir -p condor/CMSDAS
>    cd condor/CMSDAS
>    mkdir -p log
>    ```
> 2. Copy three files into this directory:
>    - modified `ParticleGun-Upsilon2MM.py`,
>    - `run_ntuple_gen.sh`,
>    - `submit_ntuple_gen.jdl`.
> 3. Edit `submit_ntuple_gen.jdl`:
>    - set `OUTPUT_DIR`,
>    - set `x509userproxy`.
> 4. Create proxy file and store it at the path used in JDL:
>    ```bash
>    voms-proxy-init -voms cms -valid 168:00 -out /afs/cern.ch/user/<u>/<username>/condor/CMSDAS/x509up
>    ```
> 5. Submit and monitor:
>    ```bash
>    condor_submit submit_ntuple_gen.jdl
>    condor_q
>    ```
>
> Typical runtime is about 1-2 hours.
> When `condor_q` no longer shows your jobs, check `OUTPUT_DIR` for produced acceptance ntuples.

## Acceptance Re-calculation with Polarized MC
> #### **Task**
> Recompute acceptance with your polarized ntuples.
>
> Hints:
> 1. In `acceptance/acceptance.C`, update the input ROOT pattern (wildcards are allowed).
> 2. Change output CSV/PDF names to avoid overwriting nominal outputs.
> 3. Run the acceptance program and store results for comparison with the nominal map.

### Polarization Correction Factors
> #### **Task**
> Write a small program that reads:
> - nominal acceptance CSV,
> - polarized acceptance CSV,
> and outputs per-bin:
> - acceptance correction factor $R_A = A_{\mathrm{pol}}/A_{\mathrm{nom}}$,
> - cross-section correction factor $R_\sigma = 1/R_A$.
>
> Hints:
> 1. Match rows by `(pt_min, pt_max, y_abs_min, y_abs_max)`.
> 2. Keep both correction factors in the output table for later propagation.
> 3. Inspect the strongest deviations in low-$p_T$ and forward-rapidity bins first.

> #### **Question**
> 1. Which kinematic regions are most sensitive to polarization effects on acceptance?
> 2. For your generated scenario, are the correction factors small enough to support the unpolarized approximation, or do they indicate a potentially non-negligible bias?
