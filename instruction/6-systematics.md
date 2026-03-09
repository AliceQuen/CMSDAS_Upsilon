# Part 6 - Systematics

This part introduces systematic uncertainties and shows one concrete example in detail.

## Sources
Every term in the cross-section expression can carry systematic uncertainty:
$$
\mathcal{B}\frac{d^2\sigma}{dp_T\,d|y|}
\propto
\frac{N}{\mathcal{L}\,A\,\epsilon}.
$$

This means uncertainty can come from:
- yield extraction ($N$),
- acceptance modeling ($A$),
- efficiency modeling ($\epsilon$),
- luminosity normalization ($\mathcal{L}$).

A common per-bin combination is
$$
\delta_{\mathrm{tot}}=\sqrt{\delta_{\mathrm{yield}}^2+\delta_A^2+\delta_{\epsilon}^2+\delta_{\mathrm{lumi}}^2}.
$$

Note: luminosity uncertainty is usually provided officially by the experiment's luminosity group.

## Example: Acceptance Systematics
Due to limited time, we will use acceptance systematics as the main hands-on example.
- acceptance depends strongly on kinematic thresholds near phase-space boundaries,
- threshold variations provide an intuitive and controlled way to probe modeling sensitivity.

Threshold-variation setup:
- nominal: $p_T^{\mu,\min}=3.1$ GeV,
- down: $p_T^{\mu,\min}=3.0$ GeV,
- up: $p_T^{\mu,\min}=3.2$ GeV.

> #### **Task**
> 1. Copy and modify the acceptance program variants:
>    ```bash
>    cd /path/to/CMSDAS/acceptance
>    cp acceptance.C acceptance_pt3p0.C
>    cp acceptance.C acceptance_pt3p2.C
>    ```
> 2. In each variant, change both the threshold and output file names to avoid overwriting nominal results.
>    Threshold line:
>    ```cpp
>    const double kMuonPtMin     = 3.1;
>    ```
> 3. Run nominal/down/up and produce three CSV files.
> 4. Compute per-bin acceptance systematic:
>    $$
>    \delta_A^{\mathrm{syst}}=\frac{\max\left(|A_{\mathrm{up}}-A_{\mathrm{nom}}|,|A_{\mathrm{down}}-A_{\mathrm{nom}}|\right)}{A_{\mathrm{nom}}}.
>    $$
> 5. Build a generated-method bias estimate (continuation of `Part 4`):
>    $$
>    R_{\mathrm{closure}}=\frac{A_{\mathrm{alt}}}{A_{\mathrm{gen}}},\qquad
>    \delta_{A,\mathrm{model}}=|1-R_{\mathrm{closure}}|.
>    $$
>    Here, $A_{\mathrm{gen}}$ is the current generated-level acceptance, and $A_{\mathrm{alt}}$ is your alternative definition (for example including reconstruction-related effects).

> #### **Checkpoint**
> Verify that low-$p_T$ bins show larger variation, consistent with threshold-dominated behavior.

> #### **Question**
> 1. Which $(p_T,|y|)$ regions are expected to show the largest closure deviation, why, and does this match your expectation from trigger/reconstruction boundary effects?
