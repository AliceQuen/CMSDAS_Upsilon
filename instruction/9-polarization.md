# Part 9 - Polarization Effects on Acceptance (Optional)

Quarkonium polarization changes the decay-angle distribution of the muons, which can modify kinematic acceptance.
A common angular variable is the polar angle $\theta^*$ in the helicity frame, where the $z$ axis is defined by the $\Upsilon$ momentum direction in the laboratory frame.

<p align="center">
  <img src="figures/HX_frame.png" alt="Helicity frame definition" width="50%">
</p>

In the helicity frame, a standard parameterization is 

$\frac{dN}{d\cos\theta^*} \propto 1 + \lambda_\vartheta \cos^2\theta^*$.

The default acceptance maps in this exercise assume unpolarized production, consistent with previous $\Upsilon$ polarization measurements. If the true production is polarized, the muon angular and kinematic distributions change.

Therefore, the acceptance changes, and the extracted cross section changes as well.

A practical strategy is to evaluate two extreme scenarios for $\theta^*$ polarization:
- transverse polarization: $\lambda_\vartheta = +1$,
- longitudinal polarization: $\lambda_\vartheta = -1$.

For any polarization scenario, you can decompose it into transverse and longitudinal components and propagate acceptance changes to cross-section corrections:

$R_A(\lambda_\vartheta) = \frac{A(\lambda_\vartheta)}{A(0)}, \qquad R_\sigma(\lambda_\vartheta) = \frac{\sigma(\lambda_\vartheta)}{\sigma(0)} = \frac{1}{R_A(\lambda_\vartheta)}$.

## Acceptance Recalculation with Polarization Re-weighting
We will start from the available unpolarized MC sample and split events into effective transverse and longitudinal classes with an event-by-event probabilistic assignment.

For the transverse and longitudinal polarizations, we have
$I_T(\cos\theta^*)=\frac{3}{8}\left(1+\cos^2\theta^*\right), \qquad I_L(\cos\theta^*)=\frac{3}{4}\left(1-\cos^2\theta^*\right)$.

For each event in the unpolarized sample, compute

$f_T(\cos\theta^*)=\frac{I_T(\cos\theta^*)}{I_T(\cos\theta^*)+I_L(\cos\theta^*)}$.

Then, generate a random number $r\in[0,1)$. If $r<f_T$, classify the event as transverse; otherwise classify it as longitudinal. The polarization assignment and the acceptance counting can be done in the same event loop as acceptance calculation.

> #### **Task**
> 1. Copy `acceptance/acceptance.C` to `acceptance/acceptance_polrw.C`.
> 2. Keep the same kinematic bins and fiducial cuts, so the only difference is the polarization treatment.
> 3. Inside the event loop, compute $\cos\theta^*$ from the `gen_dimuon_p4`, `gen_muonP_p4`, and `gen_muonN_p4` using the definitions in the figure above. Calculate $I_T$, $I_L$, and $f_T$, and perform the random assignment to the transverse or longitudinal class.
> 4. Use separate counters for the two classes, for example `All_T`, `Passed_T`, `All_L`, and `Passed_L`.
> 5. Write separate outputs (`.csv` and `.pdf`) for the transverse and longitudinal cases, with distinct file names to avoid overwriting the nominal results.
> 6. Run your new program:
>    ```bash
>    root -l acceptance_polrw.C
>    ```
>
> Hints:
> 1. In the HX definition used here, $\theta^*$ is the polar angle of one muon (in the **$\Upsilon$ rest frame**) with respect to the $\Upsilon$ flight direction (in the **laboratory frame**).

### Polarization Correction Factors
> #### **Task**
> Write a small program that reads:
> - nominal acceptance,
> - transverse acceptance,
> - longitudinal acceptance,
> and outputs polarization correction factors per-bin:
>   - $R_A^T = A_T/A_{\mathrm{nom}}$ and $R_A^L = A_L/A_{\mathrm{nom}}$,
>   - $R_\sigma^T = 1/R_A^T$ and $R_\sigma^L = 1/R_A^L$.


> #### **Question**
> 1. Which $(p_\mathrm{T},|y|)$ bins show the largest acceptance correction factors, and why?
> 2. Conduct literature review. Clarify the trend of the polarization effects on acceptance as functions of $p_\mathrm{T}$ and rapidity.
