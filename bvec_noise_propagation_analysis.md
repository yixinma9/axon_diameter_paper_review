## Gradient direction–dependent noise propagation analysis (b-vector GNL)

### Background

Gradient nonlinearity (GNL) affects diffusion MRI measurements by spatially varying the effective gradient field. The standard scalar b-value correction accounts for changes in gradient *magnitude* (b-scaling), but GNL also rotates the effective gradient direction. For the Connectome scanners, which use high-performance gradient coils (C2: $G_\text{max} = 495$ mT/m; C1: $G_\text{max} = 290$ mT/m), these spatial deviations can be substantial. This analysis quantifies the impact of the full gradient deviation tensor—including both b-scaling and b-rotation—on AxCaliber-SMT axon radius estimates, and evaluates whether a scalar $b_\text{scale}$ correction is sufficient.

### Method

#### Gradient deviation model

For each voxel, the gradient deviation is described by a $3 \times 3$ tensor $\mathbf{L}$, obtained from FSL's `calc_grad_perc_dev` applied to the gradient coil's nonlinearity field. The effective gradient for a nominal direction $\hat{g}$ is:

$$\mathbf{g}_\text{eff} = (\mathbf{I} + \mathbf{L}) \, \hat{g} = \mathbf{M} \, \hat{g}$$

This changes both the gradient magnitude ($|\mathbf{g}_\text{eff}|$) and direction ($\hat{g}_\text{eff}$). The direction-specific b-scaling factor is $|\mathbf{g}_\text{eff}|^2$, and the scalar (direction-averaged) b-scale used for correction is:

$$b_\text{scale} = \text{tr}(\mathbf{M}^\top \mathbf{M}) \,/\, 3$$

#### Signal generation

For each protocol (C2, C1), we load the actual diffusion gradient tables (bvec/bval files) and the gradient deviation maps (grad_dev). The protocol uses 16 b-value shells with $\delta = 8$ ms, $\Delta = 19$ ms (short) or $\Delta = 49$ ms (long), spanning $b = 0.05$–$17.8$ ms/$\mu$m$^2$. Each shell has 32 or 64 gradient directions.

We sample $N = 2000$ white-matter voxels from the WM mask (excluding CSF and GM), and for each voxel draw random ground-truth parameters: axon radius $r \in [0.1, 5]$ $\mu$m, intra-axonal fraction $f \in [0.5, 1]$, CSF fraction $f_\text{csf} \in [0, 0.2]$, and extra-axonal radial diffusivity $D_e^\perp \in [0.5, 1.5]$ $\mu$m$^2$/ms.

For each gradient direction $\hat{g}_k$ in shell $s$, the single-fiber signal is:

$$S_k = (1 - f_\text{csf}) \big[ f \, S_a(\hat{g}_k) + (1-f) \, S_e(\hat{g}_k) \big] + f_\text{csf} \, e^{-b_s D_\text{csf}}$$

where:
- $S_a$ is the intra-axonal (restricted cylinder) signal using the Van Gelderen model with 10 Bessel zeros, depending on $r$, $D_0$, $D_a$, $\delta$, $\Delta$, and the perpendicular gradient component $g_\perp^2$
- $S_e$ is the extra-axonal (zeppelin) signal: $S_e = \exp\!\big(-b (D_e^\perp + (D_e^\parallel - D_e^\perp) \cos^2\theta)\big)$
- $\theta$ is the angle between the gradient direction and the fiber axis ($\hat{z}$)

The powder-averaged signal per shell is then $\bar{S}_s = \frac{1}{N_s} \sum_{k=1}^{N_s} S_k$.

Fixed diffusivities: $D_0 = 2$, $D_a = 1.7$, $D_e^\parallel = 1.7$, $D_\text{csf} = 3$ $\mu$m$^2$/ms. SNR: 38 (C2) and 18 (C1).

#### Three scenarios

For each voxel, a single Rician noise realization is shared across all three scenarios to isolate the effect of GNL:

1. **Baseline**: Signal generated with nominal gradient directions (no GNL). Fit with nominal b-values.
2. **Uncorrected**: Signal generated with GNL-perturbed gradients (per-direction b-scaling and b-rotation via $\mathbf{M}$). Fit with nominal b-values (no correction).
3. **Corrected**: Same GNL signal as (2). Fit with scalar-corrected b-values ($b_\text{corr} = b \cdot b_\text{scale}$).

The GNL signal (scenarios 2–3) applies the full $\mathbf{M}$ matrix to all gradient directions simultaneously:
- Effective gradient: $\mathbf{g}_\text{eff} = \mathbf{M} \, \hat{g}_k$
- Direction-specific b-scaling: $b_{\text{eff},k} = b_s \cdot |\mathbf{g}_\text{eff}|^2$
- Rotated angle: $\theta_k = \arccos(|\hat{g}_{\text{eff},k} \cdot \hat{z}|)$

#### Fitting

All fitting uses `AxCaliberSMT.mcmc` from the C2 protocol design library (Lee et al.), with the Van Gelderen restricted diffusion model and 5 MCMC restarts. This is identical to the fitting method used in the scalar b-scale noise propagation analysis, ensuring a fair comparison. The only difference between the two analyses is how GNL affects the synthetic signal: scalar b-scaling only (b-scale analysis) versus full directional perturbation (this analysis).

#### Filtering

Voxels where any of the three scenarios produced $r_\text{fit}$ at the optimization bounds ($r < 0.05$ or $r > 4.95$ $\mu$m) are excluded from all statistics and plots.

### Results

*(To be added)*
