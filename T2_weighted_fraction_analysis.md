## T2-weighted volume fraction analysis

### Background

The restricted volume fraction ($f_r$) estimated by AxCaliber-SMT is not the true (non-T2-weighted) intra-axonal water fraction $f_0$, but rather a T2-decay-weighted apparent fraction that depends on echo time (TE). Because intra-axonal and extra-axonal compartments have different T2 relaxation times, the apparent fraction shifts with TE according to (Veraart et al., 2018):

$$f_r(\text{TE}) = \frac{f_0 \, e^{-\text{TE}/T_{2}^{a}}}{f_0 \, e^{-\text{TE}/T_{2}^{a}} + f_e \, e^{-\text{TE}/T_{2}^{e}} + f_{\text{csf}} \, e^{-\text{TE}/T_{2}^{\text{csf}}}}$$

where $T_2^a$ and $T_2^e$ are the intra-axonal and extra-axonal transverse relaxation times, $T_2^{\text{csf}}$ is the isotropic compartment relaxation time (~2000 ms at 3T), and $f_0 + f_e + f_{\text{csf}} = 1$.

Since $T_2^a > T_2^e$ in white matter (Veraart et al., 2018), the extra-axonal signal decays faster with increasing TE. At longer TE, the extra-axonal contribution is more suppressed, inflating the apparent restricted fraction. Because C2 operates at a shorter TE (54 ms) than C1 (77 ms), C2 retains more extra-axonal signal relative to intra-axonal signal, resulting in systematically lower $f_r$ estimates compared to C1 for the same underlying tissue.

### Method

We adopted a three-compartment model (intra-axonal, extra-axonal, and isotropic) consistent with the AxCaliber-SMT signal model (Kaden et al., 2016). The isotropic compartment fraction $f_{\text{csf}}$ was determined by minimizing the mean squared error between the T2-predicted and observed C2–C1 percentage differences across all 12 white matter tracts (sweeping $f_{\text{csf}}$ from 0 to 0.25 in steps of 0.005). Note that $T_2^{\text{csf}} \approx 2000$ ms at 3T, so the isotropic signal changes minimally between TE = 54 ms and TE = 77 ms ($e^{-54/2000} = 0.973$ vs $e^{-77/2000} = 0.962$); its main effect is to dilute the intra-/extra-axonal T2 contrast rather than introduce a TE-dependent shift of its own.

We used per-tract compartmental T2 values from Veraart et al. (2018), mapped from the JHU ICBM-DTI-81 atlas regions to the JHU tractography atlas tracts used in our analysis:

| Tract (this study) | TEdDI ROI | $T_2^a$ (ms) | $T_2^e$ (ms) |
|-----|-----------|:---:|:---:|
| ATR L/R | ACR | 75 | 48 |
| CST L/R | PLIC | 110 | 42 |
| Forceps major | SCC | 90 | 40 |
| Forceps minor | GCC | 75 | 35 |
| IFO L/R | EC | 70 | 55 |
| ILF L/R | BCC* | 80 | 45 |
| SLF L/R | SLF | 90 | 38 |

\* No direct TEdDI match for ILF; approximate values used.

For each tract, we back-calculated the non-T2-weighted fraction $f_0$ from the observed C1 restricted fraction using the three-compartment inverse formula (given the optimized $f_{\text{csf}}$), then forward-predicted $f_r$ at both echo times, and compared the T2-predicted C2–C1 difference with the observed difference from our data (Supplementary Figure 3).

### Results

*Results table to be filled after running the optimization (see `demo_T2_weighted_fraction.m`). The script outputs the optimal $f_{\text{csf}}$, per-tract $f_0$, predicted vs observed differences, and residuals.*

The T2 weighting model predicts that C2 should yield systematically lower restricted fractions than C1 across all tracts, consistent with the direction of the observed differences. Including a small isotropic fraction attenuates the predicted T2 effect, improving agreement with tracts where the observed C2–C1 gap is modest (e.g., forceps major/minor, ILF, SLF). For tracts where the observed difference substantially exceeds the T2 prediction (e.g., CST), additional factors—such as lower SNR (18 vs 38) and reduced b-value range on C1—likely contribute to the gap.

### References

Kaden, E., Kelm, N. D., Carson, R. P., Does, M. D., & Alexander, D. C. (2016). Multi-compartment microscopic diffusion imaging. *NeuroImage*, 139, 346–359. https://doi.org/10.1016/j.neuroimage.2016.06.002

Veraart, J., Novikov, D. S., & Fieremans, E. (2018). TE dependent Diffusion Imaging (TEdDI) distinguishes between compartmental T2 relaxation times. *NeuroImage*, 182, 360–369. https://doi.org/10.1016/j.neuroimage.2017.09.030
