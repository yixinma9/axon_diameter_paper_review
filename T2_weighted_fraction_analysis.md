## T2-weighted volume fraction analysis

### Background

The restricted volume fraction ($f_r$) estimated by AxCaliber-SMT is not the true (non-T2-weighted) intra-axonal water fraction $f_0$, but rather a T2-decay-weighted apparent fraction that depends on echo time (TE). Because intra-axonal and extra-axonal compartments have different T2 relaxation times, the apparent fraction shifts with TE according to (Veraart et al., 2018):

$$f_r(\text{TE}) = \frac{f_0 \, e^{-\text{TE}/T_{2}^{a}}}{f_0 \, e^{-\text{TE}/T_{2}^{a}} + (1 - f_0) \, e^{-\text{TE}/T_{2}^{e}}}$$

where $T_2^a$ and $T_2^e$ are the intra-axonal and extra-axonal transverse relaxation times, respectively.

Since $T_2^a > T_2^e$ in white matter (Veraart et al., 2018), the extra-axonal signal decays faster with increasing TE. At longer TE, the extra-axonal contribution is more suppressed, inflating the apparent restricted fraction. Because C2 operates at a shorter TE (54 ms) than C1 (77 ms), C2 retains more extra-axonal signal relative to intra-axonal signal, resulting in systematically lower $f_r$ estimates compared to C1 for the same underlying tissue.

### Method

We adopted a two-compartment model (intra-axonal and extra-axonal, $f_{\text{csf}} = 0$) for this analysis, following the approach of Veraart et al. (2018) and Kaden et al. (2016). Although AxCaliber-SMT includes an isotropic compartment in its signal model, the free-water fraction in deep white matter is small and, crucially, CSF has a very long $T_2$ (~2000 ms at 3T) that produces negligible signal change between TE = 54 ms and TE = 77 ms ($e^{-54/2000} = 0.973$ vs $e^{-77/2000} = 0.962$, a ~1% difference). Including the isotropic compartment would therefore dilute the intra-/extra-axonal T2 contrast without meaningfully contributing to the C1–C2 difference in $f_r$.

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

Rather than assuming a uniform ground-truth fraction, we used a data-driven approach: for each tract, we back-calculated the non-T2-weighted fraction $f_0$ from the observed C1 restricted fraction using the inverse of the T2-weighting formula:

$$f_0 = \frac{f_r \, e^{-\text{TE}/T_{2}^{e}}}{f_r \, e^{-\text{TE}/T_{2}^{e}} + (1 - f_r) \, e^{-\text{TE}/T_{2}^{a}}}$$

We then forward-predicted $f_r$ at both echo times from this tract-specific $f_0$, and compared the T2-predicted C2–C1 difference with the observed difference from our data (Supplementary Figure 3).

### Results

| Tract | $f_0$ | $f_r^{C1}$ obs | $f_r^{C2}$ obs | $f_r^{C2}$ pred | Obs C2–C1 | T2-pred C2–C1 | Residual |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| ATR L | 0.390 | 0.57 | 0.52 | 0.50 | −8.8% | −7.5% | −1.3% |
| ATR R | 0.369 | 0.55 | 0.49 | 0.48 | −10.9% | −7.5% | −3.4% |
| CST L | 0.530 | 0.72 | 0.55 | 0.65 | −23.6% | −11.2% | −12.4% |
| CST R | 0.505 | 0.70 | 0.54 | 0.63 | −22.9% | −11.2% | −11.7% |
| Forceps major | 0.383 | 0.57 | 0.55 | 0.50 | −3.5% | −11.3% | +7.8% |
| Forceps minor | 0.359 | 0.56 | 0.55 | 0.48 | −1.8% | −11.5% | +9.7% |
| IFO L | 0.375 | 0.47 | 0.40 | 0.44 | −14.9% | −4.2% | −10.7% |
| IFO R | 0.404 | 0.50 | 0.48 | 0.47 | −4.0% | −4.2% | +0.2% |
| ILF L | 0.380 | 0.55 | 0.55 | 0.49 | 0.0% | −9.0% | +9.0% |
| ILF R | 0.380 | 0.55 | 0.55 | 0.49 | 0.0% | −9.0% | +9.0% |
| SLF L | 0.367 | 0.56 | 0.55 | 0.48 | −1.8% | −11.5% | +9.7% |
| SLF R | 0.486 | 0.67 | 0.63 | 0.59 | −6.0% | −11.5% | +5.5% |

The T2 weighting model predicts that C2 should yield 4–12% lower restricted fractions than C1 across all tracts, consistent with the direction of the observed differences. For tracts such as ATR, the predicted and observed differences closely agree (−7.5% predicted vs −9 to −11% observed). For other tracts (e.g., CST), the observed difference substantially exceeds the T2 prediction, indicating that additional factors—such as lower SNR (18 vs 38) and reduced b-value range on C1—contribute to the gap. Conversely, for tracts where the observed difference is smaller than predicted (e.g., forceps, ILF, SLF L), compensating factors may partially offset the T2 effect.

### References

Kaden, E., Kelm, N. D., Carson, R. P., Does, M. D., & Alexander, D. C. (2016). Multi-compartment microscopic diffusion imaging. *NeuroImage*, 139, 346–359. https://doi.org/10.1016/j.neuroimage.2016.06.002

Veraart, J., Novikov, D. S., & Fieremans, E. (2018). TE dependent Diffusion Imaging (TEdDI) distinguishes between compartmental T2 relaxation times. *NeuroImage*, 182, 360–369. https://doi.org/10.1016/j.neuroimage.2017.09.030
