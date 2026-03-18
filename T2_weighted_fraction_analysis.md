## T2-weighted volume fraction analysis

### Background

The restricted volume fraction ($f_r$) estimated by AxCaliber-SMT is not the true (non-T2-weighted) intra-axonal water fraction $f_0$, but rather a T2-decay-weighted apparent fraction that depends on echo time (TE). Because intra-axonal and extra-axonal compartments have different T2 relaxation times, the apparent fraction shifts with TE according to (Veraart et al., 2018):

$$f_r(\text{TE}) = \frac{f_0 \, e^{-\text{TE}/T_{2}^{a}}}{f_0 \, e^{-\text{TE}/T_{2}^{a}} + (1 - f_0) \, e^{-\text{TE}/T_{2}^{e}}}$$

where $T_2^a$ and $T_2^e$ are the intra-axonal and extra-axonal transverse relaxation times, respectively.

Since $T_2^a > T_2^e$ in white matter (Veraart et al., 2018), the extra-axonal signal decays faster with increasing TE. At longer TE, the extra-axonal contribution is more suppressed, inflating the apparent restricted fraction. Because the C1 protocol operates at a longer TE (77 ms) than C2 (54 ms), C1 is expected to yield systematically higher $f_r$ estimates than C2 for the same underlying tissue.

### Method

We computed the predicted T2-weighted restricted fraction for 10 major white matter ROIs at both echo times, using compartmental T2 values read from Figure 5 of Veraart et al. (2018) and a representative non-T2-weighted fraction of $f_0 = 0.45$. The ROIs and their T2 values are:

| ROI  | $T_2^a$ (ms) | $T_2^e$ (ms) |
|------|:---:|:---:|
| PLIC | 110 | 42 |
| ALIC |  85 | 48 |
| SLF  |  90 | 38 |
| EC   |  70 | 55 |
| GCC  |  75 | 35 |
| BCC  |  90 | 45 |
| SCC  |  90 | 40 |
| ACR  |  75 | 48 |
| SCR  | 100 | 48 |
| PCR  |  95 | 48 |

### Results

| ROI  | $f_r$ (C2, TE=54 ms) | $f_r$ (C1, TE=77 ms) | Difference |
|------|:---:|:---:|:---:|
| PLIC | 0.645 | 0.717 | +11.2% |
| ALIC | 0.572 | 0.621 | +8.6%  |
| SLF  | 0.650 | 0.725 | +11.5% |
| EC   | 0.503 | 0.524 | +4.2%  |
| GCC  | 0.651 | 0.725 | +11.5% |
| BCC  | 0.599 | 0.655 | +9.4%  |
| SCC  | 0.634 | 0.706 | +11.3% |
| ACR  | 0.551 | 0.592 | +7.5%  |
| SCR  | 0.595 | 0.652 | +9.7%  |
| PCR  | 0.588 | 0.643 | +9.4%  |
| **Mean** | **0.599** | **0.656** | **+9.4%** |

Across all ROIs, the T2 weighting effect predicts that C1 yields 4–12% higher restricted fractions than C2 (mean 9.4%). The effect is largest in ROIs with the greatest $T_2^a / T_2^e$ ratio (e.g., SLF, GCC, SCC) and smallest where the two compartmental T2 values are closer (e.g., EC).

### References

Veraart, J., Novikov, D. S., & Fieremans, E. (2018). TE dependent Diffusion Imaging (TEdDI) distinguishes between compartmental T2 relaxation times. *NeuroImage*, 182, 360–369. https://doi.org/10.1016/j.neuroimage.2017.09.030
