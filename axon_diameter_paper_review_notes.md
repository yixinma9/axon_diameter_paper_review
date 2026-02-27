# Gradient Nonlinearity Correction: b_scale for AxCaliber-SMT

Reference: QT (Qiuyun Fan), WashU HCP pipeline

---

## Step 1. Compute GNC warp field

Use Qiuyun's script to compute the gradient nonlinearity warp field (C2 scanner).

- **Input**: any reference image in native DWI space
- **Output**: unwarped image + `*_deform_grad_rel.nii.gz` (relative warp field)

```bash
/space/scheherazade/2/users/qfan/tools/preproc/hcps_diff_prep_v2/hcps_gnc_c2.sh \
    -i /autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/b0s_mean.nii.gz \
    -o /autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc/b0s_mean_gnc.nii.gz \
    -interp spline
```

## Step 2. Compute grad_dev

Compute the 3x3 gradient deviation tensor per voxel from the relative warp field.

- **Input**: `*_deform_grad_rel.nii.gz` from Step 1
- **Output**: `grad_dev.nii.gz` `[x,y,z,9]` — 9 elements of the L matrix per voxel

```bash
${FSLDIR}/bin/calc_grad_perc_dev \
    --fullwarp=/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc/b0s_mean_deform_grad_rel.nii.gz \
    -o /autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc/grad_dev
```

## Step 3. Compute b_scale map in MATLAB

`b_scale = trace(M'*M)/3` where `M = I + L`

This is the direction-averaged b-value scaling factor per voxel (unitless, ~1.0 for ideal gradients).

```matlab
g = niftiread('grad_dev.nii.gz');
dims = size(g, 1:3);
b_scale = zeros(dims, 'single');
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            L = reshape(squeeze(g(i,j,k,:)), 3, 3);
            M = eye(3) + L;
            b_scale(i,j,k) = trace(M' * M) / 3;
        end
    end
end
```

Then pass to AxCaliberSMT:

```matlab
extradata.b_scale = b_scale;
out = smt.estimate(dwi, mask, extradata, fitting);
```
