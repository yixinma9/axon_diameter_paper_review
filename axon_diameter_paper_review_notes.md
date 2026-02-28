# Gradient Nonlinearity Correction: b_scale for AxCaliber-SMT

Reference: QT (Qiuyun Fan)

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
- **Output**: `grad_dev_{x,y,z}.nii.gz` — 3 files, each `[x,y,z,3]`

```bash
${FSLDIR}/bin/calc_grad_perc_dev \
    --fullwarp=/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc/b0s_mean_deform_grad_rel.nii.gz \
    -o /autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc/grad_dev
```

Merge into a single 9-volume file (matching QT's `grad_dev.nii.gz` format):

```bash
fslmerge -t grad_dev.nii.gz grad_dev_x.nii.gz grad_dev_y.nii.gz grad_dev_z.nii.gz
```

## Step 3. Correct bvals/bvecs per voxel (QT's method)

For a given voxel `(i, j, k)`:

```matlab
g = read_avw('grad_dev.nii.gz');   % [x,y,z,9]
bvecs = load('bvecs');              % 3xN
bvals = load('bvals');              % 1xN
L = reshape(squeeze(g(i,j,k,:)), 3, 3);
I = eye(3);

% correct bvecs and calculate their norm
v = (I+L)*bvecs;
n = sqrt(sum(v.^2));

% Normalise corrected bvecs and correct bvals
new_bvecs = v ./ repmat(n,3,1);
new_bvals = n.^2 .* bvals;
```

## Step 4. Compute b_scale map for spherical mean

For spherical mean signal (direction-averaged), the scalar b_scale per voxel is the mean of `n^2` over all directions, which equals `trace((I+L)'*(I+L))/3`:

```matlab
g = read_avw('grad_dev.nii.gz');
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
