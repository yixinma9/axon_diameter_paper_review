clear; close all;
addpath(genpath('/autofs/space/linen_001/users/Yixin/C2_protocoldesign-main/lib'));

%% Step 1. Compute temporal SNR from b0.nii.gz (first 10 volumes)
data_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001';
gnc_dir  = fullfile(data_dir, 'gnc');

b0 = single(niftiread(fullfile(data_dir, 'b0s.nii.gz')));
Nvol = size(b0, 4);
fprintf('Using %d b0 volumes for tSNR\n', Nvol);

% tSNR = mean / std across repeats
tSNR_map = mean(b0, 4) ./ std(b0, 0, 4);
tSNR_map(isnan(tSNR_map) | isinf(tSNR_map)) = 0;

%% Step 2. Load b_scale map (gradient nonlinearity correction)
gx = single(niftiread(fullfile(gnc_dir, 'grad_dev_x.nii.gz')));
gy = single(niftiread(fullfile(gnc_dir, 'grad_dev_y.nii.gz')));
gz = single(niftiread(fullfile(gnc_dir, 'grad_dev_z.nii.gz')));
g  = cat(4, gx, gy, gz);  % [x,y,z,9]

dims = size(g, 1:3);
b_scale_map = zeros(dims, 'single');
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            L = reshape(squeeze(g(i,j,k,:)), 3, 3) / 100;
            M = eye(3) + L;
            b_scale_map(i,j,k) = trace(M' * M) / 3;
        end
    end
end

%% Step 3. Load white matter mask
wm_mask = niftiread(fullfile(data_dir, 'sub-001_mask_nocsfgm.nii.gz')) > 0;

%% Step 4. Compute resolution limit maps
% Following Nilsson et al., NMR in Biomed 2017;30:e3711
% and Hong-Hsi Lee's demo_resolution_limit.m

% Resolution limit of perpendicular diffusion signal (Eq. 15)
dmin_func = @(sbar, D, del, G) ...
    (768/7 * sbar .* D ./ del ./ (G/40*0.0107).^2).^(1/4);

% Resolution limit of spherical mean signal (Eq. 40)
dmin_SMT_func = @(sbar, D, Da, del, G, b) ...
    dmin_func(sbar, D, del, G) ./ (sqrt(pi/4./b./Da) .* erf(sqrt(b.*Da))).^(1/4);

% C2 parameters (from Hong-Hsi's demo_resolution_limit.m)
z_alpha = 1.64;         % one-sided 95% confidence
n_dir   = 64;           % number of gradient directions
D0      = 1.7;          % intra-axonal intrinsic diffusivity [um^2/ms]
Da      = 1.7;          % intra-axonal axial diffusivity [um^2/ms]
Gmax    = 495;          % C2 actual max gradient [mT/m]
Smax    = 600;          % C2 slew rate [T/m/s]
delta_vals = [6, 8];    % pulse durations [ms] — two shells

% b-max for each delta (using fixed bmax=6 as in Hong-Hsi's landscape plot)
bmax = 6;               % ms/um^2

% Voxel-wise sbar (noise threshold, Eq. 6)
sbar_map = z_alpha ./ tSNR_map / sqrt(n_dir);
sbar_map(isnan(sbar_map) | isinf(sbar_map)) = 0;

% --- Nominal resolution limit (no GNL correction) ---
dmin_nom = zeros([dims, numel(delta_vals)], 'single');
for id = 1:numel(delta_vals)
    dmin_nom(:,:,:,id) = dmin_SMT_func(sbar_map, D0, Da, delta_vals(id), Gmax, bmax);
end
dmin_nominal = mean(dmin_nom, 4);  % average over delta values
dmin_nominal(isnan(dmin_nominal) | isinf(dmin_nominal)) = 0;

% --- GNL-corrected resolution limit ---
% b_scale = (G_eff / G_nominal)^2, so G_eff = G_nominal * sqrt(b_scale)
% b_eff = b_nominal * b_scale
G_eff_map = Gmax * sqrt(b_scale_map);
b_eff_map = bmax * b_scale_map;

dmin_cor = zeros([dims, numel(delta_vals)], 'single');
for id = 1:numel(delta_vals)
    dmin_cor(:,:,:,id) = dmin_SMT_func(sbar_map, D0, Da, delta_vals(id), G_eff_map, b_eff_map);
end
dmin_corrected = mean(dmin_cor, 4);
dmin_corrected(isnan(dmin_corrected) | isinf(dmin_corrected)) = 0;

% Apply mask
dmin_nominal(~wm_mask)   = 0;
dmin_corrected(~wm_mask) = 0;

%% Step 5. Summary statistics
dmin_nom_wm = dmin_nominal(wm_mask);
dmin_cor_wm = dmin_corrected(wm_mask);
tsnr_wm     = tSNR_map(wm_mask);
bscale_wm   = b_scale_map(wm_mask);

fprintf('\n=== White matter summary ===\n');
fprintf('tSNR:  mean=%.1f, median=%.1f, range=[%.1f, %.1f]\n', ...
    mean(tsnr_wm), median(tsnr_wm), min(tsnr_wm), max(tsnr_wm));
fprintf('b_scale: mean=%.4f, median=%.4f, range=[%.4f, %.4f]\n', ...
    mean(bscale_wm), median(bscale_wm), min(bscale_wm), max(bscale_wm));
fprintf('\nResolution limit (nominal):   mean=%.3f, median=%.3f um\n', ...
    mean(dmin_nom_wm), median(dmin_nom_wm));
fprintf('Resolution limit (corrected): mean=%.3f, median=%.3f um\n', ...
    mean(dmin_cor_wm), median(dmin_cor_wm));
fprintf('Difference (corr - nom): mean=%.4f, median=%.4f um\n', ...
    mean(dmin_cor_wm - dmin_nom_wm), median(dmin_cor_wm - dmin_nom_wm));

%% Step 6. Plot maps
% Find axial slice with most WM voxels
wm_per_slice = squeeze(sum(sum(wm_mask, 1), 2));
[~, best_slice] = max(wm_per_slice);

figure('unit','inch','position',[0 0 16 10]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Row 1: tSNR, b_scale, tSNR in WM
nexttile;
imagesc(rot90(tSNR_map(:,:,best_slice))); axis image off; colorbar;
title(sprintf('tSNR (slice %d)', best_slice)); clim([0 60]);

nexttile;
tmp = b_scale_map(:,:,best_slice); tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(rot90(tmp)); axis image off; colorbar;
title('b-scale (WM)'); clim([0.8 1.2]);

nexttile;
tmp = tSNR_map(:,:,best_slice); tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(rot90(tmp)); axis image off; colorbar;
title('tSNR (WM)'); clim([0 60]);

% Row 2: resolution limit nominal, corrected, difference
nexttile;
tmp = dmin_nominal(:,:,best_slice); tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(rot90(tmp)); axis image off; colorbar;
title('Resolution limit — nominal (\mum)'); clim([0 5]);

nexttile;
tmp = dmin_corrected(:,:,best_slice); tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(rot90(tmp)); axis image off; colorbar;
title('Resolution limit — GNL corrected (\mum)'); clim([0 5]);

nexttile;
tmp = dmin_corrected(:,:,best_slice) - dmin_nominal(:,:,best_slice);
tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(rot90(tmp)); axis image off; cb = colorbar;
title('Difference (corrected - nominal, \mum)');
clim([-0.3 0.3]); colormap(gca, 'jet');

exportgraphics(gcf, 'fig_resolution_limit_maps.pdf', 'ContentType', 'vector');

%% Step 7. Histogram of resolution limits
figure('unit','inch','position',[0 0 12 4]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram(tsnr_wm, 50, 'FaceColor', [0.5 0.5 0.5]);
xlabel('tSNR'); ylabel('Count'); title('tSNR in WM');
xline(median(tsnr_wm), 'r--', sprintf('median=%.1f', median(tsnr_wm)), 'LineWidth', 1.5);

nexttile;
histogram(dmin_nom_wm, 50, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.6); hold on;
histogram(dmin_cor_wm, 50, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.6);
xlabel('Resolution limit (\mum)'); ylabel('Count');
title('Resolution limit in WM');
legend({'Nominal', 'GNL corrected'});

nexttile;
diff_wm = dmin_cor_wm - dmin_nom_wm;
histogram(diff_wm, 50, 'FaceColor', [0.3 0.7 0.3]);
xlabel('\Delta resolution limit (\mum)'); ylabel('Count');
title('Difference (corrected - nominal)');
xline(0, 'k--', 'HandleVisibility', 'off');
xline(mean(diff_wm), 'r--', sprintf('mean=%.4f', mean(diff_wm)), 'LineWidth', 1.5);

exportgraphics(gcf, 'fig_resolution_limit_histograms.pdf', 'ContentType', 'vector');
