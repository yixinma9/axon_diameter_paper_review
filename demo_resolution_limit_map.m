clear; close all;
addpath(genpath('/autofs/space/linen_001/users/Yixin/C2_protocoldesign-main/lib'));

%% Step 1. Compute temporal SNR from b0s.nii.gz (first 10 volumes)
data_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001';

b0_all = single(niftiread(fullfile(data_dir, 'b0s.nii.gz')));
b0 = b0_all(:,:,:,1:min(10, size(b0_all,4)));
Nvol = size(b0, 4);
fprintf('Using first %d b0 volumes (out of %d) for tSNR\n', Nvol, size(b0_all,4));
clear b0_all;

% tSNR = mean / std across repeats
tSNR_map = mean(b0, 4) ./ std(b0, 0, 4);
tSNR_map(isnan(tSNR_map) | isinf(tSNR_map)) = 0;

%% Step 2. Load white matter mask
wm_mask = niftiread(fullfile(data_dir, 'sub-001_mask_nocsfgm.nii.gz')) > 0;
dims = size(tSNR_map, 1:3);

%% Step 3. Compute resolution limit map
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
delta_vals = [6, 8];    % pulse durations [ms] — two shells
bmax = 6;               % ms/um^2

% Voxel-wise sbar (noise threshold, Eq. 6)
sbar_map = z_alpha ./ tSNR_map / sqrt(n_dir);
sbar_map(isnan(sbar_map) | isinf(sbar_map)) = 0;

% Resolution limit map (average over two delta values)
dmin_all = zeros([dims, numel(delta_vals)], 'single');
for id = 1:numel(delta_vals)
    dmin_all(:,:,:,id) = dmin_SMT_func(sbar_map, D0, Da, delta_vals(id), Gmax, bmax);
end
dmin_map = mean(dmin_all, 4);
dmin_map(isnan(dmin_map) | isinf(dmin_map)) = 0;
dmin_map(~wm_mask) = 0;

%% Step 4. Summary statistics
dmin_wm = dmin_map(wm_mask);
tsnr_wm = tSNR_map(wm_mask);

fprintf('\n=== White matter summary ===\n');
fprintf('tSNR:  mean=%.1f, median=%.1f, range=[%.1f, %.1f]\n', ...
    mean(tsnr_wm), median(tsnr_wm), min(tsnr_wm), max(tsnr_wm));
fprintf('Resolution limit: mean=%.3f, median=%.3f um\n', ...
    mean(dmin_wm), median(dmin_wm));

%% Step 5. Plot maps
% Find axial slice with most WM voxels
wm_per_slice = squeeze(sum(sum(wm_mask, 1), 2));
[~, best_slice] = max(wm_per_slice);

figure('unit','inch','position',[0 0 15 5]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(flipud(rot90(tSNR_map(:,:,best_slice),2))); axis image off; colorbar;
title(sprintf('tSNR (slice %d)', best_slice)); clim([0 60]);

nexttile;
tmp = tSNR_map(:,:,best_slice); tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(flipud(rot90(tmp,2))); axis image off; colorbar;
title('tSNR (WM)'); clim([0 60]);

nexttile;
tmp = dmin_map(:,:,best_slice); tmp(~wm_mask(:,:,best_slice)) = NaN;
imagesc(flipud(rot90(tmp,2))); axis image off; colorbar;
title('Resolution limit (\mum)'); clim([0 5]);

exportgraphics(gcf, 'fig_resolution_limit_maps.pdf', 'ContentType', 'vector');

%% Step 6. Histograms
figure('unit','inch','position',[0 0 10 4]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram(tsnr_wm, 50, 'FaceColor', [0.5 0.5 0.5]);
xlabel('tSNR'); ylabel('Count'); title('tSNR in WM');
xline(median(tsnr_wm), 'r--', sprintf('median=%.1f', median(tsnr_wm)), 'LineWidth', 1.5);

nexttile;
histogram(dmin_wm, 50, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Resolution limit (\mum)'); ylabel('Count');
title('Resolution limit in WM');
xline(median(dmin_wm), 'r--', sprintf('median=%.3f \\mum', median(dmin_wm)), 'LineWidth', 1.5);

exportgraphics(gcf, 'fig_resolution_limit_histograms.pdf', 'ContentType', 'vector');
