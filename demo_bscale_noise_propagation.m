clear; close all;
addpath(genpath('/autofs/space/linen_001/users/Yixin/C2_protocoldesign-main/lib'));

%% Step 1. Load b_scale map and white matter mask
gnc_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc';
mask_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001';

% Load grad_dev (QT's method) — 3 separate files, each [x,y,z,3]
gx = single(niftiread(fullfile(gnc_dir, 'grad_dev_x.nii.gz')));
gy = single(niftiread(fullfile(gnc_dir, 'grad_dev_y.nii.gz')));
gz = single(niftiread(fullfile(gnc_dir, 'grad_dev_z.nii.gz')));
g  = cat(4, gx, gy, gz);  % [x,y,z,9]

% Compute b_scale per voxel (direction-averaged n^2)
% Following QT: L = reshape(squeeze(g(i,j,k,:)),3,3); v = (I+L)*bvecs; n^2 = sum(v.^2)
% For spherical mean: b_scale = mean(n^2) over directions = trace((I+L)'*(I+L))/3
dims = size(g, 1:3);
b_scale_map = zeros(dims, 'single');
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            L = reshape(squeeze(g(i,j,k,:)), 3, 3) / 100;  % calc_grad_perc_dev outputs percentage
            M = eye(3) + L;
            b_scale_map(i,j,k) = trace(M' * M) / 3;
        end
    end
end

% Load white matter mask
wm_mask = niftiread(fullfile(mask_dir, 'sub-001_mask_nocsfgm.nii.gz')) > 0;

% Extract b_scale values in white matter
b_scale_wm = b_scale_map(wm_mask);
b_scale_wm = b_scale_wm(b_scale_wm > 0);  % remove zeros
fprintf('WM voxels: %d\n', numel(b_scale_wm));
fprintf('b_scale range: [%.4f, %.4f]\n', min(b_scale_wm), max(b_scale_wm));

%% Step 2. Plot b_scale distribution and split into 2 groups (< 1 and > 1)
figure('unit','inch','position',[0 0 10 4]);
histogram(b_scale_wm, 100, 'FaceColor', [0.7 0.7 0.7]); hold on;
xline(1.0, 'r--', 'LineWidth', 2);
xlabel('b-scale'); ylabel('count');
title('b-scale distribution in white matter');
legend({'b-scale', 'b-scale = 1.0'});

group1_idx = b_scale_wm < 1;
group2_idx = b_scale_wm > 1;
fprintf('Group 1 (b_scale < 1): %d voxels\n', sum(group1_idx));
fprintf('Group 2 (b_scale > 1): %d voxels\n', sum(group2_idx));

%% Step 3. Sample 2000 b_scale values per group (following the distribution)
Nsample = 2000;
rng(42);

vals1 = double(b_scale_wm(group1_idx)); vals1 = vals1(randperm(numel(vals1), min(Nsample, numel(vals1))));
vals2 = double(b_scale_wm(group2_idx)); vals2 = vals2(randperm(numel(vals2), min(Nsample, numel(vals2))));

Ngroups = 2;
b_scale_samples = {vals1, vals2};
group_names = {sprintf('b-scale < 1 (%.3f-%.3f)', min(vals1), max(vals1)), ...
               sprintf('b-scale > 1 (%.3f-%.3f)', min(vals2), max(vals2))};

%% Step 4. Setup AxCaliberSMT protocol (TractCaliber)
b = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
     0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta = [8*ones(8,1); 8*ones(8,1)];
Delta = [19*ones(8,1); 49*ones(8,1)];

D0   = 2;
Da   = 1.7;
DeL  = 1.7;
Dcsf = 3;
snr  = 30;
model = 'VanGelderen';

% Random ground truth parameters (shared across all groups)
rng(0);
Xrange = [0.1 5; 0.5 1; 0 0.2; 0.5 1.5];
Xinit = Xrange(:,1) + diff(Xrange,1,2) .* rand(4, Nsample);
r_true    = Xinit(1,:)';
f_true    = Xinit(2,:)';
fcsf_true = Xinit(3,:)';
DeR_true  = Xinit(4,:)';

%% Step 5. Noise propagation for each group
% For each group:
%   - Generate signal with b_eff = b * b_scale (ground truth)
%   - Fit with nominal b (uncorrected)
%   - Fit with b_eff (corrected)

results = struct();

for ig = 1:Ngroups
    fprintf('\n=== %s ===\n', group_names{ig});
    bscales = b_scale_samples{ig};
    N = numel(bscales);

    r_fit_uncorr = zeros(N, 1);
    f_fit_uncorr = zeros(N, 1);
    fcsf_fit_uncorr = zeros(N, 1);
    DeR_fit_uncorr = zeros(N, 1);

    r_fit_corr = zeros(N, 1);
    f_fit_corr = zeros(N, 1);
    fcsf_fit_corr = zeros(N, 1);
    DeR_fit_corr = zeros(N, 1);

    % Nominal AxCaliberSMT object (for uncorrected fitting)
    smt_nominal = AxCaliberSMT(b, delta, Delta, D0, Da, DeL, Dcsf);

    tic;
    for i = 1:N
        if mod(i, 100) == 0; fprintf('  voxel %d/%d\n', i, N); end

        % Effective b-value for this voxel
        b_eff = b * bscales(i);

        % Generate ground truth signal using effective b-values
        smt_true = AxCaliberSMT(b_eff, delta, Delta, D0, Da, DeL, Dcsf);
        S_clean = smt_true.AxonDiameterFWD([r_true(i), f_true(i), fcsf_true(i), DeR_true(i)], model);

        % Add Rician noise
        S_noisy = double(abs(S_clean + 1/snr*randn(size(S_clean)) + 1j/snr*randn(size(S_clean))));

        % Fit with nominal b (uncorrected)
        [r_fit_uncorr(i), f_fit_uncorr(i), fcsf_fit_uncorr(i), DeR_fit_uncorr(i)] = ...
            smt_nominal.mcmc(S_noisy, model, 5);

        % Fit with corrected b
        smt_corr = AxCaliberSMT(b_eff, delta, Delta, D0, Da, DeL, Dcsf);
        [r_fit_corr(i), f_fit_corr(i), fcsf_fit_corr(i), DeR_fit_corr(i)] = ...
            smt_corr.mcmc(S_noisy, model, 5);
    end
    t = toc;
    fprintf('  Elapsed: %.1f s\n', t);

    results(ig).group_name = group_names{ig};
    results(ig).bscales = bscales;
    results(ig).r_true = r_true(1:N);
    results(ig).r_fit_uncorr = r_fit_uncorr;
    results(ig).r_fit_corr = r_fit_corr;
    results(ig).f_fit_uncorr = f_fit_uncorr;
    results(ig).f_fit_corr = f_fit_corr;
    results(ig).fcsf_fit_uncorr = fcsf_fit_uncorr;
    results(ig).fcsf_fit_corr = fcsf_fit_corr;
    results(ig).DeR_fit_uncorr = DeR_fit_uncorr;
    results(ig).DeR_fit_corr = DeR_fit_corr;
end

%% Step 6. Quantitative evaluation — stratified by axon radius
r_bins = [0 1; 1 2; 2 3; 3 4];
Nbins = size(r_bins, 1);

for ig = 1:Ngroups
    fprintf('\n========== %s ==========\n', results(ig).group_name);
    fprintf('%-15s %6s %10s %10s %10s %10s\n', 'r range', 'N', 'Bias(unc)', 'Bias(cor)', 'RMSE(unc)', 'RMSE(cor)');

    for ib = 1:Nbins
        idx = results(ig).r_true >= r_bins(ib,1) & results(ig).r_true < r_bins(ib,2);
        if sum(idx) == 0; continue; end

        err_unc = results(ig).r_fit_uncorr(idx) - results(ig).r_true(idx);
        err_cor = results(ig).r_fit_corr(idx)   - results(ig).r_true(idx);

        fprintf('[%g-%g] um %6d %+10.4f %+10.4f %10.4f %10.4f\n', ...
                r_bins(ib,1), r_bins(ib,2), sum(idx), ...
                mean(err_unc), mean(err_cor), sqrt(mean(err_unc.^2)), sqrt(mean(err_cor.^2)));
    end

    % Overall
    err_unc = results(ig).r_fit_uncorr - results(ig).r_true;
    err_cor = results(ig).r_fit_corr   - results(ig).r_true;
    fprintf('%-15s %6d %+10.4f %+10.4f %10.4f %10.4f\n', 'Overall', ...
            numel(err_unc), mean(err_unc), mean(err_cor), ...
            sqrt(mean(err_unc.^2)), sqrt(mean(err_cor.^2)));

    results(ig).bias_uncorr = mean(err_unc);
    results(ig).bias_corr   = mean(err_cor);
    results(ig).rmse_uncorr = sqrt(mean(err_unc.^2));
    results(ig).rmse_corr   = sqrt(mean(err_cor.^2));
end

%% Step 7. Plot: scatter (2 groups x 2 methods) with metrics in title
figure('unit','inch','position',[0 0 10 8]);
for ig = 1:Ngroups
    err_unc = results(ig).r_fit_uncorr - results(ig).r_true;
    err_cor = results(ig).r_fit_corr   - results(ig).r_true;

    % Uncorrected
    subplot(2, Ngroups, ig);
    scatter(results(ig).r_true, results(ig).r_fit_uncorr, 10, '.');
    hold on; hr = refline(1,0); set(hr, 'color', 'k');
    xlabel('Ground truth r (\mum)'); ylabel('Fitted r (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s\nUncorrected | bias=%.3f, RMSE=%.3f', ...
          results(ig).group_name, mean(err_unc), sqrt(mean(err_unc.^2))));

    % Corrected
    subplot(2, Ngroups, ig + Ngroups);
    scatter(results(ig).r_true, results(ig).r_fit_corr, 10, '.');
    hold on; hr = refline(1,0); set(hr, 'color', 'k');
    xlabel('Ground truth r (\mum)'); ylabel('Fitted r (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s\nCorrected | bias=%.3f, RMSE=%.3f', ...
          results(ig).group_name, mean(err_cor), sqrt(mean(err_cor.^2))));
end

%% Step 8. Bar plot: bias & RMSE by r-bin for each group
figure('unit','inch','position',[0 0 14 8]);
for ig = 1:Ngroups
    bias_unc_bins = zeros(Nbins, 1);
    bias_cor_bins = zeros(Nbins, 1);
    rmse_unc_bins = zeros(Nbins, 1);
    rmse_cor_bins = zeros(Nbins, 1);

    for ib = 1:Nbins
        idx = results(ig).r_true >= r_bins(ib,1) & results(ig).r_true < r_bins(ib,2);
        if sum(idx) == 0; continue; end
        err_unc = results(ig).r_fit_uncorr(idx) - results(ig).r_true(idx);
        err_cor = results(ig).r_fit_corr(idx)   - results(ig).r_true(idx);
        bias_unc_bins(ib) = mean(err_unc);
        bias_cor_bins(ib) = mean(err_cor);
        rmse_unc_bins(ib) = sqrt(mean(err_unc.^2));
        rmse_cor_bins(ib) = sqrt(mean(err_cor.^2));
    end

    bin_labels = arrayfun(@(i) sprintf('%g-%g', r_bins(i,1), r_bins(i,2)), 1:Nbins, 'uni', 0);

    % Bias
    subplot(2, Ngroups, ig);
    bar(categorical(bin_labels, bin_labels), [bias_unc_bins, bias_cor_bins]);
    ylabel('Bias (\mum)'); xlabel('r true (\mum)');
    legend({'Uncorrected', 'Corrected'}, 'Location', 'best');
    title(sprintf('%s — Bias', results(ig).group_name));
    yline(0, 'k--'); grid on;

    % RMSE
    subplot(2, Ngroups, ig + Ngroups);
    bar(categorical(bin_labels, bin_labels), [rmse_unc_bins, rmse_cor_bins]);
    ylabel('RMSE (\mum)'); xlabel('r true (\mum)');
    legend({'Uncorrected', 'Corrected'}, 'Location', 'best');
    title(sprintf('%s — RMSE', results(ig).group_name));
    grid on;
end

%% Step 9. Plot: error distribution (corrected vs uncorrected per group)
figure('unit','inch','position',[0 0 10 4]);
for ig = 1:Ngroups
    subplot(1, Ngroups, ig);
    err_unc = results(ig).r_fit_uncorr - results(ig).r_true;
    err_cor = results(ig).r_fit_corr   - results(ig).r_true;
    histogram(err_unc, 50, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.5); hold on;
    histogram(err_cor, 50, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.5);
    xline(mean(err_unc), 'r--', 'LineWidth', 2);
    xline(mean(err_cor), 'b--', 'LineWidth', 2);
    xline(0, 'k-', 'LineWidth', 1);
    xlabel('Fitted r - True r (\mum)'); ylabel('Count');
    title(sprintf('%s\nmean b-scale=%.3f', results(ig).group_name, mean(results(ig).bscales)));
    legend({'Uncorrected', 'Corrected', ...
            sprintf('Bias=%.3f', mean(err_unc)), sprintf('Bias=%.3f', mean(err_cor))});
end

%% Save results
save('bscale_noise_propagation_results.mat', 'results', 'Xrange', 'b_scale_samples', 'group_names', ...
     'r_true', 'f_true', 'fcsf_true', 'DeR_true', 'r_bins');
