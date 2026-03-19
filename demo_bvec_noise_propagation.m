function demo_bvec_noise_propagation()
clear; close all;

%% b-vector GNL noise propagation (full L-matrix: scaling + rotation)
%
% Same structure as demo_bscale_noise_propagation.m, but signal generation
% uses the directional model with the full per-voxel L matrix, capturing
% both b-value scaling AND b-vector rotation from gradient nonlinearity.
%
% For each sampled WM voxel:
%   - Reconstruct its 3x3 L matrix from grad_dev
%   - Generate signal per gradient direction using VG + zeppelin, average per shell
%   - Add Rician noise
%   - Fit with nominal b (uncorrected) and b*b_scale (corrected)
%
% This is slower than demo_bscale_noise_propagation because signal
% generation loops over Ndirs directions per shell per voxel.

%% Step 1. Load grad_dev and white matter mask
gnc_dir  = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc';
mask_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001';

% Load grad_dev (QT's method) — 3 separate files, each [x,y,z,3]
gx = single(niftiread(fullfile(gnc_dir, 'grad_dev_x.nii.gz')));
gy = single(niftiread(fullfile(gnc_dir, 'grad_dev_y.nii.gz')));
gz = single(niftiread(fullfile(gnc_dir, 'grad_dev_z.nii.gz')));
g  = cat(4, gx, gy, gz);  % [x,y,z,9]
clear gx gy gz;

% Compute b_scale per voxel: trace((I+L)'*(I+L))/3
dims = size(g, 1:3);
b_scale_map = zeros(dims, 'single');
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            L = reshape(squeeze(g(i,j,k,:)), 3, 3) / 100;  % percentage to fraction
            M = eye(3) + L;
            b_scale_map(i,j,k) = trace(M' * M) / 3;
        end
    end
end

% Load white matter mask
wm_mask = niftiread(fullfile(mask_dir, 'sub-001_mask_nocsfgm.nii.gz')) > 0;

% Extract WM voxel data: b_scale + full L matrices
wm_idx   = find(wm_mask(:));
N_wm     = numel(wm_idx);
g_flat   = reshape(g, [], 9);               % [Nvox_total x 9]
L_wm     = double(g_flat(wm_idx, :)) / 100; % [N_wm x 9], fraction
b_scale_wm = b_scale_map(wm_idx);           % [N_wm x 1]
b_scale_wm(b_scale_wm == 0) = [];
clear g g_flat;

fprintf('WM voxels: %d\n', N_wm);
fprintf('b_scale range: [%.4f, %.4f]\n', min(b_scale_wm), max(b_scale_wm));

%% Step 2. Plot b_scale distribution and split into 2 groups (< 1 and > 1)
figure('unit','inch','position',[0 0 10 4]);
histogram(b_scale_wm, 100, 'FaceColor', [0.7 0.7 0.7]); hold on;
xline(1.0, 'r--', 'LineWidth', 2);
xlabel('b-scale'); ylabel('count');
title('b-scale distribution in white matter');
legend({'b-scale', 'b-scale = 1.0'});
exportgraphics(gcf, 'fig1_bvec_bscale_distribution.pdf', 'ContentType', 'vector');

group1_mask = b_scale_map(wm_idx) < 1;
group2_mask = b_scale_map(wm_idx) > 1;
fprintf('Group 1 (b_scale < 1): %d voxels\n', sum(group1_mask));
fprintf('Group 2 (b_scale > 1): %d voxels\n', sum(group2_mask));

%% Step 3. Sample 2000 voxels per group (with full L matrices)
Nsample = 2000;
rng(42);

idx1 = find(group1_mask); idx1 = idx1(randperm(numel(idx1), min(Nsample, numel(idx1))));
idx2 = find(group2_mask); idx2 = idx2(randperm(numel(idx2), min(Nsample, numel(idx2))));

Ngroups = 2;
L_samples       = {L_wm(idx1,:), L_wm(idx2,:)};       % each [Nsample x 9]
b_scale_samples = {b_scale_map(wm_idx(idx1)), b_scale_map(wm_idx(idx2))};

vals1 = double(b_scale_samples{1});
vals2 = double(b_scale_samples{2});
group_names = {sprintf('b-scale < 1 (%.3f-%.3f)', min(vals1), max(vals1)), ...
               sprintf('b-scale > 1 (%.3f-%.3f)', min(vals2), max(vals2))};

%% Step 4. Setup protocol and gradient directions
b = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
     0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta = [8*ones(8,1); 8*ones(8,1)];
Delta = [19*ones(8,1); 49*ones(8,1)];
Nsh   = numel(b);

D0   = 2;
Da   = 1.7;
DeL  = 1.7;
Dcsf = 3;
snr  = 30;

g_mag = sqrt(b ./ (delta.^2 .* (Delta - delta/3)));   % gradient magnitude per shell
g2    = b ./ (delta.^2 .* (Delta - delta/3));

% Van Gelderen Bessel zeros
bm  = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

% Gradient directions from sub-001.bvec (16 representative directions)
bvecs = [-0.5899516186,  0.9725294444,  0.4117928911,  0.1356963446, ...
          0.7850186249,  0.8149950102,  0.5622431322,  0.0866392126, ...
          0.4973663712,  0.8953048492, -0.0273308969,  0.2337241520, ...
         -0.5618138772, -0.7473683543,  0.3695585598, -0.2369838544; ...
         -0.6903964370, -0.0619626359, -0.3605186270,  0.7493190415, ...
         -0.1336973403, -0.3552195534,  0.7208258343,  0.0180902564, ...
          0.1228170165, -0.4451790664, -0.4351433515,  0.9562599154, ...
          0.3376072842, -0.3371839958, -0.6891298219,  0.9555028625; ...
          0.4187001880,  0.2243816202,  0.8369306629, -0.6481569841, ...
          0.6048725318, -0.4578233308,  0.4053304539,  0.9960754939, ...
          0.8588030468,  0.0156469120,  0.8999462683,  0.1758976830, ...
          0.7552393588, -0.5724923545, -0.6233189869,  0.1756500284];
Ndirs = size(bvecs, 2);

fiber_dir = [0; 0; 1];   % fiber along z-axis

% Random ground truth parameters (shared across all groups, same seed as bscale version)
rng(0);
Xrange = [0.1 5; 0.5 1; 0 0.2; 0.5 1.5];
Xinit = Xrange(:,1) + diff(Xrange,1,2) .* rand(4, Nsample);
r_true    = Xinit(1,:)';
f_true    = Xinit(2,:)';
fcsf_true = Xinit(3,:)';
DeR_true  = Xinit(4,:)';

%% Step 5. Noise propagation (directional signal model with full L matrix)
% KEY DIFFERENCE from demo_bscale_noise_propagation:
%   - Signal generation uses per-direction VG + zeppelin model with full L matrix
%   - This captures both b-value scaling AND b-vector rotation from GNL
%   - Fitting still uses the standard SMT powder-average model (lsqnonlin)

results = struct();

for ig = 1:Ngroups
    fprintf('\n=== %s ===\n', group_names{ig});
    L_grp   = L_samples{ig};
    bscales = double(b_scale_samples{ig});
    N = size(L_grp, 1);

    r_fit_uncorr = zeros(N, 1);
    f_fit_uncorr = zeros(N, 1);
    fcsf_fit_uncorr = zeros(N, 1);
    DeR_fit_uncorr = zeros(N, 1);

    r_fit_corr = zeros(N, 1);
    f_fit_corr = zeros(N, 1);
    fcsf_fit_corr = zeros(N, 1);
    DeR_fit_corr = zeros(N, 1);

    tic;
    for i = 1:N
        if mod(i, 100) == 0; fprintf('  voxel %d/%d (%.0fs)\n', i, N, toc); end

        % This voxel's full L matrix
        L = reshape(L_grp(i,:), 3, 3);
        M = eye(3) + L;
        b_scale_scalar = trace(M'*M) / 3;

        % Per-direction effective gradient quantities
        g_eff_all       = M * bvecs;                       % [3 x Ndirs]
        b_scale_per_dir = sum(g_eff_all.^2, 1);            % |M*g_hat|^2 [1 x Ndirs]
        g_eff_hat       = g_eff_all ./ vecnorm(g_eff_all); % unit effective dirs
        cos_theta       = abs(g_eff_hat' * fiber_dir);     % [Ndirs x 1]
        sin2_theta      = 1 - cos_theta.^2;

        % ---- Generate signal: directional model, averaged per shell ----
        S_clean = zeros(Nsh, 1);
        for is = 1:Nsh
            Sd = zeros(Ndirs, 1);
            for id = 1:Ndirs
                b_eff = b(is) * b_scale_per_dir(id);
                g_eff = g_mag(is) * sqrt(b_scale_per_dir(id));
                g2p   = g_eff^2 * sin2_theta(id);

                C_r = vg_atten(r_true(i), g2p, delta(is), Delta(is), D0, bm2);
                Sa  = exp(-b_eff * Da * cos_theta(id)^2 - C_r);
                Se  = exp(-b_eff * (DeR_true(i) + (DeL - DeR_true(i)) * cos_theta(id)^2));
                Sc  = exp(-b_eff * Dcsf);
                Sd(id) = (1-fcsf_true(i)) * (f_true(i)*Sa + (1-f_true(i))*Se) + fcsf_true(i)*Sc;
            end
            S_clean(is) = mean(Sd);
        end

        % Add Rician noise
        S_noisy = double(abs(S_clean + 1/snr*randn(Nsh,1) + 1j/snr*randn(Nsh,1)));

        % Fit with nominal b (uncorrected)
        p = fit_smt(b, g2, delta, Delta, Da, DeL, D0, Dcsf, bm2, S_noisy);
        r_fit_uncorr(i) = p.r; f_fit_uncorr(i) = p.f;
        fcsf_fit_uncorr(i) = p.fcsf; DeR_fit_uncorr(i) = p.DeR;

        % Fit with corrected b (scalar b_scale)
        b_corr  = b * b_scale_scalar;
        g2_corr = b_corr ./ (delta.^2 .* (Delta - delta/3));
        p = fit_smt(b_corr, g2_corr, delta, Delta, Da, DeL, D0, Dcsf, bm2, S_noisy);
        r_fit_corr(i) = p.r; f_fit_corr(i) = p.f;
        fcsf_fit_corr(i) = p.fcsf; DeR_fit_corr(i) = p.DeR;
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

%% Step 6. Filter out boundary hits and quantitative evaluation
% Fit bounds: r in [0, 5]. Exclude estimates hitting boundaries.
r_lb = 0.05;  r_ub = 4.95;

for ig = 1:Ngroups
    valid_unc = results(ig).r_fit_uncorr > r_lb & results(ig).r_fit_uncorr < r_ub;
    valid_cor = results(ig).r_fit_corr   > r_lb & results(ig).r_fit_corr   < r_ub;
    valid = valid_unc & valid_cor;
    fprintf('\n%s: %d/%d valid (excluded %d boundary hits)\n', ...
            results(ig).group_name, sum(valid), numel(valid), sum(~valid));

    results(ig).r_true_filt       = results(ig).r_true(valid);
    results(ig).r_fit_uncorr_filt = results(ig).r_fit_uncorr(valid);
    results(ig).r_fit_corr_filt   = results(ig).r_fit_corr(valid);
    results(ig).bscales_filt      = results(ig).bscales(valid);
end

%% Step 6b. Quantitative evaluation — stratified by axon radius (finer bins)
r_bins = [0 0.5; 0.5 1; 1 1.5; 1.5 2; 2 3; 3 4; 4 5];
Nbins = size(r_bins, 1);

for ig = 1:Ngroups
    rt = results(ig).r_true_filt;
    ru = results(ig).r_fit_uncorr_filt;
    rc = results(ig).r_fit_corr_filt;

    fprintf('\n========== %s (filtered) ==========\n', results(ig).group_name);
    fprintf('%-15s %6s %10s %10s %10s %10s %10s %10s\n', ...
            'r range', 'N', 'Bias(unc)', 'Bias(cor)', 'Med(unc)', 'Med(cor)', 'RMSE(unc)', 'RMSE(cor)');

    for ib = 1:Nbins
        idx = rt >= r_bins(ib,1) & rt < r_bins(ib,2);
        if sum(idx) == 0; continue; end

        err_unc = ru(idx) - rt(idx);
        err_cor = rc(idx) - rt(idx);

        fprintf('[%.1f-%.1f] um %5d %+10.4f %+10.4f %+10.4f %+10.4f %10.4f %10.4f\n', ...
                r_bins(ib,1), r_bins(ib,2), sum(idx), ...
                mean(err_unc), mean(err_cor), median(err_unc), median(err_cor), ...
                sqrt(mean(err_unc.^2)), sqrt(mean(err_cor.^2)));
    end

    % Overall
    err_unc = ru - rt;
    err_cor = rc - rt;
    fprintf('%-15s %5d %+10.4f %+10.4f %+10.4f %+10.4f %10.4f %10.4f\n', 'Overall', ...
            numel(err_unc), mean(err_unc), mean(err_cor), ...
            median(err_unc), median(err_cor), ...
            sqrt(mean(err_unc.^2)), sqrt(mean(err_cor.^2)));

    results(ig).bias_uncorr = mean(err_unc);
    results(ig).bias_corr   = mean(err_cor);
    results(ig).rmse_uncorr = sqrt(mean(err_unc.^2));
    results(ig).rmse_corr   = sqrt(mean(err_cor.^2));
end

%% Step 6c. Scatter: error vs r_true (filtered)
figure('unit','inch','position',[0 0 12 5]);
for ig = 1:Ngroups
    subplot(1, Ngroups, ig);
    rt = results(ig).r_true_filt;
    err_unc = results(ig).r_fit_uncorr_filt - rt;
    err_cor = results(ig).r_fit_corr_filt   - rt;
    scatter(rt, err_unc, 8, [0.8 0.3 0.3], '.', 'MarkerFaceAlpha', 0.3); hold on;
    scatter(rt, err_cor, 8, [0.3 0.3 0.8], '.', 'MarkerFaceAlpha', 0.3);
    yline(0, 'k-', 'LineWidth', 1);
    xlabel('r true (\mum)'); ylabel('r fit - r true (\mum)');
    title(sprintf('%s (filtered)', results(ig).group_name));
    legend({'Uncorrected', 'Corrected'}, 'Location', 'best');
    grid on;
end
exportgraphics(gcf, 'fig2_bvec_error_vs_rtrue.pdf', 'ContentType', 'vector');

%% Step 7. Plot: scatter (2 groups x 2 methods, filtered) with metrics in title
figure('unit','inch','position',[0 0 10 8]);
for ig = 1:Ngroups
    rt = results(ig).r_true_filt;
    ru = results(ig).r_fit_uncorr_filt;
    rc = results(ig).r_fit_corr_filt;
    err_unc = ru - rt;
    err_cor = rc - rt;

    % Uncorrected
    subplot(2, Ngroups, ig);
    scatter(rt, ru, 10, '.');
    hold on; hr = refline(1,0); set(hr, 'color', 'k');
    xlabel('Ground truth r (\mum)'); ylabel('Fitted r (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s\nUncorrected | bias=%.3f, RMSE=%.3f', ...
          results(ig).group_name, mean(err_unc), sqrt(mean(err_unc.^2))));

    % Corrected
    subplot(2, Ngroups, ig + Ngroups);
    scatter(rt, rc, 10, '.');
    hold on; hr = refline(1,0); set(hr, 'color', 'k');
    xlabel('Ground truth r (\mum)'); ylabel('Fitted r (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s\nCorrected | bias=%.3f, RMSE=%.3f', ...
          results(ig).group_name, mean(err_cor), sqrt(mean(err_cor.^2))));
end
exportgraphics(gcf, 'fig3_bvec_scatter_gt_vs_fit.pdf', 'ContentType', 'vector');

%% Step 8. Bar plot: bias & RMSE by r-bin for each group (filtered)
r_bins_bar = [0 1; 1 2; 2 3; 3 4; 4 5];
Nbins_bar = size(r_bins_bar, 1);

figure('unit','inch','position',[0 0 14 8]);
for ig = 1:Ngroups
    rt = results(ig).r_true_filt;
    ru = results(ig).r_fit_uncorr_filt;
    rc = results(ig).r_fit_corr_filt;

    bias_unc_bins = zeros(Nbins_bar, 1);
    bias_cor_bins = zeros(Nbins_bar, 1);
    rmse_unc_bins = zeros(Nbins_bar, 1);
    rmse_cor_bins = zeros(Nbins_bar, 1);

    for ib = 1:Nbins_bar
        idx = rt >= r_bins_bar(ib,1) & rt < r_bins_bar(ib,2);
        if sum(idx) == 0; continue; end
        err_unc = ru(idx) - rt(idx);
        err_cor = rc(idx) - rt(idx);
        bias_unc_bins(ib) = mean(err_unc);
        bias_cor_bins(ib) = mean(err_cor);
        rmse_unc_bins(ib) = sqrt(mean(err_unc.^2));
        rmse_cor_bins(ib) = sqrt(mean(err_cor.^2));
    end

    bin_labels = arrayfun(@(i) sprintf('%g-%g', r_bins_bar(i,1), r_bins_bar(i,2)), 1:Nbins_bar, 'uni', 0);

    % Bias
    subplot(2, Ngroups, ig);
    bar(categorical(bin_labels, bin_labels), [bias_unc_bins, bias_cor_bins]);
    ylabel('Bias (\mum)'); xlabel('r true (\mum)');
    legend({'Uncorrected', 'Corrected'}, 'Location', 'best');
    title(sprintf('%s — Bias (filtered)', results(ig).group_name));
    yline(0, 'k--', 'HandleVisibility', 'off'); grid on;

    % RMSE
    subplot(2, Ngroups, ig + Ngroups);
    bar(categorical(bin_labels, bin_labels), [rmse_unc_bins, rmse_cor_bins]);
    ylabel('RMSE (\mum)'); xlabel('r true (\mum)');
    legend({'Uncorrected', 'Corrected'}, 'Location', 'best');
    title(sprintf('%s — RMSE (filtered)', results(ig).group_name));
    grid on;
end
exportgraphics(gcf, 'fig4_bvec_bar_bias_rmse.pdf', 'ContentType', 'vector');

%% Step 9. Plot: error distribution (filtered)
figure('unit','inch','position',[0 0 10 4]);
for ig = 1:Ngroups
    subplot(1, Ngroups, ig);
    rt = results(ig).r_true_filt;
    err_unc = results(ig).r_fit_uncorr_filt - rt;
    err_cor = results(ig).r_fit_corr_filt   - rt;
    histogram(err_unc, 50, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.5); hold on;
    histogram(err_cor, 50, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.5);
    xline(mean(err_unc), 'r--', 'LineWidth', 2);
    xline(mean(err_cor), 'b--', 'LineWidth', 2);
    xline(0, 'k-', 'LineWidth', 1);
    xlabel('Fitted r - True r (\mum)'); ylabel('Count');
    title(sprintf('%s (filtered)\nmean b-scale=%.3f', results(ig).group_name, mean(results(ig).bscales_filt)));
    legend({'Uncorrected', 'Corrected', ...
            sprintf('Bias=%.3f', mean(err_unc)), sprintf('Bias=%.3f', mean(err_cor))});
end
exportgraphics(gcf, 'fig5_bvec_error_distribution.pdf', 'ContentType', 'vector');

%% Save results
save('bvec_noise_propagation_results.mat', 'results', 'Xrange', 'b_scale_samples', 'group_names', ...
     'r_true', 'f_true', 'fcsf_true', 'DeR_true', 'r_bins');
fprintf('\nDone. Results saved to bvec_noise_propagation_results.mat\n');

end  % main function

%% ======================================================================
%  Helper functions
%  ======================================================================

function C = vg_atten(r, g2_perp, delta, Delta, D0, bm2)
    % Van Gelderen restricted radial attenuation for a cylinder of radius r
    bm2 = [3.3899 28.4238 72.9367 137.0305 221.0266 324.5583 447.9337 591.1395 753.4878 935.6762];
    td = r^2 / D0;
    C = 0;
    for k = 1:numel(bm2)
        bk  = bm2(k);
        bkd = bk * delta / td;
        bkD = bk * Delta / td;
        C = C + (2 / (bk^3 * (bk - 1))) * ...
            (-2 + 2*bkd + 2*exp(-bkd) + 2*exp(-bkD) ...
             - exp(-bkD - bkd) - exp(-bkD + bkd));
    end
    C = C * D0 * g2_perp * td^3;
end

function pars = fit_smt(b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data)
    % Fit SMT powder-average model using lsqnonlin
    S_data = double(S_data(:));
    cost = @(x) smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data);
    opts = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',500, ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12);
    x = lsqnonlin(cost, [2 0.7 0.05 0.8], [0 0.01 0 0.01], [5 1 1 DeL], opts);
    pars = struct('r',x(1),'f',x(2),'fcsf',x(3),'DeR',x(4));
end

function res = smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data)
    % SMT powder-averaged forward model residual
    r = x(1); f = x(2); fcsf = x(3); DeR = x(4);
    N = numel(b);
    S = zeros(N,1);
    for is = 1:N
        C   = vg_atten_fit(r, g2(is), del(is), Del(is), D0, bm2);
        arg = max(b(is)*Da - C, 1e-10);
        Sa  = sqrt(pi/(4*arg)) * exp(-C) * erf(sqrt(arg));
        dDe = max(DeL - DeR, 1e-10);
        Se  = sqrt(pi/(4*dDe*b(is))) * exp(-b(is)*DeR) * erf(sqrt(b(is)*dDe));
        S(is) = (1-fcsf)*(f*Sa + (1-f)*Se) + fcsf*exp(-b(is)*Dcsf);
    end
    res = min(S,1) - S_data;
end

function C = vg_atten_fit(r, g2p, del, Del, D0, bm2)
    td = r^2/D0; C = 0;
    for k = 1:numel(bm2)
        bk = bm2(k); bkd = bk*del/td; bkD = bk*Del/td;
        C = C + (2/(bk^3*(bk-1)))*(-2+2*bkd+2*exp(-bkd)+2*exp(-bkD)-exp(-bkD-bkd)-exp(-bkD+bkd));
    end
    C = C*D0*g2p*td^3;
end
