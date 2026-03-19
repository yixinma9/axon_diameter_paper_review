function demo_bvec_noise_propagation()
clear; close all;

%% b-vector GNL noise propagation (full L-matrix: scaling + rotation)
%
% For each protocol (C2, C1):
%   - Load grad_dev, extract full 3x3 L matrix per WM voxel
%   - Sample Nsample voxels, generate signal using directional model
%     (per-direction VG + zeppelin, averaged per shell)
%   - Add Rician noise
%   - Fit with nominal b (uncorrected) and b*b_scale (corrected)
%   - Produce scatter plots, bar plots, error distributions

%% 1. Fixed model parameters
D0 = 2; Da = 1.7; DeL = 1.7; Dcsf = 3;

%% 2. Protocol definitions
b_C2     = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
            0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta_C2 = [8*ones(8,1); 8*ones(8,1)];
Delta_C2 = [19*ones(8,1); 49*ones(8,1)];

b_C1     = b_C2 * (290/495)^2;
delta_C1 = delta_C2;
Delta_C1 = Delta_C2;

protocols = struct( ...
    'name',     {'C2', 'C1'}, ...
    'gnc_dir',  {'/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc', ...
                 '/autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016/gnc'}, ...
    'mask_dir', {'/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001', ...
                 '/autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016'}, ...
    'mask_file',{'sub-001_mask_nocsfgm.nii.gz', 'sub-016_mask_nocsfgm.nii.gz'}, ...
    'b',        {b_C2, b_C1}, ...
    'delta',    {delta_C2, delta_C1}, ...
    'Delta',    {Delta_C2, Delta_C1}, ...
    'snr',      {38, 18});

%% 3. Gradient directions from sub-001.bvec (16 representative directions)
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
fiber_dir = [0; 0; 1];

%% 4. Van Gelderen Bessel zeros
bm  = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

%% 5. Random ground truth (shared across protocols)
Nsample = 2000;
rng(0);
Xrange = [0.1 5; 0.5 1; 0 0.2; 0.5 1.5];  % [r, f, fcsf, DeR]
Xinit = Xrange(:,1) + diff(Xrange,1,2) .* rand(4, Nsample);
r_true    = Xinit(1,:)';
f_true    = Xinit(2,:)';
fcsf_true = Xinit(3,:)';
DeR_true  = Xinit(4,:)';

%% 6. Main loop over protocols
all_results = struct();

for ip = 1:numel(protocols)
    pname = protocols(ip).name;
    fprintf('\n====== Protocol: %s ======\n', pname);

    % --- Load grad_dev and WM mask ---
    gnc_dir  = protocols(ip).gnc_dir;
    mask_dir = protocols(ip).mask_dir;

    gx = single(niftiread(fullfile(gnc_dir, 'grad_dev_x.nii.gz')));
    gy = single(niftiread(fullfile(gnc_dir, 'grad_dev_y.nii.gz')));
    gz = single(niftiread(fullfile(gnc_dir, 'grad_dev_z.nii.gz')));
    g  = cat(4, gx, gy, gz);  % [x,y,z,9]
    clear gx gy gz;

    dims = size(g, 1:3);
    b_scale_map = zeros(dims, 'single');
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                Lv = reshape(squeeze(g(i,j,k,:)), 3, 3) / 100;
                Mv = eye(3) + Lv;
                b_scale_map(i,j,k) = trace(Mv'*Mv) / 3;
            end
        end
    end

    wm_mask = niftiread(fullfile(mask_dir, protocols(ip).mask_file)) > 0;
    wm_idx  = find(wm_mask(:));
    N_wm    = numel(wm_idx);

    g_flat = reshape(g, [], 9);
    L_wm   = double(g_flat(wm_idx, :)) / 100;   % [N_wm x 9]
    bscale_wm = double(b_scale_map(wm_idx));
    clear g g_flat b_scale_map;

    fprintf('%s WM voxels: %d\n', pname, N_wm);
    fprintf('%s b_scale range: [%.4f, %.4f]\n', pname, min(bscale_wm), max(bscale_wm));

    % --- Sample Nsample voxels from all WM ---
    rng(42);
    sample_idx = randperm(N_wm, min(Nsample, N_wm));
    L_samp     = L_wm(sample_idx, :);        % [Nsample x 9]
    bscale_samp = bscale_wm(sample_idx);      % [Nsample x 1]
    N = numel(sample_idx);

    % --- Protocol quantities ---
    b_sh   = protocols(ip).b;
    del_sh = protocols(ip).delta;
    Del_sh = protocols(ip).Delta;
    snr    = protocols(ip).snr;
    Nsh    = numel(b_sh);
    g_mag  = sqrt(b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3)));
    g2_sh  = b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3));

    % --- Noise propagation ---
    r_fit_uncorr    = zeros(N,1); f_fit_uncorr    = zeros(N,1);
    fcsf_fit_uncorr = zeros(N,1); DeR_fit_uncorr  = zeros(N,1);
    r_fit_corr      = zeros(N,1); f_fit_corr      = zeros(N,1);
    fcsf_fit_corr   = zeros(N,1); DeR_fit_corr    = zeros(N,1);

    tic;
    for i = 1:N
        if mod(i, 100) == 0; fprintf('  %s voxel %d/%d (%.0fs)\n', pname, i, N, toc); end

        % This voxel's full L matrix
        L = reshape(L_samp(i,:), 3, 3);
        M = eye(3) + L;
        b_scale_scalar = trace(M'*M) / 3;

        % Per-direction effective quantities
        g_eff_all       = M * bvecs;                       % [3 x Ndirs]
        b_scale_per_dir = sum(g_eff_all.^2, 1);            % [1 x Ndirs]
        g_eff_hat       = g_eff_all ./ vecnorm(g_eff_all);
        cos_theta       = abs(g_eff_hat' * fiber_dir);     % [Ndirs x 1]
        sin2_theta      = 1 - cos_theta.^2;

        % Generate signal: directional model, averaged per shell
        S_clean = zeros(Nsh, 1);
        for is = 1:Nsh
            Sd = zeros(Ndirs, 1);
            for id = 1:Ndirs
                b_eff = b_sh(is) * b_scale_per_dir(id);
                g_eff = g_mag(is) * sqrt(b_scale_per_dir(id));
                g2p   = g_eff^2 * sin2_theta(id);

                C_r = vg_atten(r_true(i), g2p, del_sh(is), Del_sh(is), D0, bm2);
                Sa  = exp(-b_eff * Da * cos_theta(id)^2 - C_r);
                Se  = exp(-b_eff * (DeR_true(i) + (DeL - DeR_true(i)) * cos_theta(id)^2));
                Sc  = exp(-b_eff * Dcsf);
                Sd(id) = (1-fcsf_true(i)) * (f_true(i)*Sa + (1-f_true(i))*Se) + fcsf_true(i)*Sc;
            end
            S_clean(is) = mean(Sd);
        end

        % Add Rician noise
        S_noisy = double(abs(S_clean + 1/snr*randn(Nsh,1) + 1j/snr*randn(Nsh,1)));

        % Fit uncorrected (nominal b)
        p = fit_smt(b_sh, g2_sh, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_noisy);
        r_fit_uncorr(i) = p.r; f_fit_uncorr(i) = p.f;
        fcsf_fit_uncorr(i) = p.fcsf; DeR_fit_uncorr(i) = p.DeR;

        % Fit corrected (b * b_scale_scalar)
        b_corr  = b_sh * b_scale_scalar;
        g2_corr = b_corr ./ (del_sh.^2 .* (Del_sh - del_sh/3));
        p = fit_smt(b_corr, g2_corr, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_noisy);
        r_fit_corr(i) = p.r; f_fit_corr(i) = p.f;
        fcsf_fit_corr(i) = p.fcsf; DeR_fit_corr(i) = p.DeR;
    end
    fprintf('  %s done: %.0f s\n', pname, toc);

    % Store results
    res.name = pname;
    res.bscales = bscale_samp(1:N);
    res.r_true = r_true(1:N);
    res.r_fit_uncorr = r_fit_uncorr;
    res.r_fit_corr   = r_fit_corr;
    res.f_fit_uncorr = f_fit_uncorr;
    res.f_fit_corr   = f_fit_corr;
    res.fcsf_fit_uncorr = fcsf_fit_uncorr;
    res.fcsf_fit_corr   = fcsf_fit_corr;
    res.DeR_fit_uncorr = DeR_fit_uncorr;
    res.DeR_fit_corr   = DeR_fit_corr;
    all_results.(pname) = res;
end

%% 7. Filter boundary hits
r_lb = 0.05; r_ub = 4.95;
pnames = fieldnames(all_results);
for ip = 1:numel(pnames)
    R = all_results.(pnames{ip});
    valid = R.r_fit_uncorr > r_lb & R.r_fit_uncorr < r_ub & ...
            R.r_fit_corr   > r_lb & R.r_fit_corr   < r_ub;
    all_results.(pnames{ip}).valid = valid;
    all_results.(pnames{ip}).r_true_filt       = R.r_true(valid);
    all_results.(pnames{ip}).r_fit_uncorr_filt = R.r_fit_uncorr(valid);
    all_results.(pnames{ip}).r_fit_corr_filt   = R.r_fit_corr(valid);
    all_results.(pnames{ip}).bscales_filt      = R.bscales(valid);
    fprintf('%s: %d/%d valid\n', pnames{ip}, sum(valid), numel(valid));
end

%% 8. Quantitative evaluation — stratified by r
r_bins = [0 0.5; 0.5 1; 1 1.5; 1.5 2; 2 3; 3 4; 4 5];
Nbins = size(r_bins, 1);

for ip = 1:numel(pnames)
    R  = all_results.(pnames{ip});
    rt = R.r_true_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;

    fprintf('\n========== %s (filtered) ==========\n', pnames{ip});
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

    err_unc = ru - rt;
    err_cor = rc - rt;
    fprintf('%-15s %5d %+10.4f %+10.4f %+10.4f %+10.4f %10.4f %10.4f\n', 'Overall', ...
            numel(err_unc), mean(err_unc), mean(err_cor), ...
            median(err_unc), median(err_cor), ...
            sqrt(mean(err_unc.^2)), sqrt(mean(err_cor.^2)));

    all_results.(pnames{ip}).bias_uncorr = mean(err_unc);
    all_results.(pnames{ip}).bias_corr   = mean(err_cor);
    all_results.(pnames{ip}).rmse_uncorr = sqrt(mean(err_unc.^2));
    all_results.(pnames{ip}).rmse_corr   = sqrt(mean(err_cor.^2));
end

%% 9. Fig 2 — Error vs r_true
Nprot = numel(pnames);
figure('unit','inch','position',[0 0 6*Nprot 5]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    subplot(1, Nprot, ip);
    scatter(rt, R.r_fit_uncorr_filt - rt, 8, [0.8 0.3 0.3], '.', 'MarkerFaceAlpha', 0.3); hold on;
    scatter(rt, R.r_fit_corr_filt   - rt, 8, [0.3 0.3 0.8], '.', 'MarkerFaceAlpha', 0.3);
    yline(0, 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} - r_{true} (\mum)');
    title(sprintf('%s (N=%d, mean b_{sc}=%.3f)', pnames{ip}, sum(R.valid), mean(R.bscales_filt)));
    legend({'Uncorrected', 'Corrected'}, 'Location', 'best');
    grid on; box on;
end
exportgraphics(gcf, 'fig2_bvec_error_vs_rtrue.pdf', 'ContentType', 'vector');

%% 10. Fig 3 — Scatter: r_true vs r_fitted (2 rows x Nprot cols)
figure('unit','inch','position',[0 0 5*Nprot 8]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;
    err_unc = ru - rt; err_cor = rc - rt;

    % Uncorrected
    subplot(2, Nprot, ip);
    scatter(rt, ru, 10, '.'); hold on;
    plot(Xrange(1,:), Xrange(1,:), 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s — Uncorrected\nbias=%.3f, RMSE=%.3f', pnames{ip}, mean(err_unc), sqrt(mean(err_unc.^2))));

    % Corrected
    subplot(2, Nprot, ip + Nprot);
    scatter(rt, rc, 10, '.'); hold on;
    plot(Xrange(1,:), Xrange(1,:), 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s — Corrected\nbias=%.3f, RMSE=%.3f', pnames{ip}, mean(err_cor), sqrt(mean(err_cor.^2))));
end
exportgraphics(gcf, 'fig3_bvec_scatter_gt_vs_fit.pdf', 'ContentType', 'vector');

%% 11. Fig 4 — Bar: bias & RMSE by r-bin
r_bins_bar = [0 1; 1 2; 2 3; 3 4; 4 5];
Nbins_bar  = size(r_bins_bar, 1);
bin_labels = arrayfun(@(i) sprintf('%g-%g', r_bins_bar(i,1), r_bins_bar(i,2)), 1:Nbins_bar, 'uni', 0);

figure('unit','inch','position',[0 0 7*Nprot 8]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;

    bias_unc = zeros(Nbins_bar,1); bias_cor = zeros(Nbins_bar,1);
    rmse_unc = zeros(Nbins_bar,1); rmse_cor = zeros(Nbins_bar,1);
    for ib = 1:Nbins_bar
        idx = rt >= r_bins_bar(ib,1) & rt < r_bins_bar(ib,2);
        if sum(idx) == 0; continue; end
        eu = ru(idx)-rt(idx); ec = rc(idx)-rt(idx);
        bias_unc(ib) = mean(eu); bias_cor(ib) = mean(ec);
        rmse_unc(ib) = sqrt(mean(eu.^2)); rmse_cor(ib) = sqrt(mean(ec.^2));
    end

    subplot(2, Nprot, ip);
    bar(categorical(bin_labels, bin_labels), [bias_unc, bias_cor]);
    ylabel('Bias (\mum)'); xlabel('r_{true} (\mum)');
    legend({'Uncorrected','Corrected'}, 'Location','best');
    title(sprintf('%s — Bias', pnames{ip}));
    yline(0, 'k--', 'HandleVisibility','off'); grid on;

    subplot(2, Nprot, ip + Nprot);
    bar(categorical(bin_labels, bin_labels), [rmse_unc, rmse_cor]);
    ylabel('RMSE (\mum)'); xlabel('r_{true} (\mum)');
    legend({'Uncorrected','Corrected'}, 'Location','best');
    title(sprintf('%s — RMSE', pnames{ip}));
    grid on;
end
exportgraphics(gcf, 'fig4_bvec_bar_bias_rmse.pdf', 'ContentType', 'vector');

%% 12. Fig 5 — Error distribution
figure('unit','inch','position',[0 0 6*Nprot 4]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    err_unc = R.r_fit_uncorr_filt - rt;
    err_cor = R.r_fit_corr_filt   - rt;

    subplot(1, Nprot, ip);
    histogram(err_unc, 50, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.5); hold on;
    histogram(err_cor, 50, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.5);
    xline(mean(err_unc), 'r--', 'LineWidth', 2);
    xline(mean(err_cor), 'b--', 'LineWidth', 2);
    xline(0, 'k-', 'LineWidth', 1);
    xlabel('r_{fit} - r_{true} (\mum)'); ylabel('Count');
    title(sprintf('%s (mean b_{sc}=%.3f)', pnames{ip}, mean(R.bscales_filt)));
    legend({'Uncorrected','Corrected', ...
            sprintf('Bias=%.3f', mean(err_unc)), sprintf('Bias=%.3f', mean(err_cor))});
end
exportgraphics(gcf, 'fig5_bvec_error_distribution.pdf', 'ContentType', 'vector');

%% 13. Save
save('bvec_noise_propagation_results.mat', 'all_results', 'Xrange', ...
     'r_true', 'f_true', 'fcsf_true', 'DeR_true', 'Nsample');
fprintf('\nDone. Saved bvec_noise_propagation_results.mat\n');

end  % main function

%% ======================================================================
%  Helper functions
%  ======================================================================

function C = vg_atten(r, g2_perp, delta, Delta, D0, ~)
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
    S_data = double(S_data(:));
    cost = @(x) smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data);
    opts = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',500, ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12);
    x = lsqnonlin(cost, [2 0.7 0.05 0.8], [0 0.01 0 0.01], [5 1 1 DeL], opts);
    pars = struct('r',x(1),'f',x(2),'fcsf',x(3),'DeR',x(4));
end

function res = smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, ~, S_data)
    bm2 = [3.3899 28.4238 72.9367 137.0305 221.0266 324.5583 447.9337 591.1395 753.4878 935.6762];
    r = x(1); f = x(2); fcsf = x(3); DeR = x(4);
    N = numel(b);
    S = zeros(N,1);
    for is = 1:N
        td = r^2/D0; C = 0;
        for k = 1:numel(bm2)
            bk = bm2(k); bkd = bk*del(is)/td; bkD = bk*Del(is)/td;
            C = C + (2/(bk^3*(bk-1)))*(-2+2*bkd+2*exp(-bkd)+2*exp(-bkD)-exp(-bkD-bkd)-exp(-bkD+bkd));
        end
        C = C*D0*g2(is)*td^3;
        arg = max(b(is)*Da - C, 1e-10);
        Sa  = sqrt(pi/(4*arg)) * exp(-C) * erf(sqrt(arg));
        dDe = max(DeL - DeR, 1e-10);
        Se  = sqrt(pi/(4*dDe*b(is))) * exp(-b(is)*DeR) * erf(sqrt(b(is)*dDe));
        S(is) = (1-fcsf)*(f*Sa + (1-f)*Se) + fcsf*exp(-b(is)*Dcsf);
    end
    res = min(S,1) - S_data;
end
