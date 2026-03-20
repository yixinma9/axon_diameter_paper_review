function demo_bvec_noise_propagation()
clear; close all;

%% b-vector GNL noise propagation (full L-matrix: scaling + rotation)
%
% For each protocol (C2, C1) with protocol-specific bvec/bval:
%   - Load bvec/bval, parse per-shell direction tables (varying Ndirs)
%   - Load grad_dev, extract full 3x3 L matrix per WM voxel
%   - Sample Nsample voxels, generate signal using directional model
%     (per-direction VG + zeppelin, averaged per shell)
%   - Add Rician noise (shared noise realization across scenarios)
%   - Three scenarios:
%       1) Baseline: no GNL signal, fit with nominal b
%       2) Uncorrected: GNL signal, fit with nominal b
%       3) Corrected: GNL signal, fit with b*b_scale
%   - Produce scatter plots, bar plots, error distributions

%% 1. Fixed model parameters
D0 = 2; Da = 1.7; DeL = 1.7; Dcsf = 3;

%% 2. Protocol definitions (b-values in ms/um^2)
b_shells = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
            0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta_sh = [8*ones(8,1); 8*ones(8,1)];
Delta_sh = [19*ones(8,1); 49*ones(8,1)];
Nsh = numel(b_shells);
bval_targets = round(b_shells * 1000);  % convert to s/mm^2 for matching

g_mag = sqrt(b_shells ./ (delta_sh.^2 .* (Delta_sh - delta_sh/3)));
g2_sh = b_shells ./ (delta_sh.^2 .* (Delta_sh - delta_sh/3));

%% 3. Protocol paths
protocols = struct( ...
    'name',     {'C2', 'C1'}, ...
    'gnc_dir',  {'/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc', ...
                 '/autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016/gnc'}, ...
    'mask_dir', {'/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001', ...
                 '/autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016'}, ...
    'mask_file',{'sub-001_mask_nocsfgm.nii.gz', 'sub-016_mask_nocsfgm.nii.gz'}, ...
    'bvec_file',{'/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/sub-001.bvec', ...
                 '/autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016/sub-016.bvec'}, ...
    'bval_file',{'/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/sub-001.bval', ...
                 '/autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016/sub-016.bval'}, ...
    'snr',      {38, 18});

%% 4. Parse bvec/bval for each protocol → per-shell direction tables
fiber_dir = [0; 0; 1];

% Pre-compute per-shell direction data for each protocol
proto_dirs = cell(numel(protocols), 1);  % each entry is a struct with per-shell info

for ip = 1:numel(protocols)
    fprintf('Loading bvec/bval for %s...\n', protocols(ip).name);

    [bv3, bvals] = load_bvec_bval(protocols(ip).bvec_file, protocols(ip).bval_file);

    % Build per-shell direction tables (tolerance-based bval matching)
    dirs_per_shell = cell(Nsh, 1);
    ndirs = zeros(Nsh, 1);
    for is = 1:Nsh
        % Match within 2% or 25 s/mm^2, whichever is larger
        tol = max(25, 0.02 * bval_targets(is));
        idx = find(abs(bvals - bval_targets(is)) < tol);
        dirs = bv3(:, idx);
        % Remove zero-norm directions (b=0 volumes)
        norms = vecnorm(dirs);
        dirs = dirs(:, norms > 0.5);
        dirs_per_shell{is} = dirs;
        ndirs(is) = size(dirs, 2);
    end

    fprintf('  %s: shells with directions: %s\n', protocols(ip).name, ...
        mat2str(ndirs'));

    % Concatenate all directions and record shell assignment
    all_dirs = cell2mat(dirs_per_shell');     % [3 x total_dirs]
    total_dirs = size(all_dirs, 2);
    shell_assign = zeros(total_dirs, 1);     % which shell each dir belongs to
    offset = 0;
    for is = 1:Nsh
        shell_assign(offset+1 : offset+ndirs(is)) = is;
        offset = offset + ndirs(is);
    end

    % Pre-compute nominal cos_theta per direction
    cos_theta_nom = abs(all_dirs' * fiber_dir);   % [total_dirs x 1]
    sin2_theta_nom = 1 - cos_theta_nom.^2;

    pd.dirs_per_shell = dirs_per_shell;
    pd.ndirs = ndirs;
    pd.all_dirs = all_dirs;
    pd.total_dirs = total_dirs;
    pd.shell_assign = shell_assign;
    pd.cos_theta_nom = cos_theta_nom;
    pd.sin2_theta_nom = sin2_theta_nom;
    proto_dirs{ip} = pd;
end

%% 5. Van Gelderen Bessel zeros
bm  = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

%% 6. Random ground truth (shared across protocols)
Nsample = 2000;
rng(0);
Xrange = [0.1 5; 0.5 1; 0 0.2; 0.5 1.5];  % [r, f, fcsf, DeR]
Xinit = Xrange(:,1) + diff(Xrange,1,2) .* rand(4, Nsample);
r_true    = Xinit(1,:)';
f_true    = Xinit(2,:)';
fcsf_true = Xinit(3,:)';
DeR_true  = Xinit(4,:)';

%% 7. Main loop over protocols
all_results = struct();

for ip = 1:numel(protocols)
    pname = protocols(ip).name;
    fprintf('\n====== Protocol: %s ======\n', pname);

    % --- Check grad_dev existence ---
    gnc_dir  = protocols(ip).gnc_dir;
    mask_dir = protocols(ip).mask_dir;
    grad_dev_file = fullfile(gnc_dir, 'grad_dev_x.nii.gz');
    if ~isfile(grad_dev_file)
        fprintf('  SKIPPING %s: grad_dev not found at %s\n', pname, gnc_dir);
        continue;
    end

    % --- Per-shell direction data ---
    pd = proto_dirs{ip};

    % --- Load grad_dev and WM mask ---
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

    snr = protocols(ip).snr;

    % --- Allocate results for 3 scenarios ---
    r_fit_base   = zeros(N,1); f_fit_base   = zeros(N,1);
    fcsf_fit_base = zeros(N,1); DeR_fit_base = zeros(N,1);
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

        % === Apply GNL to ALL directions at once ===
        g_eff_all = M * pd.all_dirs;                        % [3 x total_dirs]
        bscale_dir = sum(g_eff_all.^2, 1);                  % [1 x total_dirs]
        g_eff_hat = g_eff_all ./ vecnorm(g_eff_all);
        cos_theta_gnl = abs(g_eff_hat' * fiber_dir);        % [total_dirs x 1]
        sin2_theta_gnl = 1 - cos_theta_gnl.^2;

        % === Generate shared noise ===
        n1 = randn(Nsh, 1);
        n2 = randn(Nsh, 1);

        % === Scenario 1: Baseline (no GNL signal, nominal fit) ===
        S_base = zeros(Nsh, 1);
        for is = 1:Nsh
            idx = (pd.shell_assign == is);
            ct = pd.cos_theta_nom(idx);
            st2 = pd.sin2_theta_nom(idx);
            nd = pd.ndirs(is);
            Sd = zeros(nd, 1);
            for id = 1:nd
                g2p = g2_sh(is) * st2(id);
                C_r = vg_atten(r_true(i), g2p, delta_sh(is), Delta_sh(is), D0, bm2);
                Sa  = exp(-b_shells(is) * Da * ct(id)^2 - C_r);
                Se  = exp(-b_shells(is) * (DeR_true(i) + (DeL - DeR_true(i)) * ct(id)^2));
                Sc  = exp(-b_shells(is) * Dcsf);
                Sd(id) = (1-fcsf_true(i)) * (f_true(i)*Sa + (1-f_true(i))*Se) + fcsf_true(i)*Sc;
            end
            S_base(is) = mean(Sd);
        end
        S_noisy_base = double(abs(S_base + 1/snr*n1 + 1j/snr*n2));

        p = fit_smt(b_shells, g2_sh, delta_sh, Delta_sh, Da, DeL, D0, Dcsf, bm2, S_noisy_base);
        r_fit_base(i) = p.r; f_fit_base(i) = p.f;
        fcsf_fit_base(i) = p.fcsf; DeR_fit_base(i) = p.DeR;

        % === Scenario 2 & 3: GNL signal ===
        S_gnl = zeros(Nsh, 1);
        for is = 1:Nsh
            idx = (pd.shell_assign == is);
            ct = cos_theta_gnl(idx);
            st2 = sin2_theta_gnl(idx);
            bs = bscale_dir(idx);
            nd = pd.ndirs(is);
            Sd = zeros(nd, 1);
            for id = 1:nd
                b_eff = b_shells(is) * bs(id);
                g2p   = g2_sh(is) * bs(id) * st2(id);

                C_r = vg_atten(r_true(i), g2p, delta_sh(is), Delta_sh(is), D0, bm2);
                Sa  = exp(-b_eff * Da * ct(id)^2 - C_r);
                Se  = exp(-b_eff * (DeR_true(i) + (DeL - DeR_true(i)) * ct(id)^2));
                Sc  = exp(-b_eff * Dcsf);
                Sd(id) = (1-fcsf_true(i)) * (f_true(i)*Sa + (1-f_true(i))*Se) + fcsf_true(i)*Sc;
            end
            S_gnl(is) = mean(Sd);
        end
        S_noisy_gnl = double(abs(S_gnl + 1/snr*n1 + 1j/snr*n2));

        % Fit uncorrected (nominal b)
        p = fit_smt(b_shells, g2_sh, delta_sh, Delta_sh, Da, DeL, D0, Dcsf, bm2, S_noisy_gnl);
        r_fit_uncorr(i) = p.r; f_fit_uncorr(i) = p.f;
        fcsf_fit_uncorr(i) = p.fcsf; DeR_fit_uncorr(i) = p.DeR;

        % Fit corrected (b * b_scale_scalar)
        b_corr  = b_shells * b_scale_scalar;
        g2_corr = b_corr ./ (delta_sh.^2 .* (Delta_sh - delta_sh/3));
        p = fit_smt(b_corr, g2_corr, delta_sh, Delta_sh, Da, DeL, D0, Dcsf, bm2, S_noisy_gnl);
        r_fit_corr(i) = p.r; f_fit_corr(i) = p.f;
        fcsf_fit_corr(i) = p.fcsf; DeR_fit_corr(i) = p.DeR;
    end
    fprintf('  %s done: %.0f s\n', pname, toc);

    % Store results
    res.name = pname;
    res.bscales = bscale_samp(1:N);
    res.r_true = r_true(1:N);
    res.r_fit_base   = r_fit_base;
    res.r_fit_uncorr = r_fit_uncorr;
    res.r_fit_corr   = r_fit_corr;
    res.f_fit_base   = f_fit_base;
    res.f_fit_uncorr = f_fit_uncorr;
    res.f_fit_corr   = f_fit_corr;
    res.fcsf_fit_base   = fcsf_fit_base;
    res.fcsf_fit_uncorr = fcsf_fit_uncorr;
    res.fcsf_fit_corr   = fcsf_fit_corr;
    res.DeR_fit_base   = DeR_fit_base;
    res.DeR_fit_uncorr = DeR_fit_uncorr;
    res.DeR_fit_corr   = DeR_fit_corr;
    all_results.(pname) = res;
end

%% 8. Filter boundary hits (require all 3 scenarios within bounds)
r_lb = 0.05; r_ub = 4.95;
pnames = fieldnames(all_results);
if isempty(pnames); fprintf('No protocols completed. Exiting.\n'); return; end

for ip = 1:numel(pnames)
    R = all_results.(pnames{ip});
    valid = R.r_fit_base   > r_lb & R.r_fit_base   < r_ub & ...
            R.r_fit_uncorr > r_lb & R.r_fit_uncorr < r_ub & ...
            R.r_fit_corr   > r_lb & R.r_fit_corr   < r_ub;
    all_results.(pnames{ip}).valid = valid;
    all_results.(pnames{ip}).r_true_filt       = R.r_true(valid);
    all_results.(pnames{ip}).r_fit_base_filt   = R.r_fit_base(valid);
    all_results.(pnames{ip}).r_fit_uncorr_filt = R.r_fit_uncorr(valid);
    all_results.(pnames{ip}).r_fit_corr_filt   = R.r_fit_corr(valid);
    all_results.(pnames{ip}).bscales_filt      = R.bscales(valid);
    fprintf('%s: %d/%d valid\n', pnames{ip}, sum(valid), numel(valid));
end

%% 9. Quantitative evaluation — stratified by r
r_bins = [0 0.5; 0.5 1; 1 1.5; 1.5 2; 2 3; 3 4; 4 5];
Nbins = size(r_bins, 1);

for ip = 1:numel(pnames)
    R  = all_results.(pnames{ip});
    rt = R.r_true_filt;
    rb = R.r_fit_base_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;

    fprintf('\n========== %s (filtered) ==========\n', pnames{ip});
    fprintf('%-15s %6s %10s %10s %10s\n', ...
            'r range', 'N', 'Bias(base)', 'Bias(unc)', 'Bias(cor)');

    for ib = 1:Nbins
        idx = rt >= r_bins(ib,1) & rt < r_bins(ib,2);
        if sum(idx) == 0; continue; end
        err_b = rb(idx) - rt(idx);
        err_u = ru(idx) - rt(idx);
        err_c = rc(idx) - rt(idx);
        fprintf('[%.1f-%.1f] um %5d %+10.4f %+10.4f %+10.4f\n', ...
                r_bins(ib,1), r_bins(ib,2), sum(idx), ...
                mean(err_b), mean(err_u), mean(err_c));
    end

    err_b = rb - rt; err_u = ru - rt; err_c = rc - rt;
    fprintf('%-15s %5d %+10.4f %+10.4f %+10.4f\n', 'Overall', ...
            numel(err_b), mean(err_b), mean(err_u), mean(err_c));
    fprintf('%-15s %5s %10.4f %10.4f %10.4f\n', 'RMSE', '', ...
            sqrt(mean(err_b.^2)), sqrt(mean(err_u.^2)), sqrt(mean(err_c.^2)));

    all_results.(pnames{ip}).bias_base   = mean(err_b);
    all_results.(pnames{ip}).bias_uncorr = mean(err_u);
    all_results.(pnames{ip}).bias_corr   = mean(err_c);
    all_results.(pnames{ip}).rmse_base   = sqrt(mean(err_b.^2));
    all_results.(pnames{ip}).rmse_uncorr = sqrt(mean(err_u.^2));
    all_results.(pnames{ip}).rmse_corr   = sqrt(mean(err_c.^2));
end

%% 10. Fig 2 — Error vs r_true
Nprot = numel(pnames);
colors = {[0.4 0.7 0.4], [0.8 0.3 0.3], [0.3 0.3 0.8]};  % base, uncorr, corr

figure('unit','inch','position',[0 0 6*Nprot 5]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    subplot(1, Nprot, ip);
    scatter(rt, R.r_fit_base_filt   - rt, 8, colors{1}, '.', 'MarkerFaceAlpha', 0.3); hold on;
    scatter(rt, R.r_fit_uncorr_filt - rt, 8, colors{2}, '.', 'MarkerFaceAlpha', 0.3);
    scatter(rt, R.r_fit_corr_filt   - rt, 8, colors{3}, '.', 'MarkerFaceAlpha', 0.3);
    yline(0, 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} - r_{true} (\mum)');
    title(sprintf('%s (N=%d, mean b_{sc}=%.3f)', pnames{ip}, sum(R.valid), mean(R.bscales_filt)));
    legend({'Baseline','Uncorrected','Corrected'}, 'Location', 'best');
    grid on; box on;
end
exportgraphics(gcf, 'fig2_bvec_error_vs_rtrue.pdf', 'ContentType', 'vector');

%% 11. Fig 3 — Scatter: r_true vs r_fitted (3 rows x Nprot cols)
figure('unit','inch','position',[0 0 5*Nprot 12]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    rb = R.r_fit_base_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;
    err_b = rb - rt; err_u = ru - rt; err_c = rc - rt;

    % Baseline
    subplot(3, Nprot, ip);
    scatter(rt, rb, 10, '.'); hold on;
    plot(Xrange(1,:), Xrange(1,:), 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s — Baseline\nbias=%.3f, RMSE=%.3f', pnames{ip}, mean(err_b), sqrt(mean(err_b.^2))));

    % Uncorrected
    subplot(3, Nprot, ip + Nprot);
    scatter(rt, ru, 10, '.'); hold on;
    plot(Xrange(1,:), Xrange(1,:), 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s — Uncorrected\nbias=%.3f, RMSE=%.3f', pnames{ip}, mean(err_u), sqrt(mean(err_u.^2))));

    % Corrected
    subplot(3, Nprot, ip + 2*Nprot);
    scatter(rt, rc, 10, '.'); hold on;
    plot(Xrange(1,:), Xrange(1,:), 'k-', 'LineWidth', 1);
    xlabel('r_{true} (\mum)'); ylabel('r_{fit} (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s — Corrected\nbias=%.3f, RMSE=%.3f', pnames{ip}, mean(err_c), sqrt(mean(err_c.^2))));
end
exportgraphics(gcf, 'fig3_bvec_scatter_gt_vs_fit.pdf', 'ContentType', 'vector');

%% 12. Fig 4 — Bar: bias & RMSE by r-bin (3 scenarios)
r_bins_bar = [0 1; 1 2; 2 3; 3 4; 4 5];
Nbins_bar  = size(r_bins_bar, 1);
bin_labels = arrayfun(@(i) sprintf('%g-%g', r_bins_bar(i,1), r_bins_bar(i,2)), 1:Nbins_bar, 'uni', 0);

figure('unit','inch','position',[0 0 7*Nprot 8]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    rb = R.r_fit_base_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;

    bias_b = zeros(Nbins_bar,1); bias_u = zeros(Nbins_bar,1); bias_c = zeros(Nbins_bar,1);
    rmse_b = zeros(Nbins_bar,1); rmse_u = zeros(Nbins_bar,1); rmse_c = zeros(Nbins_bar,1);
    for ib = 1:Nbins_bar
        idx = rt >= r_bins_bar(ib,1) & rt < r_bins_bar(ib,2);
        if sum(idx) == 0; continue; end
        eb = rb(idx)-rt(idx); eu = ru(idx)-rt(idx); ec = rc(idx)-rt(idx);
        bias_b(ib) = mean(eb); bias_u(ib) = mean(eu); bias_c(ib) = mean(ec);
        rmse_b(ib) = sqrt(mean(eb.^2)); rmse_u(ib) = sqrt(mean(eu.^2)); rmse_c(ib) = sqrt(mean(ec.^2));
    end

    subplot(2, Nprot, ip);
    bar(categorical(bin_labels, bin_labels), [bias_b, bias_u, bias_c]);
    ylabel('Bias (\mum)'); xlabel('r_{true} (\mum)');
    legend({'Baseline','Uncorrected','Corrected'}, 'Location','best');
    title(sprintf('%s — Bias', pnames{ip}));
    yline(0, 'k--', 'HandleVisibility','off'); grid on;

    subplot(2, Nprot, ip + Nprot);
    bar(categorical(bin_labels, bin_labels), [rmse_b, rmse_u, rmse_c]);
    ylabel('RMSE (\mum)'); xlabel('r_{true} (\mum)');
    legend({'Baseline','Uncorrected','Corrected'}, 'Location','best');
    title(sprintf('%s — RMSE', pnames{ip}));
    grid on;
end
exportgraphics(gcf, 'fig4_bvec_bar_bias_rmse.pdf', 'ContentType', 'vector');

%% 13. Fig 5 — Error distribution (3 scenarios)
figure('unit','inch','position',[0 0 6*Nprot 4]);
for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    err_b = R.r_fit_base_filt   - rt;
    err_u = R.r_fit_uncorr_filt - rt;
    err_c = R.r_fit_corr_filt   - rt;

    subplot(1, Nprot, ip);
    histogram(err_b, 50, 'FaceColor', colors{1}, 'FaceAlpha', 0.4); hold on;
    histogram(err_u, 50, 'FaceColor', colors{2}, 'FaceAlpha', 0.4);
    histogram(err_c, 50, 'FaceColor', colors{3}, 'FaceAlpha', 0.4);
    xline(mean(err_b), '-', 'Color', colors{1}, 'LineWidth', 2);
    xline(mean(err_u), '-', 'Color', colors{2}, 'LineWidth', 2);
    xline(mean(err_c), '-', 'Color', colors{3}, 'LineWidth', 2);
    xline(0, 'k-', 'LineWidth', 1);
    xlabel('r_{fit} - r_{true} (\mum)'); ylabel('Count');
    title(sprintf('%s (mean b_{sc}=%.3f)', pnames{ip}, mean(R.bscales_filt)));
    legend({sprintf('Base (%.3f)', mean(err_b)), ...
            sprintf('Uncorr (%.3f)', mean(err_u)), ...
            sprintf('Corr (%.3f)', mean(err_c))});
end
exportgraphics(gcf, 'fig5_bvec_error_distribution.pdf', 'ContentType', 'vector');

%% 14. Save
save('bvec_noise_propagation_results.mat', 'all_results', 'Xrange', ...
     'r_true', 'f_true', 'fcsf_true', 'DeR_true', 'Nsample');
fprintf('\nDone. Saved bvec_noise_propagation_results.mat\n');

end  % main function

%% ======================================================================
%  Helper functions
%  ======================================================================

function [bv3, bvals] = load_bvec_bval(bvec_file, bval_file)
    % Load bvec file (3 rows x Nvol columns)
    bv3 = dlmread(bvec_file);
    if size(bv3, 1) ~= 3
        bv3 = bv3';  % transpose if needed
    end
    % Load bval file (1 row x Nvol)
    bvals = dlmread(bval_file);
    bvals = bvals(:)';  % ensure row vector
    % Verify
    assert(size(bv3, 2) == numel(bvals), ...
        'bvec (%d) and bval (%d) volume count mismatch', size(bv3, 2), numel(bvals));
    fprintf('  Loaded %d volumes from bvec/bval\n', numel(bvals));
end

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
