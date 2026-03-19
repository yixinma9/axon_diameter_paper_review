clear; close all;

%% Simulate AxCaliber-SMT signal with full GNL (b-scaling + b-rotation)
%
% Three scenarios compared for BOTH C1 and C2:
%   Baseline  : No GNL            → fit with standard SMT
%   Scenario 1: GNL (scale+rot)   → fit with NO correction
%   Scenario 2: GNL (scale+rot)   → fit with scalar b_scale correction
%
% Signal generation uses the *directional* model (angle-dependent),
% then averaged per shell to obtain the spherical mean.
% Fitting uses the standard SMT (orientation-free) forward model.

addpath(genpath('/autofs/space/linen_001/users/Yixin/C2_protocoldesign-main/lib'));

%% 1. Ground truth tissue parameters
a_true    = 4;        % axon diameter [um]
f_true    = 0.5;      % restricted fraction
fcsf_true = 0.05;     % CSF fraction
DeR_true  = 0.5;      % extra-cellular radial diffusivity [um2/ms]

% Fixed model parameters
D0   = 2.0;   % intra-cellular intrinsic diffusivity [um2/ms]
Da   = 1.7;   % intra-cellular axial diffusivity [um2/ms]
DeL  = 1.7;   % extra-cellular axial diffusivity [um2/ms]
Dcsf = 3.0;   % CSF diffusivity [um2/ms]

%% 2. Protocol definitions for C1 and C2
% C2
b_C2     = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
            0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';   % [ms/um2]
delta_C2 = 8 * ones(16,1);                                % [ms]
Delta_C2 = [19*ones(8,1); 49*ones(8,1)];                  % [ms]

% C1 — same timing, b-values scaled by (Gmax_C1/Gmax_C2)^2
Gmax_C1 = 290; Gmax_C2 = 495;
b_C1     = b_C2 * (Gmax_C1/Gmax_C2)^2;
delta_C1 = delta_C2;
Delta_C1 = Delta_C2;

%% 3. Data paths (for loading actual b-tables and grad_dev)
% C2: /autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001
% C1: /autofs/cluster/connectome2/Bay8_C1/bids/derivatives/processed_dwi_old_way/sub-016
% C2 has gnc/ subfolder with grad_dev_{x,y,z}.nii.gz
% C1 does not have gnc yet
%
% To load actual b-table:
%   bvals = load('<path>/sub-001.bval');
%   bvecs = load('<path>/sub-001.bvec');  % [3 x Nvol]
%
% To load grad_dev L-matrix from gnc/:
%   Lx = niftiread('<path>/gnc/grad_dev_x.nii.gz');  % [x,y,z,3]
%   Ly = niftiread('<path>/gnc/grad_dev_y.nii.gz');
%   Lz = niftiread('<path>/gnc/grad_dev_z.nii.gz');
%   L = [Lx(vox,1) Lx(vox,2) Lx(vox,3);
%        Ly(vox,1) Ly(vox,2) Ly(vox,3);
%        Lz(vox,1) Lz(vox,2) Lz(vox,3)];

%% 4. Simulation parameters
% Fiber direction (arbitrary, for directional signal generation)
fiber_dir = [0; 0; 1];   % along z

% Gradient directions (uniform sampling; replace with actual bvecs if desired)
Ndirs = 16;
rng(42);
bvecs = randn(3, Ndirs);
bvecs = bvecs ./ vecnorm(bvecs);  % [3 x Ndirs] unit vectors

% GNL L-matrix: gradient deviation tensor, g_eff = (I + L) * g_nominal
% Base L from realistic gradient coil anisotropy:
L_base = [ 0.03  0.01  0.005;
           0.01 -0.02  0.008;
           0.005 0.008  0.04];   % ~3-4% diagonal, off-diagonal ~1%

gnl_scales = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0];  % 0 = no GNL
Ngnl = numel(gnl_scales);

% Van Gelderen Bessel function zeros
bm = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

%% 5. Define protocols to loop over
protocols = struct();
protocols(1).name  = 'C2';  protocols(1).b = b_C2;
protocols(1).delta = delta_C2; protocols(1).Delta = Delta_C2;
protocols(2).name  = 'C1';  protocols(2).b = b_C1;
protocols(2).delta = delta_C1; protocols(2).Delta = Delta_C1;

%% 6. Run simulation for each protocol and GNL level
all_results = struct();

for ip = 1:2
    pname = protocols(ip).name;
    b_shells     = protocols(ip).b;
    delta_shells = protocols(ip).delta;
    Delta_shells = protocols(ip).Delta;
    Nshells      = numel(b_shells);
    g_shells     = sqrt(b_shells ./ (delta_shells.^2 .* (Delta_shells - delta_shells/3)));

    res = struct();
    res.gnl_scales   = gnl_scales;
    res.b_scale_true = zeros(Ngnl, 1);
    res.a_baseline   = zeros(Ngnl, 1);  res.f_baseline   = zeros(Ngnl, 1);
    res.a_nocorr     = zeros(Ngnl, 1);  res.f_nocorr     = zeros(Ngnl, 1);
    res.a_bscale     = zeros(Ngnl, 1);  res.f_bscale     = zeros(Ngnl, 1);

    for ig = 1:Ngnl
        L = gnl_scales(ig) * L_base;
        M = eye(3) + L;

        % Per-direction effective quantities
        g_eff_all = M * bvecs;                       % [3 x Ndirs]
        b_scale_per_dir = sum(g_eff_all.^2, 1);      % |M*g_hat|^2
        g_eff_hat = g_eff_all ./ vecnorm(g_eff_all);

        % Scalar b_scale = trace(M'*M)/3
        b_scale_scalar = trace(M' * M) / 3;
        res.b_scale_true(ig) = b_scale_scalar;

        % Angles with fiber — GNL case
        cos_theta = abs(g_eff_hat' * fiber_dir);
        sin2_theta = 1 - cos_theta.^2;
        cos2_theta = cos_theta.^2;

        % Angles with fiber — nominal case
        cos_theta_nom = abs(bvecs' * fiber_dir);
        sin2_nom = 1 - cos_theta_nom.^2;
        cos2_nom = cos_theta_nom.^2;

        % Generate spherical mean signals
        S_gnl = zeros(Nshells, 1);
        S_baseline = zeros(Nshells, 1);

        for is = 1:Nshells
            b_nom = b_shells(is);
            g_nom = g_shells(is);
            del   = delta_shells(is);
            Del   = Delta_shells(is);

            % --- Baseline: no GNL ---
            S_dirs = zeros(Ndirs, 1);
            for id = 1:Ndirs
                g2p = g_nom^2 * sin2_nom(id);
                C_r = vangelderen_atten(a_true/2, g2p, del, Del, D0, bm2);
                Sa  = exp(-b_nom * Da * cos2_nom(id) - C_r);
                Se  = exp(-b_nom * (DeR_true + (DeL - DeR_true) * cos2_nom(id)));
                Sc  = exp(-b_nom * Dcsf);
                S_dirs(id) = (1-fcsf_true)*(f_true*Sa + (1-f_true)*Se) + fcsf_true*Sc;
            end
            S_baseline(is) = mean(S_dirs);

            % --- GNL: effective b and direction ---
            S_dirs = zeros(Ndirs, 1);
            for id = 1:Ndirs
                b_eff = b_nom * b_scale_per_dir(id);
                g_eff = g_nom * sqrt(b_scale_per_dir(id));
                g2p   = g_eff^2 * sin2_theta(id);

                C_r = vangelderen_atten(a_true/2, g2p, del, Del, D0, bm2);
                Sa  = exp(-b_eff * Da * cos2_theta(id) - C_r);
                Se  = exp(-b_eff * (DeR_true + (DeL - DeR_true) * cos2_theta(id)));
                Sc  = exp(-b_eff * Dcsf);
                S_dirs(id) = (1-fcsf_true)*(f_true*Sa + (1-f_true)*Se) + fcsf_true*Sc;
            end
            S_gnl(is) = mean(S_dirs);
        end

        % --- Fit ---
        smt = gpuAxCaliberSMT(b_shells, delta_shells, Delta_shells, D0, Da, DeL, Dcsf);

        pars_bl = fit_smt_lsq(smt, S_baseline, bm2);
        res.a_baseline(ig) = pars_bl.a;
        res.f_baseline(ig) = pars_bl.f;

        pars_s1 = fit_smt_lsq(smt, S_gnl, bm2);
        res.a_nocorr(ig) = pars_s1.a;
        res.f_nocorr(ig) = pars_s1.f;

        pars_s2 = fit_smt_lsq_bscale(smt, S_gnl, bm2, b_scale_scalar);
        res.a_bscale(ig) = pars_s2.a;
        res.f_bscale(ig) = pars_s2.f;

        fprintf('%s | GNL=%.1f  b_scale=%.4f | a: bl=%.2f nc=%.2f bs=%.2f | f: bl=%.3f nc=%.3f bs=%.3f\n', ...
            pname, gnl_scales(ig), b_scale_scalar, ...
            res.a_baseline(ig), res.a_nocorr(ig), res.a_bscale(ig), ...
            res.f_baseline(ig), res.f_nocorr(ig), res.f_bscale(ig));
    end

    all_results.(pname) = res;
end

%% 7. Plot results — C1 vs C2 side by side
figure('Units','inches','Position',[0 0 14 10]);

pnames = {'C2', 'C1'};
colors = struct('baseline', [0 0 0], 'nocorr', [0.8 0.2 0.2], 'bscale', [0.2 0.4 0.8]);

for ip = 1:2
    res = all_results.(pnames{ip});

    % Axon diameter
    subplot(2,2,(ip-1)*2+1);
    plot(res.gnl_scales, res.a_baseline, '-o', 'Color', colors.baseline, 'LineWidth', 2, 'DisplayName', 'Baseline (no GNL)'); hold on;
    plot(res.gnl_scales, res.a_nocorr,   '-s', 'Color', colors.nocorr,  'LineWidth', 2, 'DisplayName', 'GNL, no correction');
    plot(res.gnl_scales, res.a_bscale,   '-^', 'Color', colors.bscale,  'LineWidth', 2, 'DisplayName', 'GNL, b-scale correction');
    yline(a_true, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('GNL severity (\times L_{base})');
    ylabel('Estimated axon diameter [\mum]');
    title(sprintf('%s — Axon diameter', pnames{ip}));
    legend('Location', 'best'); grid on; box on;
    set(gca, 'FontSize', 12);

    % Volume fraction
    subplot(2,2,(ip-1)*2+2);
    plot(res.gnl_scales, res.f_baseline, '-o', 'Color', colors.baseline, 'LineWidth', 2, 'DisplayName', 'Baseline (no GNL)'); hold on;
    plot(res.gnl_scales, res.f_nocorr,   '-s', 'Color', colors.nocorr,  'LineWidth', 2, 'DisplayName', 'GNL, no correction');
    plot(res.gnl_scales, res.f_bscale,   '-^', 'Color', colors.bscale,  'LineWidth', 2, 'DisplayName', 'GNL, b-scale correction');
    yline(f_true, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('GNL severity (\times L_{base})');
    ylabel('Estimated restricted fraction');
    title(sprintf('%s — Restricted fraction', pnames{ip}));
    legend('Location', 'best'); grid on; box on;
    set(gca, 'FontSize', 12);
end

exportgraphics(gcf, 'fig_bvec_gnl_simulation.pdf', 'ContentType', 'vector');

%% ======================================================================
%  Helper functions
%  ======================================================================

function C = vangelderen_atten(r, g2_perp, delta, Delta, D0, bm2)
    % Van Gelderen restricted radial attenuation for a cylinder
    td = r^2 / D0;
    C = 0;
    for k = 1:numel(bm2)
        bk = bm2(k);
        bardelta = delta / td;
        barDelta = Delta / td;
        bkd = bk * bardelta;
        bkD = bk * barDelta;
        term = (2 / (bk^3 * (bk - 1))) * ...
               (-2 + 2*bkd + 2*exp(-bkd) + 2*exp(-bkD) ...
                - exp(-bkD - bkd) - exp(-bkD + bkd));
        C = C + term;
    end
    C = C * D0 * g2_perp * td^3;
end

function pars = fit_smt_lsq(smt, S_data, bm2)
    % Fit standard SMT model (no GNL correction)
    b     = double(smt.b);
    delta = double(smt.delta);
    Delta = double(smt.Delta);
    Da    = double(smt.Da);
    DeL   = double(smt.DeL);
    D0    = double(smt.D0);
    Dcsf  = double(smt.Dcsf);
    S_data = double(S_data(:));
    g2 = b ./ (delta.^2 .* (Delta - delta/3));

    cost = @(x) smt_residual(x, b, g2, delta, Delta, Da, DeL, D0, Dcsf, bm2, S_data);
    x0 = [3, 0.5, 0.05, 0.5];
    lb = [0.1, 0, 0, 0.01];
    ub = [20, 1, 1, DeL];
    opts = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 500);
    x_fit = lsqnonlin(cost, x0, lb, ub, opts);

    pars.a = x_fit(1); pars.f = x_fit(2);
    pars.fcsf = x_fit(3); pars.DeR = x_fit(4);
end

function pars = fit_smt_lsq_bscale(smt, S_data, bm2, b_scale)
    % Fit SMT model with scalar b_scale correction
    b     = double(smt.b) * b_scale;
    delta = double(smt.delta);
    Delta = double(smt.Delta);
    Da    = double(smt.Da);
    DeL   = double(smt.DeL);
    D0    = double(smt.D0);
    Dcsf  = double(smt.Dcsf);
    S_data = double(S_data(:));
    g2 = b ./ (delta.^2 .* (Delta - delta/3));

    cost = @(x) smt_residual(x, b, g2, delta, Delta, Da, DeL, D0, Dcsf, bm2, S_data);
    x0 = [3, 0.5, 0.05, 0.5];
    lb = [0.1, 0, 0, 0.01];
    ub = [20, 1, 1, DeL];
    opts = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 500);
    x_fit = lsqnonlin(cost, x0, lb, ub, opts);

    pars.a = x_fit(1); pars.f = x_fit(2);
    pars.fcsf = x_fit(3); pars.DeR = x_fit(4);
end

function res = smt_residual(x, b, g2, delta, Delta, Da, DeL, D0, Dcsf, bm2, S_data)
    % SMT forward model residual (powder-averaged signal)
    a = x(1); f = x(2); fcsf = x(3); DeR = x(4);
    r = a / 2;
    Nshells = numel(b);
    S_pred = zeros(Nshells, 1);

    for is = 1:Nshells
        C = vangelderen_atten_static(r, g2(is), delta(is), Delta(is), D0, bm2);
        arg = max(b(is)*Da - C, 1e-10);
        Sa = sqrt(pi / (4*arg)) * exp(-C) * erf(sqrt(arg));
        dDe = max(DeL - DeR, 1e-10);
        Se = sqrt(pi / (4*dDe*b(is))) * exp(-b(is)*DeR) * erf(sqrt(b(is)*dDe));
        Sc = exp(-b(is) * Dcsf);
        S_pred(is) = (1-fcsf) * (f*Sa + (1-f)*Se) + fcsf*Sc;
    end

    S_pred = min(S_pred, 1);
    res = S_pred - S_data;
end

function C = vangelderen_atten_static(r, g2_perp, delta, Delta, D0, bm2)
    td = r^2 / D0;
    C = 0;
    for k = 1:numel(bm2)
        bk = bm2(k);
        bkd = bk * delta / td;
        bkD = bk * Delta / td;
        term = (2 / (bk^3 * (bk - 1))) * ...
               (-2 + 2*bkd + 2*exp(-bkd) + 2*exp(-bkD) ...
                - exp(-bkD - bkd) - exp(-bkD + bkd));
        C = C + term;
    end
    C = C * D0 * g2_perp * td^3;
end
