clear; close all;

%% Simulate AxCaliber-SMT signal with full GNL (b-scaling + b-rotation)
%
% Three scenarios compared for BOTH C1 and C2:
%   Baseline  : No GNL            → fit with standard SMT
%   Scenario 1: GNL (scale+rot)   → fit with NO correction
%   Scenario 2: GNL (scale+rot)   → fit with scalar b_scale correction
%
% Signal generation uses the directional model (angle-dependent),
% then averaged per shell to obtain the spherical mean.
% Fitting uses the standard SMT (orientation-free) forward model.

addpath(genpath('/Users/quantum/gacelle/AxCaliberSMT'));
addpath(genpath('/Users/quantum/gacelle/utils'));

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

%% 2. Protocol definitions
% C2 protocol
b_C2     = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
            0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';   % [ms/um2]
delta_C2 = 8 * ones(16,1);                                % [ms]
Delta_C2 = [19*ones(8,1); 49*ones(8,1)];                  % [ms]

% C1 protocol — same timing, b-values scaled by (Gmax_C1/Gmax_C2)^2
Gmax_C1 = 290; Gmax_C2 = 495;
b_C1     = b_C2 * (Gmax_C1/Gmax_C2)^2;
delta_C1 = delta_C2;
Delta_C1 = Delta_C2;

%% 3. Actual gradient directions from C2 b-table (sub-001.bvec, first b=50 shell)
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
% bvecs: [3 x 16], each column is a unit gradient direction
Ndirs = size(bvecs, 2);

%% 4. Fiber direction
fiber_dir = [0; 0; 1];   % along z-axis

%% 5. GNL L-matrix sweep
% L is the gradient deviation tensor: g_eff = (I + L) * g_nominal
% Realistic anisotropic L (~3-4% diagonal, ~1% off-diagonal)
L_base = [ 0.03   0.01   0.005;
           0.01  -0.02   0.008;
           0.005  0.008  0.04];

gnl_scales = 0:0.5:3.0;   % 0 = no GNL, 3 = severe (~12% deviation)
Ngnl = numel(gnl_scales);

%% 6. Van Gelderen Bessel function zeros
bm = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

%% 7. Run simulation for both protocols
protocols = struct( ...
    'name',  {'C2', 'C1'}, ...
    'b',     {b_C2, b_C1}, ...
    'delta', {delta_C2, delta_C1}, ...
    'Delta', {Delta_C2, Delta_C1});

all_results = struct();

for ip = 1:2
    pname  = protocols(ip).name;
    b_sh   = protocols(ip).b;
    del_sh = protocols(ip).delta;
    Del_sh = protocols(ip).Delta;
    Nsh    = numel(b_sh);
    g_sh   = sqrt(b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3)));

    res.gnl_scales   = gnl_scales;
    res.b_scale_true = zeros(Ngnl, 1);
    res.a_baseline   = zeros(Ngnl, 1);  res.f_baseline   = zeros(Ngnl, 1);
    res.a_nocorr     = zeros(Ngnl, 1);  res.f_nocorr     = zeros(Ngnl, 1);
    res.a_bscale     = zeros(Ngnl, 1);  res.f_bscale     = zeros(Ngnl, 1);
    res.fcsf_baseline = zeros(Ngnl,1);  res.DeR_baseline  = zeros(Ngnl,1);
    res.fcsf_nocorr   = zeros(Ngnl,1);  res.DeR_nocorr    = zeros(Ngnl,1);
    res.fcsf_bscale   = zeros(Ngnl,1);  res.DeR_bscale    = zeros(Ngnl,1);

    for ig = 1:Ngnl
        L = gnl_scales(ig) * L_base;
        M = eye(3) + L;

        % Per-direction effective quantities
        g_eff_all       = M * bvecs;                      % [3 x Ndirs]
        b_scale_per_dir = sum(g_eff_all.^2, 1);           % |M*g_hat|^2, [1 x Ndirs]
        g_eff_hat       = g_eff_all ./ vecnorm(g_eff_all); % unit effective dirs

        % Scalar b_scale = trace(M'*M)/3
        b_scale_scalar = trace(M' * M) / 3;
        res.b_scale_true(ig) = b_scale_scalar;

        % Angles with GNL
        cos_theta  = abs(g_eff_hat' * fiber_dir);   % [Ndirs x 1]
        sin2_theta = 1 - cos_theta.^2;
        cos2_theta = cos_theta.^2;

        % Angles without GNL (nominal)
        cos_theta_nom = abs(bvecs' * fiber_dir);
        sin2_nom = 1 - cos_theta_nom.^2;
        cos2_nom = cos_theta_nom.^2;

        % ---- Generate shell-averaged signals ----
        S_baseline = zeros(Nsh, 1);
        S_gnl      = zeros(Nsh, 1);

        for is = 1:Nsh
            b0  = b_sh(is);
            g0  = g_sh(is);
            del = del_sh(is);
            Del = Del_sh(is);

            % Baseline (no GNL)
            Sd = zeros(Ndirs, 1);
            for id = 1:Ndirs
                g2p = g0^2 * sin2_nom(id);
                C_r = vg_atten(a_true/2, g2p, del, Del, D0, bm2);
                Sa  = exp(-b0*Da*cos2_nom(id) - C_r);
                Se  = exp(-b0*(DeR_true + (DeL-DeR_true)*cos2_nom(id)));
                Sc  = exp(-b0*Dcsf);
                Sd(id) = (1-fcsf_true)*(f_true*Sa + (1-f_true)*Se) + fcsf_true*Sc;
            end
            S_baseline(is) = mean(Sd);

            % With GNL (scale + rotation)
            Sd = zeros(Ndirs, 1);
            for id = 1:Ndirs
                b_eff = b0 * b_scale_per_dir(id);
                g_eff = g0 * sqrt(b_scale_per_dir(id));
                g2p   = g_eff^2 * sin2_theta(id);

                C_r = vg_atten(a_true/2, g2p, del, Del, D0, bm2);
                Sa  = exp(-b_eff*Da*cos2_theta(id) - C_r);
                Se  = exp(-b_eff*(DeR_true + (DeL-DeR_true)*cos2_theta(id)));
                Sc  = exp(-b_eff*Dcsf);
                Sd(id) = (1-fcsf_true)*(f_true*Sa + (1-f_true)*Se) + fcsf_true*Sc;
            end
            S_gnl(is) = mean(Sd);
        end

        % ---- Fit ----
        smt = gpuAxCaliberSMT(b_sh, del_sh, Del_sh, D0, Da, DeL, Dcsf);

        % Baseline
        p = fit_smt(smt, S_baseline, bm2);
        res.a_baseline(ig) = p.a;  res.f_baseline(ig) = p.f;
        res.fcsf_baseline(ig) = p.fcsf; res.DeR_baseline(ig) = p.DeR;

        % Scenario 1: GNL, no correction
        p = fit_smt(smt, S_gnl, bm2);
        res.a_nocorr(ig) = p.a;  res.f_nocorr(ig) = p.f;
        res.fcsf_nocorr(ig) = p.fcsf; res.DeR_nocorr(ig) = p.DeR;

        % Scenario 2: GNL, scalar b_scale correction
        p = fit_smt_bscale(smt, S_gnl, bm2, b_scale_scalar);
        res.a_bscale(ig) = p.a;  res.f_bscale(ig) = p.f;
        res.fcsf_bscale(ig) = p.fcsf; res.DeR_bscale(ig) = p.DeR;

        fprintf('%s | GNL=%.1f  b_scale=%.4f | a: bl=%.2f nc=%.2f bs=%.2f | f: bl=%.3f nc=%.3f bs=%.3f\n', ...
            pname, gnl_scales(ig), b_scale_scalar, ...
            res.a_baseline(ig), res.a_nocorr(ig), res.a_bscale(ig), ...
            res.f_baseline(ig), res.f_nocorr(ig), res.f_bscale(ig));
    end

    all_results.(pname) = res;
end

%% 8. Plot
figure('Units','inches','Position',[0 0 14 10]);
pnames = {'C2', 'C1'};
col_bl = [0 0 0]; col_nc = [0.8 0.2 0.2]; col_bs = [0.2 0.4 0.8];

for ip = 1:2
    r = all_results.(pnames{ip});

    % Axon diameter
    subplot(2,2,(ip-1)*2+1);
    plot(r.gnl_scales, r.a_baseline, '-o', 'Color', col_bl, 'LineWidth', 2, 'DisplayName', 'Baseline (no GNL)'); hold on;
    plot(r.gnl_scales, r.a_nocorr,   '-s', 'Color', col_nc, 'LineWidth', 2, 'DisplayName', 'GNL, no correction');
    plot(r.gnl_scales, r.a_bscale,   '-^', 'Color', col_bs, 'LineWidth', 2, 'DisplayName', 'GNL, b-scale corr.');
    yline(a_true, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('GNL severity (\times L_{base})'); ylabel('Axon diameter [\mum]');
    title(sprintf('%s — Axon diameter (truth = %.1f \\mum)', pnames{ip}, a_true));
    legend('Location', 'best'); grid on; box on; set(gca, 'FontSize', 12);

    % Volume fraction
    subplot(2,2,(ip-1)*2+2);
    plot(r.gnl_scales, r.f_baseline, '-o', 'Color', col_bl, 'LineWidth', 2, 'DisplayName', 'Baseline (no GNL)'); hold on;
    plot(r.gnl_scales, r.f_nocorr,   '-s', 'Color', col_nc, 'LineWidth', 2, 'DisplayName', 'GNL, no correction');
    plot(r.gnl_scales, r.f_bscale,   '-^', 'Color', col_bs, 'LineWidth', 2, 'DisplayName', 'GNL, b-scale corr.');
    yline(f_true, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('GNL severity (\times L_{base})'); ylabel('Restricted fraction');
    title(sprintf('%s — Restricted fraction (truth = %.2f)', pnames{ip}, f_true));
    legend('Location', 'best'); grid on; box on; set(gca, 'FontSize', 12);
end

sgtitle(sprintf('GNL simulation: a_{true}=%.1f\\mum, f=%.2f, f_{csf}=%.2f, D_{eR}=%.1f', ...
    a_true, f_true, fcsf_true, DeR_true), 'FontSize', 14);
exportgraphics(gcf, 'fig_bvec_gnl_simulation.pdf', 'ContentType', 'vector');

%% Print summary table
fprintf('\n===== Summary =====\n');
fprintf('%-4s  %-5s  %-7s  %-12s  %-12s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Prot','GNL','b_scale','a_baseline','a_nocorr','a_bscale','f_baseline','f_nocorr','f_bscale');
for ip = 1:2
    r = all_results.(pnames{ip});
    for ig = 1:Ngnl
        fprintf('%-4s  %4.1f  %7.4f  %10.3f  %10.3f  %10.3f  %10.4f  %10.4f  %10.4f\n', ...
            pnames{ip}, r.gnl_scales(ig), r.b_scale_true(ig), ...
            r.a_baseline(ig), r.a_nocorr(ig), r.a_bscale(ig), ...
            r.f_baseline(ig), r.f_nocorr(ig), r.f_bscale(ig));
    end
end

%% ======================================================================
%  Helper functions
%  ======================================================================

function C = vg_atten(r, g2_perp, delta, Delta, D0, bm2)
    % Van Gelderen restricted radial attenuation for a cylinder
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

function pars = fit_smt(smt, S_data, bm2)
    % Fit standard SMT model using lsqnonlin
    b  = double(smt.b);  del = double(smt.delta);  Del = double(smt.Delta);
    Da = double(smt.Da); DeL = double(smt.DeL);
    D0 = double(smt.D0); Dcsf = double(smt.Dcsf);
    g2 = b ./ (del.^2 .* (Del - del/3));
    S_data = double(S_data(:));

    cost = @(x) smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data);
    opts = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',500, ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12);
    x = lsqnonlin(cost, [3 0.5 0.05 0.5], [0.1 0 0 0.01], [20 1 1 DeL], opts);
    pars = struct('a',x(1),'f',x(2),'fcsf',x(3),'DeR',x(4));
end

function pars = fit_smt_bscale(smt, S_data, bm2, b_scale)
    % Fit SMT model with scalar b_scale correction
    b  = double(smt.b) * b_scale;
    del = double(smt.delta);  Del = double(smt.Delta);
    Da = double(smt.Da); DeL = double(smt.DeL);
    D0 = double(smt.D0); Dcsf = double(smt.Dcsf);
    g2 = b ./ (del.^2 .* (Del - del/3));
    S_data = double(S_data(:));

    cost = @(x) smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data);
    opts = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',500, ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12);
    x = lsqnonlin(cost, [3 0.5 0.05 0.5], [0.1 0 0 0.01], [20 1 1 DeL], opts);
    pars = struct('a',x(1),'f',x(2),'fcsf',x(3),'DeR',x(4));
end

function res = smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data)
    % SMT powder-averaged forward model residual
    r = x(1)/2; f = x(2); fcsf = x(3); DeR = x(4);
    N = numel(b);
    S = zeros(N,1);
    for is = 1:N
        C   = vg_atten_s(r, g2(is), del(is), Del(is), D0, bm2);
        arg = max(b(is)*Da - C, 1e-10);
        Sa  = sqrt(pi/(4*arg)) * exp(-C) * erf(sqrt(arg));
        dDe = max(DeL - DeR, 1e-10);
        Se  = sqrt(pi/(4*dDe*b(is))) * exp(-b(is)*DeR) * erf(sqrt(b(is)*dDe));
        S(is) = (1-fcsf)*(f*Sa + (1-f)*Se) + fcsf*exp(-b(is)*Dcsf);
    end
    res = min(S,1) - S_data;
end

function C = vg_atten_s(r, g2p, del, Del, D0, bm2)
    td = r^2/D0; C = 0;
    for k = 1:numel(bm2)
        bk = bm2(k); bkd = bk*del/td; bkD = bk*Del/td;
        C = C + (2/(bk^3*(bk-1)))*(-2+2*bkd+2*exp(-bkd)+2*exp(-bkD)-exp(-bkD-bkd)-exp(-bkD+bkd));
    end
    C = C*D0*g2p*td^3;
end
