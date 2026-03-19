function demo_bvec_gnl_simulation()
clear; close all;

%% Simulate AxCaliber-SMT signal with full GNL (b-scaling + b-rotation)
%
% Three scenarios compared for BOTH C1 and C2:
%   Baseline  : No GNL            -> fit with standard SMT
%   Scenario 1: GNL (scale+rot)   -> fit with NO correction
%   Scenario 2: GNL (scale+rot)   -> fit with scalar b_scale correction
%
% Two GNL regimes: b_scale > 1 (L_base) and b_scale < 1 (-L_base)
% Output: bar plot of bias for [a, f] x [C1, C2] x [b>1, b<1]

% No external dependencies needed — all helpers are local functions

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

%% 3. Actual gradient directions from C2 b-table (sub-001.bvec)
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

%% 4. Fiber direction and Van Gelderen Bessel zeros
fiber_dir = [0; 0; 1];
bm  = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

%% 5. GNL L-matrix
% L is the gradient deviation tensor: g_eff = (I + L) * g_nominal
% Realistic anisotropic L (~3-4% diagonal, ~1% off-diagonal)
L_base = [ 0.03   0.01   0.005;
           0.01  -0.02   0.008;
           0.005  0.008  0.04];

gnl_scale = 1.0;   % representative GNL severity

% Two regimes: +L gives b_scale > 1, -L gives b_scale < 1
regimes = struct('label', {'b_{scale} > 1', 'b_{scale} < 1'}, ...
                 'L',     {gnl_scale * L_base, -gnl_scale * L_base});

%% 6. Protocols
protocols = struct( ...
    'name',  {'C2', 'C1'}, ...
    'b',     {b_C2, b_C1}, ...
    'delta', {delta_C2, delta_C1}, ...
    'Delta', {Delta_C2, Delta_C1});

%% 7. Run simulation
% Store bias: bias_a(ip, ir, 1:3) and bias_f(ip, ir, 1:3)
%   ip=1:C2, ip=2:C1;  ir=1:b>1, ir=2:b<1;  dim3 = [baseline, nocorr, bscale]
bias_a = zeros(2, 2, 3);
bias_f = zeros(2, 2, 3);
bscale_values = zeros(2, 2);  % for reporting

for ip = 1:2
    pname  = protocols(ip).name;
    b_sh   = protocols(ip).b;
    del_sh = protocols(ip).delta;
    Del_sh = protocols(ip).Delta;
    Nsh    = numel(b_sh);
    g_sh   = sqrt(b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3)));

    for ir = 1:2
        L = regimes(ir).L;
        M = eye(3) + L;

        % Per-direction effective quantities
        g_eff_all       = M * bvecs;                       % [3 x Ndirs]
        b_scale_per_dir = sum(g_eff_all.^2, 1);            % |M*g_hat|^2
        g_eff_hat       = g_eff_all ./ vecnorm(g_eff_all); % unit effective dirs

        % Scalar b_scale = trace(M'*M)/3
        b_scale_scalar = trace(M' * M) / 3;
        bscale_values(ip, ir) = b_scale_scalar;

        % Angles with GNL
        cos_theta  = abs(g_eff_hat' * fiber_dir);
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

            % Baseline (no GNL): directional signal, averaged over directions
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
        g2_sh = b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3));

        % Baseline: no GNL, standard fit
        p_bl = fit_smt(b_sh, g2_sh, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_baseline);

        % Scenario 1: GNL, no correction
        p_nc = fit_smt(b_sh, g2_sh, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_gnl);

        % Scenario 2: GNL, scalar b_scale correction
        b_corr  = b_sh * b_scale_scalar;
        g2_corr = b_corr ./ (del_sh.^2 .* (Del_sh - del_sh/3));
        p_bs = fit_smt(b_corr, g2_corr, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_gnl);

        % Store bias
        bias_a(ip, ir, :) = [p_bl.a - a_true, p_nc.a - a_true, p_bs.a - a_true];
        bias_f(ip, ir, :) = [p_bl.f - f_true, p_nc.f - f_true, p_bs.f - f_true];

        fprintf('%s | %s | b_scale=%.4f | a bias: bl=%+.3f  nc=%+.3f  bs=%+.3f | f bias: bl=%+.4f  nc=%+.4f  bs=%+.4f\n', ...
            pname, regimes(ir).label, b_scale_scalar, ...
            p_bl.a-a_true, p_nc.a-a_true, p_bs.a-a_true, ...
            p_bl.f-f_true, p_nc.f-f_true, p_bs.f-f_true);
    end
end

%% 8. Bar plot
figure('Units','inches','Position',[1 1 13 5]);

scenarios = {'Baseline', 'No correction', 'b-scale correction'};
colors    = [0.3 0.3 0.3; 0.85 0.25 0.25; 0.25 0.45 0.85];

% Group order: C2 b>1, C2 b<1, C1 b>1, C1 b<1
group_labels = {sprintf('C2\n(b_{sc}=%.3f)', bscale_values(1,1)), ...
                sprintf('C2\n(b_{sc}=%.3f)', bscale_values(1,2)), ...
                sprintf('C1\n(b_{sc}=%.3f)', bscale_values(2,1)), ...
                sprintf('C1\n(b_{sc}=%.3f)', bscale_values(2,2))};

% Reshape: [4 groups x 3 scenarios]
a_data = [squeeze(bias_a(1,1,:))'; squeeze(bias_a(1,2,:))'; ...
          squeeze(bias_a(2,1,:))'; squeeze(bias_a(2,2,:))'];
f_data = [squeeze(bias_f(1,1,:))'; squeeze(bias_f(1,2,:))'; ...
          squeeze(bias_f(2,1,:))'; squeeze(bias_f(2,2,:))'];

% --- Axon diameter bias ---
subplot(1,2,1);
bh = bar(a_data, 'grouped');
for k = 1:3, bh(k).FaceColor = colors(k,:); end
set(gca, 'XTickLabel', group_labels, 'FontSize', 11);
ylabel('Axon diameter bias [\mum]');
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
legend(scenarios, 'Location', 'best', 'FontSize', 9);
title(sprintf('Axon diameter (truth = %.1f \\mum)', a_true));
box on; grid on;

% --- Restricted fraction bias ---
subplot(1,2,2);
bh = bar(f_data, 'grouped');
for k = 1:3, bh(k).FaceColor = colors(k,:); end
set(gca, 'XTickLabel', group_labels, 'FontSize', 11);
ylabel('Restricted fraction bias');
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
legend(scenarios, 'Location', 'best', 'FontSize', 9);
title(sprintf('Restricted fraction (truth = %.2f)', f_true));
box on; grid on;

sgtitle(sprintf('GNL bias simulation (gnl\\_scale = %.1f)', gnl_scale), 'FontSize', 13);
exportgraphics(gcf, 'fig4_bar_bias_rmse.pdf', 'ContentType', 'vector');
fprintf('\nSaved fig4_bar_bias_rmse.pdf\n');

end  % main function

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

function pars = fit_smt(b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data)
    S_data = double(S_data(:));
    cost = @(x) smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data);
    opts = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',500, ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12);
    x = lsqnonlin(cost, [3 0.5 0.05 0.5], [0.1 0 0 0.01], [20 1 1 DeL], opts);
    pars = struct('a',x(1),'f',x(2),'fcsf',x(3),'DeR',x(4));
end

function res = smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data)
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
