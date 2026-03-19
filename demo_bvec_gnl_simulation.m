function demo_bvec_gnl_simulation()
clear; close all;

%% GNL simulation with full b-vector effects (scaling + rotation)
%
% Three scenarios × two protocols (C2, C1) × two GNL regimes (b>1, b<1):
%   Baseline  : No GNL in signal  -> fit with nominal b
%   Scenario 1: GNL in signal     -> fit with NO correction
%   Scenario 2: GNL in signal     -> fit with scalar b_scale correction
%
% Signal generation: directional model (per-gradient VG + zeppelin), averaged per shell
% Fitting: SMT powder-average forward model (lsqnonlin)
% No external dependencies.

%% 1. Fixed model parameters
D0 = 2.0; Da = 1.7; DeL = 1.7; Dcsf = 3.0;

%% 2. Protocol definitions
b_C2     = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
            0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta_C2 = 8 * ones(16,1);
Delta_C2 = [19*ones(8,1); 49*ones(8,1)];

b_C1     = b_C2 * (290/495)^2;
delta_C1 = delta_C2;
Delta_C1 = Delta_C2;

%% 3. Gradient directions from C2 b-table (sub-001.bvec)
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

%% 4. Fiber direction and Bessel zeros
fiber_dir = [0; 0; 1];
bm  = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 24.3113 27.4571 30.6019];
bm2 = bm.^2;

%% 5. GNL L-matrix
L_base = [0.03 0.01 0.005; 0.01 -0.02 0.008; 0.005 0.008 0.04];
gnl_scale = 1.0;

regimes = struct('label', {'b_{sc}>1', 'b_{sc}<1'}, ...
                 'L', {gnl_scale*L_base, -gnl_scale*L_base});

%% 6. Protocols
protocols = struct('name', {'C2','C1'}, ...
    'b',     {b_C2, b_C1}, ...
    'delta', {delta_C2, delta_C1}, ...
    'Delta', {Delta_C2, Delta_C1}, ...
    'snr',   {38, 18});

%% 7. Random ground truth (matching demo_bscale_noise_propagation ranges)
Nsample = 1000;
rng(42);
Xrange = [0.1 5; 0.5 1; 0 0.2; 0.5 1.5];  % [r, f, fcsf, DeR]
Xinit  = Xrange(:,1) + diff(Xrange,1,2) .* rand(4, Nsample);
r_gt    = Xinit(1,:)';   % radius [um]
f_gt    = Xinit(2,:)';
fcsf_gt = Xinit(3,:)';
DeR_gt  = Xinit(4,:)';

%% 8. Pre-compute angular quantities
cos_theta_nom = abs(bvecs' * fiber_dir);     % [Ndirs x 1]
sin2_nom      = 1 - cos_theta_nom.^2;

for ir = 1:2
    M = eye(3) + regimes(ir).L;
    regimes(ir).b_scale_scalar = trace(M'*M) / 3;
    g_eff = M * bvecs;
    regimes(ir).b_scale_per_dir = sum(g_eff.^2, 1);     % [1 x Ndirs]
    g_hat = g_eff ./ vecnorm(g_eff);
    regimes(ir).cos_theta = abs(g_hat' * fiber_dir);     % [Ndirs x 1]
    regimes(ir).sin2      = 1 - regimes(ir).cos_theta.^2;
end

%% 9. Main simulation loop
results = struct();
for ip = 1:2
    for ir = 1:2
        results(ip,ir).r_bl = zeros(Nsample,1);
        results(ip,ir).r_nc = zeros(Nsample,1);
        results(ip,ir).r_bs = zeros(Nsample,1);
        results(ip,ir).f_bl = zeros(Nsample,1);
        results(ip,ir).f_nc = zeros(Nsample,1);
        results(ip,ir).f_bs = zeros(Nsample,1);
    end
end

for ip = 1:2
    b_sh   = protocols(ip).b;
    del_sh = protocols(ip).delta;
    Del_sh = protocols(ip).Delta;
    snr    = protocols(ip).snr;
    Nsh    = numel(b_sh);
    g_sh   = sqrt(b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3)));
    g2_sh  = b_sh ./ (del_sh.^2 .* (Del_sh - del_sh/3));

    for ir = 1:2
        bsc    = regimes(ir).b_scale_scalar;
        bsd    = regimes(ir).b_scale_per_dir;
        ct_gnl = regimes(ir).cos_theta;
        s2_gnl = regimes(ir).sin2;

        b_corr  = b_sh * bsc;
        g2_corr = b_corr ./ (del_sh.^2 .* (Del_sh - del_sh/3));

        fprintf('\n%s | %s | b_scale=%.4f | SNR=%d\n', ...
                protocols(ip).name, regimes(ir).label, bsc, snr);
        tic;

        for i = 1:Nsample
            if mod(i, 200) == 0
                fprintf('  %d/%d (%.0fs)\n', i, Nsample, toc);
            end

            r = r_gt(i);

            % --- Baseline signal (no GNL) ---
            S_bl = zeros(Nsh, 1);
            for is = 1:Nsh
                Sd = zeros(Ndirs, 1);
                for id = 1:Ndirs
                    g2p = g_sh(is)^2 * sin2_nom(id);
                    C_r = vg_atten(r, g2p, del_sh(is), Del_sh(is), D0, bm2);
                    Sa  = exp(-b_sh(is)*Da*cos_theta_nom(id)^2 - C_r);
                    Se  = exp(-b_sh(is)*(DeR_gt(i) + (DeL-DeR_gt(i))*cos_theta_nom(id)^2));
                    Sc  = exp(-b_sh(is)*Dcsf);
                    Sd(id) = (1-fcsf_gt(i))*(f_gt(i)*Sa + (1-f_gt(i))*Se) + fcsf_gt(i)*Sc;
                end
                S_bl(is) = mean(Sd);
            end
            S_bl_noisy = abs(S_bl + (1/snr)*randn(Nsh,1) + 1j*(1/snr)*randn(Nsh,1));

            % --- GNL signal (scale + rotation) ---
            S_gnl = zeros(Nsh, 1);
            for is = 1:Nsh
                Sd = zeros(Ndirs, 1);
                for id = 1:Ndirs
                    b_eff = b_sh(is) * bsd(id);
                    g_eff = g_sh(is) * sqrt(bsd(id));
                    g2p   = g_eff^2 * s2_gnl(id);
                    C_r = vg_atten(r, g2p, del_sh(is), Del_sh(is), D0, bm2);
                    Sa  = exp(-b_eff*Da*ct_gnl(id)^2 - C_r);
                    Se  = exp(-b_eff*(DeR_gt(i) + (DeL-DeR_gt(i))*ct_gnl(id)^2));
                    Sc  = exp(-b_eff*Dcsf);
                    Sd(id) = (1-fcsf_gt(i))*(f_gt(i)*Sa + (1-f_gt(i))*Se) + fcsf_gt(i)*Sc;
                end
                S_gnl(is) = mean(Sd);
            end
            S_gnl_noisy = abs(S_gnl + (1/snr)*randn(Nsh,1) + 1j*(1/snr)*randn(Nsh,1));

            % --- Fit ---
            p = fit_smt(b_sh, g2_sh, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_bl_noisy);
            results(ip,ir).r_bl(i) = p.r; results(ip,ir).f_bl(i) = p.f;

            p = fit_smt(b_sh, g2_sh, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_gnl_noisy);
            results(ip,ir).r_nc(i) = p.r; results(ip,ir).f_nc(i) = p.f;

            p = fit_smt(b_corr, g2_corr, del_sh, Del_sh, Da, DeL, D0, Dcsf, bm2, S_gnl_noisy);
            results(ip,ir).r_bs(i) = p.r; results(ip,ir).f_bs(i) = p.f;
        end
        fprintf('  Done: %.0f s\n', toc);
    end
end

%% 10. Filter boundary hits
r_lb = 0.05; r_ub = 4.95;
for ip = 1:2
    for ir = 1:2
        R = results(ip,ir);
        valid = R.r_bl > r_lb & R.r_bl < r_ub & ...
                R.r_nc > r_lb & R.r_nc < r_ub & ...
                R.r_bs > r_lb & R.r_bs < r_ub;
        results(ip,ir).valid = valid;
        fprintf('%s | %s: %d/%d valid\n', protocols(ip).name, regimes(ir).label, ...
                sum(valid), Nsample);
    end
end

%% 11. Fig 2 — Error vs r_true (scatter)
figure('Units','inches','Position',[0 0 12 8]);
col_bl = [0.3 0.3 0.3]; col_nc = [0.85 0.25 0.25]; col_bs = [0.25 0.45 0.85];

for ip = 1:2
    for ir = 1:2
        subplot(2, 2, (ip-1)*2 + ir);
        v = results(ip,ir).valid;
        rt = r_gt(v);
        err_bl = results(ip,ir).r_bl(v) - rt;
        err_nc = results(ip,ir).r_nc(v) - rt;
        err_bs = results(ip,ir).r_bs(v) - rt;

        scatter(rt, err_bl, 6, col_bl, '.', 'MarkerFaceAlpha', 0.3); hold on;
        scatter(rt, err_nc, 6, col_nc, '.', 'MarkerFaceAlpha', 0.3);
        scatter(rt, err_bs, 6, col_bs, '.', 'MarkerFaceAlpha', 0.3);
        yline(0, 'k-', 'LineWidth', 1);
        xlabel('r_{true} (\mum)'); ylabel('r_{fit} - r_{true} (\mum)');
        title(sprintf('%s | %s (b_{sc}=%.3f)', protocols(ip).name, regimes(ir).label, ...
              regimes(ir).b_scale_scalar), 'FontSize', 10);
        grid on; box on;
        if ip == 1 && ir == 1
            legend({'Baseline','No correction','b-scale corr.'}, 'Location','best','FontSize',8);
        end
    end
end
exportgraphics(gcf, 'fig2_error_vs_rtrue.pdf', 'ContentType', 'vector');

%% 12. Fig 3 — Scatter: r_true vs r_fitted
% 3 rows (baseline, no corr, b-scale) x 4 cols (C2 b>1, C2 b<1, C1 b>1, C1 b<1)
figure('Units','inches','Position',[0 0 16 10]);
scenario_names = {'Baseline (no GNL)', 'GNL, no correction', 'GNL, b-scale corr.'};
fields_r = {'r_bl', 'r_nc', 'r_bs'};
rlim = Xrange(1,:);

for iscen = 1:3
    for ip = 1:2
        for ir = 1:2
            col_idx = (ip-1)*2 + ir;
            subplot(3, 4, (iscen-1)*4 + col_idx);
            v  = results(ip,ir).valid;
            rt = r_gt(v);
            rf = results(ip,ir).(fields_r{iscen})(v);
            err = rf - rt;

            scatter(rt, rf, 8, '.'); hold on;
            plot(rlim, rlim, 'k-', 'LineWidth', 1);
            xlim(rlim); ylim(rlim); pbaspect([1 1 1]); box on; grid on;
            xlabel('r_{true} (\mum)'); ylabel('r_{fit} (\mum)');
            title(sprintf('%s %s\n%s\nbias=%.3f RMSE=%.3f', ...
                  protocols(ip).name, regimes(ir).label, scenario_names{iscen}, ...
                  mean(err), sqrt(mean(err.^2))), 'FontSize', 8);
        end
    end
end
exportgraphics(gcf, 'fig3_scatter_gt_vs_fit.pdf', 'ContentType', 'vector');

%% 13. Fig 4 — Bar plot: bias by r-bin
r_bins = [0 1; 1 2; 2 3; 3 4; 4 5];
Nbins  = size(r_bins, 1);
bin_labels = arrayfun(@(i) sprintf('%g-%g', r_bins(i,1), r_bins(i,2)), 1:Nbins, 'uni', 0);

figure('Units','inches','Position',[0 0 14 10]);
colors = [col_bl; col_nc; col_bs];
scenario_short = {'Baseline', 'No corr.', 'b-scale corr.'};

for ip = 1:2
    for ir = 1:2
        subplot(2, 2, (ip-1)*2 + ir);
        v  = results(ip,ir).valid;
        rt = r_gt(v);

        bias_mat = zeros(Nbins, 3);  % [bins x scenarios]
        for iscen = 1:3
            rf = results(ip,ir).(fields_r{iscen})(v);
            for ib = 1:Nbins
                idx = rt >= r_bins(ib,1) & rt < r_bins(ib,2);
                if sum(idx) > 0
                    bias_mat(ib, iscen) = mean(rf(idx) - rt(idx));
                end
            end
        end

        bh = bar(categorical(bin_labels, bin_labels), bias_mat, 'grouped');
        for k = 1:3, bh(k).FaceColor = colors(k,:); end
        ylabel('Bias (\mum)'); xlabel('r_{true} (\mum)');
        yline(0, 'k--', 'HandleVisibility', 'off');
        title(sprintf('%s | %s (b_{sc}=%.3f)', protocols(ip).name, regimes(ir).label, ...
              regimes(ir).b_scale_scalar), 'FontSize', 10);
        legend(scenario_short, 'Location', 'best', 'FontSize', 8);
        grid on; box on;
    end
end
exportgraphics(gcf, 'fig4_bar_bias_rmse.pdf', 'ContentType', 'vector');

%% 14. Print summary
fprintf('\n===== Summary =====\n');
fprintf('%-4s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Prot','Regime','b_scale','bias_bl','bias_nc','bias_bs','rmse_bl','rmse_nc','rmse_bs');
for ip = 1:2
    for ir = 1:2
        v = results(ip,ir).valid;
        rt = r_gt(v);
        e_bl = results(ip,ir).r_bl(v) - rt;
        e_nc = results(ip,ir).r_nc(v) - rt;
        e_bs = results(ip,ir).r_bs(v) - rt;
        fprintf('%-4s  %-8s  %7.4f  %+9.4f  %+9.4f  %+9.4f  %9.4f  %9.4f  %9.4f\n', ...
            protocols(ip).name, regimes(ir).label, regimes(ir).b_scale_scalar, ...
            mean(e_bl), mean(e_nc), mean(e_bs), ...
            sqrt(mean(e_bl.^2)), sqrt(mean(e_nc.^2)), sqrt(mean(e_bs.^2)));
    end
end

%% 15. Save results
save('bvec_gnl_simulation_results.mat', 'results', 'protocols', 'regimes', ...
     'r_gt', 'f_gt', 'fcsf_gt', 'DeR_gt', 'Xrange', 'Nsample', 'gnl_scale');
fprintf('\nSaved results to bvec_gnl_simulation_results.mat\n');

end  % main function

%% ======================================================================
%  Helper functions
%  ======================================================================

function C = vg_atten(r, g2_perp, delta, Delta, D0, bm2)
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
    x = lsqnonlin(cost, [2 0.7 0.05 0.8], [0 0 0 0.01], [5 1 1 DeL], opts);
    pars = struct('r',x(1),'f',x(2),'fcsf',x(3),'DeR',x(4));
end

function res = smt_res(x, b, g2, del, Del, Da, DeL, D0, Dcsf, bm2, S_data)
    r = x(1); f = x(2); fcsf = x(3); DeR = x(4);
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
