clear; close all;

%% T2 decay weighted vs original restricted fraction (data-driven)
% Reference: Veraart, Novikov, Fieremans (2018) — TEdDI, NeuroImage
%
% Primary analysis: 2-compartment model (fCSF = 0)
%   f_r(TE) = f0*exp(-TE/T2a) / [f0*exp(-TE/T2a) + (1-f0)*exp(-TE/T2e)]
%
% Sensitivity analysis: 3-compartment with fCSF = 0.02, 0.05, 0.10

% Protocol TEs
TE_C2 = 54;   % ms
TE_C1 = 77;   % ms
T2csf = 2000; % ms (at 3T, used only in sensitivity analysis)

%% 1. Per-tract data: T2 values (mapped from TEdDI ROIs) + observed f_r
% T2 values from TEdDI paper Figure 5 (Veraart et al. 2018)
% Observed f_r read from Supplementary Figure 3(b) (mean across 20 subjects)
%
% Mapping: tractography atlas ROI → TEdDI ICBM-DTI-81 ROI
%   ATR  → ACR (anterior corona radiata)
%   CST  → PLIC (posterior limb of internal capsule)
%   Fmaj → SCC (splenium of corpus callosum)
%   Fmin → GCC (genu of corpus callosum)
%   IFO  → EC  (external capsule)
%   SLF  → SLF (superior longitudinal fasciculus)

%            Tract name       T2a  T2e   f_obs_C1  f_obs_C2
roi_data = {
    'ATR L',                  75,  48,   0.57,     0.52;
    'ATR R',                  75,  48,   0.55,     0.49;
    'CST L',                 110,  42,   0.72,     0.55;
    'CST R',                 110,  42,   0.70,     0.54;
    'Forceps major',          90,  40,   0.57,     0.55;
    'Forceps minor',          75,  35,   0.56,     0.55;
    'IFO L',                  70,  55,   0.47,     0.40;
    'IFO R',                  70,  55,   0.50,     0.48;
    'ILF L',                  80,  45,   0.55,     0.55;
    'ILF R',                  80,  45,   0.55,     0.55;
    'SLF L',                  90,  38,   0.56,     0.55;
    'SLF R',                  90,  38,   0.67,     0.63;
};

roi_names  = roi_data(:,1);
T2a_all    = cell2mat(roi_data(:,2));
T2e_all    = cell2mat(roi_data(:,3));
f_obs_C1   = cell2mat(roi_data(:,4));
f_obs_C2   = cell2mat(roi_data(:,5));
Nroi       = numel(roi_names);

% Observed % difference (C2 - C1, negative = C2 lower)
obs_diff = (f_obs_C2 - f_obs_C1) ./ f_obs_C1 * 100;

%% ===================================================================
%%  PRIMARY ANALYSIS: 2-compartment (fCSF = 0)
%% ===================================================================

%% 2. Back-calculate f0 from observed C1 f_r
% Inverse: f0 = f*exp(-TE/T2e) / [f*exp(-TE/T2e) + (1-f)*exp(-TE/T2a)]
f0_from_C1 = zeros(Nroi, 1);
for k = 1:Nroi
    f  = f_obs_C1(k);
    ea = exp(-TE_C1 / T2a_all(k));
    ee = exp(-TE_C1 / T2e_all(k));
    f0_from_C1(k) = f * ee / (f * ee + (1-f) * ea);
end

%% 3. Forward predict f_r at both TEs
f_pred_C1 = zeros(Nroi, 1);
f_pred_C2 = zeros(Nroi, 1);
for k = 1:Nroi
    f0  = f0_from_C1(k);
    ea1 = exp(-TE_C1 / T2a_all(k));
    ee1 = exp(-TE_C1 / T2e_all(k));
    f_pred_C1(k) = f0*ea1 / (f0*ea1 + (1-f0)*ee1);

    ea2 = exp(-TE_C2 / T2a_all(k));
    ee2 = exp(-TE_C2 / T2e_all(k));
    f_pred_C2(k) = f0*ea2 / (f0*ea2 + (1-f0)*ee2);
end

pred_diff   = (f_pred_C2 - f_pred_C1) ./ f_pred_C1 * 100;
unexplained = obs_diff - pred_diff;

%% 4. Print summary table
fprintf('=== Primary analysis: fCSF = 0 (2-compartment) ===\n');
fprintf('%-16s  f0     fC1obs  fC2obs  fC2pred  obs%%    pred%%   resid%%\n', 'Tract');
fprintf('%-16s  -----  ------  ------  -------  ------  ------  ------\n', '');
for k = 1:Nroi
    fprintf('%-16s  %.3f  %.3f   %.3f   %.3f    %+.1f%%  %+.1f%%  %+.1f%%\n', ...
        roi_names{k}, f0_from_C1(k), f_obs_C1(k), f_obs_C2(k), ...
        f_pred_C2(k), obs_diff(k), pred_diff(k), unexplained(k));
end
fprintf('%-16s  %.3f                           %+.1f%%  %+.1f%%  %+.1f%%\n', ...
    'Mean', mean(f0_from_C1), mean(obs_diff), mean(pred_diff), mean(unexplained));

%% 5. Figure 1: Observed vs T2-predicted % difference (primary)
figure('unit','inch','position',[0 0 12 5]);
bar_data_diff = [obs_diff, pred_diff];
b3 = bar(categorical(roi_names, roi_names), bar_data_diff);
b3(1).FaceColor = [0.3 0.3 0.3];
b3(2).FaceColor = [0.8 0.5 0.2];
ylabel('(f_{C2} - f_{C1}) / f_{C1}  (%)');
title('C2 vs C1 restricted fraction difference: observed vs T2-predicted (f_{CSF} = 0)');
legend({'Observed', 'T2-predicted'}, 'Location', 'southwest');
yline(0, 'k--', 'HandleVisibility', 'off');
grid on; box on;
set(gca, 'XTickLabelRotation', 45);

exportgraphics(gcf, 'fig_T2_fr_diff_obs_vs_pred.pdf', 'ContentType', 'vector');

%% 6. Figure 2: Scatter — observed diff vs predicted diff (primary)
figure('unit','inch','position',[0 0 6 5]);
scatter(pred_diff, obs_diff, 80, 'filled', 'MarkerFaceColor', [0.2 0.5 0.7]);
hold on;

lims = [min([pred_diff; obs_diff])-2, max([pred_diff; obs_diff])+2];
plot(lims, lims, 'k--', 'LineWidth', 1.5);

for k = 1:Nroi
    text(pred_diff(k)+0.3, obs_diff(k)+0.3, roi_names{k}, 'FontSize', 7);
end

xlabel('T2-predicted difference (%)');
ylabel('Observed difference (%)');
title('T2 weighting explains partial C2–C1 gap (f_{CSF} = 0)');
pbaspect([1 1 1]); grid on; box on;
xlim(lims); ylim(lims);

[r, p] = corr(pred_diff, obs_diff);
text(lims(1)+1, lims(2)-2, sprintf('r = %.2f, p = %.3f', r, p), 'FontSize', 10);

exportgraphics(gcf, 'fig_T2_fr_scatter_pred_vs_obs.pdf', 'ContentType', 'vector');

%% 7. Figure 3: f(TE) vs f0 sweep per ROI (3x4 grid, 2-compartment)
fr_sweep = linspace(0.3, 0.6, 50);

figure('unit','inch','position',[0 0 16 12]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:Nroi
    nexttile;
    T2a = T2a_all(k);
    T2e = T2e_all(k);

    for it = 1:2
        if it == 1; TE = TE_C1; col = [0.2 0.4 0.8]; else; TE = TE_C2; col = [0.8 0.2 0.2]; end
        fr_app = zeros(size(fr_sweep));
        for i = 1:numel(fr_sweep)
            f = fr_sweep(i);
            fr_app(i) = f*exp(-TE/T2a) / (f*exp(-TE/T2a) + (1-f)*exp(-TE/T2e));
        end
        plot(fr_sweep, fr_app, '-', 'Color', col, 'LineWidth', 1.5); hold on;
    end

    plot([0.3 0.6], [0.3 0.6], 'k--', 'LineWidth', 0.8);

    xlim([0.3 0.6]); ylim([0.3 0.85]);
    pbaspect([1 1 1]); grid on; box on;
    title(sprintf('%s (T2a=%d, T2e=%d)', roi_names{k}, T2a, T2e), 'FontSize', 8);

    if k == 1
        legend({'C1 (TE=77)', 'C2 (TE=54)', 'Identity'}, ...
            'FontSize', 6, 'Location', 'northwest');
    end
    if mod(k-1,4) == 0; ylabel('T2-weighted f_r'); end
    if k > 8; xlabel('f_0'); end
end

exportgraphics(gcf, 'fig_T2_fr_sweep_per_ROI.pdf', 'ContentType', 'vector');

%% ===================================================================
%%  SENSITIVITY ANALYSIS: effect of fCSF
%% ===================================================================

fcsf_values = [0.02, 0.05, 0.10];
colors = {[0.2 0.6 0.3], [0.8 0.5 0.2], [0.6 0.2 0.6]};

%% 8. Figure 4: 3-panel bar — sensitivity to fCSF
figure('unit','inch','position',[0 0 18 5]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for ip = 1:3
    fcsf = fcsf_values(ip);
    pd = compute_pred_diff(fcsf, f_obs_C1, T2a_all, T2e_all, TE_C1, TE_C2, T2csf, Nroi);

    nexttile;
    b = bar(categorical(roi_names, roi_names), [obs_diff, pd]);
    b(1).FaceColor = [0.3 0.3 0.3];
    b(2).FaceColor = colors{ip};
    ylabel('(f_{C2} - f_{C1}) / f_{C1}  (%)');
    title(sprintf('f_{CSF} = %.2f', fcsf));
    legend({'Observed', 'T2-predicted'}, 'Location', 'southwest');
    yline(0, 'k--', 'HandleVisibility', 'off');
    grid on; box on;
    set(gca, 'XTickLabelRotation', 45);
    ylim([min(obs_diff)-5, 15]);
end

exportgraphics(gcf, 'fig_T2_fr_sensitivity_fcsf.pdf', 'ContentType', 'vector');

%% 9. Figure 5: Scatter overlay — all fCSF values
figure('unit','inch','position',[0 0 7 6]);
hold on;

% fCSF = 0 (primary)
scatter(pred_diff, obs_diff, 70, 'filled', ...
    'MarkerFaceColor', [0.2 0.5 0.7], 'MarkerFaceAlpha', 0.8);

for ip = 1:3
    fcsf = fcsf_values(ip);
    pd = compute_pred_diff(fcsf, f_obs_C1, T2a_all, T2e_all, TE_C1, TE_C2, T2csf, Nroi);
    scatter(pd, obs_diff, 50, 'filled', ...
        'MarkerFaceColor', colors{ip}, 'MarkerFaceAlpha', 0.6);
end

lims2 = [min(obs_diff)-3, max(obs_diff)+3];
plot(lims2, lims2, 'k--', 'LineWidth', 1.5);
xlabel('T2-predicted difference (%)');
ylabel('Observed difference (%)');
title('Sensitivity of T2-predicted diff to f_{CSF}');
legend([{'f_{CSF}=0'}, arrayfun(@(x) sprintf('f_{CSF}=%.2f',x), fcsf_values, 'UniformOutput', false)], ...
    'Location', 'northwest');
pbaspect([1 1 1]); grid on; box on;
xlim(lims2); ylim(lims2);

exportgraphics(gcf, 'fig_T2_fr_scatter_sensitivity.pdf', 'ContentType', 'vector');

%% ---- Local function ----
function [pred_diff, f0_all, f_pred_C2] = compute_pred_diff(fcsf, f_obs_C1, T2a, T2e, TE_C1, TE_C2, T2csf, Nroi)
    f0_all    = zeros(Nroi, 1);
    f_pred_C1 = zeros(Nroi, 1);
    f_pred_C2 = zeros(Nroi, 1);

    for k = 1:Nroi
        ea1 = exp(-TE_C1 / T2a(k));
        ee1 = exp(-TE_C1 / T2e(k));
        ec1 = exp(-TE_C1 / T2csf);

        fr = f_obs_C1(k);
        f0 = fr * (ee1 - fcsf*(ee1 - ec1)) / (ea1*(1-fr) + fr*ee1);
        f0_all(k) = f0;
        fe = 1 - f0 - fcsf;

        f_pred_C1(k) = f0*ea1 / (f0*ea1 + fe*ee1 + fcsf*ec1);

        ea2 = exp(-TE_C2 / T2a(k));
        ee2 = exp(-TE_C2 / T2e(k));
        ec2 = exp(-TE_C2 / T2csf);
        f_pred_C2(k) = f0*ea2 / (f0*ea2 + fe*ee2 + fcsf*ec2);
    end

    pred_diff = (f_pred_C2 - f_pred_C1) ./ f_pred_C1 * 100;
end
