clear; close all;

%% T2 decay weighted vs original restricted fraction (data-driven)
% Reference: Veraart, Novikov, Fieremans (2018) — TEdDI, NeuroImage
% 2-compartment model (no CSF):
%   f(TE) = f0*exp(-TE/T2a) / [f0*exp(-TE/T2a) + (1-f0)*exp(-TE/T2e)]
% Inverse (solve for f0 from observed f):
%   f0 = f*exp(-TE/T2e) / [f*exp(-TE/T2e) + (1-f)*exp(-TE/T2a)]

% Protocol TEs
TE_C2 = 54;   % ms
TE_C1 = 77;   % ms

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

%% 2. Back-calculate f0 from observed C1 f_r (longer TE → more T2 weighting)
% f0 = f*exp(-TE/T2e) / [f*exp(-TE/T2e) + (1-f)*exp(-TE/T2a)]
f0_from_C1 = zeros(Nroi, 1);
for k = 1:Nroi
    f  = f_obs_C1(k);
    ea = exp(-TE_C1 / T2a_all(k));
    ee = exp(-TE_C1 / T2e_all(k));
    f0_from_C1(k) = f * ee / (f * ee + (1-f) * ea);
end

%% 3. Forward predict f_r at both TEs using back-calculated f0
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

% Differences (C2 - C1, negative = C2 lower)
obs_diff      = (f_obs_C2  - f_obs_C1)  ./ f_obs_C1  * 100;
pred_diff     = (f_pred_C2 - f_pred_C1) ./ f_pred_C1 * 100;
unexplained   = obs_diff - pred_diff;  % residual not explained by T2

%% 4. Print summary table
fprintf('%-16s  f0     f_C1_obs  f_C1_pred  f_C2_obs  f_C2_pred  obs%%    pred%%   residual%%\n', 'Tract');
fprintf('%-16s  -----  --------  ---------  --------  ---------  ------  ------  ---------\n', '');
for k = 1:Nroi
    fprintf('%-16s  %.3f  %.3f     %.3f      %.3f     %.3f      %+.1f%%  %+.1f%%  %+.1f%%\n', ...
        roi_names{k}, f0_from_C1(k), f_obs_C1(k), f_pred_C1(k), ...
        f_obs_C2(k), f_pred_C2(k), obs_diff(k), pred_diff(k), unexplained(k));
end
fprintf('%-16s  %.3f                                             %+.1f%%  %+.1f%%  %+.1f%%\n', ...
    'Mean', mean(f0_from_C1), mean(obs_diff), mean(pred_diff), mean(unexplained));

%% 5. Figure 1: Observed vs predicted f_r for C1 and C2
figure('unit','inch','position',[0 0 14 5]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% C1
nexttile;
bar_data_C1 = [f_obs_C1, f_pred_C1];
b1 = bar(categorical(roi_names, roi_names), bar_data_C1);
b1(1).FaceColor = [0.2 0.4 0.8];
b1(2).FaceColor = [0.6 0.7 0.9];
ylabel('Restricted volume fraction');
title('C1 (TE = 77 ms)');
legend({'Observed', 'Predicted (T2 model)'}, 'Location', 'northwest');
ylim([0 0.85]); grid on; box on;
set(gca, 'XTickLabelRotation', 45);

% C2
nexttile;
bar_data_C2 = [f_obs_C2, f_pred_C2];
b2 = bar(categorical(roi_names, roi_names), bar_data_C2);
b2(1).FaceColor = [0.8 0.2 0.2];
b2(2).FaceColor = [0.9 0.6 0.6];
ylabel('Restricted volume fraction');
title('C2 (TE = 54 ms)');
legend({'Observed', 'Predicted (T2 model)'}, 'Location', 'northwest');
ylim([0 0.85]); grid on; box on;
set(gca, 'XTickLabelRotation', 45);

exportgraphics(gcf, 'fig_T2_fr_obs_vs_pred.pdf', 'ContentType', 'vector');

%% 6. Figure 2: Observed vs T2-predicted % difference (C2 - C1)
figure('unit','inch','position',[0 0 12 5]);
bar_data_diff = [obs_diff, pred_diff];
b3 = bar(categorical(roi_names, roi_names), bar_data_diff);
b3(1).FaceColor = [0.3 0.3 0.3];
b3(2).FaceColor = [0.8 0.5 0.2];
ylabel('(f_{C2} - f_{C1}) / f_{C1}  (%)');
title('C2 vs C1 restricted fraction difference: observed vs T2-predicted');
legend({'Observed', 'T2-predicted'}, 'Location', 'southwest');
yline(0, 'k--', 'HandleVisibility', 'off');
grid on; box on;
set(gca, 'XTickLabelRotation', 45);

exportgraphics(gcf, 'fig_T2_fr_diff_obs_vs_pred.pdf', 'ContentType', 'vector');

%% 7. Figure 3: Scatter — observed diff vs predicted diff
figure('unit','inch','position',[0 0 6 5]);
scatter(pred_diff, obs_diff, 80, 'filled', 'MarkerFaceColor', [0.2 0.5 0.7]);
hold on;

% Identity line
lims = [min([pred_diff; obs_diff])-2, max([pred_diff; obs_diff])+2];
plot(lims, lims, 'k--', 'LineWidth', 1.5);

% Label each point
for k = 1:Nroi
    text(pred_diff(k)+0.3, obs_diff(k)+0.3, roi_names{k}, 'FontSize', 7);
end

xlabel('T2-predicted difference (%)');
ylabel('Observed difference (%)');
title('T2 weighting explains partial C2–C1 gap');
pbaspect([1 1 1]); grid on; box on;
xlim(lims); ylim(lims);

% Correlation
[r, p] = corr(pred_diff, obs_diff);
text(lims(1)+1, lims(2)-2, sprintf('r = %.2f, p = %.3f', r, p), 'FontSize', 10);

exportgraphics(gcf, 'fig_T2_fr_scatter_pred_vs_obs.pdf', 'ContentType', 'vector');

%% 8. Figure 4: f(TE) vs original f sweep for each ROI (3x4 grid)
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

    % Identity line
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
