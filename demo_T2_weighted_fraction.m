clear; close all;

%% T2 decay weighted vs original restricted fraction (data-driven, 3-compartment)
% Reference: Veraart, Novikov, Fieremans (2018) — TEdDI, NeuroImage
% 3-compartment model:
%   f_r(TE) = f0*exp(-TE/T2a) / [f0*exp(-TE/T2a) + fe*exp(-TE/T2e) + fcsf*exp(-TE/T2csf)]
%   where fe = 1 - f0 - fcsf

% Protocol TEs
TE_C2 = 54;   % ms
TE_C1 = 77;   % ms

% CSF T2
T2csf = 2000;  % ms (at 3T)

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

% Observed % difference
obs_diff = (f_obs_C2 - f_obs_C1) ./ f_obs_C1 * 100;

%% 2. Sweep fCSF to find optimal value (minimize MSE of predicted vs observed diff)
fcsf_sweep = 0:0.005:0.25;
mse_sweep  = zeros(size(fcsf_sweep));

for ifc = 1:numel(fcsf_sweep)
    fcsf = fcsf_sweep(ifc);
    pred = zeros(Nroi, 1);

    for k = 1:Nroi
        ea1 = exp(-TE_C1 / T2a_all(k));
        ee1 = exp(-TE_C1 / T2e_all(k));
        ec1 = exp(-TE_C1 / T2csf);

        % Back-calculate f0 from observed C1 f_r
        fr = f_obs_C1(k);
        f0 = fr * (ee1 - fcsf*(ee1 - ec1)) / (ea1*(1 - fr) + fr*ee1);

        % Check validity
        if f0 <= 0 || f0 + fcsf >= 1
            pred(k) = NaN;
            continue;
        end

        % Forward predict C2
        ea2 = exp(-TE_C2 / T2a_all(k));
        ee2 = exp(-TE_C2 / T2e_all(k));
        ec2 = exp(-TE_C2 / T2csf);
        fe  = 1 - f0 - fcsf;
        f_C1 = f0*ea1 / (f0*ea1 + fe*ee1 + fcsf*ec1);
        f_C2 = f0*ea2 / (f0*ea2 + fe*ee2 + fcsf*ec2);
        pred(k) = (f_C2 - f_C1) / f_C1 * 100;
    end

    mse_sweep(ifc) = nanmean((pred - obs_diff).^2);
end

[~, best_idx] = min(mse_sweep);
fcsf_opt = fcsf_sweep(best_idx);
fprintf('Optimal fCSF = %.3f (MSE = %.2f)\n', fcsf_opt, mse_sweep(best_idx));

%% 3. Compute final results with optimal fCSF
fcsf = fcsf_opt;
f0_from_C1 = zeros(Nroi, 1);
f_pred_C1  = zeros(Nroi, 1);
f_pred_C2  = zeros(Nroi, 1);

for k = 1:Nroi
    ea1 = exp(-TE_C1 / T2a_all(k));
    ee1 = exp(-TE_C1 / T2e_all(k));
    ec1 = exp(-TE_C1 / T2csf);

    fr = f_obs_C1(k);
    f0 = fr * (ee1 - fcsf*(ee1 - ec1)) / (ea1*(1 - fr) + fr*ee1);
    f0_from_C1(k) = f0;
    fe = 1 - f0 - fcsf;

    f_pred_C1(k) = f0*ea1 / (f0*ea1 + fe*ee1 + fcsf*ec1);

    ea2 = exp(-TE_C2 / T2a_all(k));
    ee2 = exp(-TE_C2 / T2e_all(k));
    ec2 = exp(-TE_C2 / T2csf);
    f_pred_C2(k) = f0*ea2 / (f0*ea2 + fe*ee2 + fcsf*ec2);
end

pred_diff   = (f_pred_C2 - f_pred_C1) ./ f_pred_C1 * 100;
unexplained = obs_diff - pred_diff;

%% 4. Print summary table
fprintf('\n=== Results with fCSF = %.3f, T2csf = %d ms ===\n', fcsf, T2csf);
fprintf('%-16s  f0     f_C1_obs  f_C2_obs  f_C2_pred  obs%%    pred%%   residual%%\n', 'Tract');
fprintf('%-16s  -----  --------  --------  ---------  ------  ------  ---------\n', '');
for k = 1:Nroi
    fprintf('%-16s  %.3f  %.3f     %.3f     %.3f      %+.1f%%  %+.1f%%  %+.1f%%\n', ...
        roi_names{k}, f0_from_C1(k), f_obs_C1(k), ...
        f_obs_C2(k), f_pred_C2(k), obs_diff(k), pred_diff(k), unexplained(k));
end
fprintf('%-16s  %.3f                                   %+.1f%%  %+.1f%%  %+.1f%%\n', ...
    'Mean', mean(f0_from_C1), mean(obs_diff), mean(pred_diff), mean(unexplained));

%% 5. Figure 0: MSE vs fCSF (justification for chosen fCSF)
figure('unit','inch','position',[0 0 6 4]);
plot(fcsf_sweep, mse_sweep, 'b-', 'LineWidth', 1.5); hold on;
plot(fcsf_opt, mse_sweep(best_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('f_{CSF}'); ylabel('MSE (predicted vs observed diff)');
title(sprintf('Optimal f_{CSF} = %.3f', fcsf_opt));
grid on; box on;
exportgraphics(gcf, 'fig_T2_fcsf_optimization.pdf', 'ContentType', 'vector');

%% 6. Figure 1: Observed vs predicted f_r for C1 and C2
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

%% 7. Figure 2: Observed vs T2-predicted % difference (C2 - C1)
figure('unit','inch','position',[0 0 12 5]);
bar_data_diff = [obs_diff, pred_diff];
b3 = bar(categorical(roi_names, roi_names), bar_data_diff);
b3(1).FaceColor = [0.3 0.3 0.3];
b3(2).FaceColor = [0.8 0.5 0.2];
ylabel('(f_{C2} - f_{C1}) / f_{C1}  (%)');
title(sprintf('C2 vs C1 restricted fraction difference (f_{CSF} = %.3f)', fcsf));
legend({'Observed', 'T2-predicted'}, 'Location', 'southwest');
yline(0, 'k--', 'HandleVisibility', 'off');
grid on; box on;
set(gca, 'XTickLabelRotation', 45);

exportgraphics(gcf, 'fig_T2_fr_diff_obs_vs_pred.pdf', 'ContentType', 'vector');

%% 8. Figure 3: Scatter — observed diff vs predicted diff
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
title(sprintf('T2 weighting explains partial C2–C1 gap (f_{CSF} = %.3f)', fcsf));
pbaspect([1 1 1]); grid on; box on;
xlim(lims); ylim(lims);

% Correlation
[r, p] = corr(pred_diff, obs_diff);
text(lims(1)+1, lims(2)-2, sprintf('r = %.2f, p = %.3f', r, p), 'FontSize', 10);

exportgraphics(gcf, 'fig_T2_fr_scatter_pred_vs_obs.pdf', 'ContentType', 'vector');

%% 9. Figure 4: f(TE) vs original f sweep for each ROI (3x4 grid)
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
            fe = 1 - f - fcsf;
            ea = exp(-TE/T2a); ee = exp(-TE/T2e); ec = exp(-TE/T2csf);
            fr_app(i) = f*ea / (f*ea + fe*ee + fcsf*ec);
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
