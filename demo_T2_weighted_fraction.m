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

%% 2. Helper: compute predicted diff for a given fCSF
compute_pred = @(fcsf_val) compute_pred_diff(fcsf_val, f_obs_C1, T2a_all, T2e_all, ...
    TE_C1, TE_C2, T2csf, Nroi);

%% 3. Three fCSF values to compare
fcsf_values = [0.02, 0.05, 0.10];
colors = {[0.2 0.6 0.3], [0.8 0.5 0.2], [0.6 0.2 0.6]};

%% 4. Figure 1: Observed vs T2-predicted % difference — 3 panels
figure('unit','inch','position',[0 0 18 5]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for ip = 1:3
    fcsf = fcsf_values(ip);
    [pred_diff, f0_all, f_pred_C2] = compute_pred(fcsf);

    nexttile;
    bar_data = [obs_diff, pred_diff];
    b = bar(categorical(roi_names, roi_names), bar_data);
    b(1).FaceColor = [0.3 0.3 0.3];
    b(2).FaceColor = colors{ip};
    ylabel('(f_{C2} - f_{C1}) / f_{C1}  (%)');
    title(sprintf('f_{CSF} = %.2f', fcsf));
    legend({'Observed', 'T2-predicted'}, 'Location', 'southwest');
    yline(0, 'k--', 'HandleVisibility', 'off');
    grid on; box on;
    set(gca, 'XTickLabelRotation', 45);
    ylim([min(obs_diff)-5, max([obs_diff; pred_diff])+5]);
end

exportgraphics(gcf, 'fig_T2_fr_diff_3panels.pdf', 'ContentType', 'vector');

%% 5. Figure 2: Scatter — observed vs predicted for all 3 fCSF values (overlay)
figure('unit','inch','position',[0 0 7 6]);
hold on;

for ip = 1:3
    fcsf = fcsf_values(ip);
    pred_diff = compute_pred(fcsf);
    scatter(pred_diff, obs_diff, 60, 'filled', ...
        'MarkerFaceColor', colors{ip}, 'MarkerFaceAlpha', 0.7);
end

lims = [min(obs_diff)-3, max(obs_diff)+3];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
xlabel('T2-predicted difference (%)');
ylabel('Observed difference (%)');
title('Predicted vs observed C2–C1 difference');
legend(arrayfun(@(x) sprintf('f_{CSF} = %.2f', x), fcsf_values, 'UniformOutput', false), ...
    'Location', 'northwest');
pbaspect([1 1 1]); grid on; box on;
xlim(lims); ylim(lims);

exportgraphics(gcf, 'fig_T2_fr_scatter_3fcsf.pdf', 'ContentType', 'vector');

%% 6. Figure 3: f(TE) vs f0 sweep per ROI — using fCSF = 0.05 (middle value)
fcsf_plot = 0.05;
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
            fe = 1 - f - fcsf_plot;
            ea = exp(-TE/T2a); ee = exp(-TE/T2e); ec = exp(-TE/T2csf);
            fr_app(i) = f*ea / (f*ea + fe*ee + fcsf_plot*ec);
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

%% 7. Print summary tables for all 3 fCSF values
for ip = 1:3
    fcsf = fcsf_values(ip);
    [pred_diff, f0_all, f_pred_C2] = compute_pred(fcsf);
    unexplained = obs_diff - pred_diff;

    fprintf('\n=== fCSF = %.2f ===\n', fcsf);
    fprintf('%-16s  f0     fC1obs  fC2obs  fC2pred  obs%%    pred%%   resid%%\n', 'Tract');
    fprintf('%-16s  -----  ------  ------  -------  ------  ------  ------\n', '');
    for k = 1:Nroi
        fprintf('%-16s  %.3f  %.3f   %.3f   %.3f    %+.1f%%  %+.1f%%  %+.1f%%\n', ...
            roi_names{k}, f0_all(k), f_obs_C1(k), f_obs_C2(k), ...
            f_pred_C2(k), obs_diff(k), pred_diff(k), unexplained(k));
    end
    fprintf('%-16s  %.3f                           %+.1f%%  %+.1f%%  %+.1f%%\n', ...
        'Mean', mean(f0_all), mean(obs_diff), mean(pred_diff), mean(unexplained));
end

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
