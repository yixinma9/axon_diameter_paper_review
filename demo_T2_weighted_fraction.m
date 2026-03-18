clear; close all;

%% T2 decay weighted vs original restricted fraction
% Reference: Veraart, Novikov, Fieremans (2017) — TEdDI, Neuroimage
% 2-compartment model (no CSF):
% f(TE) = f0*exp(-TE/T2a) / [f0*exp(-TE/T2a) + (1-f0)*exp(-TE/T2e)]

% Protocol TEs
TE_C2 = 54;   % ms
TE_C1 = 77;   % ms

% Per-ROI T2 values from TEdDI paper Figure 5 (Veraart et al. 2017)
%           ROI name    T2a(ms)  T2e(ms)
roi_data = {
    'PLIC',             110,     42;
    'ALIC',              85,     48;
    'SLF',               90,     38;
    'EC',                70,     55;
    'GCC',               75,     35;
    'BCC',               90,     45;
    'SCC',               90,     40;
    'ACR',               75,     48;
    'SCR',              100,     48;
    'PCR',               95,     48;
};
roi_names = roi_data(:,1);
T2a_all   = cell2mat(roi_data(:,2));
T2e_all   = cell2mat(roi_data(:,3));
Nroi = numel(roi_names);

% Assumed ground truth fraction
f0 = 0.45;

%% 1. Compute T2-weighted f for each ROI at both TEs (2-compartment, no CSF)
f_C2 = zeros(Nroi, 1);
f_C1 = zeros(Nroi, 1);

for k = 1:Nroi
    T2a = T2a_all(k);
    T2e = T2e_all(k);
    fe  = 1 - f0;

    % C2
    S0_C2 = f0*exp(-TE_C2/T2a) + fe*exp(-TE_C2/T2e);
    f_C2(k) = f0*exp(-TE_C2/T2a) / S0_C2;

    % C1
    S0_C1 = f0*exp(-TE_C1/T2a) + fe*exp(-TE_C1/T2e);
    f_C1(k) = f0*exp(-TE_C1/T2a) / S0_C1;
end

% Percentage difference: (C2 - C1) / C1 * 100  (negative = C2 lower)
pct_diff = (f_C2 - f_C1) ./ f_C1 * 100;

%% 2. Print table
fprintf('%-6s  T2a  T2e   f_C2   f_C1   diff(%%)  \n', 'ROI');
fprintf('------  ---  ---  -----  -----  --------\n');
for k = 1:Nroi
    fprintf('%-6s  %3d  %3d  %.3f  %.3f  %+.1f%%\n', ...
        roi_names{k}, T2a_all(k), T2e_all(k), f_C2(k), f_C1(k), pct_diff(k));
end
fprintf('------  ---  ---  -----  -----  --------\n');
fprintf('%-6s            %.3f  %.3f  %+.1f%%\n', 'Mean', mean(f_C2), mean(f_C1), mean(pct_diff));

%% 3. Figure 1: Bar plot — predicted f at C1 vs C2 per ROI
figure('unit','inch','position',[0 0 12 5]);
bar_data = [f_C1, f_C2];
b = bar(categorical(roi_names, roi_names), bar_data);
b(1).FaceColor = [0.2 0.4 0.8];  % C1 blue
b(2).FaceColor = [0.8 0.2 0.2];  % C2 red
ylabel('T2-weighted restricted fraction');
title('Predicted T2-weighted f_r per ROI (C1 vs C2)');
legend({'C1 (TE=77ms)', 'C2 (TE=54ms)'}, 'Location', 'northwest');
ylim([0 0.8]);
grid on; box on;
exportgraphics(gcf, 'fig_T2_fr_per_ROI.pdf', 'ContentType', 'vector');

%% 4. Figure 2: Percentage difference per ROI
figure('unit','inch','position',[0 0 12 4]);
bar(categorical(roi_names, roi_names), pct_diff, 'FaceColor', [0.5 0.7 0.3]);
ylabel('(f_{C2} - f_{C1}) / f_{C1}  (%)');
title('Predicted % difference in f_r due to T2 weighting (C2 vs C1)');
yline(0, 'k--', 'HandleVisibility', 'off');
grid on; box on;
exportgraphics(gcf, 'fig_T2_fr_pct_diff.pdf', 'ContentType', 'vector');

%% 5. Figure 3: f(TE) vs original f sweep — per ROI (3x4 grid)
fcsf = 0.05;
T2csf = 2000;
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
        plot(fr_sweep, fr_app, '-', 'Color', col, 'LineWidth', 2); hold on;
    end

    plot([0.3 0.8], [0.3 0.8], 'k--', 'LineWidth', 1);

    xlim([0.3 0.6]); ylim([0.3 0.8]);
    grid on; box on;
    set(gca, 'FontSize', 11);
    title(sprintf('%s (T2a=%d, T2e=%d)', roi_names{k}, T2a, T2e), 'FontSize', 13);

    if k == 1
        legend({'C1 (TE=77)', 'C2 (TE=54)', 'Identity'}, ...
            'FontSize', 10, 'Location', 'northwest');
    end
    if mod(k-1,4) == 0; ylabel('T2-weighted f_r', 'FontSize', 12); end
    if k > 8; xlabel('f_0', 'FontSize', 12); end
end

% Use last tiles (3x4=12, only 10 ROIs) for annotation
nexttile; axis off;
text(0.1, 0.5, sprintf('f_0 = %.2f\nf_{CSF} = %.2f\nT2_{CSF} = %d ms', f0, fcsf, T2csf), ...
    'FontSize', 14, 'VerticalAlignment', 'middle');

exportgraphics(gcf, 'fig_T2_weighted_fraction_per_ROI.pdf', 'ContentType', 'vector');

%% 6. Figure 4: f_r vs TE (30–60 ms) per ROI
TE_sweep = linspace(30, 60, 100);

figure('unit','inch','position',[0 0 8 5]);
cmap = lines(Nroi);
hold on;
for k = 1:Nroi
    T2a = T2a_all(k);
    T2e = T2e_all(k);
    fr_te = zeros(size(TE_sweep));
    for i = 1:numel(TE_sweep)
        TE = TE_sweep(i);
        fr_te(i) = f0*exp(-TE/T2a) / (f0*exp(-TE/T2a) + (1-f0)*exp(-TE/T2e));
    end
    plot(TE_sweep, fr_te, '-', 'Color', cmap(k,:), 'LineWidth', 1.5);
end
xline(TE_C2, 'r--', 'C2 (54ms)', 'LineWidth', 1.2, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xlabel('TE (ms)');
ylabel('T2-weighted restricted fraction');
title(sprintf('f_r vs TE per ROI (f_0 = %.2f)', f0));
legend(roi_names, 'Location', 'eastoutside', 'FontSize', 7);
grid on; box on;
exportgraphics(gcf, 'fig_T2_fr_vs_TE_per_ROI.pdf', 'ContentType', 'vector');
