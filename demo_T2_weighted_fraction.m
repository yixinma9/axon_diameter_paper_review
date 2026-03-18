clear; close all;

%% T2 decay weighted vs original restricted fraction
% Reference: Veraart, Novikov, Fieremans (2017) — TEdDI, Neuroimage
% f(TE) = f0*exp(-TE/T2a) / [f0*exp(-TE/T2a) + (1-f0-fcsf)*exp(-TE/T2e) + fcsf*exp(-TE/T2csf)]

% Protocol TEs
TE_C2 = 54;   % ms
TE_C1 = 77;   % ms
T2csf = 2000; % ms

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

% Assumed ground truth fractions (same for all ROIs for simplicity)
f0   = 0.45;
fcsf = 0.1;

%% 1. Compute T2-weighted f for each ROI at both TEs
f_C2 = zeros(Nroi, 1);
f_C1 = zeros(Nroi, 1);

for k = 1:Nroi
    T2a = T2a_all(k);
    T2e = T2e_all(k);
    fe  = 1 - f0 - fcsf;

    % C2
    S0_C2 = f0*exp(-TE_C2/T2a) + fe*exp(-TE_C2/T2e) + fcsf*exp(-TE_C2/T2csf);
    f_C2(k) = f0*exp(-TE_C2/T2a) / S0_C2;

    % C1
    S0_C1 = f0*exp(-TE_C1/T2a) + fe*exp(-TE_C1/T2e) + fcsf*exp(-TE_C1/T2csf);
    f_C1(k) = f0*exp(-TE_C1/T2a) / S0_C1;
end

% Percentage difference: (C1 - C2) / C2 * 100
pct_diff = (f_C1 - f_C2) ./ f_C2 * 100;

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
ylabel('(f_{C1} - f_{C2}) / f_{C2}  (%)');
title('Predicted % difference in f_r due to T2 weighting (C1 vs C2)');
yline(0, 'k--', 'HandleVisibility', 'off');
grid on; box on;
exportgraphics(gcf, 'fig_T2_fr_pct_diff.pdf', 'ContentType', 'vector');

%% 5. Figure 3: f(TE) vs original f sweep (for mean T2 values)
fr_sweep = linspace(0.3, 0.6, 50);
T2a_mean = mean(T2a_all);
T2e_mean = mean(T2e_all);

figure('unit','inch','position',[0 0 6 5]);
TEs = [TE_C2, TE_C1];
markers = {'-o', '-s'};
labels = {sprintf('f_r (C2: TE=%dms)', TE_C2), sprintf('f_r (C1: TE=%dms)', TE_C1)};

for it = 1:2
    TE = TEs(it);
    fr_app = zeros(size(fr_sweep));
    for i = 1:numel(fr_sweep)
        fr = fr_sweep(i);
        fe = 1 - fr - fcsf;
        S0 = fr*exp(-TE/T2a_mean) + fe*exp(-TE/T2e_mean) + fcsf*exp(-TE/T2csf);
        fr_app(i) = fr*exp(-TE/T2a_mean) / S0;
    end
    plot(fr_sweep, fr_app, markers{it}, 'MarkerSize', 5, 'LineWidth', 1.5); hold on;
end
xlabel('Original Restricted Fraction');
ylabel('T2 decay weighted Restricted Fraction');
title({'T2 decay weighted and original restricted fraction', ...
       sprintf('mean T2a=%.0fms, T2e=%.0fms', T2a_mean, T2e_mean)});
legend(labels, 'Location', 'northwest');
pbaspect([1 1 1]); box on; grid on;
exportgraphics(gcf, 'fig_T2_weighted_fraction.pdf', 'ContentType', 'vector');
