clear; close all;

%% T2 decay weighted vs original restricted fraction
% Reference: Veraart, Novikov, Fieremans (2017) — TEdDI, Neuroimage
% f(TE) = f0*exp(-TE/T2a) / [f0*exp(-TE/T2a) + (1-f0-fcsf)*exp(-TE/T2e) + fcsf*exp(-TE/T2csf)]

% T2 values [ms]
T2r   = 90;     % restricted (intra-axonal)
T2e   = 50;     % extra-axonal
T2csf = 2000;   % CSF

% Protocol TEs
TE_C2 = 54;   % ms
TE_C1 = 77;   % ms

% Sweep original restricted fraction
fr_sweep = linspace(0.3, 0.6, 50);
fcsf_val = 0.1;

figure('unit','inch','position',[0 0 6 5]);
TEs = [TE_C2, TE_C1];
markers = {'-o', '-s'};
labels = {sprintf('f_r (C2: TE=%dms)', TE_C2), sprintf('f_r (C1: TE=%dms)', TE_C1)};

for it = 1:2
    TE = TEs(it);
    fr_app = zeros(size(fr_sweep));
    for i = 1:numel(fr_sweep)
        fr = fr_sweep(i);
        fe = 1 - fr - fcsf_val;
        S0 = fr*exp(-TE/T2r) + fe*exp(-TE/T2e) + fcsf_val*exp(-TE/T2csf);
        fr_app(i) = fr*exp(-TE/T2r) / S0;
    end
    plot(fr_sweep, fr_app, markers{it}, 'MarkerSize', 5, 'LineWidth', 1.5); hold on;
end
xlabel('Original Restricted Fraction');
ylabel('T2 decay weighted Restricted Fraction');
title({'T2 decay weighted and original restricted fraction', 'at two TEs'});
legend(labels, 'Location', 'northwest');
pbaspect([1 1 1]); box on; grid on;
exportgraphics(gcf, 'fig_T2_weighted_fraction.pdf', 'ContentType', 'vector');
