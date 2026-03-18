clear; close all;
addpath(genpath('/autofs/space/linen_001/users/Yixin/C2_protocoldesign-main/lib'));

%% 1. Parameters
Nsample = 1000;
rng(42);

% T2 values [ms]
T2r   = 90;     % restricted (intra-axonal)
T2e   = 50;     % extra-axonal
T2csf = 2000;   % CSF

% Tissue parameters
D0   = 2;
Da   = 1.7;
DeL  = 1.7;
Dcsf = 3;
model = 'VanGelderen';

% Ground truth
r_true    = rand(Nsample,1) * (5 - 0.1) + 0.1;
f_true    = rand(Nsample,1) * (0.6 - 0.3) + 0.3;
fcsf_true = 0.1 * ones(Nsample,1);
DeR_true  = rand(Nsample,1) * (1.5 - 0.5) + 0.5;

%% 2. Protocol definitions
% C2
b_C2 = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
        0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta_C2 = 8 * ones(16,1);
Delta_C2 = [19*ones(8,1); 49*ones(8,1)];
SNR_C2 = 38;
TE_C2  = 54;  % ms

% C1 — same timing, b-values scaled by (Gmax_C1/Gmax_C2)^2
Gmax_C1 = 290; Gmax_C2_val = 495;
b_C1 = b_C2 * (Gmax_C1/Gmax_C2_val)^2;
delta_C1 = delta_C2;
Delta_C1 = Delta_C2;
SNR_C1 = 18;
TE_C1  = 77;  % ms

protocols = struct();
protocols(1).name = 'C2'; protocols(1).b = b_C2; protocols(1).delta = delta_C2;
protocols(1).Delta = Delta_C2; protocols(1).SNR = SNR_C2; protocols(1).TE = TE_C2;
protocols(2).name = 'C1'; protocols(2).b = b_C1; protocols(2).delta = delta_C1;
protocols(2).Delta = Delta_C1; protocols(2).SNR = SNR_C1; protocols(2).TE = TE_C1;

%% 3. Simulate and fit
for ip = 1:2
    fprintf('\n=== Protocol %s (TE=%dms, SNR=%d) ===\n', ...
        protocols(ip).name, protocols(ip).TE, protocols(ip).SNR);

    smt = AxCaliberSMT(protocols(ip).b, protocols(ip).delta, protocols(ip).Delta, ...
                        D0, Da, DeL, Dcsf);
    TE  = protocols(ip).TE;
    snr = protocols(ip).SNR;

    r_fit = zeros(Nsample,1); f_fit = zeros(Nsample,1);
    fcsf_fit = zeros(Nsample,1); DeR_fit = zeros(Nsample,1);
    f_app = zeros(Nsample,1); fcsf_app = zeros(Nsample,1);

    tic;
    for i = 1:Nsample
        if mod(i,100)==0; fprintf('  %d/%d\n', i, Nsample); end

        % T2-weighted apparent fractions
        w_r   = f_true(i) * exp(-TE/T2r);
        w_e   = (1 - f_true(i) - fcsf_true(i)) * exp(-TE/T2e);
        w_csf = fcsf_true(i) * exp(-TE/T2csf);
        S0 = w_r + w_e + w_csf;
        f_app(i)    = w_r / S0;
        fcsf_app(i) = w_csf / S0;

        % Generate signal with T2-weighted fractions
        S = smt.AxonDiameterFWD([r_true(i), f_app(i), fcsf_app(i), DeR_true(i)], model);

        % Add Gaussian noise
        S_noisy = double(S + randn(size(S))/snr);

        % Fit (5 random restarts + kmeans)
        [r_fit(i), f_fit(i), fcsf_fit(i), DeR_fit(i)] = smt.mcmc(S_noisy, model, 5);
    end
    fprintf('  Done in %.1f s\n', toc);

    protocols(ip).r_fit = r_fit; protocols(ip).f_fit = f_fit;
    protocols(ip).fcsf_fit = fcsf_fit; protocols(ip).DeR_fit = DeR_fit;
    protocols(ip).f_app = f_app; protocols(ip).fcsf_app = fcsf_app;
end

save('C1_C2_comparison_results.mat', 'protocols', 'r_true', 'f_true', ...
     'fcsf_true', 'DeR_true');

%% 4. Figure 1: 2x4 scatter plots with density coloring
figure('unit','inch','position',[0 0 16 8]);
tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

param_labels = {'Axon diameter (\mum)', 'f', 'f_{csf}', 'D_e (\mum^2/ms)'};

for ip = 1:2
    % GT: T2-weighted for f/fcsf; true for r/DeR
    gt  = {2*r_true, protocols(ip).f_app, protocols(ip).fcsf_app, DeR_true};
    est = {2*protocols(ip).r_fit, protocols(ip).f_fit, ...
           protocols(ip).fcsf_fit, protocols(ip).DeR_fit};

    for jp = 1:4
        nexttile;
        x = gt{jp}; y = est{jp};

        % Density coloring via 2D histogram
        nbins = 30;
        [N, xe, ye] = histcounts2(x, y, nbins);
        bx = discretize(x, xe); by = discretize(y, ye);
        bx = min(bx, size(N,1)); by = min(by, size(N,2));
        valid = ~isnan(bx) & ~isnan(by);
        c = zeros(size(x));
        c(valid) = arrayfun(@(k) N(bx(k),by(k)), find(valid));
        [~, ord] = sort(c);  % plot low density first

        scatter(x(ord), y(ord), 15, c(ord), 'filled', 'MarkerFaceAlpha', 0.7);
        hold on;
        mn = min([x;y]); mx = max([x;y]);
        plot([mn mx], [mn mx], 'k-', 'LineWidth', 1.5);
        xlim([mn mx]); ylim([mn mx]);
        pbaspect([1 1 1]); box on; grid on;
        colormap(gca, 'parula');

        if ip == 1; title(param_labels{jp}); end
        xlabel('Ground truth');
        if jp == 1
            ylabel(sprintf('Protocol %s — Estimated', protocols(ip).name));
        else
            ylabel('Estimated');
        end
    end
end
exportgraphics(gcf, 'fig_C1_C2_scatter.pdf', 'ContentType', 'vector');

%% 5. Figure 2: T2 decay weighted fraction relationship
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
