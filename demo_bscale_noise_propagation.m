clear; close all;
addpath(genpath('/autofs/space/linen_001/users/Yixin/C2_protocoldesign-main/lib'));

%% Step 1. Load b_scale map and white matter mask
gnc_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001/gnc';
mask_dir = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-001';

% Load grad_dev (QT's method)
% First run: fslmerge -t grad_dev.nii.gz grad_dev_x.nii.gz grad_dev_y.nii.gz grad_dev_z.nii.gz
g = read_avw(fullfile(gnc_dir, 'grad_dev.nii.gz'));  % [x,y,z,9]

% Compute b_scale per voxel (direction-averaged n^2)
% Following QT: L = reshape(squeeze(g(i,j,k,:)),3,3); v = (I+L)*bvecs; n^2 = sum(v.^2)
% For spherical mean: b_scale = mean(n^2) over directions = trace((I+L)'*(I+L))/3
dims = size(g, 1:3);
b_scale_map = zeros(dims, 'single');
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            L = reshape(squeeze(g(i,j,k,:)), 3, 3);
            M = eye(3) + L;
            b_scale_map(i,j,k) = trace(M' * M) / 3;
        end
    end
end

% Load white matter mask
wm_mask = niftiread(fullfile(mask_dir, 'sub-001_mask_nocsfgm.nii.gz')) > 0;

% Extract b_scale values in white matter
b_scale_wm = b_scale_map(wm_mask);
b_scale_wm = b_scale_wm(b_scale_wm > 0);  % remove zeros
fprintf('WM voxels: %d\n', numel(b_scale_wm));
fprintf('b_scale range: [%.4f, %.4f]\n', min(b_scale_wm), max(b_scale_wm));

%% Step 2. Plot b_scale distribution and split into 3 groups by percentile
prc = prctile(b_scale_wm, [33.3, 66.7]);
fprintf('Percentile boundaries: %.4f, %.4f\n', prc(1), prc(2));

group1_idx = b_scale_wm <= prc(1);
group2_idx = b_scale_wm > prc(1) & b_scale_wm <= prc(2);
group3_idx = b_scale_wm > prc(2);

figure('unit','inch','position',[0 0 10 4]);
histogram(b_scale_wm, 100, 'FaceColor', [0.7 0.7 0.7]); hold on;
xline(prc(1), 'r--', 'LineWidth', 2);
xline(prc(2), 'r--', 'LineWidth', 2);
xlabel('b-scale'); ylabel('count');
title('b-scale distribution in white matter');
legend({'b-scale', sprintf('33rd prc = %.3f', prc(1)), sprintf('67th prc = %.3f', prc(2))});

%% Step 3. Sample 500 b_scale values per group
Nsample = 500;
rng(42);

vals1 = b_scale_wm(group1_idx); vals1 = vals1(randperm(numel(vals1), min(Nsample, numel(vals1))));
vals2 = b_scale_wm(group2_idx); vals2 = vals2(randperm(numel(vals2), min(Nsample, numel(vals2))));
vals3 = b_scale_wm(group3_idx); vals3 = vals3(randperm(numel(vals3), min(Nsample, numel(vals3))));

b_scale_samples = {vals1, vals2, vals3};
group_names = {sprintf('Group 1 (%.3f-%.3f)', min(vals1), max(vals1)), ...
               sprintf('Group 2 (%.3f-%.3f)', min(vals2), max(vals2)), ...
               sprintf('Group 3 (%.3f-%.3f)', min(vals3), max(vals3))};

%% Step 4. Setup AxCaliberSMT protocol (TractCaliber)
b = [0.05 0.35 0.8 1.5 2.4 3.45 4.75 6 ...
     0.2  0.95 2.3 4.25 6.75 9.85 13.5 17.8]';
delta = [8*ones(8,1); 8*ones(8,1)];
Delta = [19*ones(8,1); 49*ones(8,1)];

D0   = 2;
Da   = 1.7;
DeL  = 1.7;
Dcsf = 3;
snr  = 30;
model = 'VanGelderen';

% Random ground truth parameters (shared across all groups)
rng(0);
Xrange = [0.1 5; 0.5 1; 0 0.2; 0.5 1.5];
Xinit = Xrange(:,1) + diff(Xrange,1,2) .* rand(4, Nsample);
r_true    = Xinit(1,:)';
f_true    = Xinit(2,:)';
fcsf_true = Xinit(3,:)';
DeR_true  = Xinit(4,:)';

%% Step 5. Noise propagation for each group
% For each group:
%   - Generate signal with b_eff = b * b_scale (ground truth)
%   - Fit with nominal b (uncorrected)
%   - Fit with b_eff (corrected)

results = struct();

for ig = 1:3
    fprintf('\n=== %s ===\n', group_names{ig});
    bscales = b_scale_samples{ig};
    N = numel(bscales);

    r_fit_uncorr = zeros(N, 1);
    f_fit_uncorr = zeros(N, 1);
    fcsf_fit_uncorr = zeros(N, 1);
    DeR_fit_uncorr = zeros(N, 1);

    r_fit_corr = zeros(N, 1);
    f_fit_corr = zeros(N, 1);
    fcsf_fit_corr = zeros(N, 1);
    DeR_fit_corr = zeros(N, 1);

    % Nominal AxCaliberSMT object (for uncorrected fitting)
    smt_nominal = AxCaliberSMT(b, delta, Delta, D0, Da, DeL, Dcsf);

    tic;
    for i = 1:N
        if mod(i, 50) == 0; fprintf('  voxel %d/%d\n', i, N); end

        % Effective b-value for this voxel
        b_eff = b * bscales(i);

        % Generate ground truth signal using effective b-values
        smt_true = AxCaliberSMT(b_eff, delta, Delta, D0, Da, DeL, Dcsf);
        S_clean = smt_true.AxonDiameterFWD([r_true(i), f_true(i), fcsf_true(i), DeR_true(i)], model);

        % Add Rician noise
        S_noisy = abs(S_clean + 1/snr*randn(size(S_clean)) + 1j/snr*randn(size(S_clean)));

        % Fit with nominal b (uncorrected)
        [r_fit_uncorr(i), f_fit_uncorr(i), fcsf_fit_uncorr(i), DeR_fit_uncorr(i)] = ...
            smt_nominal.mcmc(S_noisy, model, 5);

        % Fit with corrected b
        smt_corr = AxCaliberSMT(b_eff, delta, Delta, D0, Da, DeL, Dcsf);
        [r_fit_corr(i), f_fit_corr(i), fcsf_fit_corr(i), DeR_fit_corr(i)] = ...
            smt_corr.mcmc(S_noisy, model, 5);
    end
    t = toc;
    fprintf('  Elapsed: %.1f s\n', t);

    results(ig).group_name = group_names{ig};
    results(ig).bscales = bscales;
    results(ig).r_true = r_true(1:N);
    results(ig).r_fit_uncorr = r_fit_uncorr;
    results(ig).r_fit_corr = r_fit_corr;
    results(ig).f_fit_uncorr = f_fit_uncorr;
    results(ig).f_fit_corr = f_fit_corr;
    results(ig).fcsf_fit_uncorr = fcsf_fit_uncorr;
    results(ig).fcsf_fit_corr = fcsf_fit_corr;
    results(ig).DeR_fit_uncorr = DeR_fit_uncorr;
    results(ig).DeR_fit_corr = DeR_fit_corr;
end

%% Step 6. Plot: actual vs predicted axon radius (3 groups x 2 methods)
figure('unit','inch','position',[0 0 15 8]);
for ig = 1:3
    % Uncorrected
    subplot(2, 3, ig);
    scatter(results(ig).r_true, results(ig).r_fit_uncorr, 10, '.');
    hold on; hr = refline(1,0); set(hr, 'color', 'k');
    xlabel('Ground truth r (\mum)'); ylabel('Fitted r (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s\nUncorrected', results(ig).group_name));

    % Corrected
    subplot(2, 3, ig + 3);
    scatter(results(ig).r_true, results(ig).r_fit_corr, 10, '.');
    hold on; hr = refline(1,0); set(hr, 'color', 'k');
    xlabel('Ground truth r (\mum)'); ylabel('Fitted r (\mum)');
    xlim(Xrange(1,:)); ylim(Xrange(1,:));
    pbaspect([1 1 1]); box on; grid on;
    title(sprintf('%s\nCorrected', results(ig).group_name));
end

%% Save results
save('bscale_noise_propagation_results.mat', 'results', 'Xrange', 'b_scale_samples', 'group_names', ...
     'r_true', 'f_true', 'fcsf_true', 'DeR_true', 'prc');
