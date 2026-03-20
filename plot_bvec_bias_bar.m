%% Bar plot of bias from bvec noise propagation results

load('bvec_noise_propagation_results.mat', 'all_results');

r_bins = [0 1; 1 2; 2 3; 3 4; 4 5];
Nbins = size(r_bins, 1);
bin_labels = arrayfun(@(i) sprintf('%g-%g', r_bins(i,1), r_bins(i,2)), 1:Nbins, 'uni', 0);

pnames = fieldnames(all_results);
Nprot = numel(pnames);
snr_map = struct('C2', 38, 'C1', 18);

colors_base   = [0.4 0.7 0.4];
colors_uncorr = [0.8 0.3 0.3];
colors_corr   = [0.3 0.3 0.8];

figure('unit','inch','position',[0 0 6*Nprot 4.5]);

for ip = 1:Nprot
    R = all_results.(pnames{ip});
    rt = R.r_true_filt;
    rb = R.r_fit_base_filt;
    ru = R.r_fit_uncorr_filt;
    rc = R.r_fit_corr_filt;

    bias_b = zeros(Nbins,1);
    bias_u = zeros(Nbins,1);
    bias_c = zeros(Nbins,1);
    for ib = 1:Nbins
        idx = rt >= r_bins(ib,1) & rt < r_bins(ib,2);
        if sum(idx) == 0; continue; end
        bias_b(ib) = mean(rb(idx) - rt(idx));
        bias_u(ib) = mean(ru(idx) - rt(idx));
        bias_c(ib) = mean(rc(idx) - rt(idx));
    end

    subplot(1, Nprot, ip);
    h = bar(categorical(bin_labels, bin_labels), [bias_b, bias_u, bias_c]);
    h(1).FaceColor = colors_base;
    h(2).FaceColor = colors_uncorr;
    h(3).FaceColor = colors_corr;
    ylabel('Bias (\mum)');
    xlabel('r_{true} (\mum)');
    legend({'Baseline','Uncorrected','Corrected'}, 'Location','northeast');
    title(sprintf('%s (SNR=%d, mean b_{scale}=%.3f)', ...
        pnames{ip}, snr_map.(pnames{ip}), mean(R.bscales_filt)));
    yline(0, 'k--', 'HandleVisibility','off');
    grid on; box on;
end

exportgraphics(gcf, 'fig_bvec_bias_bar.pdf', 'ContentType', 'vector');
fprintf('Saved fig_bvec_bias_bar.pdf\n');
