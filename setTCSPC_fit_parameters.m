function params = setTCSPC_fit_parameters()

% function to set parameters for TCSPC fitting regimes

params = struct();

params.dt = 0.025; % size of time bin in ns
params.bin_size_xy = 5;    % xy bin size in pixels, must be odd
params.bin_size_t = 2;     % bin size in time direction
params.cost_type = 'MLE';
params.fit_bg = true;
params.fit_shift = false;
params.error_type = '95CI';
params.show_irf_estimate = true;

params.thr_snr = 4; % count distribution cutoff based on assumed SNR

params.x0 = [1 4 7]; % initial life time guess for fit of overall counts that also fits the irf
% fit parameters for the irf
params.x0_irf = [1.2, 0.15, 1.2, 0.3, 0.05, 13.5, 0.2, 0.02, 13.5, 0.3, 0.02];
params.n_exp_im_fit = 1; % number of exponents fitted for image