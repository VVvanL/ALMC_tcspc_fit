function params = setTCSPC_fit_parameters()

% function to set parameters for TCSPC fitting regimes

params = struct();

params.dt = 0.025; % size of time bin in ns
params.bin_size_xy = 7;    % xy bin size in pixels, must be odd
params.bin_size_t = 4;     % bin size in time direction
params.cost_type = 'MLE';
params.fit_bg = true;
params.fit_shift = false;
params.error_type = '95CI';
params.show_irf_estimate = true;

params.thr_snr = 3; % count distribution cutoff based on assumed SNR

params.x0 = [1 4 7]; % initial life time guess for fit of overall counts that also fits the irf
params.lb = [0.1];
params.ub = [10];
% fit parameters for the irf
params.x0_irf = [1.2, 0.15, 1.2, 0.3, 0.05, 13.5, 0.2, 0.02, 13.5, 0.3, 0.02];
params.lb_irf = [0.5, 0.1, 0.5, 0.1, 0.005, 12.8, 0.2, 0.002, 12.7, 0.2, 0.002];
params.ub_irf = [1.9, 1.5, 1.9, 3.0, 0.500, 14.2, 2.0, 0.20, 14.2, 3.0, 0.20];

params.n_exp_im_fit = 1; % number of exponents fitted for image