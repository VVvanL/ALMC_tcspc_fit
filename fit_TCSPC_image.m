function fit_TCSPC_image(data, dt, bin_size_xy, bin_size_t)
% takes a full TCSPC image, extracts the IRF, bins it an fits it
% data: 3d-dataset with axis x,y,t
% n_IRF_peaks: number of IRF peaks

% get summed
data_xysum = squeeze(sum(data,[1,2]));
data_xysum = data_xysum(:);
irf = extract_IRF_2(data_xysum, 0.005, 0.7, true);
irf = irf(:);

% apply bin_size_t
irf_tbin = bin_array(irf, bin_size_t, 1);
data_xysum_tbin = bin_array(data_xysum, bin_size_t,1);

% get t axis
t = (0:length(irf_tbin)-1)*dt*bin_size_t;

% fit summed data with LS
x0 = [1,4,0.15];
lb = [0.2,1,-1];
ub = [3,6,1];
cost_type = 'MLE';
fit_bg = true;
error_type = '95CI';
[results, final_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data_xysum_tbin, irf_tbin, x0, lb, ub, cost_type, fit_bg, error_type);
figure;
semilogy(t,data_xysum_tbin,'.');
hold on;
semilogy(t,final_fit);
hold off
results


% 