% scratch code for prototyping TCSPC image fitting

% load file using standard script


%% Mask generation
% sum up all time bins to generate normal 2D image
data_t_sum = squeeze(sum(data,3));
figure; 
imagesc(data_t_sum);
colormap('parula');
colorbar;
axis equal

t = (0:size(data,3)-1)*params.dt;

% plot histogram of binned counts
figure; histogram(data_t_sum , 'Normalization', 'percentage')
figure; histogram(data_t_sum)

% set threshold based on count distribution
params.max_count = max(data_t_sum, [], 'all');
params.threshold = params.max_count / params.thr_snr;
params.mask = data_t_sum > params.threshold;

figure; imagesc(params.mask); axis equal

%% Sum up all pixels in mask
n_t_bins = size(data,3);
data_xy_sum = zeros(1,1,n_t_bins);
for i = 1 : n_t_bins
    dmy = data(:,:,i);
    data_xy_sum(1,1,i) = sum(dmy(params.mask),'all');
end

%% Fit IRF and data with biexponential fit
params.x0 = [1,4]; params.lb = [0.1, 2]; params.ub = [5, 10]; % parameters specific for biexponential fit
[r_fitirf, r_fitirf_fit, irf_fit] = ...
    fit_tcspc_gauss_irf_varpro(t, data_xy_sum, params);

% plot decay with IRFs and fit
figure;
semilogy(t,squeeze(data_xy_sum),'.','DisplayName','data');
hold on;
semilogy(t,r_fitirf_fit,'DisplayName','fit');
semilogy(t,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
hold off
grid on;
legend;
xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels in mask with biexponential fit');
% ylim([min([min(r_fitirf_fit) min(mask_data_xy_sum)]) max(r_fitirf_fit)*1.05]);
ylim([min(r_fitirf_fit)*0.9 max(r_fitirf_fit)*1.1])



