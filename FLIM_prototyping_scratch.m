% scratch code to customize FLIM fitting and tau mapping
clearvars; close all

params = setTCSPC_fit_parameters();

%% select file and load
[file,path] = uigetfile('*.*','select TCSPC file');
% get get dimensions of dataseries in file
[~,XYZTC] = evalc('bf_file_info(strcat(path,file))');
dataseries = 1;
[warnings, raw_data] = evalc('bf_load_parts_v7(strcat(path,file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
data = squeeze(raw_data);
disp('Done')


%% Mask generation
% sum up all time bins to generate normal 2D image
data_t_sum = squeeze(sum(data,3));
% threshold = prctile(data_t_sum,75,'all'); % set threshold based on count distribution
% mask = data_t_sum > threshold;
% plot image intensity map 
figure; 
imagesc(data_t_sum);
colormap('parula');
colorbar;
axis equal
% hold on;
% % overlay mask
% red_mask = cat(3, ones(size(mask)), zeros(size(mask)), zeros(size(mask)));
% h_mask = imshow(red_mask);
% opacity = 0.3;
% set(h_mask, 'AlphaData', mask * opacity);
% hold off;

%% calculate bin_t / bin_xy image

im_data_tbin = bin_array(data, params.bin_size_t, 3);
n_layers = size(im_data_tbin, 3);
im_data_tbin_xybin = im_data_tbin;
for i = 1:n_layers
    im_data_tbin_xybin(:,:,i) = conv2(im_data_tbin(:,:,i), ones(params.bin_size_xy, params.bin_size_xy), 'same');
end
total_count_im = sum(im_data_tbin_xybin, 3); % not really useful, just convolve sum image
figure; 
imagesc(total_count_im);
colormap('parula');
colorbar;
axis equal

% binned image histogram (all pixels)
tcspc_trace = squeeze(sum(im_data_tbin_xybin,[1,2]));
t_bin = (0:n_layers - 1) * (params.dt * params.bin_size_t);
figure; semilogy(t_bin,tcspc_trace);
grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of all pixels (tbin,xybin)')

% plot histogram of binned counts
figure; histogram(total_count_im, 'Normalization', 'percentage')
% set threshold based on count distribution
params.max_count = max(total_count_im, [], 'all');
params.threshold = params.max_count / params.thr_snr;
params.mask = total_count_im > params.threshold;

% calculate mask TCSPC image (xy and t binned)
mask_data_xy_sum = zeros(1,1,n_layers);
for i = 1 : n_layers
    dmy = im_data_tbin_xybin(:,:,i);
    mask_data_xy_sum(1,1,i) = sum(dmy(params.mask),'all');
end

% plot in semilogarithmic plot
figure;
semilogy(t_bin,squeeze(mask_data_xy_sum));
grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of pixels in mask')

% plot in semilogarithmic plot
figure;
% t = (0:size(data,3)-1)*dt;
t = (0:n_t_bins- 1) * (params.dt * params.bin_size_t);
semilogy(t,squeeze(data_xy_sum));
grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of pixels in mask')

%% Fit IRF and data with biexponential fit
params.x0 = [1,4]; params.lb = [0.1, 2]; params.ub = [5, 10]; % parameters specific for biexponential fit
[r_fitirf, r_fitirf_fit, irf_fit] = ...
    fit_tcspc_gauss_irf_varpro(t_bin, mask_data_xy_sum, params);

% plot decay with IRFs and fit
figure;
semilogy(t_bin,squeeze(mask_data_xy_sum),'.','DisplayName','data');
hold on;
semilogy(t_bin,r_fitirf_fit,'DisplayName','fit');
semilogy(t_bin,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
hold off
grid on;
legend;
xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels in mask with biexponential fit');
% ylim([min(r_fitirf_fit) min(mask_data_xy_sum)]) max(r_fitirf_fit)*1.05]);

ant_str = ([char(964), '_1: ',num2str(r_fitirf.taus(1)),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1)),' ns', ...
    newline, char(964), '_2: ',num2str(r_fitirf.taus(2)),' ns', char(177),' ', num2str(r_fitirf.err_vals.taus(2)),' ns', ...
    newline, 'background: ',num2str(r_fitirf.background), char(177),' ', num2str(r_fitirf.err_vals.background)]);
dim = [.66 .7, .1 .1];
annotation('textbox', 'Position',dim, String = ant_str, FontSize = 13)

%% Fit IRF and data with monoexponential fit
params.x0 = 3; params.lb = 0.1; params.ub = 10;
[r_fitirf, r_fitirf_fit, irf_fit] = fit_tcspc_gauss_irf_varpro(t_bin, mask_data_xy_sum, params);

figure; 
semilogy(t_bin,squeeze(mask_data_xy_sum),'.','DisplayName','data');
hold on;
semilogy(t_bin,r_fitirf_fit,'DisplayName','fit');
semilogy(t_bin,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
hold off
grid on;
legend;
xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels in mask with monoexponential fit');
ylim([min([min(r_fitirf_fit) min(mask_data_xy_sum)]) max(r_fitirf_fit)*1.05]);
% print fit parameters
ant_str = ([char(964), '_1: ',num2str(r_fitirf.taus(1)),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1)),' ns', ...
    newline, 'background: ',num2str(r_fitirf.background), char(177),' ', num2str(r_fitirf.err_vals.background)]);
dim = [.66 .7, .1 .1];
annotation('textbox', 'Position',dim, String = ant_str, FontSize = 13)






% plot in semilogarithmic plot
figure;
t = (0:size(data,3)-1)*params.dt;
semilogy(t,squeeze(data_xy_sum));
grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of pixels in image')

% turn on parallel pool if available
if canUseParallelPool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
end

%% run the fit
% r_im = ...
%     fit_irf_and_tcspc_image(data, dt, x0, x0_irf, bin_size_xy,  ...
%     bin_size_t, n_exp_im_fit, threshold, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate);
r_im = fit_irf_and_tcspc_image(data, params);

%% plot false colour image
figure;
imagesc(r_im.taus(:,:,2));
% colormap(parula(7))
colorbar;
axis equal

%% plot fit result in pixel
figure
i = 37;
j = 121;
use_log = true;
if r_im.mask(i,j)
    pixel_results = get_pixel_results(r_im, i, j);
    plot_tcspc_results(r_im.t, r_im.raw_data(i,j,:), r_im.irf, pixel_results, use_log);
else
    disp('Pixel was not fitted.')
end