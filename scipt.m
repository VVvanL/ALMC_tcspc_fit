% test script 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set fit parameters
% data = raw_data;
dt = 0.025;
bin_size_xy = 5;
bin_size_t = 4;
n_exp = 2;
cost_type = 'MLE';
fit_bg = true;
fit_shift = false;
error_type = '95CI';
threshold = 10000;
show_irf_estimate = true;

% turn on parallel pool if available
if canUseParallelPool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
end

r_im = autofit_tcspc_image(data, dt, bin_size_xy, bin_size_t, threshold, n_exp, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate);


%% data output to be saved in ome.tiff files using file
addpath('bfmatlab');
javaaddpath('bfmatlab\bioformats_package.jar');

save_fit_data_as_tiff(r_im, data, 'test.ome.tiff');

%% plot ratios
numerator_ind = 1;
denominator_ind = 2;
parameter = 'total_counts';

mask = r_im.mask;
if strcmp(parameter,'total_counts')
    num = r_im.amplitudes(:,:,numerator_ind).*r_im.taus(:,:,numerator_ind);
    den = r_im.amplitudes(:,:,denominator_ind).*r_im.taus(:,:,denominator_ind);
else 
    num = r_im.(parameter);
    num = num(:,:,numerator_ind);
    den = r_im.(parameter);
    den = den(:,:,denominator_ind);
end

ratio = zeros(size(den));
ratio(mask) = num(mask)./den(mask);
figure;
imagesc(ratio);

cmap = parula(256);
cmap(1,:) = [0 0 0]; 
colormap(cmap);

imagesc(ratio);
clim([min(ratio(mask)), max(ratio(mask))]); % This clips all values <= 0.3 to the first color (black)
colorbar;

%% plot fit result in pixel
i = 26;
j = 211;
use_log = true;

if r_im.mask(i,j)
    pixel_results = get_pixel_results(r_im, i, j);
    plot_tcspc_results(r_im.t, r_im.raw_data(i,j,:), r_im.irf, pixel_results, use_log)
else
    disp('Pixel was not fitted.')
end


