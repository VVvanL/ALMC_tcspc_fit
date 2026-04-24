%% script for fitting a TCSPC image using fit_irf_and_tcspc_image
%
%% set parameters and set up Matlab environment:
% parameters
dt = 0.025; % size of time bin in ns
bin_size_xy = 5; % xy bin size in pixels, must be odd
bin_size_t = 2; % bin size in time direction
threshold = 5000; % minimum number of photon in TCSPC trace of a pixel so that the pixel will be fitted
x0 = [4]; % initial life time guess 
x0_irf = [1.3, 0.15, 1.3, 0.3, 0.05, 13.5, 0.2, 0.02, 13.5, 0.3, 0.01];
fit_bg = true;
fit_shift = false;
cost_type = 'MLE';
error_type = '95CI';
show_irf_estimate = true;

% turn on parallel pool if available
if canUseParallelPool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    endsudo 
end

% add path to bioformats
addpath('bfmatlab');
javaaddpath('bfmatlab\bioformats_package.jar');

%% select TCSPC data and show data dimensions of dataseries in file
% select file
[file,path] = uigetfile('*.*','select TCSPC file');

% get get dimensions of dataseries in file
XYZTC = bf_file_info(strcat(path,file));
% display dimensions of dateseries
disp(XYZTC);

%% load data
dataseries = 1;
raw_data = bf_load_parts_v7(strcat(path,file),dataseries,-1,-1,-1,-1,-1);
data = squeeze(raw_data);

%% display TCSPC trace, helpful to see time positions of laser pulse, i.e. t0 for x0_irf
t = (0:size(data,3)-1)*dt;
tcspc_trace = squeeze(sum(data,[1,2]));
semilogy(t,tcspc_trace);
grid on;
xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels')

%% run the fit
% adjust x0_irf
x0_irf = [2.5, 0.15, 2.5, 0.3, 0.05, 15, 0.2, 0.02, 15, 0.3, 0.01];
r_im = fit_irf_and_tcspc_image(data, dt, x0, x0_irf, bin_size_xy, bin_size_t, threshold, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate);

%% plot false colour image
imagesc(r_im.taus(:,:,1));
colorbar;
axis equal

%% plot ratios
numerator_ind = 1;
denominator_ind = 2;
parameter = 'taus';

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
i = 40;
j = 70;
use_log = true;
if r_im.mask(i,j)
    pixel_results = get_pixel_results(r_im, i, j);
    plot_tcspc_results(r_im.t, r_im.raw_data(i,j,:), r_im.irf, pixel_results, use_log);
else
    disp('Pixel was not fitted.')
end