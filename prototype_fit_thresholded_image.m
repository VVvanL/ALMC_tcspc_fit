% script for prototyping fitting masked datasets from TCSPC image

%% fit parameters
% select fit option
% 1: uses autofit_tcspc_image
% 2: uses fit_irf_and_tcspc_image
fit_option = 2; 
% fit_type = 1; % 1 for monoexponential; 2 for biexponential
n_exp_im_fit = 1; % number of exponents fitted for image 
% monoexponential fit parameters
x0 = [1]; lb = [0.1]; ub = [10];
x0_irf = [1.2, 0.15, 1.2, 0.3, 0.05, 13.5, 0.2, 0.02, 13.5, 0.3, 0.02];
lb_irf = [0.5, 0.1, 0.5, 0.1, 0.005, 12.8, 0.2, 0.002, 12.7, 0.2, 0.002];
ub_irf = [1.9, 1.5, 1.9, 3.0, 0.500, 14.2, 2.0, 0.20, 14.2, 3.0, 0.20];
cost_type = 'MLE';
fit_bg = true;
fit_shift = false;
error_type = '95CI';

% parameters necessary for both options
dt = 0.025; % size of time bin in ns
bin_size_xy = 7;    % xy bin size in pixels, must be odd
bin_size_t = 4;     % bin size in time direction
threshold = 0;  % minimum number of photon in TCSPC trace of a pixel so that the pixel will be fitted
show_irf_estimate = true;



%%
% select file
[file,path] = uigetfile('*.*','select TCSPC file');
% get get dimensions of dataseries in file
[warnings,XYZTC] = evalc('bf_file_info(strcat(path,file))');
dataseries = 1;
[warnings, raw_data] = evalc('bf_load_parts_v7(strcat(path,file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
data = squeeze(raw_data);
disp('Done')

% generate time axis t
figure;
t = (0:size(data,3)-1)*dt;
% sum up all the pixels in the image
tcspc_trace = squeeze(sum(data,[1,2]));
% plot in semilogarithmic plot
semilogy(t,tcspc_trace);
grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of all pixels')



% turn on parallel pool if available
if canUseParallelPool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
end


r_im = ...
    fit_irf_and_tcspc_image(data, dt, x0, x0_irf, bin_size_xy, bin_size_t, n_exp_im_fit, threshold, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate);