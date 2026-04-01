function r_im = autofit_tcspc_image(data, dt, bin_size_xy, bin_size_t, threshold, n_exp, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate)
% automatic fit function for TCSPC images
% set parameters
% data = raw_data;
% dt = 0.025;
% bin_size_xy = 5;
% bin_size_t = 4;
% n_exp = 2;
% cost_type = 'MLE';
% fit_bg = true;
% fit_shift = false;
% error_type = '95CI';
% threshold = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some parameters that usually do not need to be set
fit_bg_init = true;
fit_shift_init = true;
fit_bg_fitirf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get irf estimate
data_xysum = squeeze(sum(data,[1,2]));
data_xysum = data_xysum(:);
[irf_extract,fp] = extract_IRF_2(data_xysum, 0.005, 0.7, show_irf_estimate);
irf_extract = irf_extract(:);

% find fit parameters with IRF estimate and full time resolution
% get t axis
t = (0:length(irf_extract)-1)*dt;
t = t(:);

% fit summed data for initial lifetime using extracted PSF
switch n_exp
    case 1
        x0 = [2,0.15];    lb = [0.1,-1];    ub = [10,1];
    case 2
        x0 = [1,4,0.15];    lb = [0.1,1,-1];    ub = [3,10,1];
    case 3
        x0 = [1,3,7,0.15];    lb = [0.1,1,4,-1];    ub = [3,10,20,1];
    otherwise
        disp('number of exponential components must be 1, 2, or 3')
end

[r_init, r_init_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data_xysum, irf_extract, x0, lb, ub, cost_type, fit_bg_init, fit_shift_init, error_type);

if show_irf_estimate
    figure;
    semilogy(t,data_xysum,'.');
    hold on;
    semilogy(t,r_init_fit);
    hold off;
    drawnow;
end

%% run fit to get irf with t_bin applied
% bin TCSPC data and get appropriate t axis
data_xysum_tbin = bin_array(data_xysum, bin_size_t,1);
t_binned = (0:length(data_xysum_tbin)-1)*dt*bin_size_t;
t_binned = t_binned(:);

% make x0_irf, lb_irf, and ub_irf
[pks,locs,w]=findpeaks(irf_extract);
[pks,I] = sort(pks,'descend');
locs = locs(I);
w = w(I);
max_t_shift = 0.25;
x0_irf = [(locs(1)-1)*dt + r_init.shift, w(1)/2.35482*dt];
lb_irf = [x0_irf(1)- max_t_shift, 0.5 * x0_irf(2)];
ub_irf = [x0_irf(1)+ max_t_shift, 2 * x0_irf(2)];
% lb_irf = x0_irf;
% ub_irf = x0_irf;
if length(pks)>1
    for i=2:length(pks)
        x0_irf = [x0_irf, (locs(i)-1)*dt + r_init.shift, w(i)/2.35482*dt, pks(i)/pks(1)];
        lb_irf = [lb_irf, x0_irf(3 + (i-2)*3) - max_t_shift, 0.5 * x0_irf(4 + (i-2)*3), x0_irf(5 + (i-2)*3)/2];
        ub_irf = [ub_irf, x0_irf(3 + (i-2)*3) + max_t_shift, 2 * x0_irf(4 + (i-2)*3), x0_irf(5 + (i-2)*3)*2];
    end
end

% make x0, lb and ub
x0 = r_init.taus;
lb = r_init.taus/3;
ub = r_init.taus*3;
[r_fitirf, r_fitirf_fit, irf_fit] = fit_tcspc_gauss_irf_varpro(t_binned, data_xysum_tbin, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg_fitirf, error_type);
if show_irf_estimate
figure;
    semilogy(t_binned,data_xysum_tbin,'.');
    hold on;
    semilogy(t_binned,r_fitirf_fit);
    hold off
drawnow;
end

%% TCSPC image fit with outer parallel loop and progress update in the command windo
% Pre-processing & Setup
im_data_tbin = bin_array(data, bin_size_t, 3);
n_layers = size(im_data_tbin, 3);
im_data_tbin_xybin = im_data_tbin;
for i = 1:n_layers
    im_data_tbin_xybin(:,:,i) = conv2(im_data_tbin(:,:,i), ones(bin_size_xy, bin_size_xy), 'same');
end

switch n_exp
    case 1
        x0 = [r_fitirf.taus',0];   lb = [0.1, x0(3)-0.25];     ub = [10, x0(3)+0.25];
    case 2
        x0 = [r_fitirf.taus',0];   lb = [0.1, 1, x0(3)-0.25];  ub = [5, 10, x0(3)+0.25];
    case 3
        x0 = [r_fitirf.taus',0];   lb = [0.1, 1, 3, x0(3)-0.25]; ub = [5, 10, 20, x0(3)+0.25];
end

total_count_im = sum(im_data_tbin_xybin, 3);
mask = total_count_im > threshold;
im_size = size(im_data_tbin_xybin);

% Pre-allocate
taus = zeros([im_size(1:2), n_exp]);
amplitudes = zeros([im_size(1:2), n_exp]);
shift = zeros([im_size(1:2), 1]);
err_taus = zeros([im_size(1:2), n_exp]);
err_amps = zeros([im_size(1:2), n_exp]);
err_shift = zeros([im_size(1:2), 1]);
bg = zeros([im_size(1:2), 1]);
err_bg = zeros([im_size(1:2), 1]);

% Parallel Processing with Timers
q = parallel.pool.DataQueue;
total_rows = im_size(1);
start_time = tic;

% Attach progress handler
afterEach(q, @(x) updateProgress(total_rows, start_time));

fprintf('Starting Parallel Fit...\n');
fprintf('Progress: '); 

if canUseParallelPool
    parfor i = 1:total_rows
        % --- Local row initialization ---
        row_taus      = zeros(im_size(2), n_exp);
        row_amps      = zeros(im_size(2), n_exp);
        row_shift     = zeros(im_size(2), 1);
        row_err_taus  = zeros(im_size(2), n_exp);
        row_err_amps  = zeros(im_size(2), n_exp);
        row_err_shift = zeros(im_size(2), 1);
        row_bg        = zeros(im_size(2), 1);
        row_err_bg    = zeros(im_size(2), 1);
    
        for j = 1:im_size(2)
            if mask(i,j)
                [r_dmy, ~] = fit_tcspc_varpro_bgoption_werrorbars(...
                    t_binned, squeeze(im_data_tbin_xybin(i,j,:)), irf_fit, x0, lb, ub, ...
                    cost_type, fit_bg, fit_shift, error_type);
        
                row_taus(j,:)      = r_dmy.taus;
                row_amps(j,:)      = r_dmy.amplitudes;
                row_shift(j)       = r_dmy.shift;
                row_err_taus(j,:)  = r_dmy.err_vals.taus;
                row_err_amps(j,:)  = r_dmy.err_vals.amplitudes;
                row_err_shift(j)   = r_dmy.err_vals.shift;
                if fit_bg
                    row_bg(j)      = r_dmy.background;
                    row_err_bg(j)  = r_dmy.err_vals.background;
                end
            end
        end
        
        % --- Assign Row Back ---
        taus(i,:,:)      = row_taus;
        amplitudes(i,:,:) = row_amps;
        shift(i,:)       = row_shift;
        err_taus(i,:,:)  = row_err_taus;
        err_amps(i,:,:)  = row_err_amps;
        err_shift(i,:)   = row_err_shift;
        bg(i,:)          = row_bg;
        err_bg(i,:)      = row_err_bg;
    
        send(q, i); 
    end
else
    for i = 1:total_rows
        % --- Local row initialization ---
        row_taus      = zeros(im_size(2), n_exp);
        row_amps      = zeros(im_size(2), n_exp);
        row_shift     = zeros(im_size(2), 1);
        row_err_taus  = zeros(im_size(2), n_exp);
        row_err_amps  = zeros(im_size(2), n_exp);
        row_err_shift = zeros(im_size(2), 1);
        row_bg        = zeros(im_size(2), 1);
        row_err_bg    = zeros(im_size(2), 1);
    
        for j = 1:im_size(2)
            if mask(i,j)
                [r_dmy, ~] = fit_tcspc_varpro_bgoption_werrorbars(...
                    t_binned, squeeze(im_data_tbin_xybin(i,j,:)), irf_fit, x0, lb, ub, ...
                    cost_type, fit_bg, fit_shift, error_type);
        
                row_taus(j,:)      = r_dmy.taus;
                row_amps(j,:)      = r_dmy.amplitudes;
                row_shift(j)       = r_dmy.shift;
                row_err_taus(j,:)  = r_dmy.err_vals.taus;
                row_err_amps(j,:)  = r_dmy.err_vals.amplitudes;
                row_err_shift(j)   = r_dmy.err_vals.shift;
                if fit_bg
                    row_bg(j)      = r_dmy.background;
                    row_err_bg(j)  = r_dmy.err_vals.background;
                end
            end
        end
        
        % --- Assign Row Back ---
        taus(i,:,:)      = row_taus;
        amplitudes(i,:,:) = row_amps;
        shift(i,:)       = row_shift;
        err_taus(i,:,:)  = row_err_taus;
        err_amps(i,:,:)  = row_err_amps;
        err_shift(i,:)   = row_err_shift;
        bg(i,:)          = row_bg;
        err_bg(i,:)      = row_err_bg;
    
        send(q, i); 
    end
end

% Final Total Time Calculation
total_duration = toc(start_time);
n_fits = sum(mask,'all');
fprintf('\nProcessing Complete. Total processing time was %s. %d pixels were fitted. %d fits per second. \n', formatTime(total_duration), n_fits, round(n_fits/total_duration));

% Re-assemble Results
r_im.taus = taus;
r_im.amplitudes = amplitudes;
r_im.shift = shift;
r_im.err_vals.taus = err_taus;
r_im.err_vals.amplitudes = err_amps;
r_im.err_vals.shift = err_shift;
if fit_bg
    r_im.background = bg;
    r_im.err_vals.background = err_bg;
end

r_im.irf = irf_fit;
r_im.t = t_binned;
r_im.cost_type = cost_type;
r_im.error_type = error_type;
r_im.mask = mask;
r_im.total = total_count_im;
r_im.raw_data = im_data_tbin_xybin;