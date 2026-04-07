function r_im = fit_irf_and_tcspc_image(data, dt, x0, x0_irf, bin_size_xy, bin_size_t, threshold, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate)
% first fit IRF and TCSPC of total counts together, then fit pixel by pixel
% in TCSPC image using the optimized IRF

% data: 3d array of TCSPC with x/y/t dimension
% dt:   size of timebin in t dimension in ns
% x0:   initial lifetime fit parameter(s). Is either 1x1, 1x2 or 1x3 array,
%       depending if mono-, bi-, or triexponential fit is wanted
% x0_irf:   initial fit parameters for irf. The IRF is simulated by a 1, 2, 
%           or three gaussians. The first gaussian is defined by the time 
%           bin position t0 (in ns) and the stdev of the gaussian, the 
%           gaussians also have a scale factor to set their amplitude 
%           relative to the first gaussian. Example array:
%           [t0_gauss1, stdev_gauss1, t0_gauss2, stdev_gauss2,
%           rel_amp_gauss2, t0_gauss3, stdev_gauss3, rel_amp_gauss3]
% bin_size_xy:  number of pixels that will be binned together, must be odd,
%               the binned image still has the same number of pixels
% bin_size_t:   number of time bins that are binned together
% threshold:    minimum counts in of entire TCSPC trace in a pixel after xy 
%               binning to do TCSPC fit of that pixel
% fit_bg:   true/false, indicates whether background should be used as fit
%           parameter when fitting the TCSPC image data. background will be
%           set to 0 if false
% fit_shift:    true/false, indicates whether shift should be optimized in
%               when fitting TCSPC image (best to set as false as IRF 
%               fitting already takes care of time shift)
% cost_type: 'MLE' or 'LS', 'MLE' is more optimal when fitting photons
% error_type:   '1sigma' or '95CI', indicates which confidence interval is 
%               used for error bars
% show_irf_estimate:    true/false, whether to show irf fit and overall fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % some parameters that usually do not need to be set
fit_bg_fitirf = true;
max_t_shift = 1;  %maximum time shift of gaussian position

%% run fit to get irf
% bin TCSPC data in t dimension and get appropriate t axis
data_xysum = squeeze(sum(data,[1,2]));
data_xysum_tbin = bin_array(data_xysum, bin_size_t,1);
t_binned = (0:length(data_xysum_tbin)-1)*dt*bin_size_t;
t_binned = t_binned(:);

% set lb and ub
lb = x0/3;
ub = x0*3;

% set x0_lb and x0_ub
n_pks = ((length(x0_irf)+1)/3);
lb_irf = [x0_irf(1)- max_t_shift, 0.5 * x0_irf(2)];
ub_irf = [x0_irf(1)+ max_t_shift, 2 * x0_irf(2)];

if n_pks>1
    for i=2:n_pks
        lb_irf = [lb_irf, x0_irf(3 + (i-2)*3) - max_t_shift, 0.5 * x0_irf(4 + (i-2)*3), x0_irf(5 + (i-2)*3)/5];
        ub_irf = [ub_irf, x0_irf(3 + (i-2)*3) + max_t_shift, 2 * x0_irf(4 + (i-2)*3), x0_irf(5 + (i-2)*3)*5];
    end
end

[r_fitirf, r_fitirf_fit, irf_fit] = fit_tcspc_gauss_irf_varpro(t_binned, data_xysum_tbin, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg_fitirf, error_type);
r_im = r_fitirf;
if show_irf_estimate
    figure;
    semilogy(t_binned,data_xysum_tbin,'.','DisplayName','data');
    hold on;
    semilogy(t_binned,r_fitirf_fit,'DisplayName','fit');
    semilogy(t_binned,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
    hold off
    grid on;
    legend;
    xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels');
    ylim([min([min(r_fitirf_fit) min(data_xysum_tbin)]) max(r_fitirf_fit)*1.05])
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

n_exp = length(x0);
% switch n_exp
%     case 1
%         x0 = [r_fitirf.taus',0];   lb = [0.1, x0(3)-0.25];     ub = [10, x0(3)+0.25];
%     case 2
%         x0 = [r_fitirf.taus',0];   lb = [0.1, 1, x0(3)-0.25];  ub = [5, 10, x0(3)+0.25];
%     case 3
%         x0 = [r_fitirf.taus',0];   lb = [0.1, 1, 3, x0(3)-0.25]; ub = [5, 10, 20, x0(3)+0.25];
% end

x0 = [r_fitirf.taus',0];   lb = [0.1, 1, x0(3)-0.25];  ub = [5, 10, x0(3)+0.25];

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

%% Re-assemble Results and generate data output
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