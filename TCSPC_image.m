% test script
% to be turned into function later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters
data = raw_data;
dt = 0.025;
bin_size_xy = 5;
bin_size_t = 4;
n_exp = 2;
cost_type = 'MLE';
fit_bg_init = true;
fit_shift_init = true;
fit_bg_fitirf = true;
fit_bg = true;
fit_shift = false;
error_type = '95CI';
threshold = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get irf estimate
data_xysum = squeeze(sum(data,[1,2]));
data_xysum = data_xysum(:);
[irf_extract,fp] = extract_IRF_2(data_xysum, 0.005, 0.7, true);
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
figure;
semilogy(t,data_xysum,'.');
hold on;
semilogy(t,r_init_fit);
hold off;
drawnow;

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
figure;
semilogy(t_binned,data_xysum_tbin,'.');
hold on;
semilogy(t_binned,r_fitirf_fit);
hold off
drawnow;

%% run fit to get irf with t_bin applied 3 component gauss
% % bin TCSPC data and get appropriate t axis
% data_xysum_tbin = bin_array(data_xysum, bin_size_t,1);
% t = (0:length(data_xysum_tbin)-1)*dt*bin_size_t;
% t = t(:);
% 
% % make x0_irf, lb_irf, and ub_irf
% [pks,locs,w]=findpeaks(irf);
% [pks,I] = sort(pks,'descend');
% locs = locs(I);
% w = w(I);
% max_t_shift = 0.25;
% x0_irf = [(locs(1)-1)*dt + results.shift, w(1)/2.35482*dt, (locs(1)-1)*dt + 0.025 + results.shift, w(1)/2.35482*1.5*dt, 0.2];
% lb_irf = [x0_irf(1)- max_t_shift, 0.5 * x0_irf(2), x0_irf(1)- max_t_shift, 0.5 * x0_irf(2), 0];
% ub_irf = [x0_irf(1)+ max_t_shift, 2 * x0_irf(2), x0_irf(1)- max_t_shift, 2 * x0_irf(2), 1];
% % lb_irf = x0_irf;
% % ub_irf = x0_irf;
% if length(pks)>1
%     for i=2:length(pks)
%         x0_irf = [x0_irf, (locs(i)-1)*dt + results.shift, w(i)/2.35482*dt, pks(i)/pks(1)];
%         lb_irf = [lb_irf, x0_irf(6 + (i-2)*3) - max_t_shift, 0.5 * x0_irf(7 + (i-2)*3), x0_irf(8 + (i-2)*3)/2];
%         ub_irf = [ub_irf, x0_irf(6 + (i-2)*3) + max_t_shift, 2 * x0_irf(7 + (i-2)*3), x0_irf(8 + (i-2)*3)*2];
%     end
% end
% 
% % make x0, lb and ub
% x0 = results.taus*1.5;
% lb = results.taus;
% ub = results.taus;
% fit_bg = true;
% [results, final_fit, normalized_irf] = fit_tcspc_gauss_irf_varpro(t, data_xysum_tbin, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg, error_type);
% figure;
% semilogy(t,data_xysum_tbin,'.');
% hold on;
% semilogy(t,final_fit);
% hold off

%% fit data parfor in inner loop
% tic;
% % bin data in time_direction
% im_data_tbin = bin_array(data,bin_size_t,3);
% 
% % bin data in xy_direction
% n_layers = size(im_data_tbin,3);
% 
% im_data_tbin_xybin = im_data_tbin;
% for i = 1:n_layers
%     im_data_tbin_xybin(:,:,i) = conv2(im_data_tbin(:,:,i),ones(bin_size_xy,bin_size_xy),'same');
% end
% 
% im_size = size(im_data_tbin_xybin);
% taus = zeros([im_size(1:2), n_exp]);
% amplitudes = zeros([im_size(1:2), n_exp]);
% shift = zeros([im_size(1:2), 1]);
% err_taus = zeros([im_size(1:2), n_exp]);
% err_amps = zeros([im_size(1:2), n_exp]);
% err_shift = zeros([im_size(1:2), 1]);
% 
% if fit_bg
%     bg = zeros([im_size(1:2), 1]);
%     err_bg = zeros([im_size(1:2), 1]);
% end
% 
% switch n_exp
%     case 1
%         x0 = [r_fitirf.taus',0];   lb = [0.1 ,x0(3)-0.25];    ub = [10, x0(3)+0.25];
%     case 2
%         x0 = [r_fitirf.taus',0];   lb = [0.1,1,x0(3)-0.25];    ub = [5,10,x0(3)+0.25];
%     case 3
%         x0 = [r_fitirf.taus',0];   lb = [0.1, 1, 3, x0(3)-0.25];    ub = [5,10, 20, x0(3)+0.25];
%     otherwise
%         disp('number of exponential components must be 1, 2, or 3')
% end
% 
% mask = sum(im_data_tbin_xybin,3)>threshold;
% for i = 1:im_size(1)
%     fprintf('Processing row %d\n', i);
%     parfor j = 1:im_size(2)
%         if mask(i,j)
%             [r_dmy, ~] = fit_tcspc_varpro_bgoption_werrorbars(t_binned, im_data_tbin_xybin(i,j,:), irf_fit, x0, lb, ub, cost_type, fit_bg, fit_shift, error_type);
% 
%             % Assign to the temporary independent arrays
%             taus(i,j,:) = r_dmy.taus;
%             amplitudes(i,j,:) = r_dmy.amplitudes;
%             shift(i,j) = r_dmy.shift;
%             err_taus(i,j,:) = r_dmy.err_vals.taus;
%             err_amps(i,j,:) = r_dmy.err_vals.amplitudes;
%             err_shift(i,j) = r_dmy.err_vals.shift;
% 
%             if fit_bg
%                 bg(i,j) = r_dmy.background;
%                 err_bg(i,j) = r_dmy.err_vals.background;
%             end
%         end
%     end
% end
% 
% r_im.taus = taus;
% r_im.amplitudes = amplitudes;
% r_im.shift = shift;
% r_im.err_vals.taus = err_taus;
% r_im.err_vals.amplitudes = err_amps;
% r_im.err_vals.shift = err_shift;
% 
% if fit_bg
%     r_im.background = bg;
%     r_im.err_vals.background = err_bg;
% end
% toc;

%% TCSPC image fit with outer parallel loop and progress update in the command windo
% 1. Pre-processing & Setup
% ... [Keep your binning and switch-case logic as it was] ...
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

mask = sum(im_data_tbin_xybin, 3) > threshold;
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

% 2. Parallel Processing with Timers
q = parallel.pool.DataQueue;
total_rows = im_size(1);
start_time = tic;

% Attach progress handler
afterEach(q, @(x) updateProgress(total_rows, start_time));

fprintf('Starting Parallel Fit...\n');
fprintf('Progress: '); 

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

% Final Total Time Calculation
total_duration = toc(start_time);
n_fits = sum(mask,'all');
fprintf('\nProcessing Complete. Total processing time was %s. %d pixels were fitted. %d fits per second. \n', formatTime(total_duration), n_fits, round(n_fits/total_duration));

% 3. Re-assemble Results
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


%% plot fit result in pixel
i = 26;
j = 111;
use_log = true;
if mask(i,j)
    pixel_results = get_pixel_results(r_im, i, j);
    plot_tcspc_results(t_binned, im_data_tbin_xybin(i,j,:), irf_fit, pixel_results, use_log)
else
    disp('Pixel was not fitted.')
end

%% plot taus and amplitudes

%% plot ratios
numerator_ind = 1;
denominator_ind = 2;
parameter = 'total_counts';

mask = sum(im_data_tbin_xybin, 3) > threshold;
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

%% data output to be saved in ome.tiff files
% total counts, fit mask, taus, amplitudes, shift if fitted, background if fitted, error vals of
% all fitted parameters ratio taus, ratio amplitudes, ratio total_counts

switch n_exp
    case 1
        n_datalayers = 2 + n_exp*4 + fit_bg * 2 + fit_shift*2;
    case 2
        n_datalayers = 2 + n_exp*4 + fit_bg * 2 + fit_shift*2 + 3;
    case 3
        n_datalayers = 2 + n_exp*4 + fit_bg * 2 + fit_shift*2 + 3;
end
channel_names = cell(1,n_datalayers);

data_2_tiff = zeros([size(r_im.taus,[1:2]),n_datalayers]);

cl = 1; % keep track of current layer

% total counts
data_2_tiff(:,:,cl) = sum(raw_data,3);
channel_names{cl} = 'total counts';
cl = cl + 1;

% fit mask
data_2_tiff(:,:,2) = mask;
channel_names{cl} = 'fit mask';
cl = cl + 1;

% taus
for i = 1:n_exp
    data_2_tiff(:,:,cl-1+i) = r_im.taus(:,:,i);
    channel_names{cl} = strcat('lifetime ',num2str(i));
    cl = cl + 1;
end

% amplitudes
for i = 1:n_exp
    data_2_tiff(:,:,cl-1+i) = r_im.amplitudes(:,:,i);
    channel_names{cl} = strcat('amplitude ',num2str(i));
    cl = cl + 1;
end

% background if fitted
if fit_bg
    data_2_tiff(:,:,cl) = r_im.background;
    channel_names{cl} = 'background';
    cl = cl+1;
end

% shift, if fitted
if fit_shift 
    data_2_tiff(:,:,cl) = r_im.shift;
    channel_names{cl} = 't0 shift';
    cl = cl+1;
end

% tau error vals
for i = 1:n_exp
    data_2_tiff(:,:,cl-1+i) = r_im.err_vals.taus(:,:,i);
    channel_names{cl} = strcat('lifetime error ',num2str(i));
    cl = cl + 1;
end

% amplitdue error vals
for i = 1:n_exp
    data_2_tiff(:,:,cl-1+i) = r_im.err_vals.amplitudes(:,:,i);
    channel_names{cl} = strcat('amplitude error ',num2str(i));
    cl = cl + 1;
end

% ratios
if n_exp>1
    ratio = zeros(size(r_im.taus,1:2));
    num = r_im.taus(:,:,1);
    den = r_im.taus(:,:,2);
    ratio(mask) = num(mask)./den(mask);
    data_2_tiff(:,:,cl) = ratio;
    channel_names{cl} = 'ratio liftetime 1 to lifetime 2';
    cl = cl + 1;

    num = r_im.amplitudes(:,:,1);
    den = r_im.amplitudes(:,:,2);
    ratio(mask) = num(mask)./den(mask);
    data_2_tiff(:,:,cl) = ratio;
    channel_names{cl} = 'ratio amplitude 1 to amplitude 2 ';
    cl = cl + 1;

    num = r_im.amplitudes(:,:,1).*r_im.taus(:,:,2);
    den = r_im.amplitudes(:,:,1).*r_im.taus(:,:,2);
    ratio(mask) = num(mask)./den(mask);
    data_2_tiff(:,:,cl) = ratio;
    channel_names{cl} = 'ratio photon counts fit component 1 to fit component 2';
    cl = cl + 1;
end