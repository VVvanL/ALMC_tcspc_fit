% Working script for processing/analyzing TCSPC images
clearvars; close all

make_subdirectories = false;  % flag as 'true' if individual images need to be organized into sub-directories

conditions = {'CaN_mTq', 'CaNwt_ST', 'sReach_mTq2'};
cnd_n = length(conditions);

params = setTCSPC_fit_parameters();

%% select parent directory

folderN = uigetdir; folderN = [folderN,filesep];
foldparts = strsplit(folderN,filesep); dirname = foldparts{end-1}; clear foldparts

if make_subdirectories; create_subdirectories(folderN, params.im_ext); end %#ok<UNRCH>
sublist = dir(folderN); sublist = sublist([sublist.isdir]); sublist(1:2) = []; sub_n = size(sublist,1);

% set up structure for aggregated data from experiment
aggregate_data = struct();
% image_data = struct();

%% loop through image data
for s = 1:sub_n
    subname = sublist(s).name; subpath = fullfile(sublist(s).folder,subname,filesep);
    disp(['Processing directory ', subname, '...'])
    fig_dir = [subpath, subname, '_figures',filesep];
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end  
    
    % structure to store TCSPC and fit data
    fit_data = struct();

    % find TCSPC image file in sub-directory and load
    dataseries = 1;
    im_file = [subpath, subname, '.obf'];
    [~ , raw_data] = ...
        evalc('bf_load_parts_v7(strcat(im_file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
    data = squeeze(raw_data);

    % sum up all time bins to generate normal 2D image,
    data_t_sum = squeeze(sum(data,3));
    % h_sum = plot_intensity_image(data_t_sum); % plot image

    % add title, save figure
    % calculate  bin_xy image (2D image)
    data_t_sum_xy_bin= conv2(data_t_sum, ones(params.bin_size_xy, params.bin_size_xy), 'same');
    h_binsum = plot_intensity_image(data_t_sum_xy_bin); % plot image
    savefig(h_binsum, [fig_dir, subname, '_binned_sum.fig'])
    
    % determine threshold for total mask and pixel fitting  (function call)
    [params, h_hist, h_mask] = determine_count_threshold(data_t_sum_xy_bin, params);
    savefig(h_hist, [fig_dir, subname, '_count_hist.fig'])
    savefig(h_mask, [fig_dir, subname, '_mask.fig'])
    
    %% calculate bin_t / bin_xy image
    im_data_tbin = bin_array(data, params.bin_size_t, 3);
    n_layers = size(im_data_tbin, 3);
    im_data_tbin_xybin = im_data_tbin;
    for i = 1:n_layers
        im_data_tbin_xybin(:,:,i) = conv2(im_data_tbin(:,:,i), ones(params.bin_size_xy, params.bin_size_xy), 'same');
    end
    t_bin = (0:n_layers - 1) * (params.dt * params.bin_size_t); % convert bins to seconds   

    % calculate mask TCSPC from binned xy, binned t image (global mask fit)
    mask_data_xy_sum = zeros(1,1,n_layers);
    for i = 1 : n_layers
        dmy = im_data_tbin_xybin(:,:,i);
        mask_data_xy_sum(1,1,i) = sum(dmy(params.mask),'all');
    end

    TCSPC_trace = squeeze(mask_data_xy_sum);

    fit_data.t_bin = t_bin;
    fit_data.TCSPC_trace = TCSPC_trace;

    %% Fit IRF and data with monoexponential fit
    fit_type = 1;
    params.x0 = 3; params.lb = 0.1; params.ub = 10; % parameters specific for monoexponential fit
    [r_fitirf, r_fitirf_fit, irf_fit] = fit_tcspc_gauss_irf_varpro(t_bin, mask_data_xy_sum, params);
    
    % plot data trace with fit
    h_mono = plot_TCSPC_fit(t_bin, TCSPC_trace, r_fitirf, r_fitirf_fit, irf_fit, fit_type);
    title([subname, ': monoexponential fit'], 'Interpreter','none')
    savefig(h_mono, [fig_dir, subname, '_TCSPC_monoexp_fit.fig'])

    fit_data.mono.r_fitirf = r_fitirf;
    fit_data.mono.r_fitirf_fit = r_fitirf_fit;
    fit_data.mono.ifr_fit = irf_fit;

    %% Fit IRF and data with biexponential fit
    fit_type = 2;
    params.x0 = [1,4]; params.lb = [0.1, 2]; params.ub = [5, 10]; % parameters specific for biexponential fit
    [r_fitirf, r_fitirf_fit, irf_fit] = ...
        fit_tcspc_gauss_irf_varpro(t_bin, mask_data_xy_sum, params);
    
    % plot data trace with fit
    h_bit = plot_TCSPC_fit(t_bin, TCSPC_trace, r_fitirf, r_fitirf_fit, irf_fit, fit_type);
    title([subname, ': biexponential fit'], 'Interpreter','none')
    savefig(h_bit, [fig_dir, subname, '_TCSPC_biexp_fit.fig'])

    fit_data.bi.r_fitirf = r_fitirf;
    fit_data.bi.r_fitirf_fit = r_fitirf_fit;
    fit_data.bi.ifr_fit = irf_fit;

    save([subpath, subname, '_data.mat'], 'params', 'fit_data')
    close all

    acq_field = matlab.lang.makeValidName(subname);
    aggregate_data.(acq_field).fit_data = fit_data;

end

%% create aggregate data-table for tau values (monoexponential only temp)

tau_data = struct();
tau_data.condition = {};
tau_data.tau = [];

acq_names = fieldnames(aggregate_data);
acq_n = length(acq_names);

for acq = 1:acq_n
    acq_field = acq_names{acq};

    for cnd = 1:cnd_n
        if contains(acq_field,conditions{cnd})
            cnd_str = conditions{cnd}; 
        else
            continue            
        end

        tau_data.condition = vertcat(tau_data.condition, cnd_str);
        
        tau = aggregate_data.(acq_field).fit_data.mono.r_fitirf.taus;
        tau_data.tau = vertcat(tau_data.tau, tau);
    
    end

end

tau_data.condition = categorical(tau_data.condition);
figure; 
boxchart(tau_data.condition, tau_data.tau)


save([folderN, dirname, '_fitdata.mat'], 'params', 'aggregate_data', 'tau_data')