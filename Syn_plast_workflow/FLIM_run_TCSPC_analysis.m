% Working script for processing/analyzing TCSPC images
clearvars; close all

make_subdirectories = false;  % flag as 'true' if individual images need to be organized into sub-directories
params = setTCSPC_fit_parameters();

%% select parent directory

folderN = uigetdir; folderN = [folderN,filesep];
foldparts = strsplit(folderN,filesep); dirname = foldparts{end-1}; clear foldparts

if make_subdirectories; create_subdirectories(folderN, params.im_ext); end %#ok<UNRCH>
sublist = dir(folderN); sublist = sublist([sublist.isdir]); sublist(1:2) = []; sub_n = size(sublist,1);

for s = 1:sub_n
    subname = sublist(s).name; subpath = fullfile(sublist(s).folder,subname,filesep);
    disp(['Processing directory ', subname, '...'])
    fig_dir = [subpath, subname, '_figures',filesep];
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    % set up structure for fit data and image data
    fit_data = struct();
    image_data = struct(); 

    % find TCSPC image file in sub-directory and load
    dataseries = 1;
    im_file = [subpath, subname, '.obf'];
    [~ , raw_data] = ...
        evalc('bf_load_parts_v7(strcat(im_file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
    data = squeeze(raw_data);

    % sum up all time bins to generate normal 2D image,
    data_t_sum = squeeze(sum(data,3));
    h_sum = plot_intensity_image(data_t_sum); % plot image
    % add title, save figure
    % calculate  bin_xy image (2D image)
    data_t_sum_xy_bin= conv2(data_t_sum, ones(params.bin_size_xy, params.bin_size_xy), 'same');
    h_binsum = plot_intensity_image(data_t_sum_xy_bin); % plot image
    
    % determine threshold for total mask and pixel fitting
    [params, h_counthist, h_mask] = determine_count_threshold(data_t_sum_xy_bin, params);

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




end