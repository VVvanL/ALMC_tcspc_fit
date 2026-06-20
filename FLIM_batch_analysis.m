% Working script for batch, first-pass analysis of TCSPC dataset, generate TCSPC maps
clearvars; close all

multi = 0; % set to 1 if looping through multiple condition/experimental directories (each with a set of acquisitional subdirectories)
params = setTCSPC_fit_parameters();

%% select parent directory
if multi == 1
    folderP = uigetdir; foldparts = strsplit(folderP,filesep); parent_name = foldparts{end}; clear foldparts
    dirlist = dir(folderP); dirlist = dirlist([dirlist.isdir]); dirlist(1:2) = [];
    dir_n = size(dirlist,1); folderP = [folderP,filesep];

else
    folderN = uigetdir; folderN = [folderN,filesep];
    dir_n = 1;
end

for d = 1:dir_n
    if multi == 1; folderN = [folderP,filesep,dirlist(d).name,filesep]; end
    foldparts = strsplit(folderN,filesep); dirname = foldparts{end-1}; clear foldparts
    sublist = dir(folderN); sublist = sublist([sublist.isdir]); sublist(1:2) = []; sub_n = size(sublist,1);

    for s = 1:sub_n
        subname = sublist(s).name; subpath = fullfile(sublist(s).folder,subname,filesep);
        disp(['Processing directory ', subname, '...'])

        % set up structure for fit data and image data
        fit_data = struct();
        image_data = struct();

        % find TCSPC image file in sub-directory and load
        dataseries = 1;
        im_file = [subpath, subname, '.obf'];
        [~ , raw_data] = ...
            evalc('bf_load_parts_v7(strcat(im_file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
        data = squeeze(raw_data);
        % TODO: create image plotting function (include pause and save?)

        %% generate sum projection and binned projection, then threshold based on SNR parameter

        % sum up all time bins to generate normal 2D image (no time-binning in this step)
        data_t_sum = squeeze(sum(data,3));        
        % calculate  bin_xy image
        total_count_im = conv2(data_t_sum, ones(params.bin_size_xy, params.bin_size_xy), 'same');
        % plot and save binned, summed image
        figure; 
        imagesc(total_count_im);
        colormap('parula');
        colorbar;
        axis equal

        % set threshold based on count distribution
        params.max_count = max(total_count_im, [], 'all');
        params.threshold = params.max_count / params.thr_snr;
        params.mask = total_count_im > params.threshold;

        %% Calculate global lifetime for all pixels in mask region (no time binning)

        % Sum up all pixels in mask
        n_t_bins = size(data,3);
        data_xy_sum = zeros(1,1,n_t_bins);
        for i = 1 : n_t_bins
            dmy = data(:,:,i);
            data_xy_sum(1,1,i) = sum(dmy(params.mask),'all');
        end

        t = (0:size(data,3)-1)*params.dt; % convert bins to sec units

        % Fit IRF and data with monoexponential fit
        [r_fitirf, r_fitirf_fit, irf_fit] = fit_tcspc_gauss_irf_varpro(t, data_xy_sum, params);

    end


end