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

        % set up structure for fit data
        fit_data = struct();

        % find TCSPC image file in sub-directory and load
        im_file = [subpath, subname, '.obf'];
        [~ , raw_data] = ...
            evalc('bf_load_parts_v7(strcat(im_file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
        data = squeeze(raw_data);

        %% generate sum projection and binned projection, then threshold based on SNR parameter

        % sum up all time bins to generate normal 2D image
        data_t_sum = squeeze(sum(data,3));        
        % calculate  bin_xy image
        total_count_im = conv2(data_t_sum, ones(params.bin_size_xy, params.bin_size_xy), 'same');

        % plot and save binned, summed image


    end


end