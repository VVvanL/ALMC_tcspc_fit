% working script to create image mask and then fit pixels within mask
clearvars; close all

%% parameters
im_ext = '*.obf'; % extensition for image files
dt = 0.025; % size of time bin in ns
bin_size_xy = 7;    % xy bin size in pixels, must be odd
bin_size_t = 4;     % bin size in time direction
cost_type = 'MLE';
fit_bg = true;
fit_shift = false;
error_type = '95CI';


%% select experiment directory and determine number of images
folderN = uigetdir; folderN = [folderN,filesep];