% organize individual image files from FLIM experiments and sort into subdirectories (one per image)
clearvars; close all

im_ext = '*.obf'; % extensition for image files
%% select experiment directory and determine number of images
folderN = uigetdir; folderN = [folderN,filesep];

im_list = dir([folderN, im_ext]);
im_n = size(im_list, 1);

%% loop through images, create subdirectory, move image
for im = 1:im_n
    im_str = im_list(im).name(1:end-4);

    im_dir = [folderN, im_str, filesep];
    if ~exist(im_dir, 'dir'); mkdir(im_dir); end

    movefile(([folderN, im_list(im).name]), ([im_dir, im_list(im).name]))

end