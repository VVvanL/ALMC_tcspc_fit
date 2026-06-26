function create_subdirectories(folderN, im_ext)

im_list = dir([folderN, im_ext]);
im_n = size(im_list, 1);

    %% loop through images, create subdirectory, move image
for im = 1:im_n
    im_str = im_list(im).name(1:end-4);

    im_dir = [folderN, im_str, filesep];
    if ~exist(im_dir, 'dir'); mkdir(im_dir); end

    movefile(([folderN, im_list(im).name]), ([im_dir, im_list(im).name]))

end





