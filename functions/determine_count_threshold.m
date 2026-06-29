function [params, h_counthist, h_mask] = determine_count_threshold(im_data, params)

% set threshold based on count distribution
params.max_count = max(im_data, [], 'all');
params.threshold = params.max_count / params.thr_snr;
params.mask = im_data > params.threshold;

% plot histogram of binned counts
h_counthist = figure; 
histogram(im_data)
xline(params.threshold, '-', 'Count threshold', 'LineWidth', 1.5)

% plot mask image
h_mask = figure;
imagesc(params.mask); colormap gray
axis equal 

