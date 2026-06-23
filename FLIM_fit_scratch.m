% scratch code for prototyping TCSPC image fitting


%% Mask generation
% sum up all time bins to generate normal 2D image
data_t_sum = squeeze(sum(data,3));