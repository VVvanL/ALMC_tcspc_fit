%% fit all TCSPC data from experiment

%% Set all parameters for the analysis pipelin
% set fit parameters

% select fit option
% 1: uses autofit_tcspc_image
% 2: uses fit_irf_and_tcspc_image
fit_option = 2; 

% parameters necessary for both options
dt = 0.025; % size of time bin in ns
bin_size_xy = 5;    % xy bin size in pixels, must be odd
bin_size_t = 4;     % bin size in time direction
n_exp = 2;
cost_type = 'MLE';
fit_bg = true;
fit_shift = false;
error_type = '95CI';
threshold = 10000;  % minimum number of photon in TCSPC trace of a pixel so that the pixel will be fitted
show_irf_estimate = true;

% parameters that are only necesary for option 2

% parameters
x0 = [1 4 7]; % initial life time guess for fit of overall counts that also fits the irf
x0_irf = [1.3, 0.15, 1.3, 0.3, 0.05, 13.5, 0.2, 0.02, 13.5, 0.3, 0.01]; % fit parameters for the irf
n_exp_im_fit = 2; % number of exponents fitted for image


% set file parameters
fn_part = {'AO_','F_','F_Special_','FD_1_','FDF_1_','FM_1_','FM_2_','FMF_1_'};
% fn_part = {'AO_','F_'};

n_files_per_cond = 10;
fe = '.obf';
output_folder = 'analyzed';

%% set up processing
% turn on parallel pool if available
if canUseParallelPool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
end

% add path to bioformats
addpath('bfmatlab');
javaaddpath('bfmatlab\bioformats_package.jar');

%% select folder with TCSPC data
path = uigetdir(pwd,'select path with TCSPC data');

path_out = strcat(path,filesep,output_folder);

fit_data = struct;

%% run pipeline
for i=1:numel(fn_part)
    for j=1:n_files_per_cond
        fn_load = strcat(path,filesep,fn_part{i},num2str(j),fe);
        
        % load data and suppress output from bioformats 
        disp(['loading   ',fn_load,newline]);
        [T,data,~,~] = evalc('bf_load_parts_v7(fn_load,1,-1,-1,-1,-1,-1)');
        data = squeeze(data);

        % run fit
        switch fit_option
            case 1
                r_im = autofit_tcspc_image(data, dt, bin_size_xy, bin_size_t, threshold, n_exp, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate);
            case 2
                r_im = fit_irf_and_tcspc_image(data, dt, x0, x0_irf, bin_size_xy,  bin_size_t, n_exp_im_fit, threshold, fit_bg, fit_shift, cost_type, error_type, show_irf_estimate);
        end

        % save data as ome.tiff
        fn_save_tiff =  strcat(path,filesep,output_folder,filesep,'fit_results_',fn_part{i},num2str(j),'.ome.tiff');
        save_fit_data_as_tiff(r_im, data, fn_save_tiff);
        disp(['Saved fit results as ',fn_save_tiff, newline])
        fit_data.(strcat(fn_part{i},num2str(j))) = r_im;
        fit_data.(strcat(fn_part{i},num2str(j))).total_raw = squeeze(sum(data,3)); 
        fit_data.(strcat(fn_part{i},num2str(j))).filepath_in = fn_load;
        fit_data.(strcat(fn_part{i},num2str(j))).filepath_out = fn_save_tiff;

    end
end

%% Save all fits in one matlab file

fit_data.fit_param.dt = dt;
fit_data.fit_param.bin_size_xy = bin_size_xy;
fit_data.fit_param.bin_size_t = bin_size_t;
fit_data.fit_param.n_exp = n_exp;
fit_data.fit_param.cost_type = cost_type;
fit_data.fit_param.fit_bg = fit_bg;
fit_data.fit_param.fit_shift = fit_shift;
fit_data.fit_param.error_type = error_type;
fit_data.fit_param.threshold = threshold;
fit_data.fit_param.show_irf_estimate = show_irf_estimate;
fit_data.fit_param.fit_option = fit_option;
if fit_option == 2
    fit_data.fit_param.x0 = x0;
    fit_data.fit_param.x0_irf = x0_irf;
    fit_data.fit_param.n_exp_im_fit = n_exp_im_fit;    
end
fn_save_mat = strcat(path,filesep,output_folder,filesep,'fit_results.mat');
save(string(fn_save_mat),'fit_data','-v7.3');

%% save all the figures

for i = 1:numel(fn_part)*n_files_per_cond
    if ishandle(i) % Check if figure i is actually open
        filename = strcat(path,filesep,output_folder,filesep,'Figure_',num2str(i),'.jpg');
        saveas(i, filename);
    end
end