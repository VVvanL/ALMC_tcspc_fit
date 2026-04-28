clearvars; close all

%% set initial parameters
dt = 0.025; % size of time bin in ns

multi = 0; % set to 1 if looping through multiple condition/experimental directories (each with a set of acquisitional subdirectories)

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
        im_list = dir([subpath,'*.tif']); im_n = length(im_list);
        for im = 1:im_n
            im_data = tiffreadVolume([subpath, im_list(im).name]);

        % generate time axis t
        figure;
        t = (0:size(roi_data,3)-1)*dt;
        % sum up all the pixels in the image
        tcspc_trace = squeeze(sum(roi_data,[1,2]));
        % plot in semilogarithmic plot
        semilogy(t,tcspc_trace);
        grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of all pixels (ROI)')

n_t_bins = size(im_data,3);
data_xy_sum = zeros(1,1,n_t_bins);
for i = 1 : n_t_bins
    dmy = im_data(:,:,i);
    data_xy_sum(1,1,i) = sum(dmy(mask),'all');
end

% plot in semilogarithmic plot
figure;
semilogy(t,squeeze(data_xy_sum));
grid on; xlabel('time (ns)'); ylabel('counts'); title('tcspc histogram of pixels in mask')


        end % image loop



    end % sub loop



end % dir loop

%% ===================================scratch Code======================
[file,path] = uigetfile('*.*','select TCSPC file');
% get get dimensions of dataseries in file
[warnings,XYZTC] = evalc('bf_file_info(strcat(path,file))');
dataseries = 1;
[warnings, raw_data] = evalc('bf_load_parts_v7(strcat(path,file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
data = squeeze(raw_data);
disp('Done')


[roi_file, roi_path] = uigetfile('*.*','select TCSPC ROI');
roi_list = ReadImageJROI([roi_path, roi_file]);
corners = roi_list{1}.vnRectBounds; % ['nTop', 'nLeft', 'nBottom', 'nRight']

roi_data = data(corners(1):corners(3), corners(2):corners(4),:);
tcspc_data = squeeze(sum(data,[1,2]));

%% Fit IRF and data with monoexponential fit
x0 = [1]; lb = [0.1]; ub = [10];
x0_irf = [1.2, 0.15, 1.2, 0.3, 0.05, 13.5, 0.2, 0.02, 13.5, 0.3, 0.02];
lb_irf = [0.5, 0.1, 0.5, 0.1, 0.005, 12.8, 0.2, 0.002, 12.7, 0.2, 0.002];
ub_irf = [1.9, 1.5, 1.9, 3.0, 0.500, 14.2, 2.0, 0.20, 14.2, 3.0, 0.20];
cost_type = 'MLE';
fit_bg = true;
error_type = '95CI';

[r_fitirf, r_fitirf_fit, irf_fit] = ...
    fit_tcspc_gauss_irf_varpro(t, tcspc_data, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg, error_type);

figure; 
semilogy(t,squeeze(tcspc_data),'.','DisplayName','data');
hold on;
semilogy(t,r_fitirf_fit,'DisplayName','fit');
semilogy(t,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
hold off
grid on;
legend;
xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels in mask with monoexponential fit');
ylim([min([min(r_fitirf_fit) min(tcspc_data)]) max(r_fitirf_fit)*1.05]);

% print fit parameters
disp([char(964), '_1: ',num2str(r_fitirf.taus(1)),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1)),' ns', ...
    newline, 'background: ',num2str(r_fitirf.background), char(177),' ', num2str(r_fitirf.err_vals.background)]);