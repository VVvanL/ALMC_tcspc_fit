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

        % plot histogram of binned counts
        figure; histogram(total_count_im, 'Normalization', 'percentage')

        % set threshold based on count distribution
        params.max_count = max(total_count_im, [], 'all');
        params.threshold = params.max_count / params.thr_snr;
        params.mask = total_count_im > params.threshold;

        % plot mask image
        figure; imagesc(params.mask)
        axis equal 

        %% calculate bin_t / bin_xy image

        im_data_tbin = bin_array(data, params.bin_size_t, 3);
        n_layers = size(im_data_tbin, 3);
        im_data_tbin_xybin = im_data_tbin;
        for i = 1:n_layers
            im_data_tbin_xybin(:,:,i) = conv2(im_data_tbin(:,:,i), ones(params.bin_size_xy, params.bin_size_xy), 'same');
        end
        t_bin = (0:n_layers - 1) * (params.dt * params.bin_size_t); % convert bins to seconds

        % calculate mask TCSPC from summed xy and t binned image (global mask fit)
        mask_data_xy_sum = zeros(1,1,n_layers);
        for i = 1 : n_layers
            dmy = im_data_tbin_xybin(:,:,i);
            mask_data_xy_sum(1,1,i) = sum(dmy(params.mask),'all');
        end

        %% Fit summed TCSPC trace for IRF and lifetimes with biexponential fit
        params.x0 = [1,4]; params.lb = [0.1, 2]; params.ub = [5, 10]; % parameters specific for biexponential fit
        [r_fitirf, r_fitirf_fit, irf_fit] = ...
            fit_tcspc_gauss_irf_varpro(t_bin, mask_data_xy_sum, params);

        % plot decay with IRFs and fit
        figure;
        semilogy(t_bin,squeeze(mask_data_xy_sum),'.','DisplayName','data');
        hold on;
        semilogy(t_bin,r_fitirf_fit,'DisplayName','fit');
        semilogy(t_bin,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
        hold off
        grid on;
        legend;
        xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels in mask with biexponential fit');
        % ylim([min([min(r_fitirf_fit) min(mask_data_xy_sum)]) max(r_fitirf_fit)*1.05]);
        ylim([min(r_fitirf_fit)*0.9 max(r_fitirf_fit)*1.1])

        ant_str = ([char(964), '_1: ',num2str(r_fitirf.taus(1)),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1)),' ns', ...
            newline, char(964), '_2: ',num2str(r_fitirf.taus(2)),' ns', char(177),' ', num2str(r_fitirf.err_vals.taus(2)),' ns', ...
            newline, 'background: ',num2str(r_fitirf.background), char(177),' ', num2str(r_fitirf.err_vals.background)]);
        dim = [.66 .7, .1 .1];
        annotation('textbox', 'Position',dim, String = ant_str, FontSize = 13)
        % TODO: get amplitudes?

        %% Fit IRF and data with monoexponential fit
        params.x0 = 3; params.lb = 0.1; params.ub = 10;
        [r_fitirf, r_fitirf_fit, irf_fit] = fit_tcspc_gauss_irf_varpro(t_bin, mask_data_xy_sum, params);

        figure; 
        semilogy(t_bin,squeeze(mask_data_xy_sum),'.','DisplayName','data');
        hold on;
        semilogy(t_bin,r_fitirf_fit,'DisplayName','fit');
        semilogy(t_bin,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
        hold off
        grid on;
        legend;
        xlabel('time (ns)'); ylabel('counts');title('tcspc histogram of all pixels in mask with monoexponential fit');
        % ylim([min([min(r_fitirf_fit) min(mask_data_xy_sum)]) max(r_fitirf_fit)*1.05]);
        ylim([min(r_fitirf_fit) max(r_fitirf_fit)*1.05])
        % print fit parameters
        ant_str = ([char(964), '_1: ',num2str(r_fitirf.taus(1)),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1)),' ns', ...
            newline, 'background: ',num2str(r_fitirf.background), char(177),' ', num2str(r_fitirf.err_vals.background)]);
        dim = [.66 .7, .1 .1];
        annotation('textbox', 'Position',dim, String = ant_str, FontSize = 13)
        %

    end


end