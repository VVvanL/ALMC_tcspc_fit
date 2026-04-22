% Analyze fit_results data

% set file parameters
fn_part = {'AO_','F_','F_Special_','FD_1_','FDF_1_','FM_1_','FM_2_','FMF_1_'};
n_files_per_cond = 10;
n_exp = 2;

% pool all the fit parameters for valid pixels
for i = 1:numel(fn_part)
    sample_name = fn_part{i}(1:end-1);
    % disp(sample_name)
    for k=1:n_exp
        fit_data.analysis.(sample_name).(['tau',num2str(k)]) = [];
        fit_data.analysis.(sample_name).(['amp',num2str(k)]) = [];
        fit_data.analysis.(sample_name).(['photons',num2str(k)]) = [];
    end
    for j = 1:n_files_per_cond
        % disp(j)
        name = strcat(fn_part{i},num2str(j));
        mask = fit_data.(name).mask;
        for k=1:n_exp
            dmy1 = fit_data.(name).taus(:,:,k);
            fit_data.analysis.(sample_name).(['tau',num2str(k)]) = [fit_data.analysis.(sample_name).(['tau',num2str(k)]); dmy1(mask)];
            dmy2 = fit_data.(name).amplitudes(:,:,k);
            fit_data.analysis.(sample_name).(['amp',num2str(k)]) = [fit_data.analysis.(sample_name).(['amp',num2str(k)]); dmy2(mask)];


            fit_data.analysis.(sample_name).(['photons',num2str(k)]) = [fit_data.analysis.(sample_name).(['photons',num2str(k)]); dmy1(mask).*dmy2(mask)];
        end
    end
end


%% set output folder
dir_out = uigetdir();

%%
edges = 0:0.01:3;
n1 = 'photons1';
n2 = 'photons2';
h = figure;
for i = 1:numel(fn_part)
    sample_name = fn_part{i}(1:end-1);
    [heights, edges] = histcounts(fit_data.analysis.(sample_name).(n1)./fit_data.analysis.(sample_name).(n2),edges);
    plot(edges(1:end-1), heights./max(heights), 'DisplayName',replace(sample_name,'_',' '),'LineWidth',2);
    hold on
end
hold off;
legend;
% title(['ratio ', n1,' to ', n2]);
title('normalized distributions of number of photons ratio of fast to slow decay component');
axis([min(edges) max(edges) 0 1.05]);
grid on;
saveas(h, strcat(dir_out,filesep,'ratio_fast_slow_distribution.png'));
close(h);


%% plot lifetime distributions

edges = 0.25:0.005:1.25;
n1 = ['tau1'];

h_fast = figure;
for i = 1:numel(fn_part)
    sample_name = fn_part{i}(1:end-1);
    [heights, edges] = histcounts(fit_data.analysis.(sample_name).(n1),edges);
    plot(edges(1:end-1), heights./max(heights), 'DisplayName',replace(sample_name,'_',' '),'LineWidth',2);
    hold on
end
hold off;
legend;
grid on;
axis([min(edges) max(edges) 0 1.05])
% title(['distributions of ', n1]);
title('normalized distributions of the fast lifetime components across all  conditions')
xlabel('lifetime (ns)');
ylabel('normalized lifetime distribution')
saveas(h_fast, strcat(dir_out,filesep,'lifetime_distribution_fast_component.png'));

% 
edges = 1.5:0.025:4.5;
n1 = ['tau2'];

h_slow = figure;
for i = 1:numel(fn_part)
    sample_name = fn_part{i}(1:end-1);
    [heights, edges] = histcounts(fit_data.analysis.(sample_name).(n1),edges);
    plot(edges(1:end-1), heights./max(heights), 'DisplayName',replace(sample_name,'_',' '),'LineWidth',2);
    hold on
end
hold off;
legend;
grid on;
axis([min(edges) max(edges) 0 1.05])
% title(['distributions of ', n1]);
title('normalized distributions of the slow lifetime components across all  conditions')
xlabel('lifetime (ns)');
ylabel('normalized lifetime distribution')
saveas(h_slow, strcat(dir_out,filesep,'lifetime_distribution_slow_component.png'));


%% plot both taus in the same plot
edges = 0:0.01:2;
n1 = 'tau1';
n2 = 'tau2';
figure;
ax = axes;
hold(ax,'on');
for i = 1:numel(fn_part)
    ax.ColorOrderIndex = i;
    sample_name = fn_part{i}(1:end-1);
    [heights, edges] = histcounts(fit_data.analysis.(sample_name).(n1),edges);
    plot(edges(1:end-1), heights./max(heights), 'DisplayName',replace([sample_name,' ',n1],'_',' '),'LineWidth',1.2,'LineStyle','-');
end
legend;
edges = 1:0.02:4.5;
for i = 1:numel(fn_part)
    ax.ColorOrderIndex = i;
    sample_name = fn_part{i}(1:end-1);
    [heights, edges] = histcounts(fit_data.analysis.(sample_name).(n2),edges);
    plot(edges(1:end-1), heights./max(heights),'DisplayName',replace([sample_name,' ',n2],'_',' '),'LineWidth',1.5,'LineStyle',':');
end
hold(ax,'off');
grid on


%% plot images ratio photons1/photons2


h_im = figure;
i = 5;
j = 4;

name = strcat(fn_part{i},num2str(j));
dmy = zeros(size(fit_data.(name).total));
mask = fit_data.(name).mask;
dmy1t = fit_data.(name).taus(:,:,1);
dmy1a = fit_data.(name).amplitudes(:,:,1);
dmy2t = fit_data.(name).taus(:,:,2);
dmy2a = fit_data.(name).amplitudes(:,:,2);

dmy(mask) = dmy1t(mask).*dmy1a(mask)./dmy2t(mask)./dmy2a(mask);

% imagesc(dmy);
s = replace(string(name),'_',' ');
% title(['Ratio photons in fast component to photons in slow component in ',replace(s,'_',' ')]);

rgbImage = plotColorBrightness(dmy, fit_data.(name).total,'parula', 0.5, s,'ratio of number of photons in fast to slow component' , 0.5, 'ScaleBar','br', 167, 50);

saveas(h_im, strcat(dir_out,filesep,name, '_fast_slow_photon_ratio_im.png'));
close(h_im);

%% plot image ratio for all images

for i = 1:numel(fn_part)
    for j = 1:n_files_per_cond
        name = strcat(fn_part{i},num2str(j));
        
        dmy = zeros(size(fit_data.(name).total));
        mask = fit_data.(name).mask;
        dmy1t = fit_data.(name).taus(:,:,1);
        dmy1a = fit_data.(name).amplitudes(:,:,1);
        dmy2t = fit_data.(name).taus(:,:,2);
        dmy2a = fit_data.(name).amplitudes(:,:,2);
        
        dmy(mask) = dmy1t(mask).*dmy1a(mask)./dmy2t(mask)./dmy2a(mask);
        
        s = replace(string(name),'_',' ');
        
        h_im = figure;
        rgbImage = plotColorBrightness(dmy, fit_data.(name).total_raw,'parula', 0.5, s,'ratio of number of photons in fast to slow component' , 0.5, 'ScaleBar','br', 167, 50);
        saveas(h_im, strcat(dir_out,filesep,name, '_fast_slow_photon_ratio_im.png'));
        close(h_im);
    end
end