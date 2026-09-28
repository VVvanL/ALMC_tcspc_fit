% working script to aggregate data from multiple FLIM experiments (in data frame, structure format)
% currently only aggregating mono-exponential tau data
clearvars; close all

%% select parent directory and get list of experimental folders
folderP = uigetdir; foldparts = strsplit(folderP,filesep); parent_name = foldparts{end}; clear foldparts
dirlist = dir(folderP); dirlist = dirlist([dirlist.isdir]); dirlist(1:2) = [];
dir_n = size(dirlist,1); folderP = [folderP,filesep];

data_structure = struct();
data_structure.condition = {};
data_structure.tau = [];

for dr = 1:dir_n
    folderN = [folderP,filesep,dirlist(dr).name,filesep];
    foldparts = strsplit(folderN,filesep); dirname = foldparts{end-1}; clear foldparts

    disp(['Aggregating data from ',dirname,'...'])
    load([folderN, dirname, '_fitdata.mat'], 'tau_data')

    data_structure.condition = vertcat(data_structure.condition, tau_data.condition);
    data_structure.tau = vertcat(data_structure.tau, tau_data.tau);

end


h_box = figure; 
boxchart(data_structure.condition, data_structure.tau)
% savefig(h_box, [folderN, dirname,'_boxchart.fig'])

clear tau_data
tau_data = data_structure;
save([folderP, parent_name, '_tauData.mat'], 'data_structure')

T = struct2table(tau_data);
writetable(T, [folderP, parent_name, '_tau_data.csv'])