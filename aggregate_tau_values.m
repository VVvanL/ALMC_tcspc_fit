% working script to aggregate lifetime data from analysis of individual acquisitions


multi = 0; % set to 1 if looping through multiple condition/experimental directories (each with a set of acquisitional subdirectories)
conditions = {'CaNwt_ST', 'CaNh92r_ST', 'CaNf470l_ST', 'CaN_d188n_ST'};
cnd_n = length(conditions);
fit_type = 1; % mono or bi-exponentional (place holder, will improve later)

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

    % preallocate arrays for data
    condition = cell(sub_n, 1);
    tau = nan(sub_n, fit_type);

    for s = 1:sub_n
        subname = sublist(s).name; subpath = fullfile(sublist(s).folder,subname,filesep);
        disp(['aggregating data from ', subname, '...'])

        condition{s} = subname(4:end);

        data_file = [subpath, subname, '_fitdata.mat'];
        load(data_file)

        tau(s) = fit_data.fitirf.taus;

    end
    data_table = table(condition, tau);
    save([folderN, dirname, '_data_table.mat'], 'data_table')


end

%% =======================scratch ploting code=====================================
boxplot(data_table.tau,data_table.condition)
ylim([3.5 4])

figure; grid on
cnd = categorical(data_table.condition);
scatter(cnd, data_table.tau, 'filled')

clear g

figure; 
g = gramm('x',  data_table.condition, 'y', data_table.tau);
g.geom_point('dodge',0.3);
g.set_names('x','','y','Lifetime (ns)','color','Sensor');
g.draw();


websave('example_data','https://github.com/piermorel/gramm/raw/master/sample_data/example_data.mat'); %Download data from repository
load example_data;
cars

g=gramm('x',cars.Model_Year,'y',cars.MPG);