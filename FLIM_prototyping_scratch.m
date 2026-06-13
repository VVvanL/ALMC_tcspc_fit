% scratch code to customize FLIM fitting and tau mapping
clearvars; close all

%% select file and load
[file,path] = uigetfile('*.*','select TCSPC file');
% get get dimensions of dataseries in file
[warnings,XYZTC] = evalc('bf_file_info(strcat(path,file))');
dataseries = 1;
[warnings, raw_data] = evalc('bf_load_parts_v7(strcat(path,file),dataseries,-1,-1,-1,-1,-1)'); % use evalc to block annoying bioformats warnings
data = squeeze(raw_data);
disp('Done')


%% Mask generateion
% sum up all time bins to generate normal 2D image
data_t_sum = squeeze(sum(data,3));