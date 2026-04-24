% add path to bioformats
addpath('bfmatlab');
javaaddpath('bfmatlab\bioformats_package.jar');

[file,path] = uigetfile;

data = bf_load_parts_v7(strcat(path,file),1,-1,-1,-1,-1,-1);