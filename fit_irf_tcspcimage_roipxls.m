clearvars; close all

multi = 0; % set to 1 if looping through multiple condition/experimental directories (each with a set of acquisitional subdirectories)

if multi == 1
    folderP = uigetdir; foldparts = strsplit(folderP,filesep); parent_name = foldparts{end}; clear foldparts
    dirlist = dir(folderP); dirlist = dirlist([dirlist.isdir]); dirlist(1:2) = [];
    dir_n = size(dirlist,1); folderP = [folderP,filesep];

    else
        folderN = uigetdir; folderN = [folderN,filesep];
        dir_n = 1;
end