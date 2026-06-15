% Working script for batch, first-pass analysis of TCSPC dataset
clearvars; close all

multi = 0; % set to 1 if looping through multiple condition/experimental directories (each with a set of acquisitional subdirectories)
params = setTCSPC_fit_parameters();