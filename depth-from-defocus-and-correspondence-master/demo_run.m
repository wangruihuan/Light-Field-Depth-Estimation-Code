% This gives a demo run of compute_LFdepth
clear all;
% two files that need to be mex
cd('required\mex\');
mex FAST_STDFILT_mex.c
mex REMAP2REFOCUS_mex.c
cd('..');
cd('..');
% file path
file_path     =  'raw2jpeg/input/dataset/IMG_5.jpg' ;
% the main function
depth_output  = compute_LFdepth(file_path);