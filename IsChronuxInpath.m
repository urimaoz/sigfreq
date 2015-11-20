function IsChronuxInpath (chronuxPath)
% function IsChronuxInpath (chronux_path)
% 
% Checks whether the Chronux package is in the path, and if not adds it and
% all its subdirectories
if (isempty(strfind(path,'chronux')))
    if (~exist('chronuxPath','var')||isempty(chronuxPath))
%         chronux_path=fullfile(matlabroot,'chronux');
        chronuxPath='/Users/urimaoz/Documents/MATLAB/chronux/';
    end
    addpath(genpath(chronuxPath));
end
