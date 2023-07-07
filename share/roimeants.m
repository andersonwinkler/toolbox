#!/usr/local/apps/octave/4.4.0/bin/octave -q
% function roimeants(varargin)
%
% Usage:
% roimeants(mghfile,srffile,roifile,output)
%
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Feb/2012
% http://brainder.org

% Default: compute the mean. If false, compute the first PC.
domean = true;

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q')
        % to be written...
        return;
    end
end
if numel(varargin) ~= 4,
   error('Invalid number of arguments')
end

% Load the data
data    = load_mgh(varargin{1});
data    = squeeze(data);
[nX,nS] = size(data);

% Load and label file with ROIs
dpxlab = dpxlabelling(varargin{3},varargin{2});
if size(dpxlab,1) ~= nX
    error([
        'File with masks not with same size as file with data.\n' ...
        '- Size of data: %d\n' ...
        '- Size of mask: %d'], nX, size(dpxlab,1));
end

% Extract the mean and print as a table
U  = unique(dpxlab);
nU = numel(U);
T  = zeros(nS,nU);
if domean
    for u = 1:nU
        idx = dpxlab == U(u);
        T(:,u) = mean(data(idx,:),1)';
    end
else
    for u = 1:nU
        idx = dpxlab == U(u);
        T(:,u) = epca(data(idx,:)');
    end
end

% Save to the disk. First column are the vertices that don't belong to any
% ROI.
dlmwrite(varargin{4},T,'delimiter',',','precision','%0.4f');
