function mergemgh(varargin)
% Merge across space multiple MGH/MGZ files into a single file.
% To instead merge across time/subjects, use mri_concat.
%
% Usage:
% mergesrf('file1.mgh','file2.mgh',...,'mergedfile.mgh')
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jun/2016
% http://brainder.org

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        
        fprintf('Merge across space multiple MGH/MGZ files into a single file.\n');
        fprintf('To instead merge across time/subjects, use mri_concat.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('mergesrf file1.mgh file2.mgh ... mergedfile.mgh \n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / University of Oxford\n');
        fprintf('Jun/2016\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

if nargin < 2,
    error('Error: insufficient number of arguments.\n')
end

mdata = [];
for m = 1:(nargin-1),
    [data,M,mr_parms] = load_mgh(varargin{m});
    mdata = cat(1,mdata,data);
end
save_mgh(mdata,varargin{nargin},M,mr_parms);


