function dpxthr(varargin)
% Threshold a DPX file (DPV/DPF) according to a relative
% value (percentile). The output is a binary (0-1) file
% that can be used with rpncalc for other calculations.
% 
% Usage:
% dpxthr(dpxfile,thr,outfile)
% 
% dpxfile : File to be thresholded, DPV or DPF
% thr     : Percentile threshold (0-100). Only the
%           top percentile points are labelled as 1.
% outfile : DPX file to be created (0-1)
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2012

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Threshold a DPX file (DPV/DPF) according to a relative\n');
        fprintf('value (percentile). The output is a binary (0-1) file\n');
        fprintf('that can be used with rpncalc for other calculations.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('dpxthr <dpxfile> <thr> <outfile>\n');
        fprintf('\n');
        fprintf('dpxfile : File to be thresholded, DPV or DPF\n');
        fprintf('thr     : Percentile threshold (0-100). Only the\n');
        fprintf('          top percentile points are labelled as 1.\n');
        fprintf('outfile : DPX file to be created (0-1)\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / University of Oxford\n');
        fprintf('Feb/2012\n');
        return;
    end
end

% Accept arguments
if nargin ~= 3,
    error('Invalid number of arguments')
end

% Read input files
[dpx,crd,idx] = dpxread(varargin{1});
thr = str2double(varargin{2});

% Remove NaN and sort values in ascending order
dpxr = dpx;
dpxr(isnan(dpx)) = [];
nXr = numel(dpxr);
dpxs = sort(dpxr);

% Upper index
tidx = round(nXr*(1-thr));

% Binarize and add the NaN back
dpxthr = zeros(size(dpx));
dpxthr(dpx >= dpxs(tidx)) = 1;
dpxthr(isnan(dpx)) = NaN; 

% Save
dpxwrite(varargin{3},dpxthr,crd,idx);
