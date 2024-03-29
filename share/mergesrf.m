function mergesrf(varargin)
% Merge multiple surface files in *.srf format
% into a single file.
%
% Usage:
% mergesrf('file1.srf','file2.srf',...,'mergedfile.srf')
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jul/2012
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
        
        fprintf('Merge multiple surface files in *.srf format\n');
        fprintf('into a single file.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('mergesrf file1.srf file2.srf [...] mergedfile.srf\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / University of Oxford\n');
        fprintf('Jul/2012\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

if nargin < 2,
    error('Error: insufficient number of arguments.\n')
end

nVnew = 0;
nFnew = 0;
vtx = cell(nargin-1,1);
fac = vtx;
for s = 1:(nargin-1),
    [vtx{s},fac{s}] = srfread(varargin{s});
    fac{s} = fac{s} + nVnew;
    nVnew = nVnew + size(vtx{s},1);
    nFnew = nFnew + size(fac{s},1);
end
vtxnew = vertcat(vtx{:});
facnew = vertcat(fac{:});
srfwrite(vtxnew,facnew,varargin{nargin});
