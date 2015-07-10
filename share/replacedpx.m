%#!/usr/bin/octave -q
function replacedpx(varargin)
% Replace values in a DPV/DPF file. The correspondence between
% old and new values is provided by a CSV table.
% 
% Usage:
% replacedpx(olddpx,table,newdpx)
% 
% olddpx : Original DPV/DPF file.
% table  : A CSV file containing 2 columns. The first contain the old values
%          to be replaced, whereas the second contains the new values. All
%          instances of the old value in the DPV/DPF file will be replaced
%          by the corresponding new value. The remaining of the original
%          DPV/DPF file remains intact.
% newdpx : New file to be created.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Sep/2011
% http://brainder.org

% Do the OCTAVE stuff, using TRY to ensure MATLAB compatibility
try
    % Get the inputs
    varargin = argv();

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('\n');
        fprintf('Replace values in a DPV/DPF file. The correspondence between\n');
        fprintf('old and new values is provided by a CSV table.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('replacedpx <olddpx> <table> <newdpx>\n');
        fprintf('\n');
        fprintf('olddpx : Orignial DPV/DPF file.\n');
        fprintf('table  : A CSV file containing 2 columns. The first contain the old values\n');
        fprintf('         to be replaced, whereas the second contains the new values. All\n');
        fprintf('         instances of the old value in the DPV/DPF file will be replaced\n');
        fprintf('         by the corresponding new value. The remaining of the original\n');
        fprintf('         DPV/DPF file remains intact.\n');
        fprintf('newdpx : New file to be created.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Sep/2011\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Accept inputs
olddpx = varargin{1};
table  = varargin{2};
newdpx = varargin{3};

% Read table and the original DPX file
T   = dlmread(table,',');
[dpx,crd,idx] = dpxread(olddpx);

% Do the replacement
for t = 1:size(T,1),
    dpx(dpx == T(t,1)) = T(t,2);
end

% Save the result
dpxwrite(newdpx,dpx,crd,idx);

