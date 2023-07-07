function dpv2xyz(varargin)
% Takes three DPV files, one for X, one for Y and one for Z
% face indices from a surface file, and creates a new surface
% file. It does the opposite of the function xyz2dpv.
% 
% Usage 
% dpv2xyz(dpvfileX,dpvfileY,dpvfileZ,facfile,srffile)
% 
% dpvfile{X,Y,Z}: Data-per-vertex files
% facfile:        SRF file to extract the  (surface)
% srffile:        SRF file to save the results (surface)
%           
% _____________________________________
% Anderson M. Winkler
% FMRIB / Oxford University
% Jan/2012

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),

        fprintf('Takes three DPV files, one for X, one for Y and one for Z\n');
        fprintf('face indices from a surface file, and creates a new surface\n');
        fprintf('file. It does the opposite of the function xyz2dpv.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('dpv2xyz <dpvfileX> <dpvfileY> <dpvfileZ> <facfile> <srffile>\n');
        fprintf('\n');
        fprintf('dpvfile{X,Y,Z}: Data-per-vertex files\n');
        fprintf('facfile:        SRF file to extract the  (surface)\n');
        fprintf('srffile:        SRF file to save the results (surface)\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / Oxford University\n');
        fprintf('Jan/2012\n');
        fprintf('Jul/2011\n');
        
        return;
    end
end

% Read the DPF files containing the coordinates
x = dpxread(varargin{1});
y = dpxread(varargin{2});
z = dpxread(varargin{3});

% Assemble vertices
vtx = [x y z];

% Read face indices
[~,fac] = srfread(varargin{4});

% Save
srfwrite(vtx,fac,varargin{5});
