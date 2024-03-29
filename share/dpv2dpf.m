function dpv2dpf(varargin)
% Convert data-per-vertex(DPV) to data-per-face (DPF) files, redistributing the
% face quantities to their vertices. Assumes that the quantity is
% homogeneously distributed within face and that the redistribution is conceptually
% correct.
%
% Usage:
% dpv2dpf(srffile,dpvfile,dpffile,method)
% 
% srffile : Input reference surface file, in ASCII format.
% dpvfile : File with the data-per-face. Has to have the
%           same number of faces as the reference surface file.
% dpffile : File to be created, with the quantities redistributed
%           to the vertices.
% method  : How the DPV should be converted to DPF. Enter a number
%           between 1-4 here.
%           There are at least four different ways to convert:
%           (1) The face value is the MEAN of its vertices
%           (2) The face value is the MODE of its vertices
%           (3) The face value is the ROUNDED mean of its vertices
%           (4) The face value is the FLOOR of the mean of its vertices
%           (5) The face value is the CEIL of the mean of its vertices
%           (6) The face value is the MINIMUM of its vertices
%           (7) The face value is the MAXIMUM of its vertices
%           (8) The face value is a fraction taken from its vertices (distributive)
% The 8th isn't implemented... it requires loop over faces. Maybe later...
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011
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
        fprintf('Convert data-per-vertex (DPV) to data-per-face (DPF) files, redistributing the\n');
        fprintf('face quantities to their vertices. Assumes that the quantity is\n');
        fprintf('homogeneously distributed within face and that the redistribution is conceptually\n');
        fprintf('correct.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('dpv2dpf(srffile,dpvfile,dpffile)\n');
        fprintf('\n');
        fprintf('srffile : Input reference surface file, in ASCII format.\n');
        fprintf('dpvfile : File with the data-per-face. Has to have the\n');
        fprintf('          same number of faces as the reference surface file.\n');
        fprintf('dpffile : File to be created, with the quantities redistributed\n');
        fprintf('          to the vertices.\n');
        fprintf('method  : How the DPV should be converted to DPF. Enter a number\n');
        fprintf('          between 1-4 here.\n');
        fprintf('          There are at least four different ways to convert:\n');
        fprintf('          (1) The face value is the MEAN of its vertices\n');
        fprintf('          (2) The face value is the MODE of its vertices\n');
        fprintf('          (3) The face value is the ROUNDED mean of its vertices\n');
        fprintf('          (4) The face value is the FLOOR of the mean of its vertices\n');
        fprintf('          (5) The face value is the CEIL of the mean of its vertices\n');
        fprintf('          (6) The face value is the MINIMUM of its vertices\n');
        fprintf('          (7) The face value is the MAXIMUM of its vertices\n');
        fprintf('          (8) The face value is a fraction taken from its vertices (distributive)\n');
        fprintf('The 8th isn''t implemented... it requires loop over faces. Maybe later...\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Aug/2011\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Get the inputs (varargin has to be used so that it works in OCTAVE)
srffile = varargin{1};
dpvfile = varargin{2};
dpffile = varargin{3};
method  = str2double(varargin{4});

% Read the reference surface
[vtx,fac] = srfread(srffile);
nF = size(fac,1);

% Read the data per vertex file
dpv = crvread(dpvfile);

% Get the values at each vertex, for each face
dpvfac = reshape(dpv(fac(:)),size(fac));

% Convert using the appropriate way
if     method == 1,
    dpf = mean(dpvfac,2);
elseif method == 2,
    dpf = mode(dpvfac,2);
elseif method == 3,
    dpf = round(mean(dpvfac,2));
elseif method == 4,
    dpf = floor(mean(dpvfac,2));
elseif method == 5,
    dpf = ceil(mean(dpvfac,2));
elseif method == 6,
    dpf = min(dpvfac,[],2);
elseif method == 7,
    dpf = max(dpvfac,[],2);
elseif method == 8,
    error('Not implemented yet');
end

% Save the new DPF
crvwrite(dpffile,dpf,fac,(0:nF-1));
