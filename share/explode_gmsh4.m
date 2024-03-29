function explode_gmsh4(varargin)
% "Explode" a Gmesh4 (.msh) file into
% .srf and .dpv files.
% 
% Usage:
% explode_gmsh input.msh prefix
%  
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Oct/2021
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
        fprintf('"Explode" a Gmesh4 (.msh) file into\n');
        fprintf('.srf and .dpv files.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('explode_gmsh input.msh prefix\n');
        fprintf(' \n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('National Institutes of Health\n');
        fprintf('Oct/2021\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Take input arguments
if nargin ~= 2,
    error('Invalid number of arguments');
end
gmsh4file = varargin{1};
prefix    = varargin{2};

% Read mesh file
gmsh4     = load_gmsh4(gmsh4file);

% Write surface
vtx       = gmsh4.nodes;
fac       = gmsh4.triangles;
srffile   = sprintf('%s.srf',prefix);
srfwrite(vtx,fac,srffile)

% Write "curvature" files
for nd = 1:numel(gmsh4.node_data)
    dpvfile = sprintf('%s_%s.dpv',prefix,gmsh4.node_data{nd}.name);
    dpv     = gmsh4.node_data{nd}.data;
    dpxwrite(dpvfile,dpv);
end
