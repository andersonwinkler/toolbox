function midsrf(srf1,srf2,srfout)
% Given two surface files, compute the mid-surface
% between them. To be valid, both surfaces must have 
% the same number of vertices and the vertices must
% have identical connectedness properties.
%
% midsrf(srf1,srf2,srfout)
% 
% - srf1:   filename for surface 1
% - srf2:   filename for surface 2
% - srfout: filename for surface to be saved (mid surface)
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Feb/2011

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);;

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Given two ASCII surface files, compute the mid-surface\n');
        fprintf('between them. To be valid, both surfaces must have\n');
        fprintf('the same number of vertices and the vertices must\n');
        fprintf('have identical connectedness properties.\n');
        fprintf('\n');
        fprintf('midsrf <srf1.srf> <srf2.srf> <srfout.srf> \n');
        fprintf('\n');
        fprintf('- srf1:   filename for surface 1\n');
        fprintf('- srf2:   filename for surface 2\n');
        fprintf('- srfout: filename for surface to be saved (mid surface)\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Feb/2011\n');
        return;
        
    else
        srf1   = varargin{1};
        srf2   = varargin{2};
        srfout = varargin{3};
    end
end

% Read the surface files
[vtx1,fac1] = srfread(srf1);
[vtx2,fac2] = srfread(srf2);

% Do some sanity checks
nV1 = size(vtx1,1);
nV2 = size(vtx2,1);
tmp = fac1 == fac2;
if nV1 ~= nV2 || sum(tmp(:)) ~= numel(tmp),
    error('Surfaces have different structures')
end

% Compute mid surface coords
vtxout = (vtx1+vtx2)/2;

% Save
srfwrite(vtxout,fac1,srfout);
