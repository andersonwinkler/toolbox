function xyz2dpv(varargin)
% Takes the (X,Y,Z) coordinates from a surface file
% and save as a DPV file. Three DPV files are created,
% one for each dimension.
% 
% Usage 
% xyz2dpv(srffile,fprefix)
% 
% srffile : Surface file to extract coordinates
% fprefix : Prefix (can include full path) for the DPV
%           files to be created. Three files will be created,
%           named as 'fprefix.{X,Y,Z}.dpf'.
%           
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jul/2011

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin   = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),

        fprintf('Takes the (X,Y,Z) coordinates from a surface file\n');
        fprintf('and save as a DPV file. Three DPV files are created,\n');
        fprintf('one for each dimension.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('xyz2dpv <srffile> <fprefix>\n');
        fprintf('\n');
        fprintf('srffile : Surface file to extract coordinates\n');
        fprintf('fprefix : Prefix (can include full path) for the DPV\n');
        fprintf('          files to be created. Three files will be created,\n');
        fprintf('          named as ''fprefix.{X,Y,Z}.dpf''.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Jul/2011\n');

        return;
    end
end


% Letter suffix
xyz = {'X';'Y';'Z'};

% Read original surface
vtx = srfread(varargin{1});

% Vertex indices
idx = (0:size(vtx,1)-1)';

% Save for each dimension
for d = 1:3,
    fid = fopen(sprintf('%s.%s.dpv',varargin{2},xyz{d}),'w');
    fprintf(fid,'%0.3d %g %g %g %g\n',[idx vtx vtx(:,d)]');
    fclose(fid);
end
