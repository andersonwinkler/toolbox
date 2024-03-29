function srf2ply(varargin)
% Convert Surface ASCII format to PLY.
% 
% Usage:
% srf2ply input.srf output.ply
%  
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Mar/2013
% http://brainder.org

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin   = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % For use later (PNG issue in current GraphicsMagick package)
    isoctave = true;

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Convert Surface ASCII format to PLY.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('srf2ply input.srf output.ply\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / University of Oxford\n');
        fprintf('Mar/2013\n');
        fprintf('http://brainder.org\n');
        return;
    end
end
if nargin ~= 2,
    error('Invalid number of arguments');
end

% Read surface
[vtx,fac] = srfread(varargin{1});

% Write PLY!
fid = fopen(varargin{2},'w');
fprintf(fid,'ply\n');
fprintf(fid,'format ascii 1.0\n');
fprintf(fid,'element vertex %s\n',num2str(size(vtx,1)));
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');
fprintf(fid,'element face %s\n',num2str(size(fac,1)));
fprintf(fid,'property list uchar int vertex_index\n');
fprintf(fid,'end_header\n');
fprintf(fid,'%g %g %g\n',vtx');
fprintf(fid,'%d %d %d %d\n',[3*ones(size(fac,1),1) fac-1]');
fclose(fid);
