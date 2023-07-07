function applymat2srf(varargin)
% Apply a 4x4 affine matrix to a surface file.
% 
% Usage:
% applyxfm2srf(srfin,xfm,srfout)
% 
% - srfin  : Input surface file
% - xfm    : Text file with the affine matrix.
%            It can be CSV or separated by spaces or tabs.
% - srfout : Surface with the coordinates changed.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014
% http://brainder.org

% Do some OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Apply a 4x4 affine matrix to a surface file.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('applyxfm2srf <srfin> <xfm> <srfout>\n');
        fprintf('\n');
        fprintf('- srfin  : Input surface file.\n');
        fprintf('- xfm    : Text file with the affine matrix.\n');
        fprintf('           It can be CSV or separated by spaces or tabs.\n');
        fprintf('- srfout : Surface with the coordinates changed.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / University of Oxford\n');
        fprintf('Feb/2014\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Take arguments & load data
if nargin ~= 3,
    error('Incorrect number of arguments.');
end
[vtx,fac] = srfread(varargin{1});
xfm = load(varargin{2},'-ascii');

% Apply the transform
vtx = [vtx ones(size(vtx,1),1)];
vtx = vtx*xfm';

% Save
srfwrite(vtx(:,1:3),fac,varargin{3});
