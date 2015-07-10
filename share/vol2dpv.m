function vol2dpv(varargin)
% Project the values from a volume into a surface. Each vertex
% of the surface receives the value from the closest voxel from
% the volume. If more than one voxel share the same closest
% distance, the value assigned is the average of both.
% 
% Usage:
% vol2dpv(volfile,srffile,dpvfile,mskfile)
% 
% volfile : Volume file, with values to be projected into the
%           surface. It can be in MGZ or NIFTI format.
% srffile : Surface file, that will receive the values. It has
%           to be in ASCII SRF format.
% dpvfile : Data per vertex file, to be created.
% mskfile : Optional. A mask file. Voxels outside the mask are
%           ignored.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014
% http://brainder.org

% Do OCTAVE stuff
try
    % Get the inputs
    varargin = argv();

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Project the values from a volume into a surface. Each vertex\n');
        fprintf('of the surface receives the value from the closest voxel from\n');
        fprintf('the volume. If more than one voxel share the same closest\n');
        fprintf('distance, the value assigned is the average of both.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('vol2dpv volfile srffile dpvfile [mskfile]\n');
        fprintf('\n');
        fprintf('volfile : Volume file, with values to be projected into the\n');
        fprintf('          surface. It can be in MGZ or NIFTI format.\n');
        fprintf('srffile : Surface file, that will receive the values. It has\n');
        fprintf('          to be in ASCII SRF format.\n');
        fprintf('dpvfile : Data per vertex file, to be created.\n');
        fprintf('mskfile : Optional. A mask file. Voxels outside the mask are\n');
        fprintf('          ignored.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / University of Oxford\n');
        fprintf('Feb/2014\n');
        fprintf('http://brainder.org\n');
        return;
    end
end
nargin = numel(varargin);

% Check number of arguments
if nargin < 2 || nargin > 4,
    error('Incorrect number of arguments');
end

% Check FreeSurfer and add the 'matlab' dir to the path
fshome = getenv('FREESURFER_HOME');
if isempty(fshome),
    error([ ...
        'FreeSurfer not found. Make sure it is installed and that\n' ...
        'the variable FREESURFER_HOME is set correctly.']);
end
addpath(fullfile(fshome,'matlab'));

% Take the inputs
volfile = varargin{1};
srffile = varargin{2};
dpvfile = varargin{3};
if nargin == 4 && ~ isempty(varargin{4}),
    mskfile = varargin{4};
else
    mskfile = '';
end

% Load the volume file, that will to be projected
[~,~,fext] = fileparts(volfile);
switch fext,
    case '.mgz',
        [vol,M] = load_mgh(volfile);
    case {'.nii','.gz'},
        hdr = load_nifti(volfile);
        M   = hdr.vox2ras;
        vol = hdr.vol;
        clear('hdr');
end

% Load the mask, or convert the volume if none is supplied
if isempty(mskfile),
    msk = logical(vol);
else
    [~,~,fext] = fileparts(mskfile);
    switch fext,
        case '.mgz',
            [msk,Mm] = load_mgh(mskfile);
        case {'.nii','.gz'},
            hdr = load_nifti(mskfile);
            Mm  = hdr.vox2ras;
            msk = logical(hdr.vol);
            clear('hdr');
    end
    Mtest = Mm ~= M;
    siztest = size(msk) ~= size(vol);
    if any(Mtest(:)) || any(siztest(:)),
        error('Size and/or affine matrix of volume and mask don''t match.');
    end
end

% Load the surface (to be projected to)
[vtx,fac] = srfread(srffile);
nV  = size(vtx,1);

% Put the vertex coordinates into RAS coordinates
vtx = [vtx(:,[1 2 3]) ones(nV,1)];
ras = vtx*inv(M)';

% Make the volume a grid
[volx,voly,volz] = ndgrid(1:size(vol,1),1:size(vol,2),1:size(vol,3));

% Round the RAS coordinates
flo = floor(ras(:,1:3))-2;
cei = ceil(ras(:,1:3))+2;

% For each vertex
dpv = zeros(nV,1);
for v = 1:nV,
    
    % Select the voxels of interest
    vbox = vol(flo(v,1):cei(v,1),flo(v,2):cei(v,2),flo(v,3):cei(v,3));
    mbox = msk(flo(v,1):cei(v,1),flo(v,2):cei(v,2),flo(v,3):cei(v,3));
    xbox = volx(flo(v,1):cei(v,1),flo(v,2):cei(v,2),flo(v,3):cei(v,3));
    ybox = voly(flo(v,1):cei(v,1),flo(v,2):cei(v,2),flo(v,3):cei(v,3));
    zbox = volz(flo(v,1):cei(v,1),flo(v,2):cei(v,2),flo(v,3):cei(v,3));
    
    % Take their values and x,y,z coordinates 
    vv = vbox(mbox);
    xv = xbox(mbox);
    yv = ybox(mbox);
    zv = zbox(mbox);
    
    % Squared Euclidean distance
    dsq = (ras(v,1)-xv).^2 + (ras(v,2)-yv).^2 + (ras(v,3)-zv).^2;
    
    % Smallest distance. If more than one, average them
    dpv(v) = mean(vv(dsq == min(dsq)));
end

% Save to disk
dpxwrite(dpvfile,dpv);
