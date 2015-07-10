function crvwrite(varargin)
% Write a curvature file (DPV or DPF), in ASCII format.
% This function is much faster than 'dlmread' for large files,
% and works only in Linux and Mac.
%
% crvwrite(filename,crv)
% crvwrite(filename,crv,crd,idx)
%
% - fname is the file name to be created
% - crv contains the values for each vertex or face
% - crd contains the vertex coordinates or face indices
% - idx contains vertex or face sequential index
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011

% File name
fname = varargin{1};

% Get the actual data
crv = varargin{2}(:);
nX  = numel(crv);

% Check if all are integers and use appropriate formating
if all(mod(crv,1)==0),
    fstr = '%d';
else
    fstr = '%f';
end

if nargin == 2,

    % Organise the data, fill the coords with zeros amd prep to save
    crv = [(0:nX-1) ; zeros(3,nX) ; crv'];

elseif nargin == 4,

    % Organise the coords
    crd = varargin{3};
    if size(crd,1) > size(crd,2),
        crd = crd';
    end

    % Take the indices
    idx = varargin{4}(:);

    % Prepare to save
    crv = [idx' ; crd ; crv'];

else
    error('Incorrect number of arguments');
end

% Save
fid = fopen(fname,'w');
fprintf(fid,['%d %g %g %g ' fstr ' \n'],crv);
fclose(fid);
