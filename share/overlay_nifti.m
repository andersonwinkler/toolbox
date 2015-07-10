function overlay_nifti(baseimg,baserange,overimg,overange,cmax,slices,dim,msize)
% Create overlays to show statistical results. The base image and the
% overlays must all be of the same size. Coregister them manually
% first if needed.
%
% - baseimg:   NIFTI image file to be used as reference
% - baserange: Intensity range for the base image. Leave empty [] for auto
% - overimg:   NIFTI image(s) to be overlayed. Wildcards accepted
% - overange:  Intensity range for the overlays (usually [0.95 1])
% - cmax:      RGB triplet for the max intensity voxel in the overlay. The
%              other voxels will be the same color, but dimmed.
% - slices:    Vector containing the list of slices to be shown in the mosaic
% - dim:       Dimension along which the slices are cut
% - msize:     Number of slices in the mosaic along vertical and horizontal
%
% Example usage:
% overlay_nifti('MNI152_T1_2mm_brain.nii',[], ...
%               'randomise*.nii.gz',[.95 1], ...
%                [1 1 0],[8:86],3,[8 9])
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Dec/2010

% Load base image and rescale to new range
[bpth,bnam,bext] = fileparts(baseimg);
gz = false;
if strcmp(bext,'.gz'),
    gz = true;
    gunzip(baseimg,bpth);
    baseimg = fullfile(bpth,bnam);
end
B = spm_read_vols(spm_vol(baseimg));
if gz, delete(baseimg); end
Br = rerange(B,baserange,[0 .5]);

% Get the list of overlay images
Olist = dir(overimg);
[Opth,~,~] = fileparts(overimg);
for o = numel(Olist):-1:1,
    if strcmp(Olist(o).name,'.') || strcmp(Olist(o).name,'..'),
        Olist(o) = []; % remove the '.' and '..'
    end
end

% Loop over overlays
for o = 1:numel(Olist),
    fprintf('Working on %g/%g [%s]\n',o,numel(Olist),Olist(o).name)
    
    % Load overlay and rescale to new range
    [~,onam,oext] = fileparts(Olist(o).name);
    gz = false;
    if strcmp(oext,'.gz'),
        gz = true;
        gunzip(fullfile(Opth,Olist(o).name),Opth);
        Olist(o).name = onam;
        [~,onam,~] = fileparts(Olist(o).name);
    end
    O = spm_read_vols(spm_vol(fullfile(Opth,Olist(o).name)));
    if gz, delete(fullfile(Opth,Olist(o).name)); end
    Or = rerange(O,overange,[0 .5]);
    
    % Sum to produce the RGBs to be saved later
    Icolor = cell(3,1);
    for c = 1:numel(Icolor),
        Icolor{c} = Br + (Or * cmax(c));
    end
    
    % Create the mosaics and save
    M = createmosaic(Icolor,slices,dim,msize);
    imwrite(M,fullfile(Opth,sprintf('%s_%g.png',onam,dim)));
end

function Vr = rerange(V,cr,nr)
% Changes the scaling (range) of V.
% cr: desired range
% nr: new range

if isempty(cr)
    d = .05; % fraction to be discarded
    Vvec = sort(V(V>0));
    v = numel(Vvec);
    cr = [Vvec(round(d*v)) Vvec(round((1-d)*v))];
end
V(V < cr(1)) = cr(1);
V(V > cr(2)) = cr(2);
Vr = nr(1) + (V-cr(1))*(nr(2)-nr(1))/(cr(2)-cr(1));

function M = createmosaic(Icolor,slices,dim,msize)
% Creates a mosaic
% Icolor: 3x1 cell, each containing a 3D volume for RGB
%         4x1 cell, each containing a 3D volume for CMYK
% slices: vector of the slices numbers to be used
% dim: dimension along the cuts
% msize: number of panels

% Permute the dimensions so that the slices are always across the 3rd
dimorder = 1:3;
dimorder = [fliplr(dimorder(dimorder ~= dim)) dim];
for c = 1:numel(Icolor),
    Icolor{c} = permute(Icolor{c},dimorder);
    Icolor{c} = flipdim(Icolor{c},1);
end

% Create the (empty) mosaic
M0 = cell(size(Icolor));
Isize = size(Icolor{1});
for c = 1:numel(Icolor),
    M0{c} = zeros(msize .* [Isize(1) Isize(2)]);
end

% Check number of slices and do Procrustes where needed
tmp = length(slices) - prod(msize);
if tmp > 0,
    fprintf('  Info: More slices selected than slots available. The last %g slices will be omitted\n',tmp);
    slices = slices(1:prod(msize));
elseif tmp < 0,
    fprintf('  Info: Less slices selected than slots available. The last %g slots will be empty\n',abs(tmp));
    slices(prod(msize)) = 0;
end

% Define which slice goes where
slpos = reshape(slices,fliplr(msize))';
slpos = kron(slpos,ones(Isize(1),Isize(2)));

% Fill the mosaic with the slices
for s = 1:numel(slices),
    if slices(s) > 0,
        curpos = slpos == slices(s);
        for c = 1:numel(M0),
            M0{c}(curpos) = Icolor{c}(:,:,slices(s));
        end
    end
end

% Convert from cell to a 3D volume
M = zeros([size(M0{1}) numel(M0)]);
for c = 1:numel(M0),
    M(:,:,c) = M0{c};
end