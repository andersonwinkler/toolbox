function cmap = png2cmap(varargin)
% Generates a colourmap based on a PNG file. The file should contain
% nothing but a vertical or horizontal gradient. The output is an
% m-file that can be used as a Matlab/Octave colourmap.
% 
% Inputs:
% - pngfile   : Input PNG filename.
% - mapname   : Name o the map to be created (may or may not include
%               the full path or the .m extension).
% - direction : Horizontal (default) or vertical. This is specified
%               as a string.
% 
% _____________________________________
% Anderson M. Winkler
% Hospital Israelita Albert Einstein
% Jul/2017
% http://brainder.org

% Parse arguments:
narginchk(2,3);
pngfile = varargin{1};
mapname = varargin{2};
if nargin == 3,
    direction = varargin{3};
else
    direction = 'horizontal';
end

% Define the name of the colourmap:
[mpth,mnam,mext] = fileparts(mapname);
if ~ any(strcmpi(mext,{'.m',''})),
    error('Unknown file extension for a colourmap: %s',mext);
end
if isempty(mext),
    mext = '.m';
end

% Read the image file, transpose if needed:
[png,map0,alpha] = imread(pngfile);
if any(strcmpi(direction,{'v','ver','vert','vertical'})),
    png   = png';
    alpha = alpha';
end

% Separate channels:
rgb = cell(3,1);
for c = 1:3,
    rgb{c} = png(:,:,c);
end

% New map:
cmap = nan(size(png,2),3);
for c = 1:3,
    [U,Iimg,IU] = unique(rgb{c},'rows');
    idx = mode(IU);
    cmap(:,c) = rgb{c}(idx,:)';
end
cmap = cmap/255;

% Save to an m-file:
fid = fopen(fullfile(mpth,strcat(mnam,mext)),'w');
fprintf(fid,'function map = %s(n)\n',mnam);
fprintf(fid,'%% Returns the "%s" colourmap, derived from the file:\n',mnam);
fprintf(fid,'%% %s\n',pngfile);
fprintf(fid,'\n');
fprintf(fid,'if nargin < 1\n');
fprintf(fid,'    n = size(get(gcf,''Colormap''),1);\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');
fprintf(fid,'values = [\n');
fprintf(fid,'%f %f %f\n',cmap');
fprintf(fid,'    ];\n');
fprintf(fid,'\n');
fprintf(fid,'P = size(values,1);\n');
fprintf(fid,'map = interp1(1:size(values,1), values, linspace(1,P,n), ''linear'');\n');
