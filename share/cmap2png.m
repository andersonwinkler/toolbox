function cmap2png(varargin)
% Create a PNG image of a given colormap
%
% Usage:
% cmap2png(mapname,mapsize,imagesize,imagefile)
% cmap2png(map,imagesize,imagefile)
% 
% - mapname   : The mapname (string) of any colormap supported
%               by MATLAB.
% - mapsize   : The resolution of the map. For best results, use
%               numbers larger than the largest image dimension.
% - imagesize : The size of the image to be generated. The color
%               gradient will run along the largest dimension.
% - imagefile : File name to be saved (in PNG format).
% - map       : Instead of the mapname and mapsize, it's possible
%               to parse an already manipulated colormap.
%
% ==> This function won't work in Octave 3.2.4 due to a bug in
%     loading the ImageMagick libraries. The failure affects
%     both Debian and Mac. Use a newer Octave version, or Matlab.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% May/2011 (first version)
% Jul/2013 (this version)
% http://brainder.org

% Take the arguments
if nargin == 3,
    map       = varargin{1};
    imagesize = varargin{2};
    imagename = varargin{3};
elseif nargin == 4,
    imagesize = varargin{3};
    imagename = varargin{4};
    
    % Create the colormap if needed
    map = eval(sprintf('%s(%g)',varargin{1},varargin{2}));
else
    error('Incorrect number of arguments.');
end
    
% Define the indices, with the gradient running along the largest dimension
mapsize = size(map,1);
idx = round(repmat(linspace(1,mapsize,max(imagesize)),[min(imagesize) 1]));

% Transpose to vertical if needed
if imagesize(1) > imagesize(2),
    idx = idx';
end

% Create the PNG and assign the RGB triplets
png = zeros([size(idx) 3]);
for c = 1:3,
    png(:,:,c) = reshape(map(idx(:),c),size(idx));
end

% Save to the disk
imwrite(png,imagename);