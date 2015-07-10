function mat2png(M,fname,mapname,minmax,mask,cmask)
% Save a 2D array as a PNG file.
%
% Usage:
% mat2png(M,fname,mapname,minmax,mask,cmask)
% 
% Inputs:
% M       : Matrix (2D array) to be saved as PNG
% fname   : File name to be created
% mapname : Map to be used (it can also be the actual map)
% minmax  : [min max] of values to be coloured
% mask    : Mask of pixels not to be coloured
% cmask   : RGB triplet of masked out pixels.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2015
% http://brainder.org

% Colour map 
if ischar(mapname),
    map = eval(sprintf('%s(%d)',mapname,2^10));
else
    map = mapname;
end
mapsiz = size(map,1);

mi = minmax(1);
ma = minmax(2);
mask = logical(mask);
idx = (M(mask)-mi)./(ma-mi).*(mapsiz-1) + 1;
idx(idx < 1)    = 1;
idx(idx > mapsiz) = mapsiz;
rgb = zeros([size(M) 3]);
for c = 1:3,
    tmp        = zeros(size(M));
    tmp(mask)  = map(round(idx),c);
    tmp(~mask) = cmask(c);
    rgb(:,:,c) = tmp;
end
%rgb = rgb*255;
imwrite(rgb,fname,'png');

