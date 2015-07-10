function [cidx,map] = data2color(varargin)
% Convert an array into a map of indexed colours, that
% can be used, e.g., to create a PNG, PLY, OBJ or other
% formats. This is a low-level command, to be used by
% other functions.
% 
% Usage:
% [cidx,cmap] = data2color(data,cp,mapname,colourE,colourM,mapsize)
% 
% Inputs:
% - data       : Data to be coloured.
% - cp         : Critical points. There are three possibilities:
%                - Empty []: In this case all data is coloured,
%                  and scaled to fit the colourmap on its entirety.
%                - Two datapoints [x1 x2]: Scale the colourmap to
%                  fit the data between x1 and x2 only, which are
%                  then the two extremes of the viewing window.
%                  Extreme values are coloured according to the
%                  extremes of the colourmap, or using the colour
%                  indicated in 'colourextr'.
%                - Four datapoints [x1 z1 z2 x2]: This case is
%                  similar to the previous, except that the area
%                  between the datapoints z1 and z2 is ignored
%                  for the colour mapping, and painted in a solid
%                  colour given by 'colourmid'. 
%                In all cases, -Inf and +Inf are considered extremes
%                and are coloured either with the extremes of the
%                colourmap, or with 'colourE' if supplied.
%                The NaNs, if present, are painted with 'colourmid'.
% - mapname    : An Octave or MATLAB colourmap. Default: 'jet'.
% - colourE    : Colour to paint the extremes when not using the
%                extremes of the colourmap.
% - colourM    : Colour to paint the datavalues between z1 and z2
%                and to paint NaNs if present. Default is [0 0 0].
% - mapsize    : Size (resolution) of the map. Default: 2^24.
% 
% Outputs:
% - cidx       : Indices for the colours, in the same size as the
%                original data.
% - map        : Map with the colours.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2012
% http://brainder.org

% Get inputs
v = struct(                ...
    'data',    [],         ...  % arg 1
    'cp',      [],         ...  % arg 2
    'mapname', 'jet',      ...  % arg 3
    'colourE',  [],        ...  % arg 4
    'colourM',  [1 1 1],   ...  % arg 5
    'mapsize', 2^12);           % arg 6
fields = fieldnames(v);
if nargin < 2 || nargin > 6,
    error('Incorrect number of arguments.');
end
for a = 1:nargin,
    v.(fields{a}) = varargin{a};
end

% Put aside eventual NaNs. This also vectorises the data.
inan = isnan(v.data);
if any(inan),
    warning('There are NaNs in the data. These will be marked with the same colour used for mid values.'); %#ok
end
d = v.data(~inan);
datasize = size(v.data);

% Deal with infinities
tmp  = isinf(d);
mind = min(d(~tmp));
maxd = max(d(~tmp));
if any(tmp),
    d(tmp & d < 0) = mind - 1;
    d(tmp & d > 0) = maxd + 1;
end

% Critical points for cases 0 and 2.
% The case 4 is treated within the switch
v.cp = sort(v.cp);
if numel(v.cp) == 0,
    x1 = mind;
    x2 = maxd;
elseif numel(v.cp) == 2,
    x1 = v.cp(1);
    x2 = v.cp(2);
elseif numel(v.cp) == 4;
    x1 = v.cp(1);
    z1 = v.cp(2);
    z2 = v.cp(3);
    x2 = v.cp(4);
end

% Behave differently depending on the number of critical points
switch numel(v.cp),
    case { 0, 2 },
        if any(inan),
            ms = v.mapsize - 1;
        else
            ms = v.mapsize;
        end
        if isempty(v.colourE),
            map = eval(sprintf('%s(%f);',v.mapname,ms));
            b = [x1 1; x2 1]\[1; ms];
            idx = round(b(1)*d + b(2));
            idx(idx < 1)  = 1;
            idx(idx > ms) = ms;
        else
            map = eval(sprintf('%s(%f);',v.mapname,ms-1));
            map = [map ; v.colourE];
            b = [x1 1; x2 1]\[1; ms-1];
            idx = round(b(1)*d + b(2));
            idx(idx < 1)  = ms;
            idx(idx > ms) = ms;
        end
        if any(inan),
            map = [map; v.colourM];
        end
    case 4,
        ms = v.mapsize;
        idx = zeros(size(d));
        if isempty(v.colourE),
            map = eval(sprintf('%s(%f);',v.mapname,ms-1));
            map = [map ; v.colourM];
            b = [x1 1; x2-(z2-z1) 1]\[1; ms-1];
            idx(d < x1)        = 1;
            idx(d>=x1 & d<=z1) = round(b(1)*d(d>=x1 & d<=z1) + b(2));
            idx(d>z1 & d<z2)   = ms;
            idx(d>=z2 & d<=x2) = round(b(1)*d(d>=z2 & d<=x2) + b(2) - b(1)*(z2-z1));
            idx(d > x2)        = ms - 1;
            idx(idx < 1)       = 1;
            idx(idx > (ms+1))  = ms - 1;
        else
            map = eval(sprintf('%s(%f);',v.mapname,ms-2));
            map = [map ; v.colourE; v.colourM];
            b = [x1 1; x2-(z2-z1) 1]\[1; ms-2];
            idx(d < x1)        = ms - 1;
            idx(d>=x1 & d<=z1) = round(b(1)*d(d>=x1 & d<=z1) + b(2));
            idx(d>z1 & d<z2)   = ms;
            idx(d>=z2 & d<=x2) = round(b(1)*d(d>=z2 & d<=x2) + b(2) - b(1)*(z2-z1));
            idx(d > x2)        = ms - 1;
            idx(idx < 1)       = ms - 1;
            idx(idx > (ms+1))  = ms - 1;
        end
    otherwise
        error('The number of critical points must be zero, two or four');
end

% Reshape to original size, and colour the NaNs
cidx = zeros(datasize);
cidx(~inan) = idx;
cidx(inan) = size(map,1);
