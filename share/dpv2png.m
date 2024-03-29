function dpv2png(varargin)
% Converts an ASCII Curvature/Surface pair into a PNG file.
%
% Usage:
% dpv2png(DPVFILE,SRFFILE,PNGFILE,MAPNAME,DATARANGE,COLORMIN,COLORMAX,...
%           MAPSIZE,INVERTSGN,LOGSCALE,STATS,DUAL,ROTATE)
%
% - DPVFILE   = Data-per-vertex file, in ASCII format
% - SRFFILE   = Reference surface file, in ASCII format (optional)
% - PNGFILE   = PNG file to be created
% - MAPNAME   = An OCTAVE colormap (default = 'jet').
%               Custom maps can be used if defined by a function in OCTAVE path
% - DATARANGE = Interval to be used to define the colorscale [min max].
%               If not specified, used min and max of the DPV file.
% - SHOWRANGE = Interval to be shown [min max].
%               If not specified, used min and max of the DPV file.
% - COLORMIN  = Color to be used for values below the datarange, in the format [R G B],
%               where the values for R, G and B are in the interval 0-1.
% - COLORMAX  = Same as colormin, but for when the value is above the curvrange.
% - MAPSIZE   = The number of colors to be used in the colormap (default = 2^16).
% - INVERTSGN = True/False. Multiply by (-1).
% - LOGSCALE  = True/False. Convert to log-scale (base 10).
% - STATS     = True/False. prints in the screen min, max, 
%               mean, stdev, median and mode (default = false).
% - DUAL      = True/False. If true, shows positive/negative regions of no-overlap
%               between DATARANGE and SHOWRANGE.
% - ROTATE    = True/False. Rotates the sphere 90 degrees around the
%               Z-axis. Useful for the left hemisphere (not for the right). 
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% May/2011

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Generate a PNG file using the geometry of a\n');
        fprintf('Surface file and the values of a Curvature file.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('dpv2png <dpvfile.dpv> <srffile.srf> <pngfile.png> [mapname] \\\n');
        fprintf('        [datarange] [colormin] [colormax] [mapsize] [invertsgn] \\\n');
        fprintf('        [logscale] [stats] [dual] [rotate]\n');
        fprintf('\n');
        fprintf('Compulsory arguments (in this order):\n');
        fprintf(' dpvfile    = Data-per-vertex file, in ASCII format\n');
        fprintf(' srffile    = Reference surface file, in ASCII format (optional)\n');
        fprintf(' pngfile    = PNG file to be created\n');
        fprintf('\n');
        fprintf('Optional arguments (in this order; use placeholders if needed):\n');
        fprintf(' mapname    = An OCTAVE colormap (default = jet).\n');
        fprintf('              Custom maps can be used if defined by\n');
        fprintf('              a function in OCTAVE''s path.\n');
        fprintf(' datarange* = Interval to be colored, in the format [min max].\n');
        fprintf('              If not specified, used min and max of the DPV file.\n');
        fprintf(' showrange* = Interval to be shown [min max].\n');
        fprintf('              If not specified, used min and max of the DPV file.\n');
        fprintf(' colormin*  = Color to be used for values below the datarange,\n');
        fprintf('              in the format [R G B], where the values\n');
        fprintf('              for R, G and B are in the interval 0-1.\n');
        fprintf(' colormax*  = Same as colormin, but for when the\n');
        fprintf('              value is above the curvrange.\n');
        fprintf(' mapsize    = The number of colors to be used in the colormap.\n');
        fprintf('              Default = 2^16.\n');
        fprintf(' invertsgn  = True/False. Multiply by (-1).\n');
        fprintf(' logscale   = True/False. Convert to log-scale (base 10).\n');
        fprintf(' stats      = True/False. Prints in the screen min, max, mean,\n');
        fprintf('              stdev, median and mode. Default = false.\n');
        fprintf(' dual       = True/False. If true, shows positive/negative regions of\n');
        fprintf('              no-overlap between datarange and showrange.\n');
        fprintf(' rotate     = True/False. Rotates the sphere 90 degrees around the\n');
        fprintf('              Z-axis. Useful for the left hemisphere (not for the right).\n');
        fprintf('* Should be entered between single quotes, e.g. ''[2 4]''.\n');
        fprintf('\n');
        fprintf('Use [] or ''[ ]'' as placeholders for optional arguments not given.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('May/2011');
        return;
    end
end

% Defaults
defmapname   = 'jet';
defmapsize   = 2^16;
definvsignal = false;
deflogscale  = false;
defprintdiag = true;
defdual      = false;
defrotate    = false;

d.fcrv = '';
d.fsrf = '';
d.fpng = '';
d.mapname   = defmapname;
d.crvrange  = [-1.7 2];
d.shwrange  = [];
d.colormin  = [];
d.colormax  = [];
d.mapsize   = defmapsize;
d.invsgn    = definvsignal;
d.logscale  = deflogscale;
d.printdiag = defprintdiag;
d.dual      = defdual;
d.rotate    = defrotate;

% Accept user arguments
nargin = numel(varargin); % Redundant in MATLAB, but fixes an OCTAVE compatibility issue
fields = {'fcrv','fsrf','fpng','mapname','crvrange','shwrange',...
    'colormin','colormax','mapsize','invsgn','logscale','printdiag','dual','rotate'};
for a = 1:nargin,
    d.(fields{a}) = varargin{a};
    if a >= 5 && ischar(d.(fields{a})),
        d.(fields{a}) = eval(d.(fields{a}));
    end
end

% Accept some defaults if needed
if isempty(d.mapname), d.mapname = defmapname; end
if isempty(d.mapsize), d.mapsize = defmapsize; end

% Define the colormap
rgb = eval(sprintf('%s(d.mapsize)',d.mapname));
if isempty(d.colormin),
    d.colormin = rgb(1,:);
end
if isempty(d.colormax),
    d.colormax = rgb(d.mapsize,:);
end
rgb = double([ d.colormin ; rgb ; d.colormax ]);

% Read the curvature file (vtx coords and assoc vals)
[crv,crd,idx] = crvread(d.fcrv);

% Read the reference surface file (should be a sphere)
if ~ isempty(d.fsrf),
    [vtx,fac] = srfread(d.fsrf);
    fac = fac-1; % face indices start at 0 here, not 1.
    nV = size(vtx,1);
    nF = size(fac,1);
else
    % Or use the coordinates present in the CRV file
    vtx = crd;
end

% Convert to spherical coordinates
if d.rotate,
    s = -1;
else
    s = 1;
end
[Az,El,R] = cart2sph(s*vtx(:,1),s*vtx(:,2),vtx(:,3));

% Convert to degrees and make sure it's all positive
Az = 180*Az/pi;
El = 180*El/pi;

% Invert the signal
if d.invsgn,
    crv = -crv;
end

% Convert to log-scale
if d.logscale,
    crv = -log10(crv);
end

% Interpolate to a regular grid
st = 1;
[lon,lat] = meshgrid(-179:st:179,-89:st:89);
crv = griddata(Az,El,crv,lon,lat);

% Print some diagnostics in the screen
cmin = min(crv(:)); cmax = max(crv(:));
if d.printdiag,
    fprintf('Image stats:\n');
    fprintf(' [ Min Max ] = [ %g %g ]\n',cmin,cmax);
    fprintf(' [ Mean Std ] = [ %g %g ]\n',mean(crv(:)),std(crv(:)));
    fprintf(' [ Median Mode ] = [ %g %g ]\n',median(crv(:)),mode(crv(:)));
end

% Compute the color indices. Already considers the out-of-range cases:
if ~isempty(d.crvrange),
    cmin = d.crvrange(1); cmax = d.crvrange(2);
end
cidx = round((crv-cmin)*(d.mapsize+1-2)/(cmax-cmin)+2);
cidx(cidx < 1) = 1;
cidx(cidx > d.mapsize+2) = d.mapsize+2;

% Interval to display
if isempty(d.shwrange),
    shwrange = [cmin cmax];
else
    shwrange = [cmin cmax];
end

% Define the color indices that will be shown with the color scale
% (i.e., the window)
if ~ d.dual,
    sidx = round((shwrange-cmin)*(d.mapsize+1-2)/(cmax-cmin)+2);
    cidx(cidx <= sidx(1)) = 1;
    cidx(cidx >= sidx(2)) = d.mapsize+2;
else
    tmp = [shwrange(1) mean(shwrange) shwrange(2)];
    sidx = round((tmp-cmin)*(d.mapsize+1-2)/(cmax-cmin)+2);
    cidx(cidx >= sidx(2) & cidx < sidx(3)) = d.mapsize+2;
    cidx(cidx <= sidx(2) & cidx > sidx(1)) = 1;
end

% Replace NaNs for the lowest color in the map
cidx(isnan(cidx)) = 1;

% Prepare the PNG file
P = zeros([size(lon) 3]);
for c = 1:3,
    C = zeros(size(lon));
    C(1:numel(crv)) = rgb(cidx(:),c);
    P(:,:,c) = C;
end

% Save it
imwrite(P,d.fpng,'png');
