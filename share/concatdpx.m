function concatdpx(varargin)
% Read a set of DPV or DPF files and save as a single, large CSV
% table, containing one face or vertex per line and one the data
% from each file (typically, a subject) per column.
% 
% Usage:
% concatdpx('file1.dpf','file2.dpf',...,'fileout.dpf')
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jul/2011

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin   = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),

        fprintf('Read a set of DPV or DPF files and save as a single, large CSV\n');
        fprintf('table, containing one face or vertex per line and one the data\n');
        fprintf('from each file (typically, a subject) per column.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('concatdpx <file1.dpf> <file2.dpf> ... <fileout.dpf> \n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Jul/2011\n');
        return;
    end
end

% Read the first file to allocate memory
fprintf('Reading %s\n',varargin{1});
crv = crvread(varargin{1});
nP  = numel(crv);
nS  = nargin-1;

% Create empty var and add 1st subject
Y = zeros(nS,nP);
Y(1,:) = crv(:);

% Loop over subjects
for s = 2:nargin-1,
    fprintf('Reading %s\n',varargin{s});
    crv = crvread(varargin{s});
    Y(s,:) = crv(:);
end

% Save the result 4D file as CSV
% Each row is a face or vertex
% Each column is an original file (e.g. a subject)
fstr = repmat('%g ',[1 nS]);
fstr = strcat(fstr(1:numel(fstr)-1),'\n');
fid = fopen(varargin{nargin},'w');
fprintf(fid,fstr,Y);
fclose(fid);
