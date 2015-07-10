% #!/usr/bin/octave -q
function annot2dpv(varargin)
% Convert an annotation file to a DPV file.
%
% Usage:
% annot2dpv(annotfile,dpvfile)
%
% Inputs:
% annotfile : Annotation file.
% dpvfile   : Output DPV file.
%
% Before running, be sure that ${FREESURFER_HOME}/matlab is
% in the OCTAVE/MATLAB path.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011
% http://brainder.org

% Do some OCTAVE stuff, but use TRY to ensure MATLAB compatibility
try
    % Get the inputs
    varargin = argv();
    nargin = numel(varargin);
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Convert an annotation file to a DPV file.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('annot2dpv <annotfile> <dpvfile>\n');
        fprintf('\n');
        fprintf('Inputs:\n');
        fprintf('annotfile : Annotation file.\n');
        fprintf('dpvfile   : Output DPV file.\n');
        fprintf('\n');
        fprintf('Before running, be sure that ${FREESURFER_HOME}/matlab is\n');
        fprintf('in the OCTAVE/MATLAB path.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Aug/2011\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Some FS commands are needed now
fspath = getenv('FREESURFER_HOME');
if isempty(fspath),
    error('FREESURFER_HOME variable not correctly set');
else
    addpath(fullfile(fspath,'matlab'));
end

% Accept arguments
annotfile = varargin{1};
dpvfile   = varargin{2};

% Read the annotation file
[~,lab,ctab] = read_annotation(annotfile);

% For each structure, replace its coded colour by its index
for s = 1:ctab.numEntries,
    lab(lab == ctab.table(s,5)) = s;
end

% Save the result
dpxwrite(dpvfile,lab)
