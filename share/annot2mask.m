function annot2mask(varargin)
% Create a mask from an annotation file
%
% Usage:
% annot2dpv(annotfile,maskfile,masktype,roinames...)
%
% Inputs:
% annotfile : Annotation file.
% maskfile  : Output mask (use extension .mgh or .mgz).
% masktype  : Either "include" or "exclude", indicating whether the regions
%             should be added to an initial empty mask, or removed from an
%             initial full mask.
% roinames  : Names or numerical indices of the regions to be
%             included or excluded in the mask.
%
% Before running, be sure that ${FREESURFER_HOME}/matlab is
% in the OCTAVE/MATLAB path.
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Mar/2024
% http://brainder.org

% Do some OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')

    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q')
        fprintf('Create a mask from an annotation file\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('annot2dpv annotfile maskfile masktype roinames ...\n');
        fprintf('\n');
        fprintf('Inputs:\n');
        fprintf('annotfile : Annotation file.\n');
        fprintf('maskfile  : Output mask (use extension .mgh or .mgz).\n');
        fprintf('masktype  : Either "include" or "exclude", indicating whether the regions\n');
        fprintf('            should be added to an initial empty mask, or removed from an\n');
        fprintf('            initial full mask.\n');
        fprintf('roinames  : Names or numerical indices of the regions to be\n');
        fprintf('            included or excluded in the mask.\n');
        fprintf('\n');
        fprintf('Before running, be sure that ${FREESURFER_HOME}/matlab is\n');
        fprintf('in the OCTAVE/MATLAB path.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('UTRGV\n');
        fprintf('Mar/2024\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Inputs
annotfile = varargin{1};
maskfile  = varargin{2};
masktype  = varargin{3};
roilist   = varargin(4:end);

% Ensure FreeSurfer is found
fshome = getenv('FREESURFER_HOME');
if isempty(fshome) && ~ exist('read_annotation','file')
    error('FreeSurfer not found.');
else
    addpath(fullfile(fshome,'matlab'));
end

% Read the annotation file
[~,lab,ctab] = read_annotation(annotfile);

% For each structure, replace its coded colour by its index
for s = 1:ctab.numEntries
    lab(lab == ctab.table(s,5)) = s;
end

% Start an overall mask with all vertices selected, or none selected
if strcmpi(masktype,'include')
    mask   = false(size(lab));
    newval = true;
elseif strcmpi(masktype,'exclude')
    mask   = true(size(lab));
    newval = false;
else
    error('Unknown option: %s',starting);
end

% For each region, replace include or exclude vertices from the mask
% according to the regions requested
for r = 1:numel(roilist)
    if ischar(roilist{r})
        roinum = str2double(roilist{r});
        if isnan(roinum)
            namematch = strcmpi(roilist{r},ctab.struct_names);
            if any(namematch)
                idx = lab == find(namematch);
                if strcmpi(roilist{r},'unknown')
                    idx = idx | lab == 0;
                end
            else
                error('Region name not found: %s',roilist{r});
            end
        else
            idx = lab == roinum;
        end
    elseif isnumeric(roilist{r})
        roinum = roilist{r};
        idx = lab == roilist{r};
    else
        error('Unknown region: %s',roilist{r})
    end
    if any(idx)
        mask(idx) = newval;
    else
        error('Region number not found: %d',roinum)
    end
end

% Save as mgh/mgz
save_mgh(mask,maskfile,eye(4));
