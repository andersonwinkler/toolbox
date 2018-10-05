function dpx2fdr(varargin)
% For a set of DPV/DPF files containing p-values, return the
% FDR threshold (single value) and the FDR corrected and FDR-adjusted
% as new DPF/DPV files.
%
% Usage:
% dpx2fdr(qvalue,dpxfiles)
%
% qvalue   : Alowed proportion of false discoveries.
% dpxfiles : Files containing the p-values, between 0 and 1. If there are
%            values outside these limits and all have the same sign,
%            assumes that the data are organised as -log10(p) or log10(p).
%            Otherwise, returns an error.
%
% The results are saved with filenames similar to the original,
% with suffixes 'cor' and 'adj'. The FDR-threshold is printed
% to the screen.
%
% Datapoints marked as NaN or Inf will be ignored.
%
% This function calls fdr.m, which contains the actual implementation
% of FDR and references.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011

% Do some OCTAVE stuff, but use TRY to ensure MATLAB compatibility
try
    % Get the inputs
    varargin = argv();
    nargin = numel(varargin);

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('For a set of DPV/DPF files containing p-values, return the\n');
        fprintf('FDR threshold (single value) and the FDR corrected and FDR-adjusted\n');
        fprintf('as new DPF/DPV files.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('dpx2fdr <qvalue> <dpxfiles>\n');
        fprintf('\n');
        fprintf('qvalue   : Alowed proportion of false discoveries.\n');
        fprintf('dpxfiles : Files containing the p-values, between 0 and 1. If there are\n');
        fprintf('           values outside these limits and all have the same sign,\n');
        fprintf('           assumes that the data are organised as -log10(p) or log10(p).\n');
        fprintf('           Otherwise, returns an error.\n');
        fprintf('\n');
        fprintf('The results are saved with filenames similar to the original,\n');
        fprintf('with suffixes ''cor'' and ''adj''. The FDR-threshold is printed\n');
        fprintf('to the screen.\n');
        fprintf('\n');
        fprintf('Datapoints marked as NaN or Inf will be ignored.\n');
        fprintf('\n');
        fprintf('This function calls fdr.m, which contains the actual implementation\n');
        fprintf('of FDR and references.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Aug/2011\n');
        return;
    end
end

% Get the inputs
qval = varargin{1};
if ischar(qval), qval = eval(qval); end
dpx = [];
crd = cell(nargin-1,1);
idx = cell(nargin-1,1);
nX = zeros(nargin-1,1);
for a = 1:nargin-1,
    [tmp,crd{a},idx{a}] = dpxread(varargin{a+1});
    dpx = [dpx; tmp];
    nX(a) = numel(tmp);
end
didx = ~ isnan(dpx);

% Check ranges and if the data is log or not
if all(dpx(didx) >= 0 | dpx(didx) <= 1),
    % All between 0-1
    ptype = 0;
    pvals = dpx(didx);
elseif all(dpx(didx) >= 0),
    % -log10(p)
    ptype = -1;
    pvals = 10.^(-dpx(didx));
elseif all(dpx(didx) <= 0),
    %  log10(p)
    ptype = 1;
    pvals = 10.^dpx(didx);
end

% Compute the FDR threshold, corrected and adjusted
[thr,corall,adjall] = fdr(pvals,qval);

% Put back to the original scale
if ptype == -1,     % -log10(p)
    corall = -log10(corall);
    adjall = -log10(adjall);
elseif ptype == 1,  %  log10(p)
    corall = log10(corall);
    adjall = log10(adjall);
end

% Put back the NaNs
cor = nan(size(dpx));
cor(didx) = corall;
adj = nan(size(dpx));
adj(didx) = adjall;

% Save each of the files
fprintf('FDR-threshold: %f\n',thr);
fprintf('Saving FDR-corrected and FDR-adjusted files.\n');
acum = 0;
for a = 1:nargin-1,
    % Define file names
    [fpth,fnam,fext] = fileparts(varargin{a+1});
    % FDR-corrected
    tmp = cor(acum+1:acum+nX(a));
    dpxwrite(fullfile(fpth,sprintf('%s-cor%s',fnam,fext)),tmp,crd{a},idx{a});
    % FDR-adjusted
    tmp = adj(acum+1:acum+nX(a));
    dpxwrite(fullfile(fpth,sprintf('%s-adj%s',fnam,fext)),tmp,crd{a},idx{a});
    % Increment cumulant
    acum = acum+nX(a);
end

% Done!
