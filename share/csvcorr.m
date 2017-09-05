function C = csvcorr(varargin)
% Compute simple and partial correlation coefficients
% among the columns of .csv files (no headings for rows or cols,
% only numbers).
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Jun/2016
% http://brainder.org

% Do the OCTAVE stuff, with TRY to ensure MATLAB compatibility
try %#ok
    % Get the inputs
    varargin = argv();
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('\n');
        return;
    end
end

% More OCTAVE stuff
nargin = numel(varargin);

if nargin > 1,
    error('Too many arguments');
end

% Load the data
A = load(varargin{1});

% Simple correlation:
rA = corrcoef(A)

% Partial correlations
irA = inv(rA);
pA = -irA.*((diag(irA)*diag(irA)').^-.5)

% Merge correlations
idxr = ~ triu(true(size(rA)));
idxp = ~ tril(true(size(pA)));
C = rA.*idxr + pA.*idxp + eye(size(pA));
disp(C);