function varargout = holm(varargin)
% Computes the Holm threshold for a vector of p-values.
%
% Usage:
% [pthr,pcor,padj] = holm(pvals)
%                    holm(pval,alpha)
%
% Inputs:
% pvals  = Vector of p-values.
% alpha  = Allowed error rate.
%          Default = 0.05.
%
% Outputs:
% pthr   = Holm threshold.
% pcor   = Holm corrected p-values.
% padj   = Holm adjusted p-values.
%
%
% Reference:
% * Holm S. A simple sequentially rejective multiple test
%   procedure. Scand J Stat. 1979;6(2):65-70.
%
% ________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2014
% http://brainder.org

% Accept arguments
switch nargin,
    case 0,
        error('Error: Not enough arguments.');
    case 1,
        pval = varargin{1};
        alph = 0.05;
    case 2,
        pval = varargin{1};
        alph = varargin{2};
    otherwise
        error('Error: Too many arguments.')
end

% Check if pval is a vector
if numel(pval) ~= length(pval),
    error('p-values should be a row or column vector, not an array.')
end

% Check if pvals are within the interval
if min(pval) < 0 || max(pval) > 1,
    error('Values out of range (0-1).')
end

% Check if alpha is within the interval
if alph < 0 || alph > 1,
    error('Significance level (alpha) out of range (0-1).')
end

% ========[PART 1: HOLM THRESHOLD]========================================

% Sort p-values
[pval,oidx] = sort(pval);

% Number of observations
V = numel(pval);

% Order (indices), in the same size as the pvalues
idx = reshape(1:V,size(pval));

% Line to be used as cutoff
thrline = alph./(V-idx+1);

% Find the largest pval, still under the line
thr = max(pval(pval<=thrline));

% Deal with the case when all the points under the line
% are equal to zero, and other points are above the line
if thr == 0,
    thr = max(thrline(pval<=thrline));
end

% Case when it does not cross
if isempty(thr), thr = 0; end

% Returns the result
varargout{1} = thr;

% ========[PART 2: HOLM CORRECTED]========================================

if nargout == 2 || nargout == 3,
    
    % p-corrected
    pcor = pval.*(V-idx+1);

    % Sort back to the original order and output
    [~,oidxR] = sort(oidx);
    varargout{2} = pcor(oidxR);
end

% ========[PART 3: HOLM ADJUSTED ]========================================

if nargout == 3,

    % Loop over each sorted original p-value
    padj = zeros(size(pval));
    prev = 1;
    for i = V:-1:1,
        % For the adjustment, see Westfall & Wolfinger (2000). For Holm,
        % the corrected and adjusted are nearly always the same, but not
        % always. The adjusted should be preferred.
        padj(i) = min(prev,pval(i)*(V-i+1));
        prev = padj(i);
    end
    varargout{3} = padj(oidxR);
end

% That's it!