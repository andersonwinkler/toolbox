function P = fisherexact(varargin)
% Compute the p-value for the Fisher's exact test, modified
% for hypothesis testing. The alternative hypothesis is that
% one or both anti-diagonal elements are equal to zero, which
% is the best possible performance if one variable can predict
% another.
% 
% Usage:
% P = fisherexact(C);
% P = fisherexact(Y,X);
% 
% C : Contingency table (2x2).
% Y : Observed dicothomous variable. It can be a Nx1 vector
%     or a NxM array.
% X : Explanatory dicothomous variable. It can be a Nx1 vector
%     or a NxM array.
% P : P-values (uncorrected for multiple testing) for all the
%     M columns of X and/or Y.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Mar/2014
% http://brainder.org

% Take inputs
if nargin == 1,
    a = varargin{1}(1);
    c = varargin{1}(2);
    b = varargin{1}(3);
    d = varargin{1}(4);
    M = 1;
elseif nargin == 2,
    Y = logical(varargin{1});
    X = logical(varargin{2});
    M = max(size(Y,2),size(X,2));
    if size(X,1) ~= size(Y,1),
        error('Input arguments must have the same number of rows.');
    end
    if size(X,2) > 1 && size(Y,2) > 1 && size(X,2) ~= size(Y,2),
        error('X and Y must have the same number of columns, or one of these must be a single column.');
    end
    a = sum(bsxfun(@and, Y, X),1);
    b = sum(bsxfun(@and, Y,~X),1);
    c = sum(bsxfun(@and,~Y, X),1);
    d = sum(bsxfun(@and,~Y,~X),1);
else
    error('Incorrect number of arguments');
end

% Margins & total
ab = a + b;
cd = c + d;
ac = a + c;
bd = b + d;
n  = unique(a + b + c + d);

% Logs for below
lfac = lfactorial(n);

% Compute the probabilities for each state of the table towards
% the alternative, that is, the case in which one or both elements
% of the anti-diagonal are 0.
P = zeros(1,M);
for s  = 0:max(max(b(:),c(:))),
    as = a + s;
    bs = b - s;
    cs = c - s;
    ds = d + s;
    m  = bs >= 0 & cs >= 0;
    lPs = ...
        sum(lfac([ab(m) cd(m) ac(m) bd(m)]+1)) - ...
        sum(lfac([as(m) bs(m) cs(m) ds(m) n]+1));
    P(m) = P(m) + exp(lPs)';
end

% ==============================================================
function lfac = lfactorial(N)
% Compute the log(factorial(0:N)), so dealing with
% precision issues.
lfac = zeros(N+1,1);
for n = 1:N,
    lfac(n+1) = log(n) + lfac(n);
end
