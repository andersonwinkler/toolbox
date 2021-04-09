function [Exp,Obs,Pexp,Pobs] = wachter(varargin)
% Return the expected and observed eigenvalues of the covariance of
% the matrix X. The expected are based on the Marcenko-Pastur distribution.
%
% These expected and observed can be used in a QQ-plot, which allows
% a better informed choice for the intrinsic dimensionality of the data
% than the scree test.
%
% Usage:
% [Exp,Obs,Pexp,Pobs] = wachter(X,Z,nrm,tail)
%
% Inputs:
% X     : Data matrix of size n by p
% Z     : Nuisance variables. By default, Z = ones(N,1),
%         i.e., an intercept for mean-centering.
%         To omit removal of any nuisance, and thus not
%         even mean-center, use Z = [].
% nrm   : Boolean indicating whether the columns of the matrix should be
%         normalized to the same variance (true) or not (false).
%         Default is false.
% tail  : Boolean indicating whether upper-tail probabilities (p-values)
%         will be provided in Pexp and Pobs, or simply the CDF.
%         Default is true, i.e., p-values are returned.
%
% Outputs:
% Exp   : Expected eigenvalues.
% Obs   : Observed eigenvalues.
%
% Notes:
% * Dividing X by 1./sqrt(n/p) requires Exp to be scaled by n/p.
%
% References:
% * Bai ZD. Methodologies in spectral analysis of large dimensional
%   random matrices, a review. Statistica Sinica 1999;9(3):611?77.
% * Johnstone IM. On the distribution of the largest eigenvalue in
%   principal components analysis. The Annals of Statistics
%   2001;29(2):295–327.
% * Marcenko VA, Pastur LA. Distribution of eigenvalues for some
%   sets of random matrices. Math USSR Sb 1967;1(4):457?83.
% * Wachter KW. Probability plotting points for principal components.
%   In: Proceedings of the Ninth Interface Symposium on Computer Science
%   and Statistics. Harvard University and Massachussetts Institute of
%   Technology: Prindle, Weber & Schmidt; 1976.
%
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Apr/2021
% http://brainder.org

% Accept and check arguments
Z    = 1;
nrm  = false;
tail = true;
if nargin < 1 || nargin > 4
    error('Incorrect number of arguments');
elseif nargin == 2
    Z    = varargin{2};
elseif nargin == 3
    Z    = varargin{2};
    nrm  = varargin{3};
elseif nargin == 4
    Z    = varargin{2};
    nrm  = varargin{3};
    tail = varargin{4};
end
X  = varargin{1};
nR = size(X,1);

% Remove nuisance variables
if isscalar(Z)
    Z = ones(nR,1);
end
if isempty(Z)
    Q = eye(nR);
else
    [Q,D,~] = svd(null(Z'),0);
    Q = Q*D;
end
X = Q'*X;
if nrm
    X = bsxfun(@rdivide,X,std(X,0,1));
end
[n,p] = size(X);
K     = min(n,p);
y     = n/p;
if y > 1, y = 1./y; end
sigsq = var(X(:),0);

% Covariance matrix its singular values
% Save some memory by working with the smallest
% possible covariance matrix
if n >= p
    SS = svd(X'*X,0)';
else
    SS = svd(X*X',0)';
end
SS = SS./min(n,p);
Obs = SS(1:K); % Observed eigenvalues

% Quantiles
if tail
    q = ((1:K)-.5)/K;
else
    q = ((K:-1:1)-.5)/K;
end
Exp  = zeros(size(q));
%Pexp = zeros(size(q));
Pobs = zeros(size(q));
for k = 1:K
    Exp(k)  = mpinv(  q(k),   y,sigsq,tail)./y; % Expected eigenvalues
    %Pexp(k) = mpcdf(Exp(k).*y,y,sigsq,tail);    % p-values of expected (these are the same as the quantiles, uncomment if you want to confirm mpinv.m is correct)
    Pobs(k) = mpcdf(Obs(k).*y,y,sigsq,tail);    % p-values of observed
end
Pexp = q;
