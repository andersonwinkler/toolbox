function [Exp,Obs] = wachter(varargin)
% Return the expected and observed eigenvalues of the covariance of
% the matrix X. The expected are based on the Marcenko-Pastur distribution.
% These expected and observed can be used in a QQ-plot, which allows
% a better informed choice for the intrinsic dimensionality of the data
% than the scree test.
%
% Usage:
% [Exp,Obs] = wachter(X,Z,nrm)
% 
% Inputs:
% X     : Data matrix
% Z     : Nuisance variables. By default, Z = ones(N,1),
%         i.e., an intercept for mean-centering.
%         To omit removal of any nuisance, and thus not
%         even mean-center, use Z = [].
% nrm   : Boolean indicating whether the columns of the matrix should be
%         normalized to the same variance (true) or not (false).
% 
% Outputs:
% Exp   : Expected eigenvalues.
% Obs   : Observed eigenvalues.
% 
% Reference:
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
Z   = 1;
nrm = false;
if nargin < 1 || nargin > 3
    error('Incorrect number of arguments');
elseif nargin == 2
    Z   = varargin{2};
elseif nargin == 3
    Z   = varargin{2};
    nrm = varargin{3};
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
    [Q,D,~] = svd(null(Z'),'econ');
    Q = Q*D;
end
X = Q'*X;
if nrm
    X = bsxfun(@rdivide,X,std(X,0,1));
end
[n,p] = size(X);
y     = n/p;
sigsq = var(X(:),0);

% Covariance matrix its eigenvalues
C = (X'*X)/n;
s = svd(C,0)';

% Quantiles
K  = min(n,p);
q  = ((K:-1:1)-.5)/K;
Exp = zeros(size(q));
for k = 1:K
    Exp(k) = mpinv(q(k),y,sigsq);
end
Exp = Exp./y;  % Expected
Obs = s(1:K);  % Observed
