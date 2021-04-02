function [scores,coeffs,evals,extras] = npca(varargin)
% Return the first p eigenvectors and eigenvalues of the 
% covariance matrix of X.
%
% Usage:
% [scores,coeffs,evals,extras] = npca(X,p,Z)
%
% Inputs:
% X     : 2D array
% p     : Number of eigenvectors and eigenvalues to be
%         returned.
% Z     : Nuisance variables. By default, Z = ones(N,1),
%         i.e., an intercept for mean-centering.
%         To omit removal of any nuisance, and thus not
%         even mean-center, use Z = [].
%
% Outputs:
% scores : Principal components.
% coeffs : Principal coefficients (PCA loadings).
% evals  : Eigenvalues.
% extras : A struct containing a bunch of self-explanatory
%          useful outputs, including variance explained
%          by each principal component and the semi-orthogonal
%          matrix used for residualisation.
%
% Notes:
% * If the data were residualised before calling this function
%   and Z = [], then you need to rescale the eigenvalues
%   as Evals*N/(N-Nz).
% * To restore the original (residualised) data, use Scores*Loads'.
%
% This function was called "pca", but after Matlab created its
% own similarly-named "pca" function, it was renamed to "npca".
% It is a further improvement over the previous "epca".
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Nov/2010 (first version)
% Nov/2020 (this version)
% http://brainder.org

% Accept and check arguments
p  = 1;
Z  = 1;
if nargin < 1 || nargin > 3
    error('Incorrect number of arguments');
elseif nargin == 2
    p = varargin{2};
elseif nargin == 3
    p = varargin{2};
    Z = varargin{3};
end
X = varargin{1};
[nR,nC] = size(X);
if p > max(nR,nC)
    error('Cannot extract more eigenvalues than rows or columns.');
end

% Remove nuisance variables
if isscalar(Z)
    Z = ones(nR,1);
end
if isempty(Z)
    Q = eye(nR);
else
    [Q,D,~] = svd(null(Z'));
    Q = Q*D;
end
X = Q'*X;

% Save some memory by working with the
% smallest possible square of X
if nR >= nC
    [~,SS,V] = svd(X'*X,0);
    Vp  = V(:,1:p);
    SSp = SS(1:p,1:p);
    Up  = X*Vp;
else
    [U,SS,~] = svd(X*X',0);
    Up  = U(:,1:p);
    SSp = SS(1:p,1:p);
    Vp  = X'*Up;
end

% Pick one sign
s = diag(sign(Up(1,:)));

% Eigenvectors and eigenvalues
scores = Up*s;
coeffs = Vp*s;
if isscalar(Z)
    df = nR - 1;
else
    df = size(X,1);
end
evals = diag(SSp)./df;

% Some extra outputs, not normally needed
if nargout == 4
    
    % Scaled eigenvectors
    extras.scores_scaled = X*Vp;
    extras.coeffs_scaled = Up'*X;
    
    % Unit norm eigenvectors
    if nR >= nC
        extras.scores_unit = scores/sqrt(SSp);
        extras.coeffs_unit = coeffs;
    else
        extras.scores_unit = scores;
        extras.coeffs_unit = (coeffs/sqrt(SSp))';
    end
    
    % Recovered data using the p eigenvectors
    extras.recovered_with_scores = Up*extras.coeffs_scaled;
    extras.recovered_with_coeffs = extras.scores_scaled*Vp';
    
    % Variance explained
    S = diag(SS)./df;
    extras.variance_explained = S./sum(S);
    
    % Semiorthogonal matrix used for residualisation
    extras.semi_ortho = Q;
end
