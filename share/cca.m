function [cc,A,B,U,V] = cca(X,Y,df)
% Do CCA via QR & SVD.
%
% Usage:
% [cc,A,B,U,V] = cca(X,Y,df)
%
% Inputs:
% X and Y : Arrays with N rows. The number of columns
%           in each may differ. The ranks of X and Y
%           aren't checked for speed. Both are assumed
%           to have been mean-centered, be free of
%           nuisance (partial CCA) via X=Rz*X and Y=Rz*Y,
%           and have degrees of freedom df=N-rank(Z).
% df      : Degrees of freedom (N - rank(Z)), where Z is
%           the matrix with nuisance variables already
%           regressed out from both X and Y. Since both
%           are supposed to be mean-centered, df should
%           be no larger than N-1. For speed, it is not
%           tested internally, though, so enter this
%           wisely.
% Outputs:
% cc       : Canonical correlations.
% A and B  : Canonical coefficients.
% U and V  : Canonical variables.
% 
% U=X*A and V=Y*B. such that each pair of columns in U and V
% are maximally correlated, under the orthonality constraint.
% 
% Based on the algorithm proposed by:
% * Bjorck A, Golub GH. Numerical methods for
%   computing angles between linear subspaces.
%   Math Comput. 1973;27(123):579-579.
%
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Nov/2018
% http://brainder.org

[Qx,Rx,iX] = qr(X,0);
[Qy,Ry,iY] = qr(Y,0);
k  = min(size(X,2),size(Y,2));
[L,D,M] = svd(Qx'*Qy,0);
A  = Rx\L(:,1:k)*sqrt(df);
B  = Ry\M(:,1:k)*sqrt(df);
cc = min(max(diag(D(:,1:k))',0),1);
A(iX,:) = A;
B(iY,:) = B;
U  = X*A;
V  = Y*B;
