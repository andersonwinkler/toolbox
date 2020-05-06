function [U,S,L,Y] = pcamiss(X)
% Compute PCA for a matrix that has missing entries.
% 
% Inputs:
% - X : Input matrix. Missing entries are marked as NaN.
% 
% Outputs:
% - U : Eigenvectors.
% - S : Eigenvalues.
% - L : Loadings.
% - Y : Recovered version of X. If there are no missing values, Y = X.
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Jul/2019
% http://brainder.org

% Vars for later
[nR,nC] = size(X);
C       = zeros(nR);

% Location of the missing entries
idx = isnan(X);

% Mean for each column, excluding missing values
avgC  = zeros(1,nC);
for c = 1:nC
    avgC(c) = mean(X(~idx(:,c),c));
end

% Mean-center columns
for c = 1:nC
    X(~idx(:,c),c) = X(~idx(:,c),c) - avgC(c);
end

% Sum of products, omitting missing entries
for s1 = 1:nR
    for s2 = s1:nR
        idxr     = ~ (idx(s1,:) | idx(s2,:));
        x1       = X(s1,idxr);
        x2       = X(s2,idxr);
        C(s1,s2) = x1*x2';
        C(s2,s1) = C(s1,s2);
    end
end

% Nearest symmetric postive definite matrix
nearC = nearestSPD(C);

% PCA of that. Outputs have unit norm.
rX = min(nR,nC);
[U,S,~] = svd(nearC,0);
U = U(:,1:rX);
S = S(1:rX,1:rX)./(nR-1);

% Rescale so that the variances match the eigenvalues
U = U/diag(std(U))*sqrt(S);

% Pick one sign
s = diag(sign(U(1,:)));
U = U*s;

% Compute the loadings
iU = pinv(U);
L = zeros(size(iU,1),size(X,2));
for r = 1:size(L,1)
    for c = 1:size(L,2)
        L(r,c) = iU(r,~idx(:,c))*X(~idx(:,c),c);
    end
end

% Compute a recovered version of X, based on the PCA.
Y = bsxfun(@plus,U*L,avgC);
