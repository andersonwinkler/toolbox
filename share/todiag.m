function [Y,idxR,idxC] = todiag(X,order,absflag)
% Takes a matrix X and ensures it has a dominant diagonal based
% on largest or smallest values.
% 
% Usage:
% [Y,idxR,idxC] = todiag(X,order,absflag)
% 
% Inputs: 
% X         : Matrix to have its rows/cols reordered
% order     : (Optional) Whether largest ('ascend') or 
%             smallest ('descend') should be used to populate
%             the main diagonal. Default is 'descend'.
% absflag   : (Optional) Whether values are onsidered as they are
%             in X, or as abs(X). Default is false.
%
% Outputs:
% Y         : This is X reordered.
% idxR,idxC : These are the indices such as Y = X(idxR,idxC).
% S         : Sign of the diagonal of Y. If X is a correlation or
%             covariance matrix, then S can be used to change the
%             signs of one set of variables after its columns have
%             been reordered with idxR or idxC.
% 
% Reordering is based on largest or smallest values in X (using 
% 'descend' or 'ascend' respectively), and represented in the 
% by the returned indices idxR and idxC, that is, Y = X(idxR,idxC).
% The absflag indicates whether diagonalization should be based
% on X (absflag=false) or abs(X) (absflag=true). In the latter
% case, you can use the S = sign(diag(Y)) to invert the signs of
% either rows or columns, such that the diagonal will be guaranteed
% to also be positive (you may also flip the the original variables
% in case X is a correlation or covariance matrix.
% 
% _____________________________________
% Anderson M. Winkler
% National Institute of Mental Health
% May/2020
% http://brainder.org

narginchk(1,3);
if nargin < 2
    order = 'descend';
end
if nargin < 3
    absflag = true;
end
if strcmpi(order,'descend')
    ext = @max;
    opp = @min;
elseif strcmpi(order,'ascend')
    ext = @min;
    opp = @max;
else
    error('Unknown order. Use ''ascend'' or ''descend''.');
end
if absflag,
    Sx = sign(X);
    X  = abs(X);
end

[nR,nC] = size(X);
N       = min(nR,nC);
Y       = X;
idxR    = zeros(1,N);
idxC    = zeros(1,N);
oX      = opp(Y(:));
for n = 1:N
    eX      = ext(Y(:));
    [r,c]   = find(Y == eX);
    idxR(n) = r;
    idxC(n) = c;
    Y(r,:)  = oX;
    Y(:,c)  = oX;
end
idxR  = [idxR setdiff(1:nR,idxR)];
idxC  = [idxC setdiff(1:nC,idxC)];
Y     = X(idxR,idxC);
if absflag,
    Y = Y.*Sx(idxR,idxC);
end