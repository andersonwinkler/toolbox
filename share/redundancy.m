function [Ru,Rv,Run,Rvn] = redundancy(Y,X,U,V) 
% Compute the CCA redundancy index.
% 
% Inputs:
% 
% Y, X     : Original data, before CCA.
% U, V     : Respective canonical variables, after CCA.
%
% Outputs:
% 
% Ru,  Rv  : Redundancy indices, for left and right sides.
% Run, Rvn : Normalized redundancy indices (such that their sum is 1).
%
% References:
% 
% * Stewart D, Love W. A general canonical correlation index.
%   Psychological Bulletin 1968;70(3, Pt.1):160?3.
% * https://brainder.org/2019/12/27/redundancy-in-canonical-correlation-analysis/
% 
% _____________________________________
% Thomas Wassenaar and Anderson Winkler
% Univ. of Oxford / Natl. Inst. of Health
% June/2020
% http://brainder.org

% The canonical loadings
At = bcorr(Y,U);
Bt = bcorr(X,V);

% Variance explained by the corresponding canonical variables 
explU = mean(At.^2,1);
explV = mean(Bt.^2,1);

% Proportion of variance of one canonical variable explained by the 
% corresponding canonical variable in the other side
rsq = dcorr(U,V).^2;

% Compute the redundancy index for each canonical variable
Ru = explU.*rsq;
Rv = explV.*rsq;

% Global redundancy scaled to unity
Run = Ru./sum(Ru);
Rvn = Rv./sum(Rv);
end

% =================================================================
function C = dcorr(X,Y)
% Efficiently compute the diagonal of the correlation matrix
X = bsxfun(@minus,  X,mean(X,1));
Y = bsxfun(@minus,  Y,mean(Y,1));
X = bsxfun(@rdivide,X,std(X,1));
Y = bsxfun(@rdivide,Y,std(Y,1));
C = sum(X.*Y,1)/size(X,1);
end

% =================================================================
function C = bcorr(X,Y)
% Basic equivalent to "corr" that works in Octave and in Matlab
% versions that don't have broadcasting.
X = bsxfun(@minus,  X,mean(X,1));
Y = bsxfun(@minus,  Y,mean(Y,1));
X = bsxfun(@rdivide,X,std(X,1));
Y = bsxfun(@rdivide,Y,std(Y,1));
C = (X'*Y)/size(X,1);
end
