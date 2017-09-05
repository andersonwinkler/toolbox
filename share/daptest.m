 function [K2,pval,b1,Zb1,pvalS,b2,Zb2,pvalK] = daptest(X)
% Run the D'Agostino-Pearson normality test, as well as skewness and
% kurtosis tests.
% 
% Usage:
% [K2,pval,b1,Zb1,pvalS,b2,Zb2,pvalK] = daptest(X);
% 
% Inputs:
% X     : 2-D matrix containing the original data.
%         The test is applied columnwise.
% 
% Outputs:
% K2    : D'Agostino-Pearson's K^2 statistic.
% pval  : p-value for the K^2 statistic, based on a Chi^2 approximation.
% b1    : Skewness
% Zb1   : Z-transformed skewness statistic.
% pvalS : p-value for the skewness statistic, based on normal approximation.
% b2    : Kurtosis
% Zb2   : Z-transformed kurtosis statistic.
% pvalS : p-value for the kurtosis statistic, based on normal approximation.
% 
% References:
% The skewness test is based on:
% * D'Agostino RB. Transformation to normality of the null distribution
%   of g1. Biometrika. 1970; 57(3):679-81.
% The kurtosis test is based on:
% * Anscombe FJ, Glynn WJ. Distribution of the kurtosis statistic b2 for
%   normal samples. Biometrika. 1983;70(1):227-34.
% The combined K^2 statistic was introduced in:
% * D'Agostino R, Pearson ES. Tests for departure from normality. Empirical
%   results for the distributions of b2 and sqrt(b1). Biometrika.
%   1973; 60(3):613-22.
% The implementation is based on:
% * D'Agostino RB, Belanger A, D'Agostino Jr RB. A suggestion for using
%   powerful and informative tests of normality. The American Statistician.
%   1990;44(4):316-21.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jun/2011
% http://brainder.org

% The equations below and variable names follow the paper by
% D'Agostino et al (1990). Due to notation issues, what is here being
% called as b1 and beta1 are in fact, sqrt(b1) and sqrt(beta1)
% in the original paper.

% Sample size
n = size(X,1);

% === [SAMPLE MOMENTS] ====================================================

% Sample mean
me  = sum(X)/n;                                                  % [Eqn  7]
rme = repmat(me,[n 1]);

% Compute the 2nd, 3rd and 4th moments
m = cell(4,1);
for k = 2:4,
    m{k} = (sum((X-rme).^k))./n;                                 % [Eqn  6]
end

% Standardise the 3rd and 4th moments
b1    = m{3}./(m{2}.^(3/2));             % Skewness              % [Eqn  4]
b2    = m{4}./(m{2}.^2);                 % Kurtosis              % [Eqn  5]

% === [SKEWNESS TEST] =====================================================

Y     = b1.*sqrt((n+1)*(n+3)/6/(n-2));                           % [Eqn  8]
beta2 = 3*(n^2+27*n-70)*(n+1)*(n+3)/(n-2)/(n+5)/(n+7)/(n+9);     % [Eqn  9]
W2    = -1+sqrt(2*(beta2-1));                                    % [Eqn 10]
delta = 1/sqrt(log(sqrt(W2)));                                   % [Eqn 11]
alpha = sqrt(2/(W2-1));                                          % [Eqn 12]
Zb1   = delta.*log(Y./alpha+sqrt((Y./alpha).^2+1));              % [Eqn 13]
pvalS = 1-normcdf(Zb1);      % p-value for skewness

% === [KURTOSIS TEST] =====================================================

Eb2   = 3*(n-1)/(n+1);                                           % [Eqn 14]
varb2 = 24*n*(n-2)*(n-3)/(n+1)/(n+1)/(n+3)/(n+5);                % [Eqn 15]
x     = (b2-Eb2)./sqrt(varb2);                                   % [Eqn 16]
beta1 = 6*(n^2-5*n+2)/(n+7)/(n+9)* ...
    sqrt(6*(n+3)*(n+5)/n/(n-2)/(n-3));                           % [Eqn 17]
A     = 6+8/beta1*(2/beta1+sqrt(1+4/(beta1^2)));                 % [Eqn 18]
Zb2   = ((1-2/9/A)- ...
    ((1-2/A)./(1+x.*sqrt(2/(A-4)))).^(1/3))./sqrt(2/9/A);        % [Eqn 19]
pvalK = 1-normcdf(Zb2);      % p-value for kurtosis

% === [OMNIBUS TEST] ======================================================

K2   = Zb1.^2 + Zb2.^2;                                          % [Eqn 20]
df   = 2;                    % degrees of freedom is always 2
pval = 1-chi2cdf(K2,df);     % p-value global

% That's it!