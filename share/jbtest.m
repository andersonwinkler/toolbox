function [LM,pval] = jbtest(X)
% Run the Jarque-Bera normality test.
% 
% Usage:
% [LM,pval] = jbtest(X);
% 
% Inputs:
% X     : 2-D matrix containing the original data.
%         The test is applied columnwise.
% 
% Outputs:
% LM    : Jarque-Bera LM statistic.
% pval  : p-value for the LM statistic, based on a Chi^2 approximation.
% 
% References:
% * Jarque CM, Bera AK. Efficient tests for normality, homoscedasticity
%   and serial independence of regression residuals. Economics Letters.
%   1980; 6(3):255-259.
% * Jarque CM, Bera AK. A test for normality of observations and regression
%   residuals. International Statistical Review. 1987; 55(2):163-172.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jun/2011
% http://brainder.org

% Due to notation issues, what is here being
% called as b1 and beta1 are in fact, sqrt(b1) and sqrt(beta1)
% The implementation below uses the notation in D'Agostino et al (1990)
% and Jarque & Bera (1987).

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

% === [JARQUE-BERA TEST] ==================================================

LM = n*((b1.^2)/6 + ((b2-3).^2)/24);
df   = 2;           % degrees of freedom is always 2
pval = 1-chi2cdf(LM,df);

% That's it!