function [Wp,pval] = sftest(Y)
% Run the Shapiro-Francia normality test. It is valid for n<=5000.
% 
% Usage:
% [Wp,pval] = sftest(Y);
% 
% Inputs:
% Y     : 2-D matrix containing the original data.
%         The test is applied columnwise.
% 
% Outputs:
% Wp    : Shapiro-Francia W' statistic.
% pval  : p-value for the W' statistic, based on a normal approximation.
% 
% References:
% * Shapiro SS, Francia R. An approximate analysis of variance test for
%   normality. Journal of the American Statistical Association.
%   1972;67(337):215-6.
% * Royston P. A toolkit for testing for non-normality in complete and
%   censored samples. The Statistician. 1993; 42(1):37.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jun/2011
% http://brainder.org

% Sort the data
Y = sort(Y);

% Sample size
n  = size(Y,1);
if n < 4,
    error('Sample is too small.');
end

% Sample size
n  = size(Y,1);

% Sample mean
mY = mean(Y);

% Blom scores
m  = norminv(((1:n)'-3/8)./(n+1/4));

% Weights
c  = (m'*m)^(-1/2).*m;
c = repmat(c,[1 size(Y,2)]);

% Shapiro-Francia W' statistic
Wp = sum(c.*Y).^2 ./ sum((Y - repmat(mY,[n 1])).^2);

% Polynomial coefficients
mu = [ 1.0521 -1.2725];
sd = [-0.26758 1.0308];

% Parameters of the transformed data
v  = log(n);
mu = polyval(mu, log(v)-v);
sd = polyval(sd,log(v)+2/v);

% Transformation function
g = log(1-Wp);

% Nomalization
Z = (g-mu)/sd;

% Significance
pval = 1-normcdf(Z);

