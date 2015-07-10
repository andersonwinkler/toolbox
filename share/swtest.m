function [W,pval] = swtest(Y)
% Run the Shapiro-Wilk normality test. It is valid for n<=2000, possibly 
% for samples up to 5000.
% 
% Usage:
% [W,pval] = swtest(Y);
% 
% Inputs:
% Y     : 2-D matrix containing the original data.
%         The test is applied columnwise.
% 
% Outputs:
% W     : Shapiro-Wilk W statistic.
% pval  : p-value for the W statistic, based on a normal approximation.
% 
% References:
% * Shapiro SS, Wilk MB. An analysis of variance test for normality
%   (complete samples). Biometrika. 1965; 52(3-4):591-611.
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

% Sample mean
mY = mean(Y);

% Blom scores (Gunnar Blom. Statistical Estimates and Transformed Beta-Variables. New York, John Wiley & Sons, 1958)
m  = norminv(((1:n)'-3/8)./(n+1/4));

% Weights for the Shapiro-Wilk test (approximation as in Royston 1992)
u = n^(-1/2);
c = (m'*m)^(-1/2).*m; % Weights as in the Shapiro-Francia test
coef_an = [-2.706056 4.434685 -2.071190 -0.147981 0.221157 c(n)];
a    = zeros(n,1);
a(n) = polyval(coef_an,u);
a(1) = -a(n);  % antisymmetric
if n <= 5,
    phi      = (m'*m - 2*m(n)^2) / (1 - 2*a(n)^2);
    a(2:n-1) = phi^(-1/2)*m(2:n-1);
else
    coef_an1 = [-3.582633 5.682633 -1.752461 -0.293762 0.042981 c(n-1)];
    a(n-1)   = polyval(coef_an1,u);
    a(2)     = -a(n-1);   % antisymmetric
    phi      = (m'*m - 2*m(n)^2 - 2*m(n-1)^2) / (1 - 2*a(n)^2 - 2*a(n-1)^2);
    a(3:n-2) = phi^(-1/2)*m(3:n-2);
end
a = repmat(a,[1 size(Y,2)]);

% Shapiro-Wilk W statistic
W = sum(a.*Y).^2 ./ sum((Y - repmat(mY,[n 1])).^2);
W(W > 1) = 0; % Prevent issues with constant and other awkward data, which should be rejected as non-normal (i.e. low W, low p)

% Normalizing transformation  (Table 1 in Royston, 1993)
if n <= 11,

    % Polynomial coefficients
    gam = [ 0.459    -2.273];
    mu  = [-0.0006714 0.025054 -0.39978 0.5440];
    lsd = [-0.0020322 0.062767 -0.77857 1.3822];

    % Parameters of the transformed data
    gam = polyval(gam,n);
    mu  = polyval(mu, n);
    lsd = polyval(lsd,n);
    sd  = exp(lsd);

    % Transformation function
    g   = -log(gam-log(1-W));

else

    % Polynomial coefficients
    mu  = [0.0038915 -0.083751 -0.31082 -1.5851];
    lsd = [0.0030302 -0.082676 -0.4803];

    % Parameters of the transformed data
    mu  = polyval(mu, log(n));
    lsd = polyval(lsd,log(n));
    sd  = exp(lsd);

    % Transformation function
    g   = log(1-W);
end

% Nomalization
Z = (g-mu)/sd;

% Significance
pval = 1-normcdf(Z);

