function P = mpcdf(x,y,sigsq)
% Cumulative density function of the Marcenko-Pastur
% distribution of the eigenvalues of random matrices.
% 
% Consider a matrix X (never seen by this function) with random complex
% entries with zero mean and variance sigsq, and with a ratio of
% number of rows by number columns converging to y. Let x be an
% eigenvalue of the covariance matrix of X. Then P is the cdf of x.
% 
% Usage:
% P = mpcdf(x,y,sigsq)
% 
% Inputs:
% x     : An eigenvalue of the matrix X.
% y     : Ratio #rows/#cols of the random matrix from which the
%         eigenvalue is characteristic.
% sigsq : Variance of the entries of X.
% 
% Outputs:
% P     : Cumulative density (i.e., 1-p, where p is the p-value).
% 
% References:
% * Bai ZD. Methodologies in spectral analysis of large dimensional
%   random matrices, a review. Statistica Sinica 1999;9(3):611?77. 
% * Marcenko VA, Pastur LA. Distribution of eigenvalues for some 
%   sets of random matrices. Math USSR Sb 1967;1(4):457?83. 
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Mar/2021
% http://brainder.org

a = sigsq.*(1 - y^.5).^2;
b = sigsq.*(1 + y^.5).^2;
P = quad(@(x)mppdf(x,y,sigsq,a,b),a,x); % 'quad' if faster than 'integral'
P = P + (y > 1).*(1 - 1./y); % add point mass if y > 1

function p = mppdf(x,y,sigsq,a,b)
% PDF of the Marcenko-Pastur distribution.
% This is Equation 2.12 of Bai (1999)
p = sqrt((b-x).*(x-a))./x./y./sigsq/2/pi .* (a<=x) .* (x<=b);
