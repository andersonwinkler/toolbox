function [chi2,df,pval] = chi2test(T)
% Run a chi-squared test in a table with
% arbitrary number of rows or columns.
% 
% [chi2,df,pval] = chi2test(T)
%
% T     : Table
% chi2  : Test statistic
% df    : Degrees of freedom
% pval  : p-value
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Mar/2017
% http://brainder.org

% Margins:
m1 = sum(T,1);
m2 = sum(T,2);

% Expected values:
E = m2*m1/sum(T(:));

% Test statistic:
chi2 = (T - E).^2./E;
chi2 = sum(chi2(:));

% Degrees of freedom:
df = prod(size(T)-1);

% p-value
pval = gammainc(chi2/2,df/2,'upper');
