function Y = randlap(siz,m,s)
% Generate Laplace-distributed random numbers.
%
% Usage:
% Y = randlap(siz,m,s)
% 
% siz : Size of the output array.
% m   : Mean of the distribution (Default = 0).
% s   : Standard deviation (Default = 1).
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2014
% http://brainder.org

if nargin == 1,
    m = 0;
    s = 1;
elseif nargin == 2,
    s = 1;
end
U = rand(siz)-0.5;
b = s/sqrt(2);
Y = m - b*sign(U).*log(1-2*abs(U));
