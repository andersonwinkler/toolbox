function O = findoutliers(varargin)
% Given an array T, find outliers in either side. These are marked
% in the output O with +1 or -1.
% 
% Usage:
% O = findoutliers(T,dim,k)
% 
% T   : Input multidimensional array.
% dim : Dimension along which outliers are sought.
% k   : Outliers are those farther than Q1+k*IQR or Q3+k*IQR,
%       where IQR is the inter-quartile range. Default k = 1.5.
% O   : Output array indicating who are the outliers.
% 
% Reference:
% Tukey, John W. Exploratory Data Analysis. Addison-Wesley, 1977.
% 
% _____________________________________
% Anderson M. Winkler
% Hospital Israelita Albert Einstein
% Jun/2017
% http://brainder.org

% Input arguments:
narginchk(1,3);
T = varargin{1};
if nargin == 1,
    dim = 1;
    k   = 1.5;
elseif nargin == 2,
    dim = varargin{2};
    k   = 1.5;
elseif nargin == 3,
    dim = varargin{2};
    k   = varargin{3};
end

% Define the quantile indices:
siz = size(T,dim);
q1 = round(siz/4);
q2 = siz/4*2;
q3 = round(siz/4*3);

% Sort the data:
Ts = sort(T,dim);

% Permute dimensions, to put the selected dim as the first:
tmp = 1:numel(size(Ts));
tmp(tmp == dim) = [];
dimpf = [dim tmp];
Ts = permute(Ts,dimpf);
dimr = size(Ts);

% These allow reshaping and permuting dimensions backwards:
dimr(1) = 1;
[~,dimpb] = sort(dimpf);

% Inter-quantile range, upper and lower limits
Q1 = reshape(Ts(q1,:),dimr);
Q3 = reshape(Ts(q3,:),dimr);
IQR = Q3 - Q1;
U = permute(Q3 + k * IQR,dimpb);
L = permute(Q1 - k * IQR,dimpb);

% Thresholds:
idxU = bsxfun(@gt,T,U);
idxL = bsxfun(@lt,T,L);

% Outliers marked:
O = idxU - idxL;
