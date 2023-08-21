function varargout = dcog(varargin)
% Given a matrix of pairwise Euclidean distances between
% points, return a vector of distances of every point to
% the center of gravity of all points. Points can receive
% different weights.
% 
% Usage:
% dist = dcog(dij,w)
% 
% - dij  : NxN symmetric matrix with distances of every
%          point to every other point. Diagonal is assumed 0.
% - w    : (Optional) Nx1 vector of weights.
%          Larger values imply lighter points.
% - dist : Nx1 vector of distances from each point to the
%          center of gravity.
% 
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Aug/2023
% http://brainder.org

narginchk(1,2);
if nargin >= 1
    dij = varargin{1};
end
N = size(dij,1);
if nargin >= 2
    w = varargin{2};
else
    w = ones(N,1);
end
w = w./sum(w);

sij = dij.^2;
cte = sum(sij(:))/2/N/N;
sic = sum(sij,2).*w - cte;
dic = sqrt(sic);
varargout{1} = dic;