function smoothed = geosmooth(dpvin,distdir,maskfile,fwhm,dpvout)
% Smooth data per vertex (DPV) with a Gaussian kernel
% of a specified width using geodesic distances as supplied.
% 
% Usage:
% smoothed = geosmooth(dpvin,distdir,maskfile,fwhm,dpvout)
% 
% Inputs:
% dpvin    : Input data per vertex file
% distdir  : Directory containing files with distances from every
%            vertex to every other vertex. One text file per vertex
%            space-separated distances to all others (1 line by nV cols).
% maskfile : A dpv file with 1 for vertices to be included, 0 for excluded.
% fwhm     : Full-width at half-maximum, in mm, to smooth the data.
% dpvout   : Output file to be created with the smoothed data
%
% Outputs  : Output smoothed data.
% 
% _____________________________________
% Anderson M. Winkler
% UTRGV
% July/2023
% http://brainder.org

% Load data, mask
[dpv,crd,idx] = dpxread(dpvin);
nV   = size(vtx,1);
mask = dpxread(maskfile);

% Files with the distances
D = ls(sprintf('-v %s',distdir));

% New data to be created
smoothed = zeros(nV,1);

% Convert FWHM to sigma (the standard deviation of the Gaussian filter)
sigma = fwhm./sqrt(8*log(2));

% For each vertex in the data to be smoothed
for v = 1:nV
  
  % Load the distance from this vertex to all others
  d = load(fullfile(distdir,deblank(D{v})));
  d = d./mask;
  
  % Create the filter, as function of distances
  % Note it is not necessary to scale to sum=1
  % because the variance will be rescaled later anyway
  g = exp(-d.^2./2/sigma^2);
  
  % Smooth via multiplication
  smoothed(v) = g*dpv;
  
  % Restore the variance
  smoothed(v) = smoothed(v)./sqrt(sum(g.^2));
end
