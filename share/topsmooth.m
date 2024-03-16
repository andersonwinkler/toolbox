function T = topsmooth(siz,fwhm)
% Create a Toeplitz matrix that smooths with a gaussian filter
% if of a pre-specified FWHM. Padding is circular.
%
% Usage:
% T = topsmooth(siz,fwhm)
%
% Inputs:
% - siz   : Length of the signal to be smoothed
% - fwhm  : Full-width at half maximum of the filter
%           Can be a scalar or a vector of size (siz). If a vector,
%           allows for non-stationarity (each row, a different FWHM).
%
% Outputs:
% - T     : Toeplitz matrix
%
% To use this matrix to smooth, use:
% sX = T*X; % where X is a column vector (or matrix) with the signal
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Mar/2024
% http://brainder.org

if numel(siz) > 1
    error('Toeplitz matrix can only be computed for 1-D signals.');
end

% Behave slightly differently if siz is odd or even
iseven = ~ rem(siz,2);

% For each filter size, create the filter, the scaling factor needed
% to restore the variance. Then shift as needed, and store for later use
frow = cell(numel(fwhm),1);
for f = 1:numel(fwhm)
    [filt,scal] = gauss(siz+iseven,fwhm(f));
    filt        = filt(1:end-iseven)*scal;
    frow{f}     = circshift(filt,-floor(siz/2)+f-1);
end

% Now assemble the matrix
if numel(fwhm) == 1
    T = toeplitz(frow{1},frow{1});
elseif numel(fwhm) == siz
    T = vertcat(frow{:});
else
    error('Wrong size for the FWHM parameter.');
end
