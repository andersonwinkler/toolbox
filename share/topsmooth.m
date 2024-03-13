function T = topsmooth(siz,fwhm)
% Create a Toeplitz matrix that smooths with a gaussian filter
% if of a pre-specified FWHM. Padding is circular.
%
% Usage:
% T = topsmooth(siz,fwhm)
% 
% Inputs:
% - siz   : Length of the signal to be smoothed
% - fwhms : Full-width at half maximum of the filter
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

if numel(siz) == 1 && numel(fwhm) == 1

    % Behave slightly differently if siz is odd or even
    iseven = ~ rem(siz,2);

    % Compute the filter and its scaling factor to restore the variance
    [f,s] = gauss(siz+iseven,fwhm);
    f = f(1:end-iseven)*s;

    % Assemble the matrix
    frow = circshift(f,-floor(siz/2));
    T = toeplitz(frow,frow);
else
    error('Toeplitz matrix can only be computed for 1-D signals');
end
