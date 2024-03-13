function [f,s,F] = gauss(siz,fwhms)
% Create a gaussian function of up to 4 dimensions.
%
% Usage:
% [f,s,F] = gauss(siz,fwhms)
% 
% Inputs:
% - siz   : A vector with the dimensions. E.g for a 3-D filter,
%           of size 20x20x30, use [20 20 30]
% - fwhms : Respective filter widths.
%
% Outputs:
% - f     : Gaussian filter.
% - s     : Scale factor to restore the variance of the filtered
%           signal (multiply to restore).
% - F     : FFT of the filter.
% 
% To use this filter to smooth, use:
% simg = real(ifftn(fftn(img).*abs(F)))*s;
% 
% _____________________________________
% Anderson M. Winkler
% LABIEM / CPGEI / UTFPR
% Sep/2006 (first version)
% Jul/2012 (this version)
% http://brainder.org

% Take and check arguments
siz  = siz(:);
fwhms = fwhms(:);
if numel(siz) ~= numel(fwhms)
    if numel(fwhms) ~= 1
        error('The domain size and FWHM need to have the same length');
    else
        fwhms = ones(size(siz))*fwhms;
    end
end
D = numel(siz);
if D > 4
    error('Too many dimensions');
end

% Make the filter for each dimension (one axis only)
sigmas = fwhms./sqrt(8*log(2));
lims   = (siz-1)/2;
ax     = cell(D,1);
g      = ax;
for d = 1:D
    ax{d} = -lims(d):lims(d);
    g{d} = exp(-(ax{d}.^2)./(2*(sigmas(d)^2)));
end

% 1D
if D >= 1
    f = g{1};
end

% 2D
if D >= 2
    f = f'*g{2};
end

% 3D
if D >= 3
    tmp = zeros(siz(1:3)');
    for s = 1:siz(3)
        tmp(:,:,s) = g{3}(s).*f;
    end
    f = tmp;
end

% 4D
if D >= 4
    tmp = zeros(siz(1:4)');
    for v = 1:siz(4)
        tmp(:,:,:,v) = g{4}(v).*f;
    end
    f = tmp;
end

% Normalise vals to unit sum
f = f/sum(f(:));

% Some sanity check (if the filter is too narrow, for instance, this can
% happen, but hard to predict depending on the dimensions, etc). So, let's
% go with a catch-all.
if any(isnan(f(:)))
    error('The input parameters are causing NaN in the filter.');
end

% Compute the scaling factor
if nargout >= 2
    
    % The restoration is based on expectations, hence no need
    % for the actual signal.
    F   = fftn(f);         % Fourier
    pF  = F.*conj(F);      % Power spectrum
    acf = real(ifftn(pF)); % Autocorrelation function
    
    % Scale factor (inverted just so that it becomes multiplicative)
    s = 1/sqrt(acf(1));
end
