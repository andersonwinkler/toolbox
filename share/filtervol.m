function Vf = filtervol(V,F)
% Apply a filter to a N-D array. Filtering is
% performed in the frequency domain.
%
% Usage:
% Vf = filtervol(V,F)
% 
% - V  : A N-D array.
% - F  : A filter of the same size as V.
% - Vf : Filtered array.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Jul/2012

% Check arguments
if ~ all(size(V) == size(F))
    error('Inputs need to be of the same size.');
end

% Filter in the frequency domain
iV = fftn(V);
iF = fftn(F);
Vf = real(ifftn(iV.*abs(iF)));

% Restore the energy of the signal
en = iF.*conj(iF);
Vf = Vf/sqrt(sum(en(:))/numel(en));
