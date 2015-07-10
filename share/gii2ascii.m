function gii2ascii(filein,fileout)
% Converts a surface or scalar map on a surface
% from GIFTI to ASCII.
% Surfaces will be in .srf format, and scalar maps
% will be on .dpx format.
%
% Usage:
% gii2ascii(filein,fileout)
%
% The GIFTI toolbox must be installed.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2014
% http://brainder.org

g = gifti(filein);
if isfield(g,'vertices') && isfield(g,'faces') && isfield(g,'mat');
    vtx = horzcat(g.vertices,ones(size(g.vertices,1),1))*g.mat;
    srfwrite(vtx(:,1:3),g.faces-1,fileout);
elseif isfield(g,'cdata'),
    dpxwrite(fileout,g.cdata,zeros(size(g.cdata,1),3),0:size(g.cdata,1)-1);
end
