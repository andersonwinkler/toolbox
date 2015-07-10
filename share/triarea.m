function Area = triarea(T)
% Compute the area of a triangle in the 3D space.
%
% Area = triarea(T)
%
% Coordinates for the points A, B and C are given in matrix form:
% T = [ xA yA zA ;
%       xB yB zB ;
%       xC yC zC ]
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Nov/2010

Txy = T; Txy(:,3) = [1 1 1];
Tyz = T; Tyz(:,1) = [1 1 1];
Txz = T; Txz(:,2) = [1 1 1];
Area = .5 * sqrt(det(Txy)^2 + det(Tyz)^2 + det(Txz)^2);
