function Sp = interpvtx3d(T,Sabc,P)
% Given 3 distinct points (A, B and C) in the 3D Cartesian space (x,y,z),
% each one containing an assigned scalar value (Sabc), returns NaN if a 
% given 4th point P lies outside the triangle determined by ABC,
% or return the interpolated scalar if it lies inside the triangle.
% 
% Sp = interpvtx3d(T,Sabc,P)
% 
% T    = Coordinates for the 3 points, in matrix form (see below).
% Sabc = Scalar attribute for each vertex (not to be confused with
%        the vertices coordinates). Sabc = [Sa Sb Sc]' (column vector)
% P    = Coordinates for the 4th point, to be tested and interpolated.
% 
% The coordinates for the points must be in the following format:
% T = [ xA yA zA ;
%       xB yB zB ;
%       xC yC zC ];
% P = [ xP yP zP ];
% 
% If the points are in the 2D plane (x,y), interpvtx2d is faster.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jan/2011


% Compute the area for the the triangle ABC and for the 3 subtriangles,
% PAB, PAC, PBC.
aT = triarea(T);
tA = T; tA(1,:) = P; aA = triarea(tA);
tB = T; tB(2,:) = P; aB = triarea(tB);
tC = T; tC(3,:) = P; aC = triarea(tC);

% P is inside if the sum of the areas for the 3 subtriangles is the same 
% as the big one (ABC). These smaller areas are also the weights for 
% the interpolation, without the need to calculate the plane equation.
if single(abs(aA+aB+aC)) == single(aT); % Avoid errors due to tiny decimal places
    Sp = [aA aB aC]./aT * Sabc; % Sabc has to be a COLUMN vector. This is not explicitly tested in the code to allow speed.
else
    Sp = NaN;
end
