function Sabc = extrapvtx3d(T,P,Sp)
% Given 3 distinct points (A, B and C) in the 3D Cartesian space (x,y,z),
% and a 4th point P with an associated scalar Sp, if P lies inside the 
% triangle, returns the incremental, extrapolated scalars for each of the
% 3 vertices, in the form of a column vector (Sabc). If P is outside,
% returns [0 0 0] (so, still incremental).
% 
% Sabc = extrapvtx3d(T,P,Sp)
% 
% T  = Coordinates for the 3 points, in matrix form (see below).
% Sp = Scalar attribute for P.
% P  = Coordinates for the 4th point, to be tested and interpolated.
% 
% The coordinates for the points must be in the following format:
% T = [ xA yA zA ;
%       xB yB zB ;
%       xC yC zC ];
% P = [ xP yP zP ];
% 
% If the points are in the 2D plane (x,y), extrapvtx2d is faster.
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
if single(abs(aA+aB+aC)) == single(aT), % Avoid errors due to tiny decimal places
    Sabc = [aA; aB; aC].*Sp./aT;
else
    Sabc = [0; 0; 0];
end
