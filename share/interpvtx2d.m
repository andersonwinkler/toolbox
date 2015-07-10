function Sp = interpvtx2d(T,Sabc,P)
% Given 3 distinct points (A, B and C) in the 2D Cartesian space (x,y),
% each one containing an assigned scalar value (Sabc), returns NaN if a 
% given 4th point P lies outside the triangle determined by ABC
% or return the interpolated scalar if it lies inside the triangle.
% 
% Sp = interpvtx2d(T,Sabc,P)
% 
% T    = Coordinates for the 3 points, in matrix form (see below).
% Sabc = Scalar attribute for each vertex (not to be confused with
%        the vertices coordinates). Sabc = [Sa Sb Sc]' (column vector)
% P    = Coordinates for the 4th point, to be tested and interpolated.
% 
% The coordinates for the points must be in the following format:
% T = [ xA yA ;
%       xB yB ;
%       xC yC ];
% P = [ xP yP ];
% 
% If the points are in the 3D plane (x,y,z), use interpvtx3d.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jan/2011


% Complete Z with ones so the determinant method becomes valid
T(:,3) = 1; P(3) = 1;

% The formula for the area is det(T)/2 if the 3 points are on (x,y) plane
% (so, no Z). The division by 2 can be omitted for efficiency, as we
% don't need the actual area value, but the areal proportions in the
% subsequent step.
aT = det(T);
tA = T; tA(1,:) = P; aA = det(tA);
tB = T; tB(2,:) = P; aB = det(tB);
tC = T; tC(3,:) = P; aC = det(tC);
aU = abs([aA aB aC]);

% Weight appropriately if P is inside, or return NaN if outside
if single(sum(aU)) == single(aT), % Avoid errors due to tiny decimal places
    Sp = aU ./ aT * Sabc; % Sabc has to be a COLUMN vector. This is not explicitly tested in the code to allow speed.
else
    Sp = NaN;
end
