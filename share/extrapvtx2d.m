function Sabc = extrapvtx2d(T,P,Sp)
% Given 3 distinct points (A, B and C) in the 2D Cartesian space (x,y),
% and a 4th point P with an associated scalar Sp, if P lies inside the 
% triangle, returns the incremental, extrapolated scalars for each of the
% 3 vertices, in the form of a column vector (Sabc). If P is outside,
% returns [0 0 0] (so, still incremental).
% 
% Sabc = interpvtx2d(T,P,Sp)
% 
% T  = Coordinates for the 3 points, in matrix form (see below).
% Sp = Scalar attribute for P.
% P  = Coordinates for the 4th point, to be tested and interpolated.
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
aU = abs([aA; aB; aC]);

% Weight appropriately if P is inside, or return NaN if outside
if single(sum(aU)) == single(abs(aT)), % Avoid errors due to tiny decimal places
    Sabc = aU .* Sp ./ aT;
else
    Sabc = [0; 0; 0];
end
