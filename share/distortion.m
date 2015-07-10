function [G,Gr,V,A,J] = distortion(fsrf1,fsrf2)
% Compute some measures of distortion between surfaces
%
% Usage:
% [G,Gr,V,D,J] = distortion(fsrf1,fsrf2)
%
% fsrf1 : File name for surface 1
% fsrf2 : File name for surface 2. Both surfaces must have similar geometries,
%         i.e., must have the same number of vertices and faces and the same
%         vertex indices for the faces. It's not necessary that both are spheres
%         but the geodesic distance will assume they are.
% G     : Geodesic distances between vertex pairs in surfaces 1 and 2
% Gr    : Simuilar as G, but discounting a global rotation first
% V     : Axis of rotation for the global rotation
% A     : Angle of rotation around the axis for the global rotation
% J     : Ratio between area per face before and after totation. This is
%         an estimate of the determinant of the Jacobian of the transformation
%         of srf1 into srf2.
%
% The surfaces must be in ASCII format (.srf).
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011

% Read the surfaces
[vtx1,fac1] = srfread(fsrf1);
[vtx2,fac2] = srfread(fsrf2);

% Make sure the surfaces share the same geometry
if any(fac1(:)~=fac2(:)),
    error('Surfaces don''t match');
end

% Get an approx to the sphere radius
[ignore,ignore,R1] = cart2sph(vtx1(:,1),vtx1(:,2),vtx1(:,3));
[ignore,ignore,R2] = cart2sph(vtx2(:,1),vtx2(:,2),vtx2(:,3));
R = mean([R1;R2]);

% Get the matrix that rotates 1 to 2
M = vtx1\vtx2;

% Rotate srf2 to approx srf1
vtx2r = vtx2*inv(M);

% Compute geodesic distances
G  = geodist(vtx1,vtx2,R);
Gr = geodist(vtx1,vtx2r,R);

% Rotation axis & angle
[V,D] = eig(M);
V = V(:,(imag(diag(D))==0)); % eigenvector which only eigenvalue is real
A = acos((sum(diag(M))-1)/2);

% Print results in the screen
fprintf('Surfaces %s and %s:\n',fsrf1,fsrf2);
fprintf('- axis: [%g %g %g]\n',V(:));
fprintf('- angle: %g rad (%g deg)\n',A,180*A/pi);

% Compute area per face of srf1
facvtx = [vtx1(fac1(:,1),:) vtx1(fac1(:,2),:) vtx1(fac1(:,3),:)];
facvtx0(:,1:6) = facvtx(:,1:6) - [facvtx(:,7:9) facvtx(:,7:9)];  % Place 3rd vtx at origin
cp = cross(facvtx0(:,1:3),facvtx0(:,4:6),2);                     % Cross product
apf1 = sqrt(sum(cp.^2,2))./2;                                    % Half of the norm

% Compute area per face of srf2
facvtx = [vtx2(fac2(:,1),:) vtx2(fac2(:,2),:) vtx2(fac2(:,3),:)];
facvtx0(:,1:6) = facvtx(:,1:6) - [facvtx(:,7:9) facvtx(:,7:9)];  % Place 3rd vtx at origin
cp = cross(facvtx0(:,1:3),facvtx0(:,4:6),2);                     % Cross product
apf2 = sqrt(sum(cp.^2,2))./2;                                    % Half of the norm

% Estimate the Jacobian
J = apf2./apf1;
