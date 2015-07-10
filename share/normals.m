function [vnormals,fnormals] = normals(vtx,fac)
% Compute the face and vertex normals for a mesh of
% triangular faces
% 
% Usage
% [vnormals,fnormals] = normals(vtx,fac)
% 
% - vtx : Vertex coordinates
% - fac : Face indices (triangles only)
% - vnormals : Vertex normals (normalised to unity)
% - fnormals : Face normals (normalised to unity)
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Oct/2011


% Number of vertices and faces
nV = size(vtx,1);
nF = size(fac,1);

% Compute the face normals (cross product of 2 edges).
% Note the orientation...
facvtx = [vtx(fac(:,1),:) vtx(fac(:,2),:) vtx(fac(:,3),:)];
facvtx0(:,1:6) = facvtx(:,1:6) - [facvtx(:,7:9) facvtx(:,7:9)];  % Place 3rd vtx at origin
fnormals = cross(facvtx0(:,1:3),facvtx0(:,4:6),2);               % Cross product

% face-vtx correspondence
facidx = repmat((1:nF)',[1 3]);
[~,idx] = sort(fac(:));
facidx = facidx(idx);
facidx = reshape(facidx,[4 nV]);

% Use a for-loop to avoid memory consumption
% It's slower, though
vnormals = zeros(nV,3);
for v = 1:nV,
    vnormals(v,:) = mean(fnormals(facidx(:,v),:));
end

% Normalise the normals
nrm = sqrt(sum(vnormals.^2,2));
vnormals = vnormals./repmat(nrm,[1 3]);
if nargout == 2,
    nrm = sqrt(sum(fnormals.^2,2));
    fnormals = fnormals./repmat(nrm,[1 3]);
end

