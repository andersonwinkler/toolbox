function retessellate(varargin)
% Retessellate a surface.
%
% Usage:
% retessellate(srf1,srf2,srf3,srf4)
%
% - srf1 : Source sphere.
% - srf2 : Target sphere (typically a common grid).
% - srf3 : Surface to be retessellated (typically not a sphere).
% - srf4 : Retessellated surface (to be created).
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jan/2012
% http://brainder.org

% Do OCTAVE stuff
if exist('argv','builtin') && ~ exist('varargin','var')
    
    % Get the inputs
    varargin = argv();
    nargin  = numel(varargin);
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Retessellate a surface.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('retessellate srf1 srf2 srf3 srf4\n');
        fprintf('\n');
        fprintf('- srf1 : Source sphere.\n');
        fprintf('- srf2 : Target sphere (typically a common grid).\n');
        fprintf('- srf3 : Surface to be retessellated (typically not a sphere).\n');
        fprintf('- srf4 : Retessellated surface (to be created).\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Jan/2012\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Default margin
marg =  0.05;

% Check number of arguments
nargin  = numel(varargin);
if nargin ~= 4,
    error('Incorrect number of arguments.');
end

% Load the surfaces
[vtx1,fac1] = srfread(varargin{1});
nF1 = size(fac1,1);
[vtx2,fac2] = srfread(varargin{2});
nV2 = size(vtx2,1);
[vtx3,fac3] = srfread(varargin{3});

% Some sanity check
if any(fac1(:) ~= fac3(:)),
    error(['The source surface and the surface to be retessellated\n',...
        '      must have the same geometry.%s'],'')
end

% Where the result is going to be stored
vtx4 = zeros(nV2,3);

% Vertices' coords per face
facvtx1 = [vtx1(fac1(:,1),:) vtx1(fac1(:,2),:) vtx1(fac1(:,3),:)];

% Face barycenter
cbary = [mean(facvtx1(:,[1 4 7]),2) ...  % Cartesian coordinates
    mean(facvtx1(:,[2 5 8]),2) mean(facvtx1(:,[3 6 9]),2)];
[sbary(:,1),sbary(:,2),sbary(:,3)] = ... % Spherical coordinates
    cart2sph(cbary(:,1),cbary(:,2),cbary(:,3));

% Pre-calculated sines and cosines of azimuth and elevation:
sinA = sin(sbary(:,1));  sinE = sin(sbary(:,2));
cosA = cos(sbary(:,1));  cosE = cos(sbary(:,2));

% Pre-calculated rotation matrices
rotM = [cosA.*cosE sinA.*cosE sinE -sinA cosA zeros(nF1,1) -cosA.*sinE -sinA.*sinE cosE];

% Include a random angle around X, so that it prevents an issue with
% perfectly horizontal edges
rndangX = rand(1)*pi;
sinX = sin(rndangX);
cosX = cos(rndangX);
rotM = [ rotM(:,1:3) ...
    rotM(:,4)*cosX+rotM(:,7)*sinX rotM(:,5)*cosX+rotM(:,8)*sinX rotM(:,6)*cosX+rotM(:,9)*sinX ...
    rotM(:,7)*cosX-rotM(:,4)*sinX rotM(:,8)*cosX-rotM(:,5)*sinX rotM(:,9)*cosX-rotM(:,6)*sinX ];

% Pre-calculated min and max for each face and bounding box
minF = [min(facvtx1(:,[1 4 7]),[],2) ...
    min(facvtx1(:,[2 5 8]),[],2) min(facvtx1(:,[3 6 9]),[],2)];
maxF = [max(facvtx1(:,[1 4 7]),[],2) ...
    max(facvtx1(:,[2 5 8]),[],2) max(facvtx1(:,[3 6 9]),[],2)];
b = repmat(max((maxF-minF),[],2),[1 3]) * marg;  % <= marg enters here
minF = minF-b;  maxF = maxF+b;

% For each source face
for f = 1:nF1;
    
    % Face vertices and associated scalars ("weights") from Curvature 1
    vidx = fac1(f,:);
    Fvtx = vtx1(vidx,:);
    
    % Candidate vertices
    Cidx = all(vtx2 >= repmat(minF(f,:),[nV2 1]) & ...
        vtx2 <= repmat(maxF(f,:),[nV2 1]), 2);
    Cvtx  = vtx2(Cidx,:);
    Cidxi = find(Cidx); % Though slower than logical indexing, it'll be necessary below...
    
    % Rotate the face vertices and candidate vertices
    Avtx = [Fvtx; Cvtx] * reshape(rotM(f,:),[3 3]);
    
    % Here is the main difference in relation to the 'distributive' method.
    % Instead of extrapolate (split) a point in the source into the face
    % vertices of the target for increments, find in the source
    % face what are the Target vertices that lie inside it and do a
    % barycentric interpolation (not to be confused with the rotation
    % of the barycenter, used some lines above to put the face under
    % analysis near the sphere equator and the meridian zero).
    
    % Convert to azimuthal gnomonic
    Gvtx = ones(size(Avtx));           % The 3rd col will remain full of ones
    Gvtx(:,1) = Avtx(:,2)./Avtx(:,1);  % Tangent of the angle on the XY plane
    Gvtx(:,2) = Avtx(:,3)./Avtx(:,1);  % Tangent of the angle on the XZ plane
    T  = Gvtx(1:3,:);                  % Face coords for the test below
    aT = det(T);                       % Face total area (2x the area, actually)
    
    % For every candidate vertex
    for v = 1:numel(Cidxi),
        
        % Compute the areas for the subtriangles (2x area actually)
        tA = T;  tA(1,:) = Gvtx(v+3,:);  aA = abs(det(tA));  % Subtriangle A
        tB = T;  tB(2,:) = Gvtx(v+3,:);  aB = abs(det(tB));  % Subtriangle B
        tC = T;  tC(3,:) = Gvtx(v+3,:);  aC = abs(det(tC));  % Subtriangle C
        
        % Test if the point is inside the face
        if single(aT) == ... % Single to avoid errors due to precision
                single(aA + aB + aC);
            
            % Weight appropriately by the areas and interpolate the
            % value between the 3 vertices, and 3 coordinates
            vtx4(Cidxi(v),:) = [aA aB aC] * vtx3(vidx,:) ./ aT;
        end
    end
end

% Save output surface
srfwrite(vtx4,fac2,varargin{4});
