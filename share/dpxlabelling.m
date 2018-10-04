% #!/usr/bin/octave -q
function dpxlab = dpxlabelling(varargin)
% Take a binary DPX file (DPV or DPF) and a surface file (SRF)
% and label contiguous vertices (DPV) or faces (DPF).
%
% Usage:
% dpxlab = dpxlabelling(dpxfile,srffile,dpxsave)
%
% dpxfile : File to be labelled, in ASCII format (DPF/DPV)
% srffile : Surface file, from where the geometry information
%           is taken
% dpxsave : (Optional) File to be created containing the
%           labelled vertex or face clusters.
% dpxlab  : Variable with the labeled clusters.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Feb/2012
% http://brainder.org

% Do the OCTAVE stuff, using TRY to ensure MATLAB compatibility
try
    % Get the inputs
    varargin = argv();
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q')
        fprintf('Take a binary DPX file (DPV or DPF) and a surface file (SRF)\n');
        fprintf('and label contiguous vertices (DPV) or faces (DPF).\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('dpxlabelling <dpxfile> <srffile> <dpxsave>\n');
        fprintf('\n');
        fprintf('dpxfile : File to be labelled, in ASCII format (DPF/DPV)\n');
        fprintf('srffile : Surface file, from where the geometry information\n');
        fprintf('          is taken\n');
        fprintf('dpxsave : File to be created containing the labelled\n');
        fprintf('          vertex or face clusters.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('FMRIB / Univ. of Oxford\n');
        fprintf('Feb/2012\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Accept arguments
if nargin < 3
    error('Invalid number of arguments')
end

% % Read input files
[dpx,crd,idx] = dpxread(varargin{1});
[vtx,fac] = srfread(varargin{2});

% Number of datapoints
nX = size(dpx,1);
nV = size(vtx,1);
nF = size(fac,1);

% Labelled mask
dpxl = zeros(size(dpx));

if nX == nV % vertexwise data
    fprintf('Working with vertexwise data.\n');
    
    % Identify faces that are entirely (3 vertices)
    % inside the mask
    facm = reshape(dpx(fac(:)),size(fac));
    faci = all(facm,2);
    %facm(~faci,:) = false;
    
    % For each face
    k = 1; % label indices (to be incremented)
    for f = find(faci)'
        flab = sort(dpxl(fac(f,:)));
        if flab(3) == 0
            dpxl(fac(f,:)) = k;
        else
            if flab(2) ~= 0
                dpxl(dpxl == flab(2)) = flab(3);
            end
            if flab(1) ~= 0
                dpxl(dpxl == flab(1)) = flab(3);
            end
            dpxl(fac(f,:)) = flab(3);
        end
        k = k + 1;
    end
    
elseif nX == nF  % facewise data %%% NOT TESTED YET %%%
    fprintf('Working with facewise data.\n');
    
    % Ignore faces that are not in the mask
    facm = bsxfun(@times,fac,D);
    
    % For each face
    k = 1; % label indices (to be incremented)
    for f = find(dpx)'
        
        % Find other faces that share 2 vtx (1 edge)
        neifacidx = ...
            facm == fac(f,1) | ...
            facm == fac(f,2) | ...
            facm == fac(f,3);
        neifacidx = sum(neifacidx,2) >= 2; % only 2+ shared vertices, not 1
        numneigh  = sum(neifacidx); % number of neighbours
        if numneigh > 1 % if not isolated face
            flab = sort(dpxl(neifacidx));
            if flab(numneigh) == 0
                dpxl(neifacidx) = k;
            else
                for n = numneigh:-1:1
                    if flab(n) ~= 0
                        dpxl(dpxl == flab(n)) = flab(numneigh);
                    end
                end
                dpxl(f) = flab(numneigh);
            end
        else % isolated face
            dpxl(f) = k;
        end
        k = k + 1;
    end
else
    error('The data does not match the surface.');
end

% Sort according to size (largest first)
dpxu = unique(dpxl);
dpxu(dpxu == 0) = [];
nC = numel(dpxu);
sizes = zeros(nC,1);
for l = 1:nC
    sizes(l) = sum(dpxl == dpxu(l));
end
[~,idxs] = sort(sizes,'descend');
dpxu = dpxu(idxs);

% Relabel starting from 1
dpxlab = zeros(size(dpxl));
for l = 1:nC
    dpxlab(dpxl == dpxu(l)) = l;
end

% Save results
if nargin == 3
    dpxwrite(varargin{3},dpxlab,crd,idx);
end
