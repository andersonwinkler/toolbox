function icodown(varargin)
% Downsample a data-per-vertex (DPV), data-per-face (DPF) or
% surface (SRF) in ASCII format from a higher-order tessellated
% icosahedron to a lower order one. Note that the DPF
% file must have been derived from a suface file from either
% FreeSurfer or Brainder tools. It may not work if the vertices
% follow a different sequence inside the file.
%
% For facewise, the downsampling method is pycnophylactic and,
% therefore, should be used for quantities that require mass
% conservation, such as areal quantities.
% For vertexwise, the method removes the redundant vertices.
% Consider smoothing if needed before downsampling.
%
% Usage:
% icodown(filein,ntarget,fileout,surface)
%
% - filein  : File to be downsampled.
% - ntarget : Icosahedron order of the downsampled file.
% - fileout : Output file, downsampled.
% - suface  : Optional. For facewise, if the face indices aren't
%             ordered in the usual manner, entering a surface
%             is necessary to merge the faces correctly.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Apr/2012 (first version)
% Jun/2025 (this version)
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
        fprintf('Downsample a DPF file from a higher-order tessellated\n');
        fprintf('Downsample a data-per-vertex (DPV), data-per-face (DPF) or\n');
        fprintf('surface (SRF) in ASCII format from a higher-order tessellated\n');
        fprintf('icosahedron to a lower order one. Note that the DPF\n');
        fprintf('file must have been derived from a suface file from either\n');
        fprintf('FreeSurfer or Brainder tools. It may not work if the vertices\n');
        fprintf('follow a different sequence inside the file.\n');
        fprintf('\n');
        fprintf('For facewise, the downsampling method is pycnophylactic and,\n');
        fprintf('therefore, should be used for quantities that require mass\n');
        fprintf('conservation, such as areal quantities.\n');
        fprintf('For vertexwise, the method removes the redundant vertices.\n');
        fprintf('Consider smoothing if needed before downsampling.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('icodown <filein> <ntarget> <fileout> [surface]\n');
        fprintf('\n');
        fprintf('- filein  : File to be downsampled.\n');
        fprintf('- ntarget : Icosahedron order of the downsampled file.\n');
        fprintf('- fileout : Output file, downsampled.\n');
        fprintf('- suface  : Optional. For facewise, if the face indices aren''t\n');
        fprintf('            ordered in the usual manner, entering a surface\n');
        fprintf('            is necessary to merge the faces correctly.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Apr/2012 (first version)\n');
        fprintf('May/2016 (this version, now with vertexwise and surfaces geometry)\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Accept arguments
filein  = varargin{1};
ntarget = varargin{2};
fileout = varargin{3};
if ischar(ntarget)
    ntarget = eval(ntarget);
end

% Constants for the icosahedron
V0 = 12;
F0 = 20;

% Read data from disk. Try first mgz/mgh, then dpx, then srf, then give up
[~,~,fext] = fileparts(filein);
if strcmpi(fext,'.mgz') || strcmpi(fext,'.mgh')
    fs_home = getenv('FREESURFER_HOME');
    if ~ isempty('FREESURFER_HOME') && exist(fullfile(fs_home,'matlab'),'dir')
        addpath(fullfile(fs_home,'matlab'));
        [dpx,extra.M,extra.mr_parms,extra.volsz] = load_mgh(filein);
        nX = size(dpx,1);
        ftype = 'mgh';
    else
        error('To open MGH/MGZ files, FreeSurfer must be installed and configured.');
    end
else
    try
        [dpx,crd,idx] = dpxread(filein);
        nX = numel(dpx);
        ftype = 'dpx';
    catch
        try
            [vtx,fac] = srfread(filein);
            nV = size(vtx,1);
            ftype = 'srf';
        catch
            error('File %s not recognised as either DPX (curvature) or SRF (surface).',filein);
        end
    end
end

switch ftype
    
    case 'srf'
        
        % Identify icosahedron order
        n = round(log((nV-2)/(V0-2))/log(4));
        
        % Sanity check
        if nV ~= 4^n*(V0-2)+2
            error('Data not from icosahedron.');
        elseif ntarget >= n
            error('This script only downsamples data (target order: %d / input order: %d).',ntarget,n);
        else
            fprintf('Downsampling surface geometry:');
        end
        
        % Remove vertices:
        vtx = vtx(1:(4^ntarget*(V0-2)+2),:);

        % Remove face indices:
        for j = (n-1):-1:ntarget
            fprintf(' %d',j);
            nVj     = 4^j*(V0-2)+2;
            nFj     = 4^j*F0;
            nVjprev = 4^(j+1)*(V0-2)+2;
            nFjprev = 4^(j+1)*F0;
            remap   = 1:nVjprev;
            for f = 1:nFjprev
                v1 = fac(f,1);
                v2 = fac(f,2);
                v3 = fac(f,3);
                remap(v1) = min(min(remap(v1),v2),v3);
            end
            facnew = zeros(nFj,3);
            for f = 1:nFj
                for v = 1:3
                    facnew(f,v) = remap(fac(f,v));
                end
            end
            fac = facnew;
        end
        fprintf('. Done.\n');
        
        % Save to disk
        srfwrite(vtx,fac,fileout);
        
    case {'dpx','mgh'}
        
        % Detect what kind of data this is
        if mod(nX,10)
            
            % Find icosahedron order
            n = round(log((nX-2)/(V0-2))/log(4));
            
            % Sanity check
            if nX ~= 4^n*(V0-2)+2
                error('Data not from icosahedron.');
            elseif ntarget >= n
                error('This script only downsamples data.');
            else
                fprintf('Downsampling vertexwise data.\n');
            end
            
            % Downsample vertices
            dpx = dpx(1:(4^ntarget*(V0-2)+2),:,:,:);
            
        else
            
            % Find icosahedron order
            n = round(log(nX/F0)/log(4));
            
            % Sanity check
            if nX ~= 4^n*F0
                error('Data not from icosahedron.');
            elseif ntarget >= n
                error('This script only downsamples data.')
            else
                fprintf('Downsampling facewise data:');
            end
            
            % If a surface was given
            if nargin == 4
                
                % Load it
                [~,fac] = srfread(filein);
                
                % Downsample faces (general case)
                for j = (n-1):-1:ntarget
                    fprintf(' %d',j);
                    nVj    = 4^j*(V0-2)+2;
                    facnew = zeros(4^j*F0,3);
                    fout   = find(all(fac > nVj,2));
                    dpfnew = dpf(fout,:,:,:);
                    for f = 1:numel(fout)
                        vidx = fac(fout(f),:);
                        fidx = sum(...
                            fac == vidx(1) | ...
                            fac == vidx(2) | ...
                            fac == vidx(3),2) == 2;
                        ftomerge = fac(fidx,:);
                        facnew(f,:) = sum(ftomerge.*(ftomerge <= nVj),1);
                        dpfnew(f,:,:,:) = dpfnew(f,:,:,:) + sum(dpf(fidx,:,:,:),1);
                    end
                    fac = facnew;
                end
                fprintf('. Done.\n');
                
            else
                % Downsample faces (platonic only)
                for d = 1:(n-ntarget)
                    siz = size(dpx);
                    siz(1) = siz(1)/4;
                    dpxnew = zeros(siz);
                    for n = 1:size(dpx,4)
                        dpxnew(:,1,1,n) = sum(reshape(dpx(:,1,1,n),[4 nX/4]))';
                    end
                end
            end
        end
        
        % Save to disk
        if strcmpi(ftype,'mgh')
            save_mgh(dpx,fileout,extra.M,extra.mr_parms);
        else
            nXdown = size(dpx,1);
            dpxwrite(fileout,dpx,crd(1:nXdown,:),idx(1:nXdown,1));
        end
end
