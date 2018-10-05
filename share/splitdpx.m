function splitdpx(varargin)
% Split a DPV/DPF file according to the labels given by another
% DPF/DPV file
%
% Usage:
% splitdpx(dpxfile,labelfile,dpxprefix)
%
% dpxfile   : DPX file to be split (*.dpv/dpf/dpx).
% labelfile : Labels per vertex or per face, in DPF or DPF format.
%             If empty, one surface file is generated for each closed
%             set of vertices (e.g., one per hemi).
% dpxprefix : File prefix (may include path) to create the new files.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Sep/2011
% http://brainder.org

% Do the OCTAVE stuff, using TRY to ensure MATLAB compatibility
try
    % Get the inputs
    varargin = argv();
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Split a surface according to the labels given by a DPF/DPV file\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('splitdpx <dpxfile> <labelfile> <dpxprefix>\n');
        fprintf('\n');
        fprintf('srffile   : DPX file to be split (*.dpv/dpf/dpx).\n');
        fprintf('labelfile : Labels per vertex or per face, in DPF or DPF format.\n');
        fprintf('            If empty, one surface file is generated for each closed\n');
        fprintf('            set of vertices (e.g., one per hemi).\n');
        fprintf('srfprefix : File prefix (may include path) to create the new files.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Sep/2011\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Accept inputs
dpxfile   = varargin{1};
labelfile = varargin{2};
dpxprefix = varargin{3};

lab = dpxread(labelfile);
U = unique(lab);
for u = 1:numel(U),
    idx = lab == U(u);
    fname = sprintf('%s.%0.4d.dpx',dpxprefix,U(u));
    dpxwrite(fname,dpx(lab));
end
