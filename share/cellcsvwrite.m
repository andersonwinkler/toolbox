function cellcsvwrite(varargin)
% Write a 2-D cell as a CSV file.
% 
% Usage:
% cellcsvwrite(cell,filename,delimiter,end-of-line)
% 
% Inputs:
% - cell  = Cell array to be saved
% - filename = CSV file to be created
% - delim    = Field separator. Default = ','
% - eol      = Record separator. Default = '\n'
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2011

% Accept the inputs
error(nargchk(1,4,nargin));
C = varargin{1};
filename = varargin{2};
if nargin == 4,
    delimiter = sprintf(varargin{3});
    eolmark = sprintf(varargin{4});
elseif nargin == 3,
    delimiter = sprintf(varargin{3});
    eolmark = sprintf('\n');
elseif nargin == 2,
    delimiter = sprintf(',');
    eolmark = sprintf('\n');
end

% Reshape as a column
nC = size(C,2);
C = reshape(C',numel(C),1);

% Replace numbers for strings
for c = 1:numel(C);
    C{c} = num2str(C{c});
end

% Create a formatting string
fstr = repmat(['%s' delimiter],[1 nC]);
fstr = fstr(1:numel(fstr) - numel(delimiter));
fstr = [fstr eolmark];

% Write to the disk
fid = fopen(filename,'w');
fprintf(fid,fstr,C{:});
fclose(fid);