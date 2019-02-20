function table = strcsvread(varargin)
% Read a CSV file containing non-numeric fields (strings) and
% return a cell array. Numbers are converted to double and, if contiguous,
% can be easily converted to arrays with cell2mat. The remaining
% are left as strings.
%
% Usage:
% table = strcsvread(filename,delimiter,end-of-line)
% 
% - filename = CSV file to be read
% - delim    = Field separator. Default = ','
% - eol      = Record separator. Default = '\n'
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Oct/2010
% http://brainder.org

% Accept the inputs
narginchk(1,3);
filename = varargin{1};
if nargin == 3,
    delimiter = sprintf(varargin{2});
    eolmark = sprintf(varargin{3});
elseif nargin == 2,
    delimiter = sprintf(varargin{2});
    eolmark = sprintf('\n');
elseif nargin == 1,
    delimiter = sprintf(',');
    eolmark = sprintf('\n');
end

% Read the whole CSV file as a single stream of data
fid = fopen(filename,'r');
table = fread(fid,Inf,'uint8=>char')';
fclose(fid);
tic
% Add an EOL at the end of the stream if not existing
% This prevents an error later
if table(end) ~= eolmark,
    table = [table eolmark];
end

% Identify where the EOLs and delimiters are
eolpos = find(table == eolmark);
delpos = find(table == delimiter);
nR = sum(table == eolmark);
nC = floor(sum(table == delimiter)/nR)+1;
table([eolpos delpos]) = [];
table = mat2cell(table,1,diff(sort([0 eolpos delpos])')-1);
table = reshape(table,nC,nR)';

% Loop over each EOL get the content between each
for r = 1:nR,
    for c = 1:nC,
        x = strtrim(table{r,c});
        if isempty(x) || strcmpi(x,'NaN'),
            table{r,c} = NaN;
        else
            testnum = str2double(table{r,c});
            if ~ isnan(testnum),
                table{r,c} = testnum;
            end
        end
    end
end
toc
