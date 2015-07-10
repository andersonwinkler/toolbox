function ped = loadped(varargin)
% Load the pedigree
% 
% Usage:
% ped = loadped('pedindex.csv','phi2.csv')
% ped = loadped('pedindex.csv','phi2.csv','house.csv')
% 
% The input files must have been converted to CSV files.
% The output is a struct that contains the headers for
% the pedindex, the pedindex itself (as a double array)
% and the matrices Phi2 and Delta7.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jul/2010 (first version)
% Mar/2014 (this version)
% http://brainder.org

if nargin < 2 || nargin > 3,
    error('Incorrect number of arguments');
end
pedfile = varargin{1};
phifile = varargin{2};

% Read the pedindex.csv file
pedindex = strcsvread(pedfile);

% Move the headers to a separate variable
headers = pedindex(1,:);
pedindex(1,:) = [];

% Subject IDs
sidx = strcmpi('ID',headers);
id = pedindex(:,sidx);

% Family IDs (if existing)
fidx = strcmpi('FAMID',headers);
if any(fidx),
    famid = pedindex(:,fidx);
end

% Remove the ID and FAMID
pedindex(:,sidx|fidx) = [];

% Convert char to num and to an array
pedindex = cell2mat(pedindex);

% Read the phi2.csv file
fid = fopen(phifile,'r');
phitable = textscan(fid,repmat('%f',[1 4]),'Delimiter',',');
fclose(fid);

% Organise it as simple table
phitable = [phitable{:}];

% Reorganise as symmetric matrices
tmpvar = phitable(:,1:2);
N = max(tmpvar(:));
Phi2 = zeros(N);    % Kinship
Delta7 = zeros(N);  % Jacquard D7 coefficients

% For each row in the phi2 file
for r = 2:size(phitable,1),
    Phi2(phitable(r,1),phitable(r,2)) = phitable(r,3);
    Phi2(phitable(r,2),phitable(r,1)) = phitable(r,3);
    Delta7(phitable(r,1),phitable(r,2)) = phitable(r,4);
    Delta7(phitable(r,2),phitable(r,1)) = phitable(r,4);
end

% Do the same for the house.csv, if that exists.
if nargin == 3,
    % Read the house.csv file
    hhfile = varargin{3};
    fid = fopen(hhfile,'r');
    hhtable = textscan(fid,repmat('%f',[1 4]),'Delimiter',',');
    fclose(fid);
    
    % Organise it as simple table
    hhtable = [hhtable{:}];
    
    % Reorganise as a symmetric matrix
    tmpvar = hhtable(:,1:2);
    N = max(tmpvar(:));
    HH = zeros(N);    % Household
    
    % For each row in the house file
    for r = 2:size(hhtable,1),
        HH(hhtable(r,1),hhtable(r,2)) = hhtable(r,3);
        HH(hhtable(r,2),hhtable(r,1)) = hhtable(r,3);
    end
end

% Organize the output
ped.hdr  = headers;
ped.id   = id;
if any(fidx),
    ped.famid = famid;
end
ped.pidx = pedindex;
ped.phi2 = Phi2;
ped.d7   = Delta7;
if nargin == 3,
    ped.house = HH;
end