function S = randstr(varargin)
% Creates a random string from specified alphabets.
% 
% Usage:
% S = randstr(N,A);
% 
% Inputs:
% - N : Number of characters in the string (i.e., its length).
% - A : Specification of the alphabet. It can be:
%       - 'DNA'  : Generates random sequences of the four 
%                  nucleotides found in DNA.
%       - 'RNA'  : Generates random sequences of the four 
%                  nucleotides found in DNA.
%       - 'aa20' : Generates sequences of the 20 aminoacids set
%       - 'aa22' : Generates sequences of the 22 aminoacids set
%       - In addition, a code comprised by one or more of 'A',
%         'a', '1' and '@' can be used to specify, respectively,
%         latin letters in upper case, latin letters in lower case
%         arabic digits and some common symbols that are safe for
%         use as file names in Linux systems. These alphabets can
%         be specified as, e.g., 'Aa' to allow both upper and lower
%         latin letters in the random string, or 'a1' to allow
%         only lower case and digits, or '1@' to allow only digits
%         and some common symbols.
% 
% Output:
% - S : The random string to be created.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% May/2012

% Accept arguments
N = varargin{1}; % string length to be created
A = varargin{2}; % specification of the alphabet

% Prepare the alphabet
if strcmpi(A,'dna'),       % for DNA nucleotides
    alphabet = 'ATGC';
elseif strcmpi(A,'rna'),   % for RNA nucleotides
    alphabet = 'AUGC';
elseif strcmpi(A,'aa20'),  % for the 20 aminoacids
    alphabet = 'ARNDCEQGHILKMFPSTWYV';
elseif strcmpi(A,'aa22'),  % for the 22 aminoacids
    alphabet = 'ARNDCEQGHILKMFPSTWYVUO';
else                       % concatenate multiple alphabets
    alphabet = '';
    for c = 1:numel(A),
        if A(c) == 'A',       % for the latin alphabet, capitals
            alphabet = strcat(alphabet,'ABCDEFGHIJKLMNOPQRSTUVWXYZ');
        elseif A(c) == 'a',   % for the latin alphabet, lower case
            alphabet = strcat(alphabet,'abcdefghijklmnopqrstuvwxyz');
        elseif A(c) == '@',   % for some common ASCII symbols
            alphabet = strcat(alphabet,'@!$#%&_-+');
        elseif A(c) == '1',   % for arabic digits
            alphabet = strcat(alphabet,'0123456789');
        else
            error('Unknown alphabet specification.');
        end
    end
end

% Random numbers. Use rng to control.
R = rand(1,N);
R = ceil(R.*numel(alphabet));
S = alphabet(R);
S = char(S); % output to return
