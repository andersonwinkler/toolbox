% #!/usr/bin/octave -q
function rpncalc(varargin)
% Do some simple calculations using RPN notation.
%
% Accepted inputs are file names for DPV/DPF files, for
% CSV files, and for Matlab/Octave MAT files containing at most
% one variable inside.
%
% The current operators available are:
% Mathematical operators (binary):
%   Additive: +, -
%   Multiplicative (elementwise): *, /, and ^
%   Multiplicative (matrix): **, //, \\ and ^^
%   Logic: <, >, <=. >=, ==, ~=
% Mathematical operators (unary): log, ln, exp
% Stack manipulation: swap dup drop
% File operations: load save
%
% All inputs must be strings (so, delimited with quotes '')
%
% Vertex coordinates (for DPV) or the face indices (for DPF)
% will be the same as for the last DPV/DPF file loaded.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jun/2011 (first version)
% Nov/2014 (this version)
% http://brainder.org

% Do the OCTAVE stuff, with TRY to ensure MATLAB compatibility
try %#ok
    % Get the inputs
    varargin = argv();
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q'),
        fprintf('Do some simple calculations using RPN notation.\n');
        fprintf('\n');
        fprintf('Accepted inputs are file names for DPV/DPF files, for\n');
        fprintf('CSV files, and for Matlab/Octave MAT files containing at most\n');
        fprintf('one variable inside.\n');
        fprintf('\n');
        fprintf('The current operators available are:\n');
        fprintf('Mathematical operators (binary):\n');
        fprintf('  Additive: +, -\n');
        fprintf('  Multiplicative (elementwise): *, /, and ^\n');
        fprintf('  Multiplicative (matrix): **, //, \\\\ and ^^\n');
        fprintf('  Logic: <, >, <=. >=, ==, ~=\n');
        fprintf('Mathematical operators (unary): log, ln, exp\n');
        fprintf('Stack manipulation: swap dup drop\n');
        fprintf('File operations: load save\n');
        fprintf('\n');
        fprintf('All inputs must be strings (so, delimited with quotes '''')\n');
        fprintf('\n');
        fprintf('Vertex coordinates (for DPV) or the face indices (for DPF)\n');
        fprintf('will be the same as for the last DPV/DPF file loaded.\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Yale University / Institute of Living\n');
        fprintf('Jun/2011 (first version)\n');
        fprintf('Nov/2014 (this version)\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% More OCTAVE stuff
nargin = numel(varargin);

% Define the operators
opadd   = {'+','-'};
opmul   = {'*','/','\','^'};
oplogic = {'<','>','<=','>=','==','~='};
opmat   = {'**','//','\\','^^'};

stack = cell(0,0);
for a = 1:nargin,

    if strcmpi('load',varargin{a}),
        
        % LOAD
        
        if exist(stack{1},'file') == 2,
            
            % If the current argument is a file, load and put it at
            % the 1st level of the stack 
            fprintf('Loading file: %s',stack{1})
            [~,~,fext] = fileparts(stack{1});
            if any(strcmpi(fext,{'.dpv','.dpf','.dpx','.asc'})),
                fprintf(' (as DPX file)\n');
                dat = dpxread(stack{1});
            elseif any(strcmpi(fext,{'.mat'})),
                fprintf(' (as MAT file)\n');
                tmp = load(stack{1});
                fields = fieldnames(tmp);
                if numel(fields) > 1,
                    error('Too many variables stored in the file: \n%s',stack{1});
                end
                dat = tmp.(fields{1});
                clear tmp fields;
            else
                fprintf(' (as CSV file)\n');
                dat = csvread(stack{1});
            end
            stack{1} = dat;
        else
            error('File %s not found.\n',stack{1});
        end
        
    elseif strcmpi('save',varargin{a}),
        
        % SAVE
        
        % Only DPX and CSV supported at the moment. The data remains
        % in the stack. Use 'drop' to get rid of it if necessary
        fprintf('Saving file: %s',stack{1})
        [~,~,fext] = fileparts(stack{1});
        if any(strcmpi(fext,{'.dpv','.dpf','.dpx','.asc'})),
            fprintf(' (as DPX file)\n');
            dpxwrite(stack{1},double(stack{2}));
        elseif any(strcmpi(fext,{'.mat'})),
            fprintf(' (as MAT file)\n');
            tmp = stack{2}; %#ok
            save(stack{1},'tmp','-v7.3');
            clear tmp;
        else
            fprintf(' (as CSV file)\n');
            csvwrite(stack{1},double(stack{2}));
        end         
        stack(1) = [];
        
    elseif any(strcmp(opadd,varargin{a})),
        
        % ADD/SUBTRACT
        
        % If the current argument is additive, execute it
        fprintf('Adding/Subtracting\n')
        stack{1} = eval(sprintf('stack{2} %s stack{1}',varargin{a}));
        stack(2) = [];
        
    elseif any(strcmp(opmul,varargin{a})),
        
        % MULTIPLY/DIVIDE/POTENTIATE (ELEMENTWISE)
        
        % If the current argument is multiplicative, execute it
        fprintf('Multiplying/Dividing/Potentiating (elementwise)\n')
        stack{1} = eval(sprintf('stack{2} .%s stack{1}',varargin{a}));
        stack(2) = [];
        
    elseif any(strcmp(opmat,varargin{a})),
        
        % MULTIPLY/DIVIDE/POTENTIATE (MATRIX)
        
        % If the current argument is multiplicative, execute it
        fprintf('Multiplying/Dividing/Potentiating (matrix)\n')
        stack{1} = eval(sprintf('stack{2} %s stack{1}',varargin{a}(1)));
        stack(2) = [];
        
    elseif any(strcmp(oplogic,varargin{a})),
        
        % LOGICAL
        
        % If the current argument is multiplicative, execute it
        fprintf('Performing a logical operation\n')
        stack{1} = eval(sprintf('stack{2} %s stack{1}',varargin{a}));
        stack(2) = [];
        
    elseif strcmpi('log',varargin{a});
        
        % LOGARITHM
        
        % Compute the log of the 1st element in the stack
        fprintf('Taking the base-10 logarithm.\n')
        stack{1} = log10(stack{1});
        
    elseif strcmpi('ln',varargin{a});
        
        % LOGARITHM
        
        % Compute the log of the 1st element in the stack
        fprintf('Taking the natural logarithm.\n')
        stack{1} = log(stack{1});
        
    elseif strcmpi('exp',varargin{a});
        
        % EXPONENTIAL
        
        % Compute the log of the 1st element in the stack
        fprintf('Exponentiating.\n')
        stack{1} = exp(stack{1});
        
    elseif strcmpi('swap',varargin{a});
        
        % SWAP
        
        % Swap the 2 top elements of the stack
        fprintf('Swapping first and second stack elements.\n')
        tmpvar = stack{1};
        stack{1} = stack{2};
        stack{2} = tmpvar;
        
    elseif strcmpi('dup',varargin{a});
        
        % DUPLICATE
        
        % Duplicate the 1st element of the stack
        fprintf('Duplicating the first stack element.\n')
        for s = numel(stack):-1:1,
            stack{s+1} = stack{s};
        end
        stack{1} = stack{2};
        
    elseif strcmpi('drop',varargin{a});
        
        % DROP
        
        % Duplicate the 1st element of the stack
        fprintf('Dropping first stack element.\n')
        stack(1) = [];
        
    elseif isreal(str2num(varargin{a})),
        
        % NUMBER AND STRING
        
        % Roll the stack first
        for s = numel(stack):-1:1,
            stack{s+1} = stack{s};
        end
        tmp = str2num(varargin{a});
        if isempty(tmp)
            fprintf('Entering string: %s\n',varargin{a});
            stack{1} = varargin{a};
        else
            fprintf('Entering scalar: %g\n',tmp);
            stack{1} = tmp;
        end
    end
end
