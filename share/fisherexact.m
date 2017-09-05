function P = fisherexact(varargin)
% Compute the p-value for the Fisher's exact test, modified
% for hypothesis testing. The alternative hypothesis is that
% one or both anti-diagonal elements are equal to zero, which
% is the best possible performance if one variable can predict
% another.
% 
% Usage:
% P = fisherexact(C,method);
% P = fisherexact(Y,X,'method');
% 
% C    : Contingency table (2x2).
% Y    : Observed dicothomous variable. It can be a Nx1 vector
%        or a NxM array.
% X    : Explanatory dicothomous variable. It can be a Nx1 vector
%        or a NxM array.
% meth : Method. It can be one of:
%        - 'Fisher'
%        - 'Pearson'
%        - 'Yates'
% P    : P-values (uncorrected for multiple testing) for all the
%        M columns of X and/or Y.
%
% Note: All results are one-tailed.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Mar/2014 (first version)
% Aug/2015 (this version)
% http://brainder.org

% Take inputs
if nargin == 1,
    v.C = varargin{1};
elseif nargin == 2
    if ischar(varargin{2}),
        v.C = varargin{1};
        v.m = varargin{2};
    else
        v.Y = varargin{1};
        v.X = varargin{2};
    end
elseif nargin == 3,
    v.Y = varargin{1};
    v.X = varargin{2};
    v.m = varargin{3};
end

if isfield(v,'C'),
    a = v.C(1);
    c = v.C(2);
    b = v.C(3);
    d = v.C(4);
    M = 1;
elseif isfield(v,'Y'),
    Y = logical(v.Y);
    X = logical(v.X);
    M = max(size(Y,2),size(X,2));
    if size(X,1) ~= size(Y,1),
        error('Input arguments must have the same number of rows.');
    end
    if size(X,2) > 1 && size(Y,2) > 1 && size(X,2) ~= size(Y,2),
        error('X and Y must have the same number of columns, or one of these must be a single column.');
    end
    a = sum(bsxfun(@and, Y, X),1);
    b = sum(bsxfun(@and, Y,~X),1);
    c = sum(bsxfun(@and,~Y, X),1);
    d = sum(bsxfun(@and,~Y,~X),1);
end
if isfield(v,'m'),
    meth = v.m;
else
    meth = 'Fisher';
end

% Margins & total
ab = a + b;
cd = c + d;
ac = a + c;
bd = b + d;
n  = unique(a + b + c + d);

% Some prepration
switch lower(meth),
    
    case {'fisher', 'hypergeometric'},
        
        % Flip the table to make sure we'll test towards the
        % most favourable alternative. This will also make it faster.
        if (a+d) < (b+c),
            anew = b;
            cnew = d;
            b    = a;
            d    = c;
            a    = anew;
            c    = cnew;
            ac   = a + c;
            bd   = b + d;
        end
        lfac = lfactorial(n);
        P    = zeros(1,M);
    case {'pearson', 'yates'},
        
        % The expected values are the outer product of the margins,
        % divided by the total.
        Ea = ab.*ac./n;
        Eb = ab.*bd./n;
        Ec = cd.*ac./n;
        Ed = cd.*bd./n;
end

% Pick one fo the methods
switch lower(meth),
    
    case 'fisher',

        % Compute the probabilities for each state of the table towards
        % the alternative, that is, the case in which one or both elements
        % of the anti-diagonal are 0.
        for s = 0:max(max(b(:),c(:))),
            as  = a + s;
            bs  = b - s;
            cs  = c - s;
            ds  = d + s;
            m   = bs >= 0 & cs >= 0;
            lPs = ...
                sum(lfac([ab(m) cd(m) ac(m) bd(m)]+1)) - ...
                sum(lfac([as(m) bs(m) cs(m) ds(m) n]+1));
            P(m) = P(m) + exp(lPs)';
        end
        
    case 'hypergeometric',
        
        % This is the same as Fisher, but using a slightly different
        % strategy, with the cdf of the hypergeometric distribution
        % directly.
        P = zeros(1,M);
        for j = a:min((a+c),(a+b)),
            lPs = ...
                sum(lfac(int64([   a+b; c+d;   a+c;   b+d  ])+1),1) - ...
                sum(lfac(int64([n; j;   a+b-j; a+c-j; d-a+j])+1),1);
            P = P + exp(lPs);
        end
        
    case 'pearson',
        
        % Pearson Chi^2 test
        X2 = ...
            (a-Ea).^2./Ea + ...
            (b-Eb).^2./Eb + ...
            (c-Ec).^2./Ec + ...
            (d-Ed).^2./Ed;
        df = ones(size(X2));
        P  = gammainc(X2/2,df/2,'upper')/2;
        
    case 'yates',
        
        % Yates continuity correction over Pearson's Chi^2.
        X2 = ...
            (abs(a-Ea)-.5).^2./Ea + ...
            (abs(b-Eb)-.5).^2./Eb + ...
            (abs(c-Ec)-.5).^2./Ec + ...
            (abs(d-Ed)-.5).^2./Ed;
        df = ones(size(X2));
        P  = gammainc(X2/2,df/2,'upper')/2;
        
    case 'student',
        
        % Linear regression using Students' t-statistic.
        % Both are mean centered.
        if isfield(v,'C'),
            YX = [ ...
                repmat([1 1],[a 1]);...
                repmat([1 0],[b 1]);...
                repmat([0 1],[c 1]);...
                repmat([0 0],[d 1])];
            Y = YX(:,1);
            X = YX(:,2);
        end
        Y   = bsxfun(@minus,Y,mean(Y));
        X   = bsxfun(@minus,X,mean(X));
        psi = zeros(1,size(Y,2));
        for p = 1:size(Y,2),
            psi(p) = X(:,p)\Y(:,p);
        end
        res = Y - bsxfun(@times,X,psi);
        sig = sqrt(sum(res.^2,1)/(n-1));
        t   = bsxfun(@times,sum(X.*X,1).^.5,psi) ./sig;
        P   = tcdf(-abs(t),n-1);
end

% ==============================================================
function lfac = lfactorial(N)
% Compute the log(factorial(0:N)), so dealing with
% precision issues.
lfac = zeros(N+1,1);
for n = 1:N,
    lfac(n+1) = log(n) + lfac(n);
end
