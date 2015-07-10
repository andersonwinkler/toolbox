function [G,df2] = pivotal(varargin)
% Compute a generic pivotal statistic G for inference when
% testing a general linear hypothesis.
%
% Usage:
% G = pivotal(M,psi,res,C,VG,Y);
%
% Inputs for a full model specification:
% M       : A n by r design matrix.
% psi     : A r by m matrix of regression coefficients
%           (already fitted). For the univariate case, m=1.
% res     : A n by m matrix with the residuals from the
%           model fit.
% C       : A r by s contrast matrix for the full model.
% VG      : A n by 1 vector where each unique value indicate
%           a variance group.
% Y       : (Optional) A n by m matrix of observations. If Y
%           is supplied, the variance weights will be those of
%           James (1951) and Welch (1951). Otherwise, the weights
%           are those of Horn et al (1975), which are superior.
%
% Outputs:
% G       : A 1 by m vector with the statistic for testing
%           the null hypothesis that C'*psi(:,c) = zeros(s,1)
%           for the c-th of the m columns of psi, or that
%           beta(:,c) is zero.
% df2     : Degrees of freedom #2. Df1 is simply rank(C).
%
% The partitioned model, using X and Z, can be provided as:
% G = pivotal([X Z],[beta; gamm],res,...
%        [eye(size(X,2)); zeros(size(Z,2),size(C,2))],VG,Y);
%
% The statistics that are computed by default depend on
% rank(C) and on the number of unique elements of VG:
% - If numel(unique(VG)) == 1, then a traditional F or
%   Student's t statistic is computed, with the error
%   variance being pooled.
% - If numel(unique(VG)) > 1, then the James, Welch's v^2 or
%   v is computed, so that the variances for each variance
%   group are not pooled.
% - If rank(C) == 1, the Student's t or Welch's v is returned.
% - If rank(C) > 1, the F or v^2 statistic is returned.
%
%
%                  +------------------+----------------------+
%                  |  Pool variances  |  Not pool variances  |
% +----------------+------------------+----------------------+
% | rank(C) > 1    |     F-ratio      | Welch's v^2 or James |
% +----------------+------------------+----------------------+
% | rank(C) == 1   |   Student's t    |       Welch's v      |
% +----------------+------------------+----------------------+
%
% References:
% - Winkler AM, Ridgway GR, Webster MG, Smith SM, Nichols TE.
%   Permutation inference for the general linear model.
%   Neuroimage. 2014;92:381-97.
% - James G. The comparison of several groups of observations when
%   the ratios of the population variances are unknown.
%   Biometrika, 1951;38(3):324-329.
% - Welch B. On the comparison of several mean values: an
%   alternative approach. Biometrika, 1951;38(3):330-336.
% - Horn SD, Horn RA, Duncan DB. Estimating heteroscedastic
%   variances in linear models. Journal of the American
%   Statistical Association, 1975;70(350):380-385.
% - MacKinnon JG, White H. Some heteroscedasticity-consistent
%   covariance matrix estimators with improved finite sample
%   properties. Journal of Econometrics, 1985;29(3):305-325.
% - Christensen R. Plane Answers to Complex Questions: The
%   Theory of Linear Models. Springer, 2002.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jul/2013 (first version)
% Aug/2014 (this version)
% http://brainder.org

% Take the inputs
if nargin < 3 || nargin > 6.
    error('Error: Invalid number of arguments');
end
M     = varargin{1};
psi   = varargin{2};
res   = varargin{3};
[n,r] = size(M);
C     = eye(r);
VG    = ones(n,1);
Y     = [];
if nargin >= 4, C  = varargin{4}; end
if nargin >= 5, VG = varargin{5}; end
if nargin == 6, Y  = varargin{6}; end

% Some small vars for later
m    = size(res,2);
rC   = rank(C);
[uVG,~,vgidx] = unique(VG);
nVG  = numel(uVG);

% Some sanity checks
if         size(psi,1) ~= r   ...
        || size(C,1)   ~= r   ...
        || size(res,1) ~= n   ...
        || size(res,2) ~= m   ...
        || (~isempty(Y)       ...
        && any(size(Y) ~= size(res))),
    error('Error: Dimension mismatch between input arguments.\n       [%d %d %d %d]',...
        [size(psi,1) ~= r
        size(C,1)    ~= r
        size(res,1)  ~= n
        size(res,2)  ~= m]);
end

% Residual forming matrix, diagonal only
H  = M*pinv(M);
R  = diag(eye(n) - H);
rM = sum(diag(H));

% Make an estimate for V (these are actually weights, but they take the
% place of V itself depending on the formulation).
% Loop over each variance block
W = ones(n,m);
if ~ isempty(Y),
    % Variance weights derived from Y as in James (1951) and Welch
    % (1951). These weights use the variances computed fom the
    % observations and so, don't consider correctly the loss of degrees
    % of freedom due to the model when all the samples are considered
    % together under homoscedasticity assumption. It appears to perform
    % well when the data is heteroscedastic, but under homoscedasticity,
    % with all the groups pooled together, it seems less powerful.
    
    % James-Welch (1951) weights:
    for v = 1:nVG,
        vidx = vgidx == v;
        W(vidx,:) = bsxfun(@rdivide,W(vidx,:),var(Y(vidx,:)));
    end
else
    % Variance weights derived as in Horn et al (1975) and recommended
    % by MacKinnon & White (1985) as their "HC2". This uses the residuals
    % and the df computed from the residual-forming matrix. These
    % weights are much better than those used in James/Welch method, taking
    % correctly into account the degrees of freedom lost with the model
    % from the residual forming matrix. It is the adequate weights to
    % use under both homo and heteroscedasticity.
    % Under hetero, but not homoscedasticity, both weights are the same.
    
    % Horn et al (1975) weights:
    for v = 1:nVG,
        vidx = vgidx == v;
        W(vidx,:) = bsxfun(@times,W(vidx,:),sum(res(vidx,:).^2));
        W(vidx,:) = bsxfun(@rdivide,sum(R(vidx)),W(vidx,:));
    end
end

% For each test (voxel, vertex, face, etc). Unfortunately, this
% for-loop cannot be avoided without making the code too complicated.
G = zeros(1,m);
for t = 1:m,
    % The statistic below is the matrix form of James (1951) or,
    % equivalently, the numerator of Welch v^2 when using their
    % variance weights under the assumption of heteroscedasticity.
    % Under homoscedasticity, with just one variance group and the
    % weights computed as in Horn et al (1975), then it's equivalent
    % to the traditional F-test. If diag(W) is replaced for V^(-1),
    % the expression becomes the statistic for a Generalized Linear
    % Squares model, as shown in Christensen (2002) (Theorem 3.8.2, p.88).
    % In this case, there is no need to include the MSE=res'*res/(n-rM)
    % as a denominator if sigma^2 is already present in V. Likewise,
    % it's already present in the weights defined above.
    
    % F-ratio or James statistic:
    G(t) = psi(:,t)'*C*pinv(C'*pinv(M'*diag(W(:,t))*M)*C)*C'*psi(:,t)/rC;
end

% If the errors are treated as homoscedastic, the statistic above is
% just the F-ratio for a conventional ANOVA/ANCOVA test. If heteroscedastic,
% it is the James (1951) statistic. Still for heteroscedastic errors, it
% needs be further corrected to become the Welch (1951) v^2 statistic.
% This means division of the original statistic by a certain term.
% When rank(C) = 1, or when homoscedasticity is assumed, den = 1 and so,
% James and Welch's v^2 are the same, and no need to compute the
% expression below.
if nVG > 1,
    bsum = zeros(size(G));
    sW = sum(W,1);
    for v = 1:nVG,
        vidx = vgidx == v;
        bsum = bsum + bsxfun(@rdivide,(1-sum(W(vidx,:))./sW).^2,sum(R(vidx)));
    end
    bsum = bsum/rC/(rC+2);

    % Welch's v^2 statistic:
    G = G./(1 + 2*(rC-1).*bsum);
    
    % Degrees of freedom
    df2 = 1/3./bsum;
else
    df2 = n - rM;
end

% If rank(C) = 1, take the square root and the appropriate sign to allow
% bidirectional hypothesis. The sqrt(v^2) is the statistic to be used for
% the Behrens-Fisher problem. The sqrt(F) is simply the Student's t.
if rC == 1,   
    G = sign(C'*psi).*sqrt(G);
end
