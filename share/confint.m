function [L,U] = confint(n,X,alpha,meth)
% Compute the confidence interval for a set of Bernoulli
% trials using one of different methods.
%
% Usage:
% [L,U] = confint(n,X,alpha,meth)
%
% - n     : Number of trials.
% - X     : Number of successful trials.
% - alpha : Coverage of the confidence interval.
%           For 95% CI, use alpha = 0.05.
% - meth  : Method to use. Can be one of:
%           - 'Wald'
%           - 'Wilson'
%           - 'Agresti-Coull'
%           - 'Jeffreys'
%           - 'Clopper-Pearson'
%           - 'arc-sine'
%           - 'logit'
%           - 'Anscombe'
% - L     : Lower bound for the confidence interval.
% - U     : Upper bound for the confidence interval.
%
% The variables n, X and alpha can be either scalars
% or arrays. When using arrays, use consistent sizes.
%
% Reference:
% Brown LD, Cai TT, DasGupta AA. Interval estimation for a
% binomial proportion. Statistical Science. 2001 16(2):101-133.
%
% _____________________________________
% Anderson Winkler & Tom Nichols
% FMRIB / University of Oxford
% Apr/2012
% http://brainder.org

k  = norminv(1-alpha/2);
p  = X./n;          % Proportion of successes
q  = 1 - p;         % Proportion of failures
Xt = X + (k.^2)/2;  % Modified number of sucesses
nt = n + k.^2;      % Modified number of trials
pt = Xt./nt;        % Modified proportion of successes
qt = 1 - pt;        % Modified proportion of failures

switch lower(meth),
    case 'wald',
        L = p - k.*sqrt(p.*q./n);
        U = p + k.*sqrt(p.*q./n);
        
    case 'wilson',
        L = pt - k.*sqrt(n.*p.*q + (k.^2)/4)./nt;
        U = pt + k.*sqrt(n.*p.*q + (k.^2)/4)./nt;
        
    case 'agresti-coull',
        L = pt - k.*sqrt(pt.*qt./nt);
        U = pt + k.*sqrt(pt.*qt./nt);
        
    case 'jeffreys',
        L = betainv(  alpha/2, X+1/2, n-X+1/2);
        U = betainv(1-alpha/2, X+1/2, n-X+1/2);
        
    case 'clopper-pearson',
        L = betainv(  alpha/2, X,   n-X+1);
        U = betainv(1-alpha/2, X+1, n-X  );
        
    case 'arc-sine',
        pa = (X+3/8)./(n+3/4);
        as = asin(sqrt(pa));
        L  = sin(as - k./(2*sqrt(n))).^2;
        U  = sin(as + k./(2*sqrt(n))).^2;
        
    case 'logit',
        lam  = log(X./(n-X));
        Vhat = n./(X.*(n-X));
        lamL = lam - k.*sqrt(Vhat);
        lamU = lam + k.*sqrt(Vhat);
        L    = exp(lamL)./(1+exp(lamL));
        U    = exp(lamU)./(1+exp(lamU));

    case 'anscombe',
        lam  = log((X+1/2)./(n-X+1/2));
        Vhat = (n+1).*(n+2)./(n.*(X+1).*(n-X+1));
        lamL = lam - k.*sqrt(Vhat);
        lamU = lam + k.*sqrt(Vhat);
        L    = exp(lamL)./(1+exp(lamL));
        U    = exp(lamU)./(1+exp(lamU));
        
    otherwise
        error('Unknown method.')
end
