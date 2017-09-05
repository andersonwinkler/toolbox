function nIN = baplot(varargin)
% Draw Bland-Altman plots for p-values.
%
% X    : Reference variable
% Y    : Test variable
% type : [Optional] 'orig' or 'modif'. Use the original
%        Bland-Altman plot in which the horizontal axis
%        has the average of X and Y or a modified version
%        in which only X is used. Default is 'orig'.
% xlim or nP : [Optional] A scalar representing the number
%        of permutations. If specified, instead of standard
%        deviation lines, it shows an ellipsoid with the
%        confidence interval. The xlim is then automatically
%        set to [0 1].
%        If instead of a scalar, a 2-element vector is given,
%        then these are used as minimum and maximum values
%        in the horizontal axis.
% ylim : [Optional] A 2-element vector used as min and max for horizontal
%        and vertical axes respectively.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2017
% http://brainder.org

% Take in arguments
narginchk(2,5);
X    = varargin{1};
Y    = varargin{2};
type = 'orig';
nP   = 0;
xlim = [];
ylim = [];
if nargin >= 3,
    type = varargin{3};
end
if nargin >= 4,
    if isscalar(varargin{4}),
        nP = varargin{4};
        if nP <= 0,
            error('nP must be larger than 0.');
        end
    else
        xlim = varargin{4};
    end
end
if nargin >= 5,
    ylim = varargin{5};
end

% BA data
if strcmpi(type,'orig'),
    Xplot = (X+Y)/2;
else
    Xplot = X;
end
Yplot = Y-X;

% Confidence inteval ellipsoid
if nP > 0,
    [L,U] = confint(nP,0:nP,.05,'wilson');
    xCI   = (0:nP)/nP;
    if strcmpi(type,'orig'),
        Xplot = (X+Y)/2;
        yU    = U-xCI;
        yL    = L-xCI;
        d     = (U-L)/2;
        yU    = +d;
        yL    = -d;
    else
        Xplot = X;
        yU    = U-xCI;
        yL    = L-xCI;
    end
    plot([0 1],[0 0],'k'); hold on
    plot(xCI,yU,'r:'); hold on
    plot(xCI,yL,'r:'); hold on
    scatter(Xplot,Yplot,'Marker','.'); hold off
    set(gca,'xlim',[0 1]);
    if ~isempty(ylim),
        set(gca,'ylim',ylim);
    end
else
    S = std(Yplot);
    scatter(Xplot,Yplot,'Marker','.'); hold on
    if isempty(xlim),
        xlim = get(gca,'xlim');
    end
    plot(xlim,[0 0],'r:'); hold on
    plot(xlim,+[1 1]*S,'r:'); hold on
    plot(xlim,-[1 1]*S,'r:'); hold off
    set(gca,'xlim',xlim);
    if ~isempty(ylim),
        set(gca,'ylim',ylim);
    end
end

% Return percentage within CI if requested
if nargout == 1,
    if nP > 0,
        [L,U] = confint(nP,Xplot*nP,.05,'wilson');
        U = (U - Xplot);
        L = (L - Xplot);
        nIN = sum((Yplot < U) & (Yplot > L));
    else
        nIN = sum((Yplot < S) & (Yplot > -S));
    end
end
