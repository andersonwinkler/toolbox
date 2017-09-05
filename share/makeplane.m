function [tstat,pval] = makeplane(varargin)
% Show the regression plane for a GLM with that has
% one response variable (Y), one independent variable
% of interest (X), and one nuisance variable (Z). The
% model is thus Y = X*b + Z*g + 1*d + e.
% 
% Usage:
% [tstat,pval] = makeplane(Y,X,Z)
% 
% Inputs:
% - Y      : Dependent variable.
% - X      : Independent variable of interest.
% - Z      : Independent nuisance variable.
% 
% If Z is omitted, a regression line is produced with the
% model Y = X*b + 1*d + e.
% 
% Outputs:
% - tstat  : t-statistic.
% - pvals  : p-values.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Jul/2015
% http://brainder.org

do2d = false;
do3d = false;
if nargin == 2,
    Y = varargin{1};
    X = varargin{2};
    do2d = true;
elseif nargin == 3,
    Y = varargin{1};
    X = varargin{2};
    Z = varargin{3};
    do3d = true;
else
    error('Incorrect number of arguments.');
end

% Prepare design and contrasts:
N = size(X,1);
O = ones(N,1);

if do2d,
    M = [X O];
    C = [1 0]';
else
    M = [X Z O];
    C = [1 0 0]';
end

% Do the GLM:
psi = M\Y;
res = Y - M*psi;
df  = N - rank(M);
tstat = psi'*C*(C'*pinv(M'*M)*C)^-.5 / sqrt((res'*res)/df);
pval  = tcdf(-tstat,df);

% Print results:
if pval <= 0.05,
    issig = '(significant)';
else
    issig = '(not significant)';
end
fprintf('Test statistic: %g\n',tstat);
fprintf('p-value: %g %s\n',pval,issig);
clf

% Make plot:
if do2d,
    yfunc = @(x) psi(1)*x + psi(2);
    x = linspace(min(X),max(X),21);
    plot(x,yfunc(x),'k--');
    hold on;
    scatter(X,Y,'.');
    hold off;
    axis([min(X) max(X) min(Y) max(Y)]);
    daspect([1 1 1]);
    xlabel('X');
    ylabel('Y');
    
elseif do3d,
    yfunc = @(x,z) psi(1)*x + psi(2)*z + psi(3);
    zfunc = @(x,y) -(psi(1)*x - y + psi(3))/psi(2);
    xvals = linspace(min(X),max(X),21)';
    yvals = linspace(min(Y),max(Y),21)';
    zvals = linspace(min(Z),max(Z),21)';
    
    subplot(2,2,4);
    [x,z] = meshgrid(xvals,zvals);
    subplot(2,2,4);
    hm = mesh(x,z,yfunc(x,z));
    hm.FaceColor = [1 1 1];
    hm.FaceAlpha = .5;
    hold on;
    hs = scatter3(X,Z,Y,'.');
    hpXY = plot3(...
        xvals,...
        ones(size(zvals))*mean(Z),...
        yfunc(xvals,mean(Z)),'k--');
    hpZY = plot3(...
        ones(size(xvals))*mean(X),...
        zvals,...
        yfunc(mean(X),zvals),'k--');
    hold off;
    axis([min(X) max(X) min(Z) max(Z) min(Y) max(Y)]);

    h = cell(3,1);
    for s = 1:3,
        h{s} = subplot(2,2,s);
        copyobj(hm,h{s});
        copyobj(hs,h{s});
        copyobj(hpXY,h{s});
        copyobj(hpZY,h{s});
    end
    for s = 1:4,
        subplot(2,2,s);
        daspect([1 1 1]);
        xlabel('X');
        ylabel('Z');
        zlabel('Y');
        axis([min(X) max(X) min(Z) max(Z) min(Y) max(Y)]);
    end
    subplot(2,2,1); view([0 -1 0]); camproj('orthographic'); title('Left view');
    subplot(2,2,2); view([1 0 0]);  camproj('orthographic'); title('Right view');
    subplot(2,2,3); view([0 0 1]);  camproj('orthographic'); title('Top view');
    subplot(2,2,4); view([20 45]);  camproj('perspective');  title('Oblique view');
end
