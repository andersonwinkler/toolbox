function h=abline(b,m,varargin)
% FORMAT h = abline(b,m,...)
% Plots y=a*x+b in dotted line
%
% ...  Other graphics options
%
% Like Splus' abline
%
% To plot a horizontal line at y:    abline('h',y)   
% To plot a vertical line at x:      abline('v',x)   
%
% $Id: abline.m,v 1.1 2006/03/08 14:48:25 nichols Exp $

if (nargin==2) & isstr(b)
  b = lower(b);
else

  if (nargin<1)
    b = 0;
  end
  if (nargin<2)
    m = 0;
  end
  
end

XX=get(gca,'Xlim');
YY=get(gca,'Ylim');

g = [];
if isstr(b) & (b=='h')

  for i=1:length(m)
    g=[g;line(XX,[m(i) m(i)],'LineStyle',':',varargin{:})];
  end

elseif isstr(b) & (b=='v')

  for i=1:length(m)
    g=line([m(i) m(i)],YY,'LineStyle',':',varargin{:});
  end

else

  g=line(XX,m*XX+b,'LineStyle',':',varargin{:});

end

if (nargout>0)
  h=g;
end
