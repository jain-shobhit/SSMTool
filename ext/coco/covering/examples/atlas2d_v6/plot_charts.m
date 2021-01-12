function [tri, X, C] = plot_charts(charts, ix, iy, iz)
%PLOT_CHARTS   Plot atlas as individual charts.
%
% Plot routine compatible with chart structure of atlas2d_v6 atlas
% algorithm. Surface is triangularization of chart base points.
%
% [tri X C] = PLOT_CHARTS(CHARTS, IX, IY, IZ)
%
% TRI    - Triangular represenation of surface.
% X      - Chart base points.
% C      - Chart colors.
% CHARTS - Array of two-dimensional charts.
% IX     - Integer index for first coordinate.
% IY     - Integer index for second coordinate.
% IZ     - Integer index for third coordinate.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: plot_charts.m 3126 2019-06-12 02:46:20Z hdankowicz $

if nargin>=4
  f = @(x) x([ix iy iz])';
else
  f = @(x) x([1 2 3])';
end
is_boundary = @(x) any(x.nb==0);
boundary = [];
tri = [];
X   = [];
C   = [];
N   = numel(charts);
for k=1:N
  chart = charts{k};
  xx = chart.x;
  x  = f(xx);
  X  = [ X ; x ];
  xo = size(X,1);
  NN = size(chart.s,2);
  if ~isfield(chart, 'v')
    chart.v = chart.R*ones(size(chart.s,2),1);
  end
  nb = circshift(chart.nb,[0 -1]);
  for l=1:NN-1
    x   = f(xx+min(chart.R,chart.v(l))*(chart.TS*chart.G*chart.s(:,l)));
    X   = [ X ; x ];
    tri = [tri ; xo xo+l xo+l+1];
    if ~is_boundary(chart)
      C = [ C ; k   ];
    else
      C = [ C ; N+1 ];
    end
    if nb(l)==0
      boundary = [boundary ; [xo+l xo+l+1]];
    end
  end
  x   = f(xx+min(chart.R,chart.v(NN))*(chart.TS*chart.G*chart.s(:,NN)));
  X   = [ X ; x ];
  tri = [tri ; xo xo+NN xo+1];
  if ~is_boundary(chart)
    C = [ C ; k   ];
  else
    C = [ C ; N+1 ];
  end
  if nb(NN)==0
    boundary = [boundary ; [xo+NN xo+1]];
  end
end
if nargout==0
  % make sure color range is used
  M   = size(X,1);
  X   = [ X ; X(1,:) ];
  tri = [ tri ; M+1 M+1 M+1 ; M+1 M+1 M+1 ];
  C   = [ C ; 1 ; N+1 ];
  if nargin>=4
    colormap(gca, [lines(numel(charts)) ; 1 1 1]);
    trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'LineWidth', 0.5);
    xoff = 0.005;
    for k=1:size(boundary,1)
      x = X(boundary(k,:),1);
      y = X(boundary(k,:),2);
      z = X(boundary(k,:),3);
      plot3(x+xoff,y,z, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black');
    end
  end
  clear tri X C
end
end
