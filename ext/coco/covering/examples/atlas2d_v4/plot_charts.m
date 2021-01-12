function [tri, X, C] = plot_charts(charts, ix, iy, iz, scale)
%PLOT_CHARTS   Plot atlas as individual charts.
%
% Plot routine compatible with chart structure of atlas2d_v6 atlas
% algorithm. Surface is triangularization of chart base points.
%
% [tri X C] = PLOT_CHARTS(CHARTS, IX, IY, IZ, SCALE)
%
% TRI    - Triangular represenation of surface.
% X      - Chart base points.
% C      - Chart colors.
% CHARTS - Array of two-dimensional charts.
% IX     - Integer index for first coordinate. 
% IY     - Integer index for second coordinate. 
% IZ     - Integer index for third coordinate.
% SCALE  - Plot scale for individual polygons.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: plot_charts.m 3122 2019-06-11 02:48:45Z hdankowicz $

if nargin>=4
  f = @(x) x([ix iy iz])';
else
  f = @(x) x([1 2 3])';
end
if nargin<5
  scale = 1;
end
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
  if isfield(chart, 'v')
    chart.v = scale*chart.v;
  else
    chart.v = (scale*chart.R)*ones(size(chart.s,2),1);
  end
  for l=1:NN-1
    x   = f(xx+min(chart.R,chart.v(l))*(chart.TS*chart.G*chart.s(:,l)));
    X   = [ X ; x ];
    tri = [tri ; xo xo+l xo+l+1];
    if chart.ep_flag<=1
      C = [ C ; k   ];
    else
      C = [ C ; N+1 ];
    end
  end
  x   = f(xx+min(chart.R,chart.v(NN))*(chart.TS*chart.G*chart.s(:,NN)));
  X   = [ X ; x ];
  tri = [tri ; xo xo+NN xo+1];
  if chart.ep_flag<=1
    C = [ C ; k   ];
  else
    C = [ C ; N+1 ];
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
    trisurf(tri, X(:,1), X(:,2), X(:,3), C);
    if numel(charts)<10
      for k=1:numel(charts)
        text(charts{k}.x(ix), charts{k}.x(iy), charts{k}.x(iz), ...
          sprintf('%d', charts{k}.id), ...
          'color', 'white', 'FontWeight', 'bold');
      end
    end
  end
  clear tri X C
end
end
