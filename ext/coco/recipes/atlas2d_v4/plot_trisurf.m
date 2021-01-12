function [tri X C] = plot_trisurf(charts, ix, iy, iz)
%PLOT_TRISURF   Plot atlas as triangulated surface.
%
% Plot routine compatible with chart structure of atlas2d_v6 atlas
% algorithm. Surface is triangularization of chart base points.
%
% [tri X C] = PLOT_TRISURF(CHARTS, IX, IY, IZ)
%
% TRI    - Triangular represenation of surface.
% X      - Chart base points.
% C      - Chart colors.
% CHARTS - Array of two-dimensional charts.
% IX     - Integer index for first coordinate. 
% IY     - Integer index for second coordinate. 
% IZ     - Integer index for third coordinate.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: plot_trisurf.m 2933 2015-11-02 04:11:21Z hdankowicz $

if nargin>=4
  f = @(x) x([ix iy iz])';
else
  f = @(x) x';
end
tri = [];
X   = [];
C   = [];
ghostcharts = find(cellfun(@(x)x.ep_flag>1,charts));
N   = numel(charts);
for k=1:N
  chart = charts{k};
  X     = [ X ; f(chart.x) ];
  ic    = [chart.nb chart.nb(1)];
  ix    = chart.id;
  for l=1:numel(ic)-1
    face = sort([ix ic(l) ic(l+1)]);
    if all(face>0) && (isempty(tri) ||  ~ismember(face, tri, 'rows'))
      tri  = [tri ; face];
      if any(ismember(face,ghostcharts)) % Use last color
        C = [ C ; N+1 ];
      else
        C = [ C ; k   ];
      end
    end
  end
end
if nargout==0
  % make sure color range is used
  X   = [ X ; X(1,:) ];
  tri = [ tri ; N+1 N+1 N+1 ; N+1 N+1 N+1 ];
  C   = [ C ; 1 ; N+1 ];
  if nargin>=4
    colormap(gca, [lines(numel(charts)) ; 1 1 1]);
    trisurf(tri, X(:,1), X(:,2), X(:,3), C);
  end
  clear tri X C
end
end
