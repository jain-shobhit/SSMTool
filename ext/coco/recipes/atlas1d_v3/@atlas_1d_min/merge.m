function [atlas cseg] = merge(atlas, cseg)
%MERGE   Merge curve segment into atlas.
%
% Check if last chart in point list is neighbor with any atlas boundary
% chart and remove overlapping directions of continuation. Reconstitute
% atlas boundary.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge.m 2839 2015-03-05 17:09:01Z fschild $

chart = cseg.ptlist{end};
R     = atlas.cont.h;
h     = atlas.cont.Rmarg*R;
nb    = cell(2,4);
for k=1:2
  sk      = chart.s(k);
  xk      = chart.x+h*(chart.TS*sk);
  nb(k,:) = {chart, xk, sk, h};
end
for i=size(atlas.boundary,1):-1:1
  chart2 = atlas.boundary{i,1};
  if atlas.isneighbor(chart, chart2)
    x2 = atlas.boundary{i,2};
    if norm(chart.TS'*(x2-chart.x))<R
      atlas.boundary(i,:) = [];
    end
    for k=size(nb,1):-1:1
      x1 = nb{k,2};
      if norm(chart2.TS'*(x1-chart2.x))<R
        nb(k,:) = [];
      end
    end
  end
end
atlas.boundary = [nb; atlas.boundary];
if isempty(atlas.boundary)
  chart.pt_type    = 'EP';
  chart.ep_flag    = 1;
  cseg.ptlist{end} = chart;
end

end
