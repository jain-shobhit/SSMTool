function [atlas, cseg] = merge(atlas, cseg)
%MERGE   Merge curve segment into atlas
%
% Check if last chart in point list is neighbor with any atlas boundary
% chart and remove overlapping directions of continuation. Reconstitute
% atlas boundary.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge.m 3087 2019-04-04 19:54:09Z hdankowicz $

chart = cseg.ptlist{end};
R     = atlas.cont.R;
h     = atlas.cont.Rmarg*R;
nb    = cell(numel(chart.s),4); % Not all directions are allowed
for k=1:numel(chart.s)
  sk      = chart.G*chart.s(k);
  nrms    = norm(sk);
  hk      = h*nrms;
  sk      = sk/nrms;
  xk      = chart.x + hk*chart.TS*sk;
  nb(k,:) = {chart, xk, sk, hk};
end
for i=size(atlas.boundary,1):-1:1
  chart2 = atlas.boundary{i,1};
  if atlas.isneighbor(chart, chart2)
    x2 = atlas.boundary{i,2};
    if norm(chart.TSp'*(x2(atlas.cont.indcs)-chart.xp))<R
      atlas.boundary(i,:) = [];
    end
    for k=size(nb,1):-1:1
      x1 = nb{k,2};
      if norm(chart2.TSp'*(x1(atlas.cont.indcs)-chart2.xp))<R
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
