function [atlas, chart1, checked] = merge_recursive(atlas, chart1, k, checked)
%MERGE_RECURSIVE   Merge chart into chart network.
%
% Remove overlapping directions and update network. Terminate when no
% longer overlapping branch through original neighboring chart.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge_recursive.m 3120 2019-06-11 02:25:07Z hdankowicz $

checked(end+1) = k;
chartk = atlas.charts{k};
if atlas.isneighbor(chart1, chartk) % if overlapping
  R   = atlas.cont.R;
  h   = atlas.cont.Rmarg*R;
  v   = @(i,ch) ch.xp + h*(ch.TSp*ch.s(:,i));
  v1  = arrayfun(@(i) v(i,chart1), chart1.bv, 'UniformOutput', false);
  idx = cellfun(@(x) (norm(chartk.TSp'*(x-chartk.xp))<R), v1);
  chart1.bv(idx)   = []; % Remove unavailable directions
  chart1.nb        = [chart1.nb, chartk.id]; % Add to nearest neighbors
  vk  = arrayfun(@(i) v(i,chartk), chartk.bv, 'UniformOutput', false);
  idx = cellfun(@(x) (norm(chart1.TSp'*(x-chart1.xp))<R), vk);
  chartk.bv(idx)  = []; % Remove unavailable directions
  chartk.nb       = [chartk.nb, chart1.id]; % Add to nearest neighbors
  atlas.charts{k} = chartk;
  check = setdiff(chartk.nb, checked);
  while ~isempty(check) % Call recursively with neighbors of chartk
    [atlas, chart1, checked] = ...
      atlas.merge_recursive(chart1, check(1), checked);
    check = setdiff(chartk.nb, checked); % Remove checked neighbors
  end
end

end
