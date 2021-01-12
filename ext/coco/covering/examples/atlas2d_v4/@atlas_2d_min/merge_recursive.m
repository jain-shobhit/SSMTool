function [atlas, chart1, checked] = merge_recursive(atlas, chart1, k, checked)
%MERGE_RECURSIVE   Merge chart into chart network.
%
% Modify polygons and update network. Terminate when no longer overlapping
% branch through original neighboring chart.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge_recursive.m 3087 2019-04-04 19:54:09Z hdankowicz $

checked(end+1) = k;
chartk = atlas.charts{k};
if atlas.isclose(chart1, chartk) % If overlapping
  dx     = chartk.xp-chart1.xp;
  phi1   = chart1.TSp'*dx;
  phik   = chartk.TSp'*(-dx);
  test1  = chart1.v.*(chart1.s'*phi1)-norm(phi1)^2/2;
  testk  = chartk.v.*(chartk.s'*phik)-norm(phik)^2/2;
  flag1  = (test1>0);
  flagk  = (testk>0);
  chart1 = ... % Modify polygon
    atlas.subtract_half_space(chart1, test1, phi1, flag1, chartk.id);
  chartk = ... % Modify polygon
    atlas.subtract_half_space(chartk, testk, phik, flagk, chart1.id);
  atlas.charts{k} = chartk;
  check = setdiff(chartk.nb, checked);
  while ~isempty(check) % Call recursively with neighbors of chartk
    [atlas, chart1, checked] = ...
      atlas.merge_recursive(chart1, check(1), checked);
    check = setdiff(chartk.nb, checked);
  end
end

end
