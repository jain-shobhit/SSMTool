function [atlas cseg] = merge(atlas, cseg)
%MERGE   Merge curve segment into atlas
%
% Check if last chart in point list is neighbor with any atlas chart and
% remove overlapping directions of continuation. Update chart network and
% reconstitute atlas boundary.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge.m 2839 2015-03-05 17:09:01Z fschild $

chart        = cseg.ptlist{end};
R            = atlas.cont.h;
h            = atlas.cont.Rmarg*R;
v            = @(i,ch) ch.x+h*(ch.TS*ch.s(:,i));
nbfunc       = @(x) atlas.isneighbor(chart, x);
close_charts = find(cellfun(nbfunc, atlas.charts)); % Find neighboring charts
for k=close_charts % Traverse linearly over all neighboring charts
  chartk          = atlas.charts{k}; 
  vx  = arrayfun(@(i) v(i,chart), chart.bv, 'UniformOutput', false);
  idx = cellfun(@(x) (norm(chartk.TS'*(x-chartk.x))<R), vx);
  chart.bv(idx)   = []; % Remove unavailable directions
  chart.nb        = [chart.nb, chartk.id]; % Add to nearest neighbors
  vx  = arrayfun(@(i) v(i,chartk), chartk.bv, 'UniformOutput', false);
  idx = cellfun(@(x) (norm(chart.TS'*(x-chart.x))<R), vx);
  chartk.bv(idx)  = []; % Remove unavailable directions
  chartk.nb       = [chartk.nb, chart.id]; % Add to nearest neighbors
  atlas.charts{k} = chartk;
end
atlas.charts   = [atlas.charts, {chart}];
atlas.boundary = [chart.id, atlas.boundary];
bd_charts      = atlas.charts(atlas.boundary);
idx            = cellfun(@(x) ~isempty(x.bv), bd_charts);
atlas.boundary = atlas.boundary(idx); % Retain charts with available directions
if isempty(atlas.boundary)
  chart.pt_type    = 'EP';
  chart.ep_flag    = 1;
  cseg.ptlist{end} = chart;
end

end
