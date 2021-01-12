function [atlas cseg] = merge(atlas, cseg)
%MERGE   Merge curve segment into atlas
%
% Implement a recursive search to check if last chart in point list is
% neighbor with any atlas chart and remove overlapping directions of
% continuation. Update chart network and reconstitute atlas boundary.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge.m 2839 2015-03-05 17:09:01Z fschild $

chart        = cseg.ptlist{end};
nbfunc       = @(x) atlas.isneighbor(chart, x);
close_charts = find(cellfun(nbfunc, atlas.charts)); % Find neighboring charts
checked      = chart.id; % Track already checked charts
while ~isempty(close_charts) % Traverse recursively by using neighbor information
  [atlas chart checked] = ...
    atlas.merge_recursive(chart, close_charts(1), checked);
  close_charts = setdiff(close_charts, checked); % Remove already checked charts
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
