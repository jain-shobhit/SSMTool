function [atlas, cseg] = merge(atlas, cseg)
%MERGE   Merge curve segment into atlas
%
% Implement a recursive search to check if last chart in point list is
% close to any atlas chart, starting with boundary charts, and remove
% overlapping directions of continuation. Update chart network and
% reconstitute atlas boundary.
%
% Identical to atlas2d_v4.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge.m 3087 2019-04-04 19:54:09Z hdankowicz $

chart        = cseg.ptlist{end};
nbfunc       = @(x) atlas.isclose(chart, x);
bd_charts    = atlas.charts(atlas.boundary);
idx          = cellfun(nbfunc, bd_charts);
close_charts = atlas.boundary(idx); % Find close charts on the atlas boundary
checked      = [0, chart.id]; % 0 represents the outside
while ~isempty(close_charts) % Traverse recursively by using neighbor information
  [atlas, chart, checked] = ...
    atlas.merge_recursive(chart, close_charts(1), checked);
  close_charts = setdiff(close_charts, checked); % Remove already checked charts
end
atlas.charts   = [atlas.charts, {chart}];
atlas.boundary = [atlas.boundary, chart.id];
bd_charts      = atlas.charts(atlas.boundary);
idx            = cellfun(@(x) ~isempty(x.bv), bd_charts);
atlas.boundary = atlas.boundary(idx);
if isempty(atlas.boundary)
  chart.pt_type    = 'EP';
  chart.ep_flag    = 1;
  cseg.ptlist{end} = chart;
end

end
