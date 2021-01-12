function [prob, atlas, cseg, flush] = init_atlas(atlas, prob, cseg)
%INIT_ATLAS   Initialize atlas.
%
% Assign successfully located initial chart to point list and flush. Order
% of successive points on solution manifolds dependent on sign of
% direction.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_atlas.m 3087 2019-04-04 19:54:09Z hdankowicz $

chart          = cseg.curr_chart;
chart.R        = atlas.cont.R;
chart.pt_type  = 'EP';
chart.ep_flag  = 1;
if ~isempty(chart.s) && chart.s(1)<0
  prob = AtlasBase.bddat_set(prob, 'ins_mode', 'prepend');
else
  prob = AtlasBase.bddat_set(prob, 'ins_mode', 'append');
end
cseg.ptlist{1} = chart;
flush          = true;

end
