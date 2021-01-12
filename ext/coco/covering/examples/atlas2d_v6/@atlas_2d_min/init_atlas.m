function [prob, atlas, cseg, flush] = init_atlas(atlas, prob, cseg)
%INIT_ATLAS   Initialize atlas.
%
% Assign successfully located initial chart to point list and flush.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_atlas.m 3087 2019-04-04 19:54:09Z hdankowicz $

chart          = cseg.curr_chart;
chart.pt_type  = 'EP';
chart.ep_flag  = 1;
cseg.ptlist{1} = chart;
flush          = true;

end
