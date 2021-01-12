function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
%INIT_ATLAS   Initialize atlas.
%
% Assign successfully located initial chart to point list and flush.
% Two possible directions of continuation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_atlas.m 2839 2015-03-05 17:09:01Z fschild $

chart           = cseg.curr_chart;
chart.pt        = 0;
chart.R         = atlas.cont.h;
chart.s         = [1, -1]*sign(atlas.cont.PtMX);
atlas.cont.PtMX = abs(atlas.cont.PtMX);
chart.pt_type   = 'EP';
chart.ep_flag   = 1;
[prob cseg]     = cseg.add_chart(prob, chart);
flush           = true;

end
