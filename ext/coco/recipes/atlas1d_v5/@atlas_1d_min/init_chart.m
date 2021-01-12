function [prob atlas cseg] = init_chart(atlas, prob, cseg)
%INIT_CHART   Initialize chart.
%
% Initialize chart after initial correction and add to curve segment point
% list.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_chart.m 2839 2015-03-05 17:09:01Z fschild $

chart           = cseg.curr_chart;
chart.pt        = 0;
chart.s         = [-1, 1];
atlas.cont.PtMX = abs(atlas.cont.PtMX);
[prob cseg]     = cseg.add_chart(prob, chart);
cseg.curr_chart = cseg.ptlist{1};

end
