function [prob atlas cseg] = init_chart(atlas, prob, cseg)
%INIT_CHART   Initialize chart.
%
% Initialize chart after initial correction and add to curve segment point
% list.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_chart.m 2839 2015-03-05 17:09:01Z fschild $

chart         = cseg.curr_chart;
chart.pt      = atlas.next_pt;  % Initialized to 0
atlas.next_pt = atlas.next_pt+1;
chart.R       = atlas.cont.h;
chart.s       = atlas.cont.s0;  % Multiple directions of continuation
chart.id      = chart.pt+1;
chart.bv      = atlas.cont.bv0; % Index array of available vertices
chart.nb      = atlas.cont.nb0; % Index array of chart id's associated with polygonal edges
chart.v       = atlas.cont.v0;  % Distances to polygonal vertices
[prob cseg]     = cseg.add_chart(prob, chart);
cseg.curr_chart = cseg.ptlist{1};

end
