function [prob, atlas, cseg] = init_chart(atlas, prob, cseg)
%INIT_CHART   Initialize chart.
%
% Initialize chart after initial correction and add to curve segment point
% list.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_chart.m 3134 2019-07-13 15:13:42Z hdankowicz $

chart         = cseg.curr_chart;
chart.pt      = atlas.next_pt;  % Initialized to 0
atlas.next_pt = atlas.next_pt+1;
chart.R       = atlas.cont.R;
chart.s       = atlas.cont.s0;  % Multiple directions of continuation
chart.id      = chart.pt+1;
chart.bv      = atlas.cont.bv0; % Index array of available vertices
chart.nb      = atlas.cont.nb0; % Index array of chart id's associated with polygonal edges
chart.v       = atlas.cont.v0;  % Distances to polygonal vertices
[prob, cseg]  = cseg.add_chart(prob, chart);
chart = cseg.ptlist{1};

% projection
chart.ics = atlas.cont.indcs;
chart.xp  = chart.x(chart.ics);
TSp       = chart.TS(chart.ics,:);
try
  chart.TSp = cseg.orth(TSp);
catch ME
  if (strcmp(ME.identifier,'MATLAB:square'))
    msg = ['The initial chart is located at the singular point (', ...
      num2str(chart.xp'),') of the projected geometry.'];
    causeException = MException('MATLAB:init_atlas:singularity',msg);
    ME = addCause(ME,causeException);
  end
  rethrow(ME)
end
chart.G = (TSp' * TSp) \ (TSp' * chart.TSp);
cseg.ptlist{1} = chart;
cseg.curr_chart = cseg.ptlist{1};

end
