function [prob, atlas, cseg, flush] = init_atlas(atlas, prob, cseg)
%INIT_ATLAS   Initialize atlas.
%
% Assign successfully located initial chart to point list and flush.
% Multiple possible directions of continuation.
%
% Differs from atlas2d_v1 by initializing network properties.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_atlas.m 3134 2019-07-13 15:13:42Z hdankowicz $

chart         = cseg.curr_chart;
chart.pt      = atlas.next_pt;  % Initialized to 0
atlas.next_pt = atlas.next_pt+1;
chart.R       = atlas.cont.R;
chart.s       = atlas.cont.s0;  % Multiple directions of continuation
chart.id      = chart.pt+1;
chart.bv      = atlas.cont.bv0; % Index array of available directions
chart.nb      = atlas.cont.nb0; % Index array of neighboring chart ids.
chart.pt_type = 'EP';
chart.ep_flag = 1;
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

flush         = true;

end
