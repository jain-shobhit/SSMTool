function [prob, atlas, cseg, flush] = add_chart(atlas, prob, cseg)
%ADD_CHART   Add chart to point list.
%
% Add successully located chart to curve segment point list and increment
% counter. Designate end point by 'EP' point type and flush. Terminate if
% successive charts do not qualify as "neighbors".
%
% Identical to atlas2d_v4.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: add_chart.m 3134 2019-07-13 15:13:42Z hdankowicz $

chart         = cseg.curr_chart;
chart.pt      = atlas.next_pt;
atlas.next_pt = atlas.next_pt+1; % Increment chart counter
chart.id      = chart.pt+1;
chart.s       = atlas.cont.s0;   % All directions
chart.bv      = atlas.cont.bv0;  % All directions are available
chart.nb      = atlas.cont.nb0;  % Not connected to network
chart.v       = atlas.cont.v0;   % Default polygonal size 
if chart.pt>=atlas.cont.PtMX
  chart.pt_type = 'EP';
  chart.ep_flag = 1;
end
[prob, cseg] = cseg.add_chart(prob, chart);
chart = cseg.ptlist{end};

% projection
chart.xp  = chart.x(chart.ics);
TSp       = chart.TS(chart.ics,:);
try
  chart.TSp = cseg.orth(TSp);
catch ME
  if (strcmp(ME.identifier,'MATLAB:square'))
    msg = ['The new chart is located at the singular point (', ...
      num2str(chart.xp'),') of the projected geometry.'];
    causeException = MException('MATLAB:init_atlas:singularity',msg);
    ME = addCause(ME,causeException);
  end
  rethrow(ME)
end
chart.G = (TSp' * TSp) \ (TSp' * chart.TSp);
cseg.ptlist{end} = chart;

flush       = true;

if ~atlas.isneighbor(cseg.ptlist{1}, cseg.ptlist{end}) % Check for gap in coverage
  cseg.ptlist{end}.pt_type = 'GAP';
  cseg.ptlist{end}.ep_flag = 2;
  cseg.Status              = cseg.CurveSegmentCorrupted;
end

end
