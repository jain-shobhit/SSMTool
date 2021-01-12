function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
%ADD_CHART   Add chart to point list.
%
% Add successully located chart to curve segment point list and increment
% counter. Designate end point by 'EP' point type and flush. Terminate if
% successive charts do not qualify as "neighbors".
%
% Identical to atlas1d_v3.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: add_chart.m 2839 2015-03-05 17:09:01Z fschild $

chart    = cseg.curr_chart;
chart.pt = chart.pt+1;
if chart.pt>=atlas.cont.PtMX
  chart.pt_type = 'EP';
  chart.ep_flag = 1;
end
[prob cseg] = cseg.add_chart(prob, chart);
flush       = true;

if ~atlas.isneighbor(cseg.ptlist{1}, cseg.ptlist{end}) % Check for gap in coverage
  cseg.ptlist{end}.pt_type = 'GAP';
  cseg.ptlist{end}.ep_flag = 2;
  cseg.Status              = cseg.CurveSegmentCorrupted;
end

end
