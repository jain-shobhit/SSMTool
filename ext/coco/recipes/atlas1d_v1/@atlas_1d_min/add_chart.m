function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
%ADD_CHART   Add chart to point list.
%
% Add successully located chart to curve segment point list and increment
% counter. Designate end point by 'EP' point type and flush.

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

end
