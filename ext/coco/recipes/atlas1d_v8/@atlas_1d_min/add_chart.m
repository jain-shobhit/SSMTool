function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
%ADD_CHART   Add chart to point list.
%
% Add successully located chart to curve segment point list and increment
% counter. Designate end point by 'EP' point type and flush.
%
% Identical to atlas1d_v2.

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

al = subspace(cseg.ptlist{1}.TS, cseg.ptlist{end}.TS); % Angle between successive tangent vectors
R  = atlas.base_chart.R; % Current step size
if al>atlas.cont.almax
  if R>atlas.cont.hmin
    atlas.base_chart.R = max(atlas.cont.hfred*R, atlas.cont.hmin); % Reduce step size
    flush              = false; % Try again
  end
elseif al<=atlas.cont.almax/2
  cseg.ptlist{end}.R = min(atlas.cont.hfinc*R, atlas.cont.hmax); % Increase step size
end

end
