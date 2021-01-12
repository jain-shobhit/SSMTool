function [prob, atlas, cseg, flush] = add_chart(atlas, prob, cseg)
%ADD_CHART   Add chart to point list.
%
% Add successully located chart to curve segment point list and increment
% counter. Designate end point by 'EP' point type and flush.
%
% Differs from atlas1d_v1 by including support for step-size adaptation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: add_chart.m 3134 2019-07-13 15:13:42Z hdankowicz $

chart    = cseg.curr_chart;
chart.pt = chart.pt+1;
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

al = subspace(cseg.ptlist{1}.TSp, cseg.ptlist{end}.TSp); % Angle between successive tangent vectors
R  = atlas.base_chart.R; % Current step size
if al>atlas.cont.almax
  if R>atlas.cont.Rmin
    atlas.base_chart.R = max(atlas.cont.Rfred*R, atlas.cont.Rmin); % Reduce step size
    flush              = false; % Try again
  end
elseif al<=atlas.cont.almax/2
  cseg.ptlist{end}.R = min(atlas.cont.Rfinc*R, atlas.cont.Rmax); % Increase step size
end

end
