function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
%INIT_ATLAS   Initialize atlas.
%
% Assign successfully located initial chart to point list and flush. Order
% of successive points on solution manifolds dependent on sign of
% direction.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_atlas.m 2839 2015-03-05 17:09:01Z fschild $

chart          = cseg.curr_chart;
chart.R        = atlas.cont.h;
chart.pt_type  = 'EP';
chart.ep_flag  = 1;
if ~isempty(chart.s) && chart.s(1)<0
  prob = AtlasBase.bddat_set(prob, 'ins_mode', 'prepend');
  chart.t = -chart.TS; % Assign tangent vector
else
  prob = AtlasBase.bddat_set(prob, 'ins_mode', 'append');
  chart.t = chart.TS;  % Assign tangent vector
end
[prob chart]   = cseg.update_p(prob, chart); % Update continuation parameters that depend on tangent vector
cseg.ptlist{1} = chart;
flush          = true;

end
