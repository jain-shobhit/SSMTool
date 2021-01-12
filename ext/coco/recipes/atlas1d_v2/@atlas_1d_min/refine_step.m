function [prob atlas cseg predict] = refine_step(atlas, prob, cseg)
%REFINE_STEP   Refine step size.
%
% Refine step size after convergence failed and try again. Flush and
% terminate if fail to converge with minimum step size.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: refine_step.m 2839 2015-03-05 17:09:01Z fschild $

predict = false;
chart   = cseg.ptlist{1};
R       = chart.R;
if R > atlas.cont.hmin
  chart.R = max(atlas.cont.hfred*R, atlas.cont.hmin);
  atlas.base_chart = chart;
  predict = true;
end

end
