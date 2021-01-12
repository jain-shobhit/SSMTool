function [prob, atlas, cseg, predict] = refine_step(atlas, prob, cseg)
%REFINE_STEP   Refine step size.
%
% Refine step size after convergence failed and try again. Flush and
% terminate if fail to converge with minimum step size.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: refine_step.m 3087 2019-04-04 19:54:09Z hdankowicz $

predict = false;
chart   = cseg.ptlist{1};
R       = chart.R;
if R > atlas.cont.Rmin
  chart.R = max(atlas.cont.Rfred*R, atlas.cont.Rmin);
  atlas.base_chart = chart;
  predict = true;
end

end
