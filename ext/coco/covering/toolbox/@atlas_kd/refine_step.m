function [prob, atlas, cseg, predict] = refine_step(atlas, prob, cseg)
%REFINE_STEP   Refine step size.
%
% Refine step size after convergence failed and try again. Flush and
% terminate if fail to converge with minimum step size.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

predict = false;
chart   = cseg.ptlist{1};
if chart.h > atlas.cont.R_min
  chart.h = max(atlas.cont.R_fac_min*chart.h, atlas.cont.R_min);
  atlas.charts{atlas.boundary(1)} = chart;
  atlas.first = true;
  predict = true;
end

end
