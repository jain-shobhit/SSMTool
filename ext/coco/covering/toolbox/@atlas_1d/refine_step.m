function [prob, atlas, cseg, predict] = refine_step(atlas, prob, cseg)
% correction failed, refine step size
chart = cseg.ptlist{1};
if chart.R <= atlas.cont.h_min
  % no convergence with minimum step size -> stop continuation in
  % current direction
  predict = false;
else
  % retry with reduced step size
  chart.R = max(atlas.cont.h_min, atlas.cont.h_fac_min*chart.R);
  atlas.chart_list{1} = chart;
  atlas.first = true;
  predict = true;
  coco_warn(prob, 3, atlas.cont.LogLevel, ...
    'no convergence, retry with step size h = % .2e\n', chart.R);
end
end
