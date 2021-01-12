function [prob, atlas, cseg, correct] = predict(atlas, prob, cseg) 
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and theta-method predictor
% and correct.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

chart = atlas.charts{atlas.boundary(1)};

fname = fullfile(prob.run.dir, sprintf('func_data%d', chart.id));
S = load(fname, 'funcs');
prob  = restore_funcs(prob, S.funcs);

if ~atlas.first && atlas.cont.NAdapt*chart.pt > 0 ...
    && chart.remesh && mod(chart.remesh,atlas.cont.NAdapt) == 0
  atlas.first = true;
  
  x0        = chart.x;
  V0        = [ chart.t chart.TS ];
  [prob, chart, x0, V0] = coco_remesh(prob, chart, x0, V0, atlas.cont.RMMX);
  nv        = repmat(sqrt(sum(V0.^2,1)), size(V0,1), 1);
  chart.x   = x0;
  chart.t   = V0(:,1)./nv(:,1);
  chart.TS  = V0(:,2:end)./nv(:,2:end);
  chart.ics = prob.efunc.p_idx;
  prcond = struct('x', chart.x, 'TS', chart.TS, ...
    's', chart.s(:, chart.bv(1)), 'h', 0);
  [prob, cseg] = CurveSegment.create(prob, chart, prcond, chart.x);
else
  atlas.first = false;
  
  s      = chart.G*chart.s(:, chart.bv(1));
  nrms   = norm(s);
  h      = chart.h/atlas.cont.Rmarg*nrms;
  s      = s/nrms;
  prcond = struct('x', chart.x, 'TS', chart.TS, 's', s, 'h', h);
  th     = atlas.cont.theta;
  if th>=0.5 && th<=1
    xp           = chart.x + (th*h)*(chart.TS*s);
    [prob, cseg] = CurveSegment.create(prob, chart, prcond, xp);
    [prob, cht2] = cseg.update_TS(prob, cseg.curr_chart);
    s            = h*cht2.TS'*chart.TS*s;
    h            = norm(s);
    s            = s/h;
    xp           = chart.x + h*cht2.TS*s;
    prcond       = struct('x', chart.x, 'TS', cht2.TS, 's', s, 'h', h);
  else
    xp           = chart.x + h*chart.TS*s;
  end
  [prob, cseg]   = CurveSegment.create(prob, chart, prcond, xp);
  [prob, chart]  = cseg.update_t(prob, cseg.ptlist{1});
  [prob, chart]  = cseg.update_p(prob, chart);
  cseg.ptlist{1} = chart;
end
correct     = true;

end
