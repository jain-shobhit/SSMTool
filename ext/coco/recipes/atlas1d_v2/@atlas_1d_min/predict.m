function [prob atlas cseg correct] = predict(atlas, prob, cseg)
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and theta method
% predictor and correct.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: predict.m 2839 2015-03-05 17:09:01Z fschild $

chart  = atlas.base_chart;
prcond = struct('x', chart.x, 'TS', chart.TS, ...
                's', chart.s, 'h', chart.R);
th     = atlas.cont.theta;
if th>=0.5 && th<=1 % A two-step algorithm
  xp          = chart.x+(th*chart.R)*(chart.TS*chart.s);
  [prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
  [prob ch2]  = cseg.update_TS(prob, cseg.curr_chart);
  h           = chart.R*chart.TS'*ch2.TS;
  xp          = chart.x+h*(ch2.TS*chart.s); % Compute new predictor
  prcond      = struct('x', chart.x, 'TS', ch2.TS, ...
                       's', chart.s, 'h', h); % Compute new projection condition
else
  xp          = chart.x+chart.R*(chart.TS*chart.s);
end
[prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
correct     = true;

end
