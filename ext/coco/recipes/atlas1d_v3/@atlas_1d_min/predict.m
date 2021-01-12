function [prob atlas cseg correct] = predict(atlas, prob, cseg)
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and theta-method predictor
% and correct.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: predict.m 2839 2015-03-05 17:09:01Z fschild $

[chart xp s h] = atlas.boundary{1,:};
prcond         = struct('x', chart.x, 'TS', chart.TS, 's', s, 'h', h);
th             = atlas.cont.theta;
if th>=0.5 && th<=1 % A two-step algorithm
  xp          = chart.x+(th*h)*(chart.TS*s);
  [prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
  [prob ch2]  = cseg.update_TS(prob, cseg.curr_chart);
  h           = h*(ch2.TS'*chart.TS);
  xp          = chart.x+h*(ch2.TS*s); % Compute new predictor
  prcond      = struct('x', chart.x, 'TS', ch2.TS, 's', s, 'h', h); % Compute new projection condition
end
[prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
correct     = true;

end
