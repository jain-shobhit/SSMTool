function [prob atlas cseg correct] = predict(atlas, prob, cseg)
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and theta-method predictor
% and correct.
%
% Identical to atlas2d_v2.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: predict.m 2839 2015-03-05 17:09:01Z fschild $

chart  = atlas.charts{atlas.boundary(1)};
s      = chart.s(:,chart.bv(1)); % First available boundary vertex
h      = atlas.cont.Rmarg*chart.R;
prcond = struct('x', chart.x, 'TS', chart.TS, 's', s, 'h', h);
th     = atlas.cont.theta;
if th>=0.5 && th<=1 % A two-step algorithm
  xp          = chart.x+(th*h)*(chart.TS*s);
  [prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
  [prob ch2]  = cseg.update_TS(prob, cseg.curr_chart);
  s           = h*(ch2.TS'*chart.TS*s);
  h           = norm(s);
  s           = s/h;
  xp          = chart.x+h*(ch2.TS*s); % Compute new predictor
  prcond      = struct('x', chart.x, 'TS', ch2.TS, 's', s, 'h', h); % Compute new projection condition
else
  xp          = chart.x+h*(chart.TS*s);
end
[prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
correct     = true;

end
