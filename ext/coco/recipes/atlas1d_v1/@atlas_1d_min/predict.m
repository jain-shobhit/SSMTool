function [prob atlas cseg correct] = predict(atlas, prob, cseg)
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and pseudo-arclength
% predictor and correct.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: predict.m 2839 2015-03-05 17:09:01Z fschild $

chart       = atlas.base_chart;
prcond      = struct('x', chart.x, 'TS', chart.TS, ...
                     's', chart.s, 'h', chart.R);
xp          = chart.x+chart.R*(chart.TS*chart.s);
[prob cseg] = CurveSegment.create(prob, chart, prcond, xp);
correct     = true;

end
