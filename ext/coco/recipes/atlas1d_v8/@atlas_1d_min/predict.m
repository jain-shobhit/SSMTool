function [prob atlas cseg correct] = predict(atlas, prob, cseg)
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and theta method predictor
% and correct. Update tangent vector and continuation parameters for event
% detection.
%
% Differs from atlas1d_v6 by supporting the application of toolbox-specific
% remesh algorithms.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: predict.m 2839 2015-03-05 17:09:01Z fschild $

chart  = atlas.base_chart;
nad   = atlas.cont.NAdapt; % Remesh frequencey
RMMX  = atlas.cont.RMMX;   % Number of remesh iterations
if nad>0 && mod(chart.pt,nad)==0
  x0                 = chart.x;            % Current solution
  V0                 = [chart.t chart.TS]; % Current tangent vector and tangent space
  [prob chart x0 V0] = coco_remesh(prob, chart, x0, V0, RMMX); % Remesh solution, tangent vector, and tangent space
  nv                 = repmat(sqrt(sum(V0.^2,1)), [size(V0,1) 1]);
  chart.x            = x0;
  chart.t            = V0(:,1)./nv(:,1);
  chart.TS           = V0(:,2:end)./nv(:,2:end);
end
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
[prob cseg]    = CurveSegment.create(prob, chart, prcond, xp);
[prob chart]   = cseg.update_t(prob, cseg.ptlist{1}); % Update tangent vector for event detection
[prob chart]   = cseg.update_p(prob, chart);          % Update continuation parameters that depend on tangent vector
cseg.ptlist{1} = chart;
correct        = true;

end
