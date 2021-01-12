function [ode, spec] = ode_get_settings(prob, tbid, ode)
%ODE_GET_SETTINGS   Merge ODE toolbox family settings with default values.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_get_settings.m 2861 2015-07-26 22:08:48Z hdankowicz $

spec = {% Name    Type   Default  Action  Args  Description
  'vectorized', 'log|on',  true, 'switch', {}, 'enable/disable vectorized evaluation'
  'autonomous', 'log|on',  true, 'switch', {}, 'indicate whether ODE is autonomous or not'
  };

ode = coco_parse_settings(prob, spec, ode, tbid);

end
