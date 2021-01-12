function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.h     = 0.1;  % Initial step size
defaults.PtMX  = 50;   % Maximum number of continuation steps
defaults.theta = 0.5;  % Theta predictor
defaults.hmax  = 0.1;  % Maximum step size
defaults.hmin  = 0.01; % Minimum step size
defaults.hfinc = 1.1;  % Factor of step size increment
defaults.hfred = 0.5;  % Factor of step size decrement
defaults.almax = 10;   % Maximum angle between successive tangent vectors
cont           = coco_merge(defaults, cont);
cont.almax     = cont.almax*pi/180;

end
