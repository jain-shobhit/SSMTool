function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 3134 2019-07-13 15:13:42Z hdankowicz $

defaults.R     = 0.1;  % Initial step size
defaults.PtMX  = 50;   % Maximum number of continuation steps
defaults.theta = 0.5;  % Theta predictor
defaults.Rmax  = 0.1;  % Maximum step size
defaults.Rmin  = 0.01; % Minimum step size
defaults.Rfinc = 1.1;  % Factor of step size increment
defaults.Rfred = 0.5;  % Factor of step size decrement
defaults.almax = 5;   % Maximum angle between successive tangent vectors
cont           = coco_merge(defaults, cont);
cont.almax     = cont.almax*pi/180;

end
