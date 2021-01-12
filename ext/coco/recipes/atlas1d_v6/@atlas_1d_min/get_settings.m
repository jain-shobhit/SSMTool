function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.FP    = true; % Detect fold points
defaults.BP    = true; % Detect branch points
defaults.h     = 0.1;  % Step size
defaults.PtMX  = 50;   % Maximum number of continuation steps
defaults.theta = 0.5;  % Theta method
defaults.almax = 10;   % Critical angle between successive tangent vectors
defaults.Rmarg = 0.95; % Margin for merging charts into boundary
cont           = coco_merge(defaults, cont);
cont.almax     = cont.almax*pi/180;

end
