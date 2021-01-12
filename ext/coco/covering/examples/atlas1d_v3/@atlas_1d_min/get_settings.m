function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 3113 2019-06-08 20:03:34Z hdankowicz $

defaults.R     = 0.1;  % Step size
defaults.PtMX  = 50;   % Maximum number of continuation steps
defaults.theta = 0.5;  % Theta method
defaults.almax = 10;   % Critical angle between successive tangent vectors
defaults.Rmarg = 0.95; % Margin for merging charts into boundary
cont           = coco_merge(defaults, cont);
cont.almax     = cont.almax*pi/180;

end
