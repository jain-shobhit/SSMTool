function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 3120 2019-06-11 02:25:07Z hdankowicz $

defaults.R     = 0.1;  % Step size
defaults.PtMX  = 50;   % Maximum number of continuation steps
defaults.theta = 0.5;  % Theta method
defaults.almax = 10;   % Critical angle between successive tangent vectors
defaults.Rmarg = 0.95; % Margin for merging charts into boundary
defaults.Ndirs = 6;    % Number of available directions
cont           = coco_merge(defaults, cont);
cont.almax     = cont.almax*pi/180;
al             = (0:cont.Ndirs-1)*(2*pi/cont.Ndirs);
cont.s0        = [cos(al); sin(al)];
cont.bv0       = 1:cont.Ndirs; % Index array of available directions
cont.nb0       = [];           % Index array of neighboring chart ids

end
