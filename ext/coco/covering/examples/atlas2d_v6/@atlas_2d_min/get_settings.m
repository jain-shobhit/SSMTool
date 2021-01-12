function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.
%
% Identical to atlas2d_v4.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 3124 2019-06-12 01:50:23Z hdankowicz $

defaults.R     = 0.1;  % Step size
defaults.PtMX  = 50;   % Maximum number of continuation steps
defaults.theta = 0.5;  % Theta method
defaults.almax = 10;   % Critical angle between successive tangent vectors
defaults.Rmarg = 0.95; % Margin for merging charts into boundary
defaults.Ndirs = 6;    % Number of available directions
cont           = coco_merge(defaults, cont);
cont.almax     = cont.almax*pi/180;
cont.Ndirs     = max(3, ceil(cont.Ndirs)); % Minimum of three directions (triangle)
al             = (0:cont.Ndirs-1)*(2*pi/cont.Ndirs);
cont.s0        = [cos(al); sin(al)];
cont.bv0       = 1:cont.Ndirs;          % Index array of available boundary vertices
cont.nb0       = zeros(1,cont.Ndirs);   % Edges identified with neighboring charts
r1             = cont.R/cos(pi/cont.Ndirs);
cont.v0        = r1*ones(cont.Ndirs,1); % Distances to polygonal vertices

end
