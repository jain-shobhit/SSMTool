function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 3110 2019-06-08 04:25:59Z hdankowicz $

defaults.R    = 0.1; % Step size
defaults.PtMX = 50 ; % Maximum number of continuation steps
cont          = coco_merge(defaults, cont);

end
