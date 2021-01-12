function cont = get_settings(cont)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.h    = 0.1; % Step size
defaults.PtMX = 50 ; % Maximum number of continuation steps
cont          = coco_merge(defaults, cont);

end
