function coco_rm_this_path
%COCO_RM_THIS_PATH   Remove toolbox from search path.
%
% Function executed by recipes utility coco_use_recipes_toolbox.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_rm_this_path.m 2839 2015-03-05 17:09:01Z fschild $

rmpath(fileparts(mfilename('fullpath')));
end
