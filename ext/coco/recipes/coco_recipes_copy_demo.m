function coco_recipes_copy_demo(demo)
% Copy demo folder to user-defined directory.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_recipes_copy_demo.m 2839 2015-03-05 17:09:01Z fschild $

recipes = fileparts(mfilename('fullpath'));
source = fullfile(recipes, demo);
start = textscan(userpath, '%s', 'Delimiter', ':');
work = uigetdir(start{1}{1},'Select target folder');
listing = dir(work);
if ~any(strcmp({listing.name},'coco_project_opts_recipes_user.m'))
  copyfile(fullfile(recipes, 'doc/coco_project_opts_recipes_user.m'), ...
    fullfile(work, 'coco_project_opts_recipes_user.m'))
end
target = fullfile(work, demo);
copyfile(source, target)
end
