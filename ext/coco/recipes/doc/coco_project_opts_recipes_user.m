function prob = coco_project_opts_recipes_user(prob)
%COCO_PROJECT_OPTS_RECIPES_USER    Use recipes settings for user project
%
% This function defines default settings identical to the examples from
% Recipes for Continuation. These are applied to all scripts that are
% executed below the path of the location of this function.

% See installation instructions for correct set-up.

prdir   = pwd;
recipes = fileparts(mfilename('fullpath'));
if strncmp(prdir, recipes, numel(recipes))
  prob = coco_set(prob, 'cont', 'linsolve', 'recipes');  % linear solver
  prob = coco_set(prob, 'cont', 'corrector', 'recipes'); % nonlinear solver
  prob = coco_set(prob, 'cont', 'atlas_classes', ...     % atlas algorithms
    { [] 'atlas_kd' ; 0 'atlas_0d_recipes' ; 1 'atlas_1d_recipes' });
end
end
