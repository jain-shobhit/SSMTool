function prob = coco_project_opts_recipes(prob)
%COCO_PROJECT_OPTS_RECIPES    Set options for recipes.
%
% This function defines default settings specific to the examples from
% Recipes for Continuation, which are applied to all scripts that are
% executed below the path of the location of COCO_PROJECT_OPTS_RECIPES.
% Type 'which coco_project_opts_recipes' in the command window to show the
% location.
%
% See also atlas_0d_recipes, atlas_1d_recipes

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
