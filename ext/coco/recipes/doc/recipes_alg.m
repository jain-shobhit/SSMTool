% The 'alg' toolbox.
% Recipes for Continuation
%
% Constructors and instance functions for algebraic continuation problem
% objects expressed in terms of problem variables and problem parameters,
% without distinguishing information on meaning of zero functions or origin
% of problem variables.
%
% Toolbox development lines:
%   <a href="matlab: coco_recipes_doc alg_v1 alg_v1">alg_v1</a>       - Basic toolbox constructor.
%   <a href="matlab: coco_recipes_doc alg_v2 alg_v2">alg_v2</a>       - Computing the linearization.
%   <a href="matlab: coco_recipes_doc alg_v3 alg_v3">alg_v3</a>       - Signals and slot functions.
%   <a href="matlab: coco_recipes_doc alg_v4 alg_v4">alg_v4</a>       - Generalized toolbox constructors.
%   <a href="matlab: coco_recipes_doc alg_v5 alg_v5">alg_v5</a>       - Embeddable toolbox constructors.
%   <a href="matlab: coco_recipes_doc alg_v6 alg_v6">alg_v6</a>       - Optional toolbox settings.
%   <a href="matlab: coco_recipes_doc alg_v7 alg_v7">alg_v7</a>       - Event detection and continuation.
%     +- <a href="matlab: coco_recipes_doc alg_v8 alg_v8">alg_v8</a>  - Moore-Spence fold continuation.
%   <a href="matlab: coco_recipes_doc alg_v9 alg_v9">alg_v9</a>       - Basic bifurcation analysis.
%   <a href="matlab: coco_recipes_doc alg_v10 alg_v10">alg_v10</a>      - A template event handler.
%
% Toolbox demos:
%   <a href="matlab: coco_recipes_doc alg_v1_demo alg_v1_demo">alg_v1_demo</a>  - Basic calling syntax and anonymous functions.
%   <a href="matlab: coco_recipes_doc alg_v2_demo alg_v2_demo">alg_v2_demo</a>  - Optional specification of explicit Jacobians.
%   <a href="matlab: coco_recipes_doc alg_v3_demo alg_v3_demo">alg_v3_demo</a>  - Automatic storing of toolbox instance data.
%   <a href="matlab: coco_recipes_doc alg_v4_demo alg_v4_demo">alg_v4_demo</a>  - Generalized constructors and argument parsing.
%   <a href="matlab: coco_recipes_doc alg_v5_demo alg_v5_demo">alg_v5_demo</a>  - Embeddable toolboxes and instance identifiers.
%   <a href="matlab: coco_recipes_doc alg_v6_demo alg_v6_demo">alg_v6_demo</a>  - Optional personalization of bifurcation data.
%   <a href="matlab: coco_recipes_doc alg_v7_demo alg_v7_demo">alg_v7_demo</a>  - Fold detection, localization, and continuation.
%   <a href="matlab: coco_recipes_doc alg_v8_demo alg_v8_demo">alg_v8_demo</a>  - Fold continuation with Moore-Spence system.
%   <a href="matlab: coco_recipes_doc alg_v9_demo alg_v9_demo">alg_v9_demo</a>  - Detection of Hopf- and neutral saddle points.
%   <a href="matlab: coco_recipes_doc alg_v10_demo alg_v10_demo">alg_v10_demo</a> - Distinction between Hopf- and neutral saddle points.
%
% See also <a href="matlab: coco_recipes_doc doc recipes_compalg">compalg</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: recipes_alg.m 2839 2015-03-05 17:09:01Z fschild $
