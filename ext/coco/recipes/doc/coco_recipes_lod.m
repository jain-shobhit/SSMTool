% Recipes for Continuation
% Harry Dankowicz and Frank Schilder
%
% List of toolbox development lines and demos ordered by chapter and
% section. The symbol * denotes content that is not included in Recipes for
% Continuation.
%
% Part I  Design Fundamentals
%
%   Chapter 1  A Continuation Paradigm
%     1.1    Problem formulation
%     1.2    An analytical solution
%     1.3    A numerical solution
%            <a href="matlab: coco_recipes_doc calcvar_demo calcvar_demo">calcvar_demo</a>
%
%   Chapter 2  Encapsulation
%     2.1    Solution measures and constraints
%     2.2    Continuation problems
%     2.3    Algorithm development
%
%   Chapter 3  Construction
%     3.1    Problem decomposition
%     3.2    Staged construction
%     3.3    The core interface
%            <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">cmds_demo</a>
%
%   Chapter 4  Toolbox Development
%     4.1    Toolbox constructors
%            <a href="matlab: coco_recipes_doc alg_v1 alg_v1">alg_v1</a> - <a href="matlab: coco_recipes_doc alg_v1_demo alg_v1_demo">alg_v1_demo</a>
%            <a href="matlab: coco_recipes_doc alg_v2 alg_v2">alg_v2</a> - <a href="matlab: coco_recipes_doc alg_v2_demo alg_v2_demo">alg_v2_demo</a>
%            <a href="matlab: coco_recipes_doc alg_v3 alg_v3">alg_v3</a> - <a href="matlab: coco_recipes_doc alg_v3_demo alg_v3_demo">alg_v3_demo</a>
%            <a href="matlab: coco_recipes_doc alg_v4 alg_v4">alg_v4</a> - <a href="matlab: coco_recipes_doc alg_v4_demo alg_v4_demo">alg_v4_demo</a>
%     4.2    Embeddability
%            <a href="matlab: coco_recipes_doc alg_v5 alg_v5">alg_v5</a> - <a href="matlab: coco_recipes_doc alg_v5_demo alg_v5_demo">alg_v5_demo</a>
%     4.3    Object-oriented design
%            <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">cmds_demo</a>
%
%   Chapter 5  Task Embedding
%     5.1    Tree decompositions
%     5.2    A composite toolbox
%            <a href="matlab: coco_recipes_doc compalg_v1 compalg_v1">compalg_v1</a> - <a href="matlab: coco_recipes_doc compalg_v1_demo compalg_v1_demo">compalg_v1_demo</a>
%            <a href="matlab: coco_recipes_doc compalg_v2 compalg_v2">compalg_v2</a> - <a href="matlab: coco_recipes_doc compalg_v2_demo compalg_v2_demo">compalg_v2_demo</a>
%     5.3    Toolbox settings
%            <a href="matlab: coco_recipes_doc alg_v6 alg_v6">alg_v6</a> - <a href="matlab: coco_recipes_doc alg_v6_demo alg_v6_demo">alg_v6_demo</a>
%
% Part II  Toolbox Templates
%
%   Chapter 6  Discretization
%     6.1    The collocation zero problem
%     6.2    Vectorization
%     6.3    A vectorized zero problem
%
%   Chapter 7  The Collocation Continuation Problem
%     7.1    Problem definition
%     7.2    Encoding
%            <a href="matlab: coco_recipes_doc coll_v1 coll_v1">coll_v1</a>
%     7.3    Examples
%            <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a>
%
%   Chapter 8  Single-Segment Continuation Problems
%     8.1    Boundary-value problems
%            <a href="matlab: coco_recipes_doc bvp_v1 bvp_v1">bvp_v1</a> - <a href="matlab: coco_recipes_doc bvp_v1_demo bvp_v1_demo">bvp_v1_demo</a>
%     8.2    Periodic orbits
%            <a href="matlab: coco_recipes_doc po_v1 po_v1">po_v1</a> - <a href="matlab: coco_recipes_doc po_v1_demo po_v1_demo">po_v1_demo</a>
%     8.3    Alternative embeddings
%            <a href="matlab: coco_recipes_doc bvp_v2 bvp_v2">bvp_v2</a> - <a href="matlab: coco_recipes_doc bvp_v2_demo bvp_v2_demo">bvp_v2_demo</a>
%
%   Chapter 9  Multi-Segment Continuation Problems
%     9.1    Boundary-value problems
%            <a href="matlab: coco_recipes_doc msbvp_v1 msbvp_v1">msbvp_v1</a>
%     9.2    Quasiperiodic invariant tori
%            <a href="matlab: coco_recipes_doc msbvp_v1_demo msbvp_v1_demo">msbvp_v1_demo</a>
%     9.3    Multi-segment periodic orbits
%            <a href="matlab: coco_recipes_doc coll_v2 coll_v2">coll_v2</a> - <a href="matlab: coco_recipes_doc coll_v2_demo coll_v2_demo">coll_v2_demo</a>*
%            <a href="matlab: coco_recipes_doc hspo_v1 hspo_v1">hspo_v1</a> - <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">hspo_v1_demo</a>
%
%   Chapter 10 The Variational Collocation Problem
%     10.1   The first variational equation
%            <a href="matlab: coco_recipes_doc varcoll_v1 varcoll_v1">varcoll_v1</a> - <a href="matlab: coco_recipes_doc varcoll_v1_demo varcoll_v1_demo">varcoll_v1_demo</a>
%     10.2   A variational zero problem
%            <a href="matlab: coco_recipes_doc varcoll_v2 varcoll_v2">varcoll_v2</a> - <a href="matlab: coco_recipes_doc varcoll_v2_demo varcoll_v2_demo">varcoll_v2_demo</a>
%     10.3   Candidate boundary conditions
%
% Part III Atlas Algorithms
%
%   Chapter 11 Covering Manifolds
%     11.1   Theory and terminology
%     11.2   A finite-state machine
%     11.3   An object-oriented implementation
%
%   Chapter 12 Single-Dimensional Atlas Algorithms
%     12.1   An advancing local cover
%            <a href="matlab: coco_recipes_doc atlas1d_v1 atlas1d_v1">atlas1d_v1</a> - <a href="matlab: coco_recipes_doc atlas1d_v1_demo atlas1d_v1_demo">atlas1d_v1_demo</a>
%     12.2   Adaptation and accelerated convergence
%            <a href="matlab: coco_recipes_doc atlas1d_v2 atlas1d_v2">atlas1d_v2</a> - <a href="matlab: coco_recipes_doc atlas1d_v2_demo atlas1d_v2_demo">atlas1d_v2_demo</a>
%     12.3   An expanding-boundary algorithm
%            <a href="matlab: coco_recipes_doc atlas1d_v3 atlas1d_v3">atlas1d_v3</a> - <a href="matlab: coco_recipes_doc atlas1d_v3_demo atlas1d_v3_demo">atlas1d_v3_demo</a>
%
%   Chapter 13 Multi-Dimensional Manifolds
%     13.1   A point-cloud algorithm
%            <a href="matlab: coco_recipes_doc atlas2d_v1 atlas2d_v1">atlas2d_v1</a> - <a href="matlab: coco_recipes_doc atlas2d_v1_demo atlas2d_v1_demo">atlas2d_v1_demo</a>*
%     13.2   A chart network
%            <a href="matlab: coco_recipes_doc atlas2d_v2 atlas2d_v2">atlas2d_v2</a> - <a href="matlab: coco_recipes_doc atlas2d_v2_demo atlas2d_v2_demo">atlas2d_v2_demo</a>
%            <a href="matlab: coco_recipes_doc atlas2d_v3 atlas2d_v3">atlas2d_v3</a> - <a href="matlab: coco_recipes_doc atlas2d_v3_demo atlas2d_v3_demo">atlas2d_v3_demo</a>*
%     13.3   Henderson's algorithm
%            <a href="matlab: coco_recipes_doc atlas2d_v4 atlas2d_v4">atlas2d_v4</a> - <a href="matlab: coco_recipes_doc atlas2d_v4_demo atlas2d_v4_demo">atlas2d_v4_demo</a>
%
%   Chapter 14 Computational Domains
%            <a href="matlab: coco_recipes_doc atlas1d_v3_demo atlas1d_v3_demo">atlas1d_v3_demo</a>
%     14.1   A 1-dimensional atlas algorithm
%            <a href="matlab: coco_recipes_doc atlas1d_v4 atlas1d_v4">atlas1d_v4</a> - <a href="matlab: coco_recipes_doc atlas1d_v4_demo atlas1d_v4_demo">atlas1d_v4_demo</a>
%            <a href="matlab: coco_recipes_doc atlas1d_v5 atlas1d_v5">atlas1d_v5</a> - <a href="matlab: coco_recipes_doc atlas1d_v5_demo atlas1d_v5_demo">atlas1d_v5_demo</a>
%     14.2   A 2-dimensional atlas algorithm
%            <a href="matlab: coco_recipes_doc atlas2d_v5 atlas2d_v5">atlas2d_v5</a> - <a href="matlab: coco_recipes_doc atlas2d_v5_demo atlas2d_v5_demo">atlas2d_v5_demo</a>
%            <a href="matlab: coco_recipes_doc atlas2d_v6 atlas2d_v6">atlas2d_v6</a> - <a href="matlab: coco_recipes_doc atlas2d_v6_demo atlas2d_v6_demo">atlas2d_v6_demo</a>
%     14.3   Manifolds of resonant periodic orbits
%            <a href="matlab: coco_recipes_doc atlas2d_v6_demo atlas2d_v6_demo">atlas2d_v6_demo</a>
%
% Part IV  Event Handling
%
%   Chapter 15 Special Points and Events
%     15.1   Theoretical framework
%     15.2   A continuation paradigm
%     15.3   The core interface
%            <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">cmds_demo</a>, <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">hspo_v1_demo</a>
%
%   Chapter 16 Atlas Events and Toolbox Integration
%     16.1   Event detection in atlas algorithms
%            <a href="matlab: coco_recipes_doc atlas1d_v6 atlas1d_v6">atlas1d_v6</a> - <a href="matlab: coco_recipes_doc atlas1d_v6_demo atlas1d_v6_demo">atlas1d_v6_demo</a>
%     16.2   A toolbox template
%            <a href="matlab: coco_recipes_doc alg_v7 alg_v7">alg_v7</a> - <a href="matlab: coco_recipes_doc alg_v7_demo alg_v7_demo">alg_v7_demo</a>
%     16.3   An alternative constructor
%            <a href="matlab: coco_recipes_doc alg_v8 alg_v8">alg_v8</a> - <a href="matlab: coco_recipes_doc alg_v8_demo alg_v8_demo">alg_v8_demo</a>*
%
%   Chapter 17 Event Handlers and Branch Switching
%     17.1   Toolbox event handlers
%            <a href="matlab: coco_recipes_doc alg_v9 alg_v9">alg_v9</a> - <a href="matlab: coco_recipes_doc alg_v9_demo alg_v9_demo">alg_v9_demo</a>
%            <a href="matlab: coco_recipes_doc alg_v10 alg_v10">alg_v10</a> - <a href="matlab: coco_recipes_doc alg_v10_demo alg_v10_demo">alg_v10_demo</a>
%     17.2   Bifurcations of periodic orbits
%            <a href="matlab: coco_recipes_doc po_v2 po_v2">po_v2</a> - <a href="matlab: coco_recipes_doc po_v2_demo po_v2_demo">po_v2_demo</a>*
%            <a href="matlab: coco_recipes_doc hspo_v2 hspo_v2">hspo_v2</a> - <a href="matlab: coco_recipes_doc hspo_v2_demo hspo_v2_demo">hspo_v2_demo</a>
%     17.3   Branch switching
%            <a href="matlab: coco_recipes_doc atlas1d_v7 atlas1d_v7">atlas1d_v7</a> - <a href="matlab: coco_recipes_doc atlas1d_v7_demo atlas1d_v7_demo">atlas1d_v7_demo</a>
%            <a href="matlab: coco_recipes_doc po_v2 po_v2">po_v2</a> - <a href="matlab: coco_recipes_doc po_v2_demo po_v2_demo">po_v2_demo</a>*
%            <a href="matlab: coco_recipes_doc hspo_v2 hspo_v2">hspo_v2</a> - <a href="matlab: coco_recipes_doc hspo_v2_demo hspo_v2_demo">hspo_v2_demo</a>
%
% Part V   Adaptation
%
%   Chapter 18 Pointwise Adaptation and Comoving Meshes
%     18.1   A brute-force approach
%            <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a>
%            <a href="matlab: coco_recipes_doc coll_v3 coll_v3">coll_v3</a> - <a href="matlab: coco_recipes_doc coll_v3_demo coll_v3_demo">coll_v3_demo</a>
%     18.2   (Co)moving meshes
%     18.3   A comoving mesh algorithm
%            <a href="matlab: coco_recipes_doc coll_v4 coll_v4">coll_v4</a>, <a href="matlab: coco_recipes_doc po_v3 po_v3">po_v3</a>* - <a href="matlab: coco_recipes_doc coll_v4_demo coll_v4_demo">coll_v4_demo</a>
%
%   Chapter 19 A Spectral Toolbox
%     19.1   The spectral continuation
%     19.2   Encoding
%            <a href="matlab: coco_recipes_doc dft_v1 dft_v1">dft_v1</a>
%     19.3   Adaptivity
%            <a href="matlab: coco_recipes_doc dft_v1 dft_v1">dft_v1</a> - <a href="matlab: coco_recipes_doc dft_v1_demo dft_v1_demo">dft_v1_demo</a>
%
%   Chapter 20 Integrating Adaptation in Atlas Algorithms
%     20.1   Adaptive mesh refinements
%            <a href="matlab: coco_recipes_doc atlas1d_v8 ">atlas1d_v8</a> - <a href="matlab: coco_recipes_doc atlas1d_v8_demo atlas1d_v8_demo">atlas1d_v8_demo</a>
%     20.2   Boundary-value problems
%            <a href="matlab: coco_recipes_doc coll_v5 coll_v5">coll_v5</a> - <a href="matlab: coco_recipes_doc coll_v5_demo coll_v5_demo">coll_v5_demo</a>
%            <a href="matlab: coco_recipes_doc coll_v6 coll_v6">coll_v6</a> - <a href="matlab: coco_recipes_doc coll_v6_demo coll_v6_demo">coll_v6_demo</a>
%     20.3   Numerical comparisons
%            <a href="matlab: coco_recipes_doc po_v3_demo po_v3_demo">po_v3_demo</a>
%            <a href="matlab: coco_recipes_doc coll_v7 coll_v7">coll_v7</a>
%            <a href="matlab: coco_recipes_doc canard_demo canard_demo">canard_demo</a>
%
% Part VI  Epilogue
%
%   Chapter 21 Toolbox projects
%     21.1   Calculus of variations
%     21.2   Nonlinear boundary conditions
%     21.3   Connecting orbits

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_recipes_lod.m 2842 2015-05-06 16:43:38Z hdankowicz $
