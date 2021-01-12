% Recipes for Continuation
% Harry Dankowicz and Frank Schilder
%
% List of figure demos ordered by chapter and section. The symbol * denotes
% content that is not included in Recipes for Continuation.
%
% Part I  Design Fundamentals
%
%   Chapter 1  A Continuation Paradigm
%     1.1    Problem formulation
%     1.2    An analytical solution
%     1.3    A numerical solution
%            <a href="matlab: coco_recipes_doc calcvar_demo calcvar_demo">Figures 1.5-1.8</a>
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
%            <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">Figures 3.2-3.3</a>
%
%   Chapter 4  Toolbox Development
%     4.1    Toolbox constructors
%     4.2    Embeddability
%     4.3    Object-oriented design
%
%   Chapter 5  Task Embedding
%     5.1    Tree decompositions
%     5.2    A composite toolbox
%     5.3    Toolbox settings
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
%     7.3    Examples
%            <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">Figures 7.1-7.4</a>
%
%   Chapter 8  Single-Segment Continuation Problems
%     8.1    Boundary-value problems
%            <a href="matlab: coco_recipes_doc bvp_v1_demo bvp_v1_demo">Figures 8.1-8.2</a>
%     8.2    Periodic orbits
%            <a href="matlab: coco_recipes_doc po_v1_demo po_v1_demo">Figures 8.3-8.6</a>
%     8.3    Alternative embeddings
%            <a href="matlab: coco_recipes_doc bvp_v2_demo bvp_v2_demo">Figure 8.7</a>
%
%   Chapter 9  Multi-Segment Continuation Problems
%     9.1    Boundary-value problems
%     9.2    Quasiperiodic invariant tori
%            <a href="matlab: coco_recipes_doc msbvp_v1_demo msbvp_v1_demo">Figures 9.1-9.2</a>
%     9.3    Multi-segment periodic orbits
%            <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">Figures 9.3-9.5</a>
%
%   Chapter 10 The Variational Collocation Problem
%     10.1   The first variational equation
%     10.2   A variational zero problem
%            <a href="matlab: coco_recipes_doc varcoll_v2_demo varcoll_v2_demo">Figures 10.1-10.5</a>
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
%     12.2   Adaptation and accelerated convergence
%     12.3   An expanding-boundary algorithm
%
%   Chapter 13 Multi-Dimensional Manifolds
%     13.1   A point-cloud algorithm
%     13.2   A chart network
%            <a href="matlab: coco_recipes_doc atlas2d_v2_demo atlas2d_v2_demo">Figure 13.1</a>
%     13.3   Henderson's algorithm
%            <a href="matlab: coco_recipes_doc atlas2d_v4_demo atlas2d_v4_demo">Figures 13.3-13.4</a>
%
%   Chapter 14 Computational Domains
%     14.1   A 1-dimensional atlas algorithm
%     14.2   A 2-dimensional atlas algorithm
%            <a href="matlab: coco_recipes_doc atlas2d_v5_demo atlas2d_v5_demo">Figure 14.1</a>
%            <a href="matlab: coco_recipes_doc atlas2d_v6_demo atlas2d_v6_demo">Figure 14.2</a>
%     14.3   Manifolds of resonant periodic orbits
%            <a href="matlab: coco_recipes_doc atlas2d_v6_demo atlas2d_v6_demo">Figures 14.3-14.4</a>
%
% Part IV  Event Handling
%
%   Chapter 15 Special Points and Events
%     15.1   Theoretical framework
%     15.2   A continuation paradigm
%     15.3   The core interface
%            <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">Figures 15.4-15.7</a>
%
%   Chapter 16 Atlas Events and Toolbox Integration
%     16.1   Event detection in atlas algorithms
%            <a href="matlab: coco_recipes_doc atlas1d_v6_demo atlas1d_v6_demo">Figures 16.1-16.3</a>
%     16.2   A toolbox template
%            <a href="matlab: coco_recipes_doc alg_v7_demo alg_v7_demo">Figures 16.4-16.5</a>
%     16.3   An alternative constructor
%
%   Chapter 17 Event Handlers and Branch Switching
%     17.1   Toolbox event handlers
%            <a href="matlab: coco_recipes_doc alg_v9_demo alg_v9_demo">Figure 17.1</a>
%            <a href="matlab: coco_recipes_doc alg_v10_demo alg_v10_demo">Figure 17.3</a>
%     17.2   Bifurcations of periodic orbits
%            <a href="matlab: coco_recipes_doc hspo_v2_demo hspo_v2_demo">Figures 17.4-17.5</a>
%     17.3   Branch switching
%            <a href="matlab: coco_recipes_doc atlas1d_v7_demo atlas1d_v7_demo">Figures 17.6-17.11</a>
%
% Part V   Adaptation
%
%   Chapter 18 Pointwise Adaptation and Comoving Meshes
%     18.1   A brute-force approach
%            <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">Figure 18.1</a>
%            <a href="matlab: coco_recipes_doc coll_v3_demo coll_v3_demo">Figure 18.3</a>
%     18.2   (Co)moving meshes
%     18.3   A comoving mesh algorithm
%            <a href="matlab: coco_recipes_doc coll_v4_demo coll_v4_demo">Figure 18.5</a>
%
%   Chapter 19 A Spectral Toolbox
%     19.1   The spectral continuation
%     19.2   Encoding
%     19.3   Adaptivity
%            <a href="matlab: coco_recipes_doc dft_v1_demo dft_v1_demo">Figure 19.2</a>
%
%   Chapter 20 Integrating Adaptation in Atlas Algorithms
%     20.1   Adaptive mesh refinements
%            <a href="matlab: coco_recipes_doc atlas1d_v8_demo atlas1d_v8_demo">Figures 20.2-20.4</a>
%     20.2   Boundary-value problems
%            <a href="matlab: coco_recipes_doc coll_v5_demo coll_v5_demo">Figure 20.5</a>
%            <a href="matlab: coco_recipes_doc coll_v6_demo coll_v6_demo">Figure 20.6</a>
%     20.3   Numerical comparisons
%            <a href="matlab: coco_recipes_doc po_v3_demo po_v3_demo">Figures 20.7-20.13</a>
%            <a href="matlab: coco_recipes_doc canard_demo canard_demo">Figures 20.14-20.19</a>
%
% Part VI  Epilogue
%
%   Chapter 21 Toolbox projects
%     21.1   Calculus of variations
%     21.2   Nonlinear boundary conditions
%     21.3   Connecting orbits

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_recipes_dbf.m 2839 2015-03-05 17:09:01Z fschild $
