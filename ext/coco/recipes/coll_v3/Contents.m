% Toolbox: 'coll'
% Version:  3.0
% 
% Source:   Sect. 18.1.2 of Recipes for Continuation
%
% This version of the 'coll' toolbox provides for estimation of the
% approximation error and includes the assignment of a terminal event when
% this estimate exceeds a given error tolerance.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_doc coll_v3 coll_isol2seg">coll_isol2seg</a>      - Append 'coll' instance constructed from initial
%                        data.
%   <a href="matlab:coco_recipes_doc coll_v3 coll_sol2seg">coll_sol2seg</a>       - Append 'coll' instance constructed from saved data.
%   <a href="matlab:coco_recipes_doc coll_v3 coll_construct_seg">coll_construct_seg</a> - Append an instance of 'coll' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit coll_v3 coll_construct_seg>coll_F">coll_construct_seg>coll_F</a>    - Collocation zero problem.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_construct_seg>coll_DFDU">coll_construct_seg>coll_DFDU</a> - Linearization of collocation zero
%                                  problem.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_construct_seg>coll_err">coll_construct_seg>coll_err</a>  - Evaluate estimate of approximation eror.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit coll_v3 coll_arg_check">coll_arg_check</a>     - Basic argument checking for 'coll' toolbox.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_get_settings">coll_get_settings</a>  - Read 'coll' toolbox instance settings.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_init_data">coll_init_data</a>     - Initialize toolbox data for an instance of 'coll'.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_init_data>coll_nodes">coll_init_data>coll_nodes</a> - Compute collocation nodes and integration
%                               weights.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_init_data>coll_L">coll_init_data>coll_L</a>     - Vectorized evaluation of Lagrange
%                               polynomials.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_init_data>coll_Lp">coll_init_data>coll_Lp</a>    - Vectorized evaluation of derivative of
%                               Lagrange polynomials.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_init_data>coll_Lm">coll_init_data>coll_Lm</a>     - Compute highest order coefficient of Lagrange polynomial.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_init_sol">coll_init_sol</a>      - Build initial solution guess.
%   <a href="matlab:coco_recipes_edit coll_v3 coll_read_solution">coll_read_solution</a> - Read 'coll' solution and toolbox data from disk.
%
% Toolbox data (struct)
%   fhan    : Function handle to vectorized encoding of vector field.
%   dfdxhan : Optional function handle to vectorized encoding of Jacobian
%             w.r.t. problem variables.
%   dfdphan : Optional function handle to vectorized encoding of Jacobian
%             w.r.t. problem parameters.
%   pnames  : Empty cell array or cell array of string labels for
%             continuation parameters tracking problem parameters.
%   dim     : State space dimension.
%   pdim    : Number of problem parameters.
%   xbp_idx : Index set for basepoint values.
%   T_idx   : Index for interval length.
%   p_idx   : Index set for problem parameters.
%   tbp_idx : Index set for basepoint values without duplication.
%   x_shp   : Shape for collocation nodes.
%   xbp_shp : Shape for basepoints.
%   p_rep   : Shape for problem parameters.
%   tbp     : Basepoint mesh.
%   x0_idx  : Index set for trajectory end point at t=0.
%   x1_idx  : Index set for trajectory end point at t=1.
%   wts1    : Array of quadrature weights.
%   wts2    : Array of quadrature weights.
%   W       : Interpolation matrix for Lagrange polynomials.
%   Wp      : Interpolation matrix for derivatives of Lagrange polynomials.
%   Wm      : Interpolation matrix for coefficient of highest polynomial
%             power.
%   wn      : Bound on mesh products.
%   dxrows  : Row index array for construction of sparse Jacobian with
%             respect to problem variables.
%   dxcols  : Column index array for construction of sparse Jacobian with
%             respect to problem variables.
%   dprows  : Row index array for construction of sparse Jacobian with
%             respect to problem parameters.
%   dpcols  : Column index array for construction of sparse Jacobian with
%             respect to problem parameters.
%   Qnum    : Number of rows for Jacobian of continuity conditions with respect to T and p.
%
% Toolbox settings (in coll field of data)
%   NTST    : Number of mesh intervals (default: 10).
%   NCOL    : Degree of polynomial interpolants (default: 4).
%   TOL     : Tolerance for discretization error estimate (default: 2/3
%             power of corrector tolerance.
%
% Continuation parameters
%   err     : Track discretization error estimate (regular).
%   err_TF  : Track error estimate scaled by tolerance (regular).
%
% Event types
%   MXCL    : Critical error estimate.
%
% Solution (struct)
%   t       : Array of basepoint time instances.
%   x       : Array of trajectory basepoint values.
%   p       : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit coll_v3 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc coll_v3_demo coll_v3_demo">coll_v3_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
