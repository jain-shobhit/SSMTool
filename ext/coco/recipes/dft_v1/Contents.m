% Toolbox: 'dft'
% Version:  1.0
% 
% Source:   Sects. 19.2 and 19.3 of Recipes for Continuation
%
% This version of the 'dft' toolbox implements a spectral collocation
% method for approximating a periodic solution of an autonomous ordinary
% differential equation. This version demonstrates an order adaptation
% scheme that does not change the dimension of the zero problem and
% requires an update slot only.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit dft_v1 dft_isol2orb">dft_isol2orb</a>      - Append 'dft' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_sol2orb">dft_sol2orb</a>       - Append 'dft' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_construct_orb">dft_construct_orb</a> - Append an instance of 'dft' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit dft_v1 dft_construct_orb>dft_F">dft_construct_orb>dft_F</a>    - Implement 'dft' toolbox zero problem
%   <a href="matlab:coco_recipes_edit dft_v1 dft_construct_orb>dft_DFDU">dft_construct_orb>dft_DFDU</a> - Linearization of 'dft' zero problem.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_construct_orb>dft_error">dft_construct_orb>dft_error</a> - Evaluate estimate of approximation error.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_construct_orb>dft_update">dft_construct_orb>dft_update</a> - Update spectral discretization.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit dft_v1 dft_arg_check">dft_arg_check</a>     - Basic argument checking for 'dft' toolbox.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_get_settings">dft_get_settings</a>  - Read 'dft' toolbox instance settings.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_init_data">dft_init_data</a>     - Initialize toolbox data for an instance of 'dft'.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_init_modes">dft_init_modes</a> - Initialize adaptable data for an instance of 'dft'.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_init_sol">dft_init_sol</a>      - Build initial solution guess.
%   <a href="matlab:coco_recipes_edit dft_v1 dft_read_solution">dft_read_solution</a> - Read 'dft' solution and toolbox data from disk.
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
%   T_idx   : Index for period.
%   p_idx   : Index set for problem parameters.
%   x_shp   : Shape for vectorization.
%   p_rep   : Shape for problem parameters.
%   Forig   : Fourier matrix.
%   Fsorig  : Kronecker transpose.
%   Finvs   : Kronecker transpose of inverse Fourier matrix.
%   dxrows  : Row index array for construction of sparse Jacobian with
%             respect to problem variables.
%   dxcols  : Column index array for construction of sparse Jacobian with
%             respect to problem variables.
%   dprows  : Row index array for construction of sparse Jacobian with
%             respect to problem parameters.
%   dpcols  : Column index array for construction of sparse Jacobian with
%             respect to problem parameters.
%   mdim    : Number of active coefficients
%   ddim    : Number of dummy coefficients.
%   xf_idx  : Index array for active modes.
%   xd_idx  : Index array for dummy modes.
%   phs     : Reference coefficients for phase condition.
%   jac     : Jacoboian of dummy and phase conditions.
%   Wp      : Differentiation matrix.
%   F       : Fourier matrix.
%   Fs      : Kronecker transpose.
%   FinvsW  : Interpolation matrix.
%
% Toolbox settings (in coll field of data)
%   NMAX    : Maximum number of Fourier modes (default: 8).
%   NMIN    : Minimum number of Fourier modes (default: 3).
%   NMOD    : Number of Fourier modes (default: 3).
%   TOL     : Tolerance for discretization error estimate (default: 2/3
%             power of corrector tolerance.
%   TOLINC  : Upper bound on adaptation window.
%   TOLDEC  : Lower bound on adaptation window.
%
% Solution (struct)
%   t       : Array of basepoint time instances over full period.
%   x       : Periodic array of trajectory basepoint values.
%   p       : Problem parameters.
%   T       : Period.
%   c       : Active Fourier coefficients.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit dft_v1 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc dft_v1_demo dft_v1_demo">dft_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
