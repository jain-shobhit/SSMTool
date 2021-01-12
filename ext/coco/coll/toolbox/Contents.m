% COLL   Toolbox for trajectory segments and collections of constrained trajectory segments.
% Part of the ODE toolbox family.
%
% The toolbox COLL implements continuation of individual trajectory
% segments, or of collections of constrained trajectory segments, for
% evolution equations of the form x'=dx/dt=f(t,x,p), where f is some
% non-linear function.
%
% Bifurcation detection:
%   * Branch- and geometric fold points [inherited from atlas class]
%
% Usage demonstrated in:
%   examples/catenary - starting and restarting continuation of trajectory
%                       segments
%   examples/linode   - starting and restarting continuation (including the
%                       variational problem) of trajectory segments in
%                       autonomous and non-autonomous vector fields
%   examples/huxley   - starting and restarting continuation of pairs of
%                       coupled 'coll' instances, and adding boundary
%                       conditions using explicit calls to coco_add_func
%   examples/doedel   - starting and restarting continuation of multiple
%                       coupled 'coll and 'ep' instances (including the
%                       variational problem), and adding boundary
%                       conditions using explicit calls to coco_add_func
%   examples/bratu    - starting, restarting and branch switching of
%                       constrained trajectory segments 
%   examples/lienard  - starting continuation of constrained trajectory
%                       segments with and without updates to boundary
%                       conditions
%   examples/torus    - starting continuation of collections of constrained
%                       trajectory segments
%
% Toolbox user interface:
%   ode_isol2coll - Start continuation of trajectory segments from initial
%                   guess (with optional inclusion of variational problem).
%   ode_coll2coll - Start continuation of trajectory segments from saved
%                   solution point (with optional inclusion of variational
%                   problem). Also allows branch-switching at branch
%                   points.
%   ode_BP2coll   - Switch branches at branch point (with optional
%                   inclusion of variational problem).
%   ode_isol2bvp  - Start continuation of collections of constrained
%                   trajectory segments from initial guess (with optional
%                   inclusion of variational problem).
%   ode_bvp2bvp   - Start continuation of collections of constrained
%                   trajectory segments from saved solution point (with
%                   optional inclusion of variational problem). Also allows
%                   branch-switching at branch points.
%   ode_BP2bvp    - Switch branches at branch point (with optional
%                   inclusion of variational problem).
%
%   coll_read_solution - Read 'coll'solution and toolbox data from disk.
%   bvp_read_solution  - Read 'bvp' solution and toolbox data from disk.
%   coll_settings      - Show and explain settings of 'coll' instance.
%
% COLL settings:
%   NTST : {10}
%      Integer number of discretization intervals.
%   NCOL : {4}
%      Integer degree of interpolating polynomials.
%   var  : on|off|true|{false}
%      Enable/disable nonembedded computation of fundamental solution to
%      variational problem.
%
% ODE settings:
%   vectorized : on|off|{true}|false
%      The vectorized property enables a reduction in the number of
%      function evaluations required to compute all the columns of the
%      Jacobian matrix of f, and might significantly reduce computation
%      time. A vectorized evaluation of f(x,p) is possible, if f is encoded
%      such that f([x1 x2 ...], [p1 p2 ...]) returns [f(x1,p1) f(x2,p2),
%      ...]. This enables fast evaluation of finite differences.
%
%      Set to off/false to disable vectorized evaluation. See section
%      "Jacobian Matrix Properties" of the documentation of ODESET for more
%      details on vectorized evaluation of functions.
%
%   See COLL_SETTINGS for a full list of options of the toolbox COLL and
%   their respective default values.
%
% Toolbox developer interface:
%   coll_add      - Add 'coll' instance trajectory segment zero problem.
%   coll_add_var  - Add 'coll' instance variational problem.
%   coll_sol_info - Construct 'coll' instance solution information structure.
%   bvp_sol_info  - Construct 'bvp' instance solution information structure.
%   bvp_init_data - Initialize 'bvp' instance data structure.
%
%   ode_add_tb_info - Add ODE toolbox family information.
%   ode_init_data   - Initialize ODE toolbox family data structure.
%   ode_settings    - Show and explain settings of toolbox family ODE.
%
% See also: ODESET

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2894 2015-10-02 15:05:02Z hdankowicz $
