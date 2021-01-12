% EP   Toolbox for equilibrium continuation.
% Part of the ODE toolbox family.
%
% The toolbox EP implements continuation of dynamic equilibrium points of
% autonomous evolution equations of the form x' = dx/dt = f(x,p), or of
% static equilibrium points of an equation of the form f(x,p) = 0, where f
% is some non-linear function.
%
% Bifurcation detection:
%   * Saddle-node bifurcations
%   * Hopf bifurcations
%   * Neutral saddle points [disabled by default]
%   * Bogdanov-Takens points
%
%   * Branch- and geometric fold points [inherited from atlas class]
%
% Continuation of bifurcation points:
%   * Saddle-node continuation
%   * Hopf continuation
%
% Usage demonstrated in:
%   examples/cusp    - starting and saddle-node continuation
%   examples/bratu   - starting, restarting, and saddle-node continuation
%   examples/brus    - starting, branch-switching, Hopf, and saddle-node
%                     continuation
%   examples/pdeeig  - starting, branch-switching, adding 'ep' instance
%                     outside call to coco, and adding a user-defined
%                     monitor function and an event associated with this
%                     monitor function
%   examples/chemosc - starting, Hopf and saddle-node continuation, adding
%                     'ep' instance outside call to coco, and adding a
%                     user-defined monitor function and an event associated
%                     with this monitor function
%   examples/isola   - starting, restarting, composite continuation problem
%                     with multiple coupled 'ep' instances
%
% Toolbox user interface:
%   ode_isol2ep - Start continuation of equilibrium points from initial
%                 guess (with optional inclusion of variational problem).
%   ode_ep2ep   - Start continuation of equilibrium points from saved
%                 solution point (with optional inclusion of variational
%                 problem). Also allows branch-switching at branch points.
%   ode_BP2ep   - Switch branches at branch point (with optional inclusion
%                 of variational problem).
%   ode_SN2SN   - Start continuation of saddle-node bifurcation points.
%   ode_ep2SN   - Start continuation of saddle-node bifurcation points.
%   ode_HB2HB   - Start continuation of Hopf bifurcation points.
%   ode_ep2HB   - Start continuation of Hopf bifurcation points.
%
%   ep_read_solution - Read solution and toolbox data from disk.
%   ep_settings      - Show and explain settings of 'ep' instance.
%
% EP settings:
%   bifus : on|off|{true}|false
%      Enable/disable detection of bifurcations. By default, the EP toolbox
%      will attempt to locate all bifurcation points for which a test
%      function is implemented. Set to off/false to disable detection of
%      all bifurcation points. It is possible to disable the detection of
%      particular bifurcation points; see EP_SETTINGS for details.
%
%      Note that this setting for the EP toolbox will not affect
%      bifurcation detection implemented by other toolboxes, for example,
%      the atlas classes.
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
%   See EP_SETTINGS for a full list of options of the toolbox EP and their
%   respective default values.
%
% Toolbox developer interface:
%   ep_add      - Add equilibrium zero problem.
%   ep_add_SN   - Add saddle-node bifurcation zero problem.
%   ep_add_HB   - Add Hopf bifurcation zero problem.
%   ep_add_var  - Add variational problem.
%   ep_sol_info - Construct EP solution information structure.
%
%   ode_add_tb_info - Add ODE toolbox family information.
%   ode_init_data   - Initialize ODE toolbox family data structure.
%   ode_settings    - Show and explain settings of toolbox family ODE.
%
% See also: ODESET

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2930 2015-10-30 17:22:54Z hdankowicz $
