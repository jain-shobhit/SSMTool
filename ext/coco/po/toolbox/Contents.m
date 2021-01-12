% PO   Toolbox for periodic orbits in smooth and hybrid dynamical systems.
% Part of the ODE toolbox family.
%
% The toolbox PO implements continuation of single-segment periodic orbits
% in smooth dynamical systems for evolution equations of the form
% x'=dx/dt=f(t,x,p), where f is some non-linear function, or of
% multi-segment periodic orbits in autonomous hybrid dynamical systems,
% defined in terms of families of vector fields, event functions, and jump
% functions.
%
% Bifurcation detection:
%   * Saddle-node bifurcations
%   * Period-doubling bifurcations
%   * Neimark-Sacker (torus) bifurcations
%   * Neutral saddle points [disabled by default]
%
%   * Branch- and geometric fold points [inherited from atlas class]
%
% Continuation of bifurcation points:
%   * Saddle-node continuation
%   * Period-doubling continuation
%   * Neimark-Sacker continuation
%
% Usage demonstrated in:
%   examples/linode    - starting and restarting continuation of periodic
%                        orbits in a non-autonomous smooth dynamical system
%   examples/bistable  - starting and saddle-node continuation of periodic
%                        orbits in a non-autonomous smooth dynamical system
%   examples/hopf      - starting continuation of periodic orbits in
%                        autonomous smooth dynamical system from Hopf
%                        bifurcation and saddle-node continuation
%   examples/tor       - starting continuation of periodic orbits in
%                        autonomous smooth dynamical system from Hopf
%                        bifurcation, saddle-node, period-doubling, and
%                        torus continuation, as well as branch switching at
%                        branch point and period-doubling point
%   examples/marsden   - starting continuation of periodic orbits in
%                        autonomous smooth dynamical system from Hopf
%                        bifurcation, restarting with fixed period 
%   examples/canard    - starting continuation of canard family of periodic
%                        orbits in autonomous smooth dynamical system from
%                        Hopf bifurcation
%   examples/piecewise - starting and restarting continuation of
%                        multi-segment periodic orbits in piecewise-smooth
%                        dynamical system
%   examples/impact    - starting continuation of multi-segment periodic
%                        orbits in hybrid dynamical system with monitoring
%                        of grazing condition, restarting continuation of
%                        grazing bifuration curve
%   examples/bangbang  - starting continuation of multi-segment periodic
%                        orbits in hybrid dynamical system, saddle-node and
%                        period-doubling continuation, branch switching at
%                        branch points and period-doubling bifurcation
%
% Toolbox user interface:
%   ode_isol2po   - Start continuation of single-segment periodic orbits in
%                   smooth dynamical systems from initial guess.
%   ode_isol2hspo - Start continuation of multi-segment periodic orbits in
%                   hybrid dynamical systems from initial guess.
%   ode_po2po     - Start continuation of single-segment periodic orbits in
%                   smooth dynamical systems from saved solution point.
%                   Also allows branch-switching at branch points.
%   ode_hspo2hspo - Start continuation of multi-segment periodic orbits in
%                   hybrid dynamical systems from saved solution point.
%                   Also allows branch-switching at branch points.
%   ode_BP2po and ode_BP2hspo - Switch branches at branch point.
%   ode_HB2po     - Start continuation of single-segment periodic orbits in
%                   smooth dynamical systems from saved Hopf bifurcation of
%                   equilibria.
%   ode_PD2po and ode_PD2hspo - Switch brances at period-doubling point.
%   ode_po2SN and ode_hspo2SN - Start continuation of saddle-node
%                   bifurcation points.
%   ode_po2PD and ode_hspo2PD - Start continuation of period-doubling
%                   bifurcation points.
%   ode_po2TR and ode_hspo2TR - Start continuation of torus bifurcation
%                   points.
%   ode_PD2PD   - Start continuation of period-doubling bifurcation points.
%   ode_TR2TR   - Start continuation of torus bifurcation points.
%
%   po_read_solution   - Read 'po'solution and toolbox data from disk.
%   hspo_read_solution - Read 'hspo' solution and toolbox data from disk.
%   po_settings        - Show and explain settings of 'po' instance.
%   hspo_settings      - Show and explain settings of 'hspo' instance.
%   po_mult_add        - Save Floquet multipliers to bifurcation data.
%
% PO settings:
%   bifus : on|off|{true}|false
%      Enable/disable detection of bifurcations. By default, the PO toolbox
%      will attempt to locate all bifurcation points for which a test
%      function is implemented. Set to off/false to disable detection of
%      all bifurcation points. It is possible to disable the detection of
%      particular bifurcation points; see PO_SETTINGS and HSPO_SETTINGS for
%      details.
%
%      Note that this setting for the PO toolbox will not affect
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
%   See PO_SETTINGS and HSPO_SETTINGS for a full list of options of the
%   toolbox PO and their respective default values.
%
% Toolbox developer interface:
%   po_add and hspo_add       - Add periodic orbit zero problem.
%   po_add_SN and hspo_add_SN - Add saddle-node bifurcation zero problem.
%   po_add_PD and hspo_add_PD - Add period-doubling bifurcation zero
%                               problem.
%   po_add_TR and hspo_add_TR - Add torus bifurcation zero problem.
%   po_sol_info    - Construct 'po' instance solution information
%                    structure.
%   hspo_sol_info  - Construct 'hspo' instance solution information
%                    structure.
%   po_init_data   - Initialize 'po' instance data structure.
%   hspo_init_data - Initialize 'hspo' instance data structure.
%
%   ode_add_tb_info - Add ODE toolbox family information.
%   ode_init_data   - Initialize ODE toolbox family data structure.
%   ode_settings    - Show and explain settings of toolbox family ODE.
%
% See also: ODESET

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2894 2015-10-02 15:05:02Z hdankowicz $
