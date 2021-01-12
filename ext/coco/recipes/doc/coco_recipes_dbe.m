% Recipes for Continuation
% Harry Dankowicz and Frank Schilder
%
% List of demos ordered alphabetically by problem name.
%
% BANGBANG, Eqs. (17.10-17.13), in <a href="matlab: coco_recipes_doc hspo_v2_demo hspo_v2_demo">hspo_v2_demo</a> and <a href="matlab: coco_recipes_doc atlas1d_v7_demo atlas1d_v7_demo">atlas1d_v7_demo</a>:
%   Multisegment boundary-value problem whose solutions correspond to
%   periodic orbits of the Duffing single-degree-of-freedom nonlinear
%   oscillator under bang-bang excitation.
%
% BRATU, Eqs. (8.11-8.12), in <a href="matlab: coco_recipes_doc bvp_v1_demo bvp_v1_demo">bvp_v1_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to the
%   steady-state temperature distribution of an exothermic reaction in a
%   1-dimensional medium with boundaries connected to heat baths.
%
% BRUSSELATOR, Eqs. (3.90-3.91), in <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">cmds_demo</a>:
%   Single-segment problem whose solutions correspond to equilibria of a
%   1-dimensional Brusselator model.
%
% CATENARY, Eqs. (1.12-1.13), in <a href="matlab: coco_recipes_doc calcvar_demo calcvar_demo">calcvar_demo</a>:
%   Algebraic problem whose solutions parameterize extremal curves of the
%   catenary problem from the calculus of variation.
%
% CATENARY, Eqs. (1.41-1.42), in <a href="matlab: coco_recipes_doc calcvar_demo calcvar_demo">calcvar_demo</a>:
%   Ooptimization problem whose solutions correspond to extremal curves of
%   an approximation by quadrature of the catenary problem from the
%   calculus of variation.
%
% CATENARY, Eqs. (7.22-7.23), in <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a>:
%   Single-segment problem whose solutions satisfy the Euler-Lagrange
%   equations for the catenary problem from the calculus of variation.
%
% CATENARY, Eqs. (1.34) and (8.6-8.7), in <a href="matlab: coco_recipes_doc calcvar_demo calcvar_demo">calcvar_demo</a> and <a href="matlab: coco_recipes_doc bvp_v1_demo bvp_v1_demo">bvp_v1_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to extremal
%   curves of the catenary problem from the calculus of variation.
%
% CUSP, Eq. (15.19), in <a href="matlab: coco_recipes_doc alg_v7_demo alg_v7_demo">alg_v7_demo</a> and <a href="matlab: coco_recipes_doc alg_v8_demo alg_v8_demo">alg_v8_demo</a>:
%   Algebraic problem whose solutions correspond to equilibria of the cusp
%   normal form.
%
% DOEDEL, Eq. (7.26), in <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a>:
%   Constrained two-segment problem whose solutions approximate a
%   heteroclinic connection between two known equilibria. Source: Doedel,
%   E.J. and Friedman, M.J., "Numerical computation of heteroclinic
%   orbits," Journal of Computational and Applied Mathematics, 26, pp.
%   155-170, 1989.
%
% DUFFING, Eq. (16.10), in <a href="matlab: coco_recipes_doc atlas1d_v6_demo atlas1d_v6_demo">atlas1d_v6_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to periodic
%   orbits of the harmonically excited, Duffing single-degree-of-freedom
%   mechanical oscillator.
%
% HENON, Eq. (3.74), in <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">cmds_demo</a>:
%   Algebraic problem whose solutions correspond to periodic orbits of the
%   two-dimensional nonlinear Henon map.
%
% HUXLEY, Eq. (7.39), in <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a>:
%   Constrained two-segment and equilibrium problem whose solutions
%   approximate heteroclinic connections between two unknown equilibria.
%   Source: Doedel, E.J. and Friedman, M.J., "Numerical computation of
%   heteroclinic orbits," Journal of Computational and Applied Mathematics,
%   26, pp. 155-170, 1989.
%
% IMPACT, Eq. (15.68-15.72), in <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">hspo_v1_demo</a>:
%   Multisegment boundary-value problem whose solutions correspond to
%   periodic orbits of a harmonically excited, single-degree-of-freedom,
%   linear mechanical oscillator with impacts.
%
% LANGFORD, Eq. (9.4), in <a href="matlab: coco_recipes_doc msbvp_v1_demo msbvp_v1_demo">msbvp_v1_demo</a>:
%   Multisegment boundary-value problem whose solutions correspond to
%   collections of trajectory segments on a quasiperiodic invariant torus.
%   Source: Langford, W.F., "Numerical studies of torus bifurcations," in
%   Numerical Methods for Bifurcation Problems, Kupper, T., Mittlemann,
%   H.D., and Weber, H. (Eds.), Birkhaeuser Verlag, Basel, pp. 285-295,
%   1984.
%
% LANGFORD, Eq. (14.4), in <a href="matlab: coco_recipes_doc atlas2d_v6_demo atlas2d_v6_demo">atlas2d_v6_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to
%   phase-locked, resonant periodic orbits. Source: Langford, W.F.,
%   "Numerical studies of torus bifurcations," in Numerical Methods for
%   Bifurcation Problems, Kupper, T., Mittlemann, H.D., and Weber, H.
%   (Eds.), Birkhaeuser Verlag, Basel, pp. 285-295, 1984.
%
% LIENARD, Eq. (8.27), in <a href="matlab: coco_recipes_doc bvp_v2_demo bvp_v2_demo">bvp_v2_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to periodic
%   orbits of a planar Lienard dynamical system.
%
% LINODE, Eq. (10.35), in <a href="matlab: coco_recipes_doc varcoll_v1_demo varcoll_v1_demo">varcoll_v1_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to periodic
%   orbits of a harmonically excited, single-degree-of-freedom, linear
%   mechanical oscillator.
%
% LINODE, not in Recipes for Continuation, in <a href="matlab: coco_recipes_doc coll_v2_demo coll_v2_demo">coll_v2_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to
%   trajectory segments of a harmonically excited, single-degree-of-freedom
%   mechanical oscillator whose first local maximum in the displacement
%   occurs at a displacement of 1.
%
% LORENZ, Eq. (10.87), in <a href="matlab: coco_recipes_doc varcoll_v2_demo varcoll_v2_demo">varcoll_v2_demo</a>:
%   Constrained two-segment and two-point boundary-value problem whose
%   solutions correspond to heteroclinic connections between a known
%   equilibrium and an unknown periodic orbit in the Lorenz system.
%
% MARSDEN, Eq. (8.24) and (20.44), in <a href="matlab: coco_recipes_doc po_v1_demo po_v1_demo">po_v1_demo</a> and <a href="matlab: coco_recipes_doc po_v3_demo po_v3_demo">po_v3_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to periodic
%   orbits emanating from a Hopf bifurcation and limiting on a homoclinic
%   trajectory. Source: Marsden, J.E. and McCracken, M., "The Hopf
%   bifurcation and Its Applications," Springer-Verlag, New York, 1976.
%
% PNETA, Eq. (18.5), in <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a>, <a href="matlab: coco_recipes_doc coll_v3_demo coll_v3_demo">coll_v3_demo</a>, <a href="matlab: coco_recipes_doc coll_v4_demo coll_v4_demo">coll_v4_demo</a>,
%   <a href="matlab: coco_recipes_doc coll_v5_demo coll_v5_demo">coll_v5_demo</a>, <a href="matlab: coco_recipes_doc coll_v6_demo coll_v6_demo">coll_v6_demo</a>, and <a href="matlab: coco_recipes_doc dft_v1_demo dft_v1_demo">dft_v1_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to periodic
%   orbits of a slow-fast oscillator.
%
% POPUL, Eq. (17.3), in <a href="matlab: coco_recipes_doc alg_v9_demo alg_v9_demo">alg_v9_demo</a> and <a href="matlab: coco_recipes_doc alg_v10_demo alg_v10_demo">alg_v10_demo</a>:
%   Algebraic problem whose solutions correspond to equilibria in a model
%   of population dynamics. Source: Kuznetsov, Yu.A., "Elements of Applied
%   Bifurcation Theory," Springer-Verlag, New York, 1998.
%
% PWLIN, Eqs. (9.35-9.38), in <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">hspo_v1_demo</a> and <a href="matlab: coco_recipes_doc varcoll_v1_demo varcoll_v1_demo">varcoll_v1_demo</a>:
%   Two-segment boundary-value problem whose solutions correspond to
%   periodic orbits of a piecewise smooth dynamical system.
%
% STICKSLIP, Eqs. (9.52-9.63), in <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">hspo_v1_demo</a>:
%   Two- and three-segment boundary-value problem whose solutions
%   correspond to periodic orbits of a hybrid dynamical system modeling a
%   harmonically excited, two-degree-of-freedom, nonlinear mechanical
%   oscillator with impacts and friction. Source: Zhao, X., Reddy, C.K.,
%   and Nayfeh, A., "Nonlinear dynamics of an electrically driven impact
%   microactuator," Nonlinear Dynamics, 40(3), pp. 227-239, 2005.
%
% TANH, Eq. (20.3), in <a href="matlab: coco_recipes_doc atlas1d_v8_demo atlas1d_v8_demo">atlas1d_v8_demo</a>:
%   Interpolation problem whose solutions correspond to an equipartitioned
%   sample of the function tanh(p*t)/tahn(p) on a given interval.
%
% TOR, not in Recipes for Continuation, in <a href="matlab: coco_recipes_doc po_v2_demo po_v2_demo">po_v2_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to periodic
%   orbits of a three-dimensional nonlinear dynamical system.
%
% VANDERPOL, Eq. (20.45), in <a href="matlab: coco_recipes_doc canard_demo canard_demo">canard_demo</a>:
%   Two-point boundary-value problem whose solutions correspond to members
%   of the canard family of periodic orbits of the forced Van-der-Pol
%   equation.


%   <a href="matlab: coco_recipes_edit coll_v2_demo
%   demo_coll_v2">demo_coll_v2.m</a> - The grazing manifold. Initial
%   conditions for a harmonically excited linear oscillator that result in
%   a given locally maximal displacement.
%                        <a href="matlab: coco_recipes_doc manifolds2d manifolds2d">manifolds2d</a>
%            <a href="matlab: coco_recipes_doc torex torex">torex</a>


% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_recipes_dbe.m 2839 2015-03-05 17:09:01Z fschild $
