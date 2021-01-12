% The continuation problem encoded below corresponds to a family of 50 zero
% functions in terms of 51 continuation variables, a single monitor
% functions that evaluates to the problem parameter, and the corresponding
% inactive continuation parameter 'Y'. Its dimensional deficit equals 0.
% The call to the coco entry-point function indicates a desired manifold
% dimension of 1. To this end, the continuation parameter 'Y' is released
% and allowed to vary during continuation. Special points associated with
% specific values of 'Y' and denoted by 'UZ' are detected and located using
% a Newton method. The family of solutions correspond to a parameterization
% of solutions to a discretized form of the catenary optimization problem,
% in terms of a continuous, piecewise polynomial function f(t) that
% satisfies the boundary conditions f(0)=1 and f(1)=Y.

t0 = 0:0.01:1;
x0 = cosh(t0);
p0 = cosh(1);

calc.fhan   = @caty_calcvar;
calc.pnames = {'Y'};
calc.NTST   = 10;
calc.NCOL   = 4;
[calc sol]  = calcvar_system(calc, t0, x0, p0);

prob = coco_prob();
prob = coco_add_func(prob, 'calcvar', @calcvar_F, @calcvar_DFDU, calc, ...
  'zero', 'u0', sol.u);
prob = coco_add_pars(prob, 'pars', calc.p_idx, 'Y');
prob = coco_add_event(prob, 'UZ', 'Y', [2 3 4]);
prob = coco_add_slot(prob, 'calcvar', @coco_save_data, calc, 'save_full');
coco(prob, 'run_quad', [], 1, 'Y', [0 5]);
