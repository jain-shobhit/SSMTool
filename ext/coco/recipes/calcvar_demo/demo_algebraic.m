% The continuation problem encoded below corresponds to a family of two
% zero functions in terms of three continuation variables, three monitor
% functions that evaluate to each of the continuation variables, and the
% corresponding three inactive continuation parameters. Its dimensional
% deficit equals -2. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameters
% 'Y', 'a', and 'b' are released and allowed to vary during continuation.
% The family of solution points corresponds to a parameterization of
% solutions to the Euler- Lagrange equations of the catenary problem of the
% form f:x->cosh(a*(x+b))/a, for which f(0)=1 and f(1)=Y.

prob = coco_prob();
prob = coco_add_func(prob, 'fun', @caty_alg, [], 'zero', ...
  'u0', [1; 0; cosh(1)]);
prob = coco_add_pars(prob, '', 1:3, {'a', 'b', 'Y'});
prob = coco_set(prob, 'cont', 'PtMX', 100);
coco(prob, 'run_alg', [], 1, {'Y', 'a', 'b'}, [0 10]);
