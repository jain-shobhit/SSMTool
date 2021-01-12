coco_use_recipes_toolbox dft_v1 % add the dft_v1 toolbox to the search path

%% Demo with autonomous chemical oscillator [Not in Recipes for Continuation]

% The continuation problem structure encoded below initially consists of
% 1605 zero functions (404 real and imaginary parts of the spectral zero
% problem, the 1200 dummy conditions, and a phase condition) in terms of
% 1614 continuation variables (404 active Fourier coefficients, 1200
% inactive Fourier coefficients, the period, and the problem parameters), a
% family of monitor functions that evaluate to the problem parameters and
% the period, and the corresponding inactive continuation parameter 'a'
% through 'i' and active continuation parameter 'dft.period'. Its
% dimensional deficit equals 0. A one-dimensional family of periodic orbits
% is obtained by releasing any one of the inactive continuation parameters
% and allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-2) are identified by 'MCXL'. Updates to the
% discretization order are reflected in the nonembedded continuation
% parameter 'dft.NMOD'.

p0 = [0.1631021 1250 0.046875 20 1.104 0.001 3 0.6 0.1175]';
x0 = [25 1.45468 0.01524586 0.1776113]';
f  = @(t,x) chemosc(x,p0);

[~, z]  = ode15s(f, [0 800], x0);
x0      = z(end,:)';
[t0 z0] = ode15s(f, [0 13.4], x0);

prob = coco_prob();
prob = coco_set(prob, 'dft', 'TOL', 1e-2, 'NMAX', 200, 'NMOD', 50);
prob = dft_isol2orb(prob,'',@chemosc, t0, z0, ...
  {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i'} , p0);
prob = coco_set(prob, 'cont', 'ItMX', 200, 'FP', true);
coco(prob, 'run1', [], 1, {'a' 'dft.err' 'dft.NMOD'}, [0.1 0.5]);

bd   = coco_bd_read('run1');
labs = coco_bd_labs(bd, 'all');

figure(1)
clf
hold on
for lab=labs
  sol = dft_read_solution('', 'run1', lab);
  plot3(sol.x(:,1), sol.x(:,2), sol.x(:,3), 'r.-')
end
hold off

coco_use_recipes_toolbox % remove the dft_v1 toolbox from the search path
