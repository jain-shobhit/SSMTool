coco_use_recipes_toolbox dft_v1 % add the dft_v1 toolbox to the search path

%% Example 19.4

% The continuation problem structure encoded below consists of 1003 zero
% functions (82 real and imaginary parts of the spectral zero problem, the
% 920 dummy conditions, and a phase condition) in terms of 1004
% continuation variables (82 active Fourier coefficients, 920 inactive
% Fourier coefficients, the period, and the problem parameter), a pair of
% monitor functions that evaluate to the problem parameter and the period,
% and the corresponding inactive continuation parameter 'eps' and active
% continuation parameter 'dft.period'. Its dimensional deficit equals 0. A
% one-dimensional family of periodic orbits is obtained by releasing 'eps'
% and allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-3) are identified by 'MCXL'. Updates to the
% discretization order are reflected in the nonembedded continuation
% parameter 'dft.NMOD'.

eps0 = 0.1;
t0 = linspace(0, 2*pi, 100)';
x0 = [sin(t0) cos(t0)];
prob = coco_set(coco_prob(), 'cont', 'ItMX', 100);
prob = coco_set(prob, 'dft', 'TOL', 1e-3);
prob2 = coco_set(prob, 'dft', 'NMAX', 250, 'NMOD', 20);
prob2 = dft_isol2orb(prob2, '', @pneta, t0, x0, 'eps', eps0);
coco(prob2, 'run1', [], 1, {'eps' 'dft.err' 'dft.NMOD'}, [0.1 20]);

% The continuation problem structure encoded below initially consists of
% 1003 zero functions (1002 real and imaginary parts of the spectral zero
% problem and a phase condition) in terms of 1004 continuation variables
% (1002 active Fourier coefficients, no inactive Fourier coefficients, the
% period, and the problem parameters), a pair of monitor functions that
% evaluate to the problem parameter and the period, and the corresponding
% inactive continuation parameter 'eps' and active continuation parameter
% 'dft.period'. Its dimensional deficit again equals 0. A one-dimensional
% family of periodic orbits is obtained by releasing 'eps' and allowing it
% to vary during continuation. Terminal special points associated with a
% discretization error estimate exceeding a given tolerance (here 1e-3) are
% identified by 'MCXL'. Updates to the discretization order are reflected
% in the nonembedded continuation parameter 'dft.NMOD'.

bd  = coco_bd_read('run1');
lab = coco_bd_labs(bd, 'EP');
prob2 = dft_sol2orb(prob, '', 'run1', lab(end));
coco(prob2, 'run2', [], 1, {'eps' 'dft.err' 'dft.NMOD'}, [0.1 20]);

coco_use_recipes_toolbox % remove the dft_v1 toolbox from the search path
