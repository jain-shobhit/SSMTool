coco_use_recipes_toolbox dft_v1 % add the dft_v1 toolbox to the search path

%% Demo with forced, damped, harmonic oscillator [Not in Recipes for Continuation]

% The continuation problem structure encoded below initially consists of
% 805 zero functions (52 real and imaginary parts of the spectral zero
% problem, the 752 dummy conditions, and a phase condition) in terms of
% 806 continuation variables (52 active Fourier coefficients, 752 inactive
% Fourier coefficients, the period, and the problem parameter), a pair of
% monitor functions that evaluate to the problem parameter and the period,
% and the corresponding inactive continuation parameter 'k' and active
% continuation parameter 'dft.period'. Its dimensional deficit equals 0. A
% one-dimensional family of periodic orbits is obtained by releasing 'k'
% and allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance are identified by 'MCXL'. Updates to the discretization order
% are reflected in the nonembedded continuation parameter 'dft.NMOD'.

p0 = 1;
t0 = (0:2*pi/100:2*pi)';
z0 = [sin(t0) cos(t0) sin(t0) cos(t0)];
prob = coco_prob();
prob = coco_set(prob, 'dft', 'NMAX', 100, 'NMOD', 6);
prob = dft_isol2orb(prob, '', @linode, t0, z0, 'k', p0);

bd = coco(prob, 'run', [], 1, {'k' 'dft.err' 'dft.NMOD'}, [0.5 2]);

figure(1)
clf
hold on
labs = coco_bd_labs(bd, 'EP');
for lab=labs
  sol = dft_read_solution('', 'run', lab);
  plot(sol.x(:,1),sol.x(:,2),'r.-')
end
hold off

coco_use_recipes_toolbox % remove the dft_v1 toolbox from the search path
