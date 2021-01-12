% This demo starts near a Hopf bifurcation and tries to compute the
% complete canard family. The demo is executed five times, once with the
% adaptive spectral discretization, once without adaptation, once for each
% moving-mesh algorithm, and once for the comoving mesh algorithm.

% Use
%
%   plot(coco_bd_col(bd,'PT'), coco_bd_col(bd,'po.seg.coll.err'))
%
% after each run to see the discretization error.

a   = -1;
eps = 0.01;
h   = 0.001; % Initial radius

x0  = [ a^3/3-a ; a ];

om  = sqrt(eps);

xi  = [ 1 ; -om/eps*1i ];
xi  = xi/norm(xi);

p0  = @(t) x0*ones(size(t)) + h*real(xi)*sin(om*t) + h*imag(xi)*cos(om*t);

t0  = linspace(0,2*pi/om, 100);
p1  = p0(t0);

%% Run with spectral method

coco_use_recipes_toolbox dft_v1 % Add the dft_v1 toolbox to the search path

prob = coco_prob();
prob = coco_set(prob, 'dft', 'TOL', 1e-2);
prob = coco_set(prob, 'dft', 'NMAX', 100, 'NMIN', 5, 'NMOD', 5);
prob = dft_isol2orb(prob, '', @vanderpol, t0', p1', {'a' 'eps'}, [a; eps]);
prob = coco_add_event(prob, 'UZ', 'dft.period', [100:100:400 460]);
prob = coco_set(prob, 'lsol', 'cond', true); % Add condition number to chart data
prob = coco_add_func_after(prob, 'mfunc', @lsol_cond); % Monitor condition number
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'h_max', 2, 'bi_direct', false);
coco(prob, '1', [], 1, ...
  {'dft.period' 'a' 'dft.NMOD' 'dft.err' 'lsol.cond'}, {[0 1000] [-1 -0.9]});


%% Run with uniform mesh

coco_use_recipes_toolbox po_v1 coll_v3 % Add the po_v1 and coll_v3 toolboxes to the search path

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 100, 'TOL', 1.0e-2);
prob = po_isol2orb(prob, '', @vanderpol, t0', p1', {'a' 'eps'}, [a; eps]);
prob = coco_add_event(prob, 'UZ', 'po.period', [100:100:400 460]);
prob = coco_set(prob, 'lsol', 'cond', true); % Add condition number to chart data
prob = coco_add_func_after(prob, 'mfunc', @lsol_cond); % Monitor condition number
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'h_max', 2, 'bi_direct', false);
coco(prob, '2', [], 1, ...
  {'po.period' 'a' 'po.seg.coll.NTST' 'po.seg.coll.err' 'lsol.cond'}, ...
  {[0 1000] [-1 -0.9]});


%% Run with fixed-order moving mesh algorithm

coco_use_recipes_toolbox po_v3 coll_v5 % Add the po_v3 and coll_v5 toolboxes to the search path

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 70);
prob = po_isol2orb(prob, '', @vanderpol, t0', p1', {'a' 'eps'}, [a;eps]);
prob = coco_add_event(prob, 'UZ', 'po.period', [100:100:400 460]);
prob = coco_set(prob, 'lsol', 'cond', true); % Add condition number to chart data
prob = coco_add_func_after(prob, 'mfunc', @lsol_cond); % Monitor condition number
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 2, 'bi_direct', false);
coco(prob, '3', [], 1, ...
  {'po.period' 'a' 'po.seg.coll.NTST' 'po.seg.coll.err' 'lsol.cond'}, ...
  {[0 1000] [-1 -0.9]});


%% Run with variable-order moving mesh algorithm

coco_use_recipes_toolbox po_v3 coll_v6 % Add the po_v3 and coll_v6 toolboxes to the search path

prob = coco_prob();
prob = po_isol2orb(prob, '', @vanderpol, t0', p1', {'a' 'eps'}, [a;eps]);
prob = coco_add_event(prob, 'UZ', 'po.period', [100:100:400 460]);
prob = coco_set(prob, 'lsol', 'cond', true); % Add condition number to chart data
prob = coco_add_func_after(prob, 'mfunc', @lsol_cond); % Monitor condition number
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 2, 'bi_direct', false);
coco(prob, '4', [], 1, ...
  {'po.period' 'a' 'po.seg.coll.NTST' 'po.seg.coll.err' 'lsol.cond'}, ...
  {[0 1000] [-1 -0.9]});


%% Run with comoving-mesh algorithm

coco_use_recipes_toolbox po_v3 coll_v7 % Add the po_v3 and coll_v7 toolboxes to the search path

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 100, 'TOL', 1.0e-2, 'hfac', 5);
prob = po_isol2orb(prob, '', @vanderpol, t0', p1', {'a' 'eps'}, [a;eps]);
prob = coco_add_event(prob, 'UZ', 'po.period', [100:100:400 460]);
prob = coco_set(prob, 'lsol', 'cond', true); % Add condition number to chart data
prob = coco_add_func_after(prob, 'mfunc', @lsol_cond); % Monitor condition number
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'h_max', 2, 'bi_direct', false);
coco(prob, '5', [], 1, ...
  {'po.period' 'a' 'po.seg.coll.NTST' 'po.seg.coll.err' 'lsol.cond'}, ...
  {[0 1000] [-1 -0.9]});

coco_use_recipes_toolbox % Remove the po_v1, po_v3, dft_v1, coll_v3, coll_v5, coll_v6, and coll_v7 toolboxes from the search path
