%% Slow-fast dynamics from Section 18.1.2 in Recipes for Continuation
%
% The continuation problem structure encoded below includes two monitor
% functions that evaluate to the interval length and the problem parameter,
% respectively, with the corresponding active continuation parameter
% 'po.period' and inactive continuation parameter 'eps'. Its dimensional
% deficit equals 0. A family of approximate periodic orbits is obtained by
% releasing 'eps' and allowing it to vary during continuation.

% Comparison of approximate periodic orbits of the dynamical system with
% vector field given by Eq. (18.5) in Recipes for Continuation at
% epsilon=20, obtained using a uniform, non-adapative mesh with N=250 and
% m=5, and using an adaptive mesh with variable discretization order. The
% high resolution for the uniform mesh results in a solution curve that
% seems acceptable, at least by visual inspection. A large subset of the
% additional mesh points contributes little to the improvement, however;
% only a small number is allocated along the problematic vertical part of
% the orbit. This is evident from the graph of the mesh points used by the
% adaptive mesh, for which N=49 and m=4.

%% Initial data

eps0 = 0.1;
t0   = linspace(0, 2*pi, 100)';
x0   = [sin(t0) cos(t0)];

pneta = @(x,p) [x(2,:); 0.5*p(1,:).*x(2,:)-p(1,:).*x(2,:).^3-x(1,:)];

%% Uniform mesh without adaptation

prob = coco_prob();
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = coco_set(prob, 'coll', 'NTST', 250, 'NCOL', 5);
prob = ode_isol2po(prob, '', pneta, t0, x0, 'eps', eps0);
prob = coco_set(prob, 'cont', 'NAdapt', 0, 'PtMX', [0 200]);
coco(prob, 'pneta1', [], 1, {'eps' 'po.orb.coll.err_TF'}, [0.1 20]);

%% Adaptive mesh

prob = coco_prob();
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_isol2po(prob, '', pneta, t0, x0, 'eps', eps0);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [0 200]);
coco(prob, 'pneta2', [], 1, {'eps' 'po.orb.coll.err_TF'}, [0.1 20]);

% Use the following command to extract the discretization order

% bd   = coco_bd_read('pneta2');
% labs = coco_bd_labs(bd, 'EP');
% 
% [~, data] = coll_read_solution('po.orb', 'pneta2', labs(end));
% [data.coll_seg.maps.NTST data.coll_seg.int.NCOL]

%% Graphical representation of stored solutions

bd  = coco_bd_read('pneta1');
lab = coco_bd_labs(bd, 'EP');

% Plot data
thm = struct();
thm.sol.RO = {'LineStyle', ':', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12};
figure(1); clf; hold on; box on; grid on; 
coco_plot_sol(thm, 'pneta1', lab(end), '', 'x', 'x')
axis([-3.0833 3.0833 -1.0794 1.0794]*1.05)


% Extract data
bd  = coco_bd_read('pneta2');
lab = coco_bd_labs(bd, 'EP');

% Plot data
thm.sol.RO = {'LineStyle', ':', 'LineWidth', 1.5, ...
  'Color', [0.8,0.8,0.8], 'Marker', '.', 'MarkerSize', 8};
coco_plot_sol(thm, 'pneta2', lab(end), '', 'x', 'x')

hold off
