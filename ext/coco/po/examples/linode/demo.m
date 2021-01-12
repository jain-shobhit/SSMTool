%% Continuation of periodic orbits from Example 10.2 in Recipes for Continuation

% Periodic solutions of the differential equations x''+x'+px=cos(t) are
% given by x(t)=(sin(t)+(p-1)*cos(t))/(p^2-2*p+2). A sample solution and
% the approximant obtained by linear interpolation of the base points found
% using continuation is shown in panel (a). Panel (b) shows the computed L2
% norm and the exact value sqrt(2*pi)/sqrt(p^2-2*p+2).

%% Initial encoding

% The continuation problem structure encoded below includes a single
% monitor function that evaluates to the problem parameter, and the
% corresponding inactive continuation parameter 'p'. Its dimensional
% deficit equals 0. A one-dimensional family of periodic orbits is obtained
% by releasing 'p' and allowing it to vary during continuation. Note that
% the vector field is non-autonomous.

pnames = 'p';
p0     = 1;

% Generate an initial solution guess
[t0, x0]   = ode45(@(t,x) linode(t,x,p0), [0 2*pi], [0; 1]);

% Construct periodic orbit zero problem for non-autonomous vector field
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'coll', 'NTST', 15);
coll_func = { @linode, @linode_DFDX, @linode_DFDP, @linode_DFDT };
coll_args = [ coll_func, { t0, x0, pnames, p0 }];
prob = ode_isol2po(prob, '', coll_args{:});

% Adaptive discretization before each continuation step
prob = coco_set(prob, 'cont', 'NAdapt', 1);

fprintf('\n Run=''%s'': Continue family of periodic orbits.\n', ...
  'run');

coco(prob, 'run', [], 1, 'p', [0.2 2]);

%% Restart continuation from stored solution

% Reconstruct periodic orbit zero problem
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = ode_po2po(prob, '', 'run', 3);

% Adaptive discretization before each continuation step
prob = coco_set(prob, 'cont', 'NAdapt', 1);

fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  'run_again', 3, 'run');

coco(prob, 'run_again', [], 1, 'p', [0.05 3]);

%% Graphical representation of stored solutions

% Plot data: panel (a)

figure(1); clf; grid on; box on; hold on; axis([-1 1 -1 1])
coco_plot_sol('run_again', 5, '', 'x', 'x')
thm = struct();
thm.sol.RO = {'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'w'};
thm.xlab = 'x_1';
thm.ylab = 'x_2';
coco_plot_sol(thm, 'run_again', 5, 'po.orb', ...
{'t', 'p'}, @(t,p) (sin(t)+(p-1)*cos(t))/(p^2-2*p+2), ...
{'t', 'p'}, @(t,p) (cos(t)-(p-1)*sin(t))/(p^2-2*p+2));
hold off

% Alternative formulation using the solution extractor
% sol = po_read_solution('run_again',5);
% hor = (sin(sol.tbp)+(sol.p-1)*cos(sol.tbp))/(sol.p^2-2*sol.p+2);
% ver = (cos(sol.tbp)-(sol.p-1)*sin(sol.tbp))/(sol.p^2-2*sol.p+2);
% plot(sol.xbp(:,1), sol.xbp(:,2), 'r.', hor, ver, 'b')

% Plot data: panel (b)

figure(2); clf; hold on; grid on; box on
thm = struct('lspec', {{'ro', 'MarkerSize', 5}});
coco_plot_bd(thm, 'run_again', 'p', 'p', @(p) sqrt(2*pi)./sqrt(p.^2-2*p+2));
coco_plot_bd('run_again', 'p', '||po.orb.x||_{L_2[0,T]}')
hold off
