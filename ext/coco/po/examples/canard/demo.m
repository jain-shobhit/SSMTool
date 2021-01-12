%% A canard explosion from Section 20.3.2 of Recipes for Continuation
%
% A family of periodic orbits connects two Hopf bifurcations of the
% dynamical system
%
% x1' = eps * (a - x2), x2' = x1 + x2 - x2^3
%
% occurring at a = +/-1. Along this family, a dramatic increase of period
% over an exceedingly small variation in the problem parameter a shows up
% as two virtually vertical branches in a bifurcation diagram.

% A bifurcation diagram is shown in panel (a). Members of both canard
% families are shown in the subsequent panels. The last two panels include
% a discretization of the canard orbit with period 400. The fat dots in (a)
% and the open circles in (b) mark the end points of mesh intervals.

%% Initial encoding

% The equilibrium continuation problem encoded below includes two monitor
% functions that evaluate to each of the problem parameters, and the
% corresponding inactive continuation parameters 'a' and 'eps'. The
% dimensional deficit is 0. A one-dimensional solution manifold results
% from releasing 'a' and allowing it to vary during continuation.

% Continue a family of equilibria through Hopf bifurcations at p1=-1 and p1=1.

vanderpol = @(x,p) [p(2,:).*(p(1,:)-x(2,:)); x(1,:)+x(2,:)-x(2,:).^3/3];
pnames = {'a' 'eps' };
p0     = [ 0;  0.01 ];
x0     = [0; 0];
prob = coco_prob();
prob = ode_isol2ep(prob, '', vanderpol, x0, pnames, p0);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');

bd1  = coco(prob, 'ep_run', [], 1, 'a', [-1.1 1.1]);

%% Start continuation along family of canard periodic orbits

% The periodic orbit continuation problems encoded below include three
% monitor functions that evaluate to each of the problem parameters and the
% orbit period, respectively, and the corresponding inactive continuation
% parameters 'a' and 'eps' and active continuation parameter 'po.period'.
% The dimensional deficit equals 0. A one-dimensional solution manifold
% results by releasing 'a' and allowing it to vary during continuation.

HBlabs = coco_bd_labs(bd1, 'HB');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 60);
prob = coco_set(prob, 'po', 'bifus', false);
prob = ode_HB2po(prob, '', 'ep_run', HBlabs(1));
prob = coco_add_event(prob, 'UZ', 'po.period', [100:100:400 460]);
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 4, 'bi_direct', false);
cont_args = { 1, {'po.period' 'a' 'po.orb.coll.err'}, {[1 1000] [-1.1 0]} };

fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  'po_run1', HBlabs(1), 'ep_run');

coco(prob, 'po_run1', [], cont_args{:});

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 60);
prob = coco_set(prob, 'po', 'bifus', false);
prob = ode_HB2po(prob, '', 'ep_run', HBlabs(2));
prob = coco_add_event(prob, 'UZ', 'po.period', [100:100:400 460]);
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'NPR', 100);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 4, 'bi_direct', false);
cont_args = { 1, {'po.period' 'a' 'po.orb.coll.err'}, {[1 1000] [0 1.1]} };

fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  'po_run2', HBlabs(2), 'ep_run');

coco(prob, 'po_run2', [], cont_args{:});

%% Graphical representation of stored solutions

% Plot data: panel (a)

figure(1); clf; hold on
thm = struct('ustab', '', 'ylab', 'period');
coco_plot_bd(thm, 'po_run1', 'a', 'po.period')
coco_plot_bd(thm, 'po_run2', 'a', 'po.period')
grid on; box on; axis([-1.1 1.1 85 500]); hold off

% Plot data: panels (b)-(i)

bd = coco_bd_read('po_run1');
labs1 = coco_bd_labs(bd, 'UZ');
for i=1:numel(labs1)
  figure(i+1); clf;
  coco_plot_sol('po_run1', labs1(i), '', 'x', 'x')
  grid on; box on; axis([-1 1 -3 3])
end

bd = coco_bd_read('po_run2');
labs2 = coco_bd_labs(bd, 'UZ');
for i=1:numel(labs2)
  figure(i+1+numel(labs1)); clf
  coco_plot_sol('po_run2', labs2(numel(labs2)-i+1), '', 'x', 'x')
  grid on; box on; axis([-1 1 -3 3])
end

% Extract data
lab  = labs1(6);
sol  = po_read_solution('po_run1', lab);

N   = (numel(sol.tbp)-1)/4;
idx = 1:4:N*4+1;

% Plot data: panel (j)

figure(numel(labs1)+numel(labs2)+2); clf; hold on

plot(sol.xbp(idx,1), sol.xbp(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
coco_plot_sol('po_run1', lab, '', 'x', 'x')

grid on; box on; axis([-1 0.8 -3 3]); hold off

% Plot data: panel (k)

figure(numel(labs1)+numel(labs2)+3); clf; hold on

tt  = sol.tbp/sol.tbp(end);
coco_plot_sol('po_run1', lab, '', 'tn', 'x')
coco_plot_sol('po_run1', lab, '', 'tn', 'x', 2)
plot(tt(idx), -2.3, 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
ylabel('x_1, x_2')

grid on; box on; axis([0 1 -2.5 2]); hold off
