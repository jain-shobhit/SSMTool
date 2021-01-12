%% The Lienard vector field from Section 8.3 of Recipes for Continuation
%
% We study periodic solutions of the Lienard dynamical system x1' = x2, x2'
% = p*x_2 - x2^3 - x1 under variations in the problem parameter p. A unique
% solution is obtained by relying on a Poincare section defined by a point
% on a previously located periodic orbit and the corresponding vector
% field. In the absence of regular updates to the Poincare section, an
% artificial fold occurs during continuation as the orbit becomes tangent
% to the section. This is eliminated by updating the Poincare section
% before each new continuation step.

% Figure 1 shows a family of periodic orbits obtained with fixed Poincare
% section. Figure 2 shows a family of periodic orbits obtained with moving
% Poincare section. In each case, the dots denote initial point on the
% orbit. An artifical fold is observed in Figure 1.

%% Encoding

% The continuation problem encoded below includes a single monitor function
% that evaluates to the problem parameter, and a corresponding inactive
% continuation parameter 'p'. Its dimensional deficit equals 0. The call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'p' is released and allowed to
% vary during continuation.

% Initialize boundary condition data for Poincare section

p0 = 1;
x0 = [0.4; -1.2];
data = struct();
data.fhan = @lienard;
data = per_bc_update(data, [], x0, [], p0);

% Initialize solution guess
f  = @(t,x) lienard(x, p0);
[t0, x0]   = ode45(f, [0 6.7], x0);

% Construct boundary-value problem continuation problem.

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30);
coll_args = {@lienard, t0, x0, p0};
bvp_args  = {@per_bc, data, @per_bc_update, 'F+DF'};
prob = ode_isol2bvp(prob, '', coll_args{:}, 'p', bvp_args{:});

% Start continuation of periodic orbits with updates to the Poincare
% section before each continuation step.

prob = coco_set(prob, 'cont', 'NPR', 2, 'NAdapt', 10);
cont_args = { 1, {'p' 'bvp.seg.coll.err' 'bvp.seg.coll.err_TF'}, [-1 1] };

fprintf('\n Run=''%s'': Continue periodic orbits with boundary conditions updates.\n', ...
  'run_moving');

coco(prob, 'run_moving', [], cont_args{:});

% Repeat continuation of periodic orbits without updates to the Poincare
% section.

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30);
bvp_args = {@per_bc, data, 'F+DF'};
prob = ode_isol2bvp(prob, '', coll_args{:}, 'p', bvp_args{:});
prob = coco_set(prob, 'cont', 'NPR', 2, 'NAdapt', 10);

fprintf('\n Run=''%s'': Continue periodic orbits without boundary conditions updates.\n', ...
  'run_fixed');

coco(prob, 'run_fixed', [], 1, 'p', [-1 1]);

%% Graphical representation of stored solutions

% Plot data: Figure 1

bd = coco_bd_read('run_fixed'); % Extract bifurcation data
labs = coco_bd_labs(bd);        % Extract labels
x0 = [];

figure(1); clf; hold on; grid on; box on; axis([-1.75 1.75 -1.4 1.23])

for lab=labs
  sol = bvp_read_solution('', 'run_fixed', lab);
  plot(sol{1}.xbp(:,1), sol{1}.xbp(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7], 'Marker', '.', 'MarkerSize', 12)
  x0 = [x0 ; sol{1}.xbp(1,:)]; %#ok<AGROW>
end
plot(x0(:,1), x0(:,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15)

hold off

% Plot data: Figure 2

bd = coco_bd_read('run_moving'); % Extract bifurcation data
labs = coco_bd_labs(bd);         % Extract labels
x0 = [];

figure(2); clf; hold on; grid on; box on; axis([-1.75 1.75 -1.4 1.23])

for lab=labs
  sol = bvp_read_solution('', 'run_moving', lab);
  plot(sol{1}.xbp(:,1), sol{1}.xbp(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7], 'Marker', '.', 'MarkerSize', 12)
  x0 = [x0 ; sol{1}.xbp(1,:)]; %#ok<AGROW>
end
plot(x0(:,1), x0(:,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15)

hold off
