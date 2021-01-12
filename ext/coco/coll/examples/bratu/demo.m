%% The Bratu boundary-value problem from Section 8.1 of Recipes for Continuation
%
% We study solutions to the two-point boundary-value problem given by the
% ordinary differential equation u_xx + p*u = 0 and u vanishing at x=0 and
% x=1 under variations in the problem parameter p.

% Figure 1 shows an initial solution guess obtained with an Euler step at
% t=0 with step size h=0.04, or by single point, and compares this
% approximation with the exact solution. The grown solution is shown in
% Figure 2. From this solution one can, finally, start a continuation in
% y1e. Figure 3 shows the labeled solutions obtained in this final
% continuation run.

%% Encoding

% The continuation problem encoded below includes a single monitor function
% that evaluates to the problem parameter, and a corresponding inactive
% continuation parameter 'p'. Its dimensional deficit equals 0. The call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'p' is released and allowed to
% vary during continuation.

% Define vector field using inline, anonymous function.

brat    = @(x,p) [x(2,:); -p(1,:).*exp(x(1,:))];
brat_bc = @(~,T,x0,x1,p) [T-1; x0(1); x1(1)];
brat_bc_DFDX = @(~,T,x0,x1,p) [1,0,0,0,0,0; 0,1,0,0,0,0; 0,0,0,1,0,0];

% Initialize continuation problem and settings associated with toolbox
% constructor.

prob = coco_prob();
prob = coco_set(prob, 'cont', 'PtMX', 50);

% Compute family of solutions using nested call to ode_isol2bvp. Initial
% solution guess is the trivial solution u = 0 for p = 0.

coll_args = {brat, [0;1], zeros(2), 0};
bvp_args  = [coll_args, 'p', {brat_bc, brat_bc_DFDX}];

fprintf('\n Run=''%s'': Continue family of constrained trajectory segments.\n', ...
  'brat1');

bd = coco(prob, 'brat1', @ode_isol2bvp, bvp_args{:}, 1, 'p', [0 4]);

%% Restart continuation from stored solution

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5); % Remesh every 5 steps

fprintf(...
  '\n Run=''%s'': Continue constrained trajectory segments from point %d in run ''%s''.\n', ...
  'brat2', 6, 'brat1');

coco(prob, 'brat2', @ode_bvp2bvp, 'brat1', 6, 1, 'p', [1 5]);

%% Graphical representation of stored solutions

% Plot data: Figure 1
figure(1); clf; hold on; grid on; box on; axis([0 1 -0.1 1.5])
coco_plot_sol('brat1', '')
hold off

% Plot data: Figure 2
figure(2); clf; hold on; grid on; box on; axis([0 1 -5 5])
coco_plot_sol('brat1', '', 't', 'x', 2)
hold off

drawnow

%% Continuation in alternative parameterization with branch point

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'NPR', 100);
prob = ode_isol2bvp(prob, '', bvp_args{:});
[data, uidx] = coco_get_func_data(prob, 'bvp.seg1.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'trans', @trans, [], 'zero', 'uidx', ...
  uidx(data.coll_seg.maps.p_idx), 'u0', 0);
uidx = coco_get_func_data(prob, 'trans', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(end), 'C');

fprintf('\n Run=''%s'': Continue primary family of constrained trajectory segments.\n', ...
  'brat3');

bd3 = coco(prob, 'brat3', [], 1, {'C' 'p'}, [0 5]);

BPlab = coco_bd_labs(bd3, 'BP');
prob  = coco_prob();
prob  = coco_set(prob, 'cont', 'NAdapt', 5,'NPR', 100);
prob  = ode_BP2bvp(prob, '', 'brat3', BPlab);
chart = coco_read_solution('', 'brat3', BPlab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
[data, uidx] = coco_get_func_data(prob, 'bvp.seg1.coll', 'data', 'uidx');
pidx  = uidx(data.coll_seg.maps.p_idx);
prob  = coco_add_func(prob, 'trans', @trans, [], 'zero', 'uidx', pidx, ...
  'u0', 2.3994, 't0', cdata.v(pidx+1));
uidx  = coco_get_func_data(prob, 'trans', 'uidx');
prob  = coco_add_pars(prob, 'pars', uidx(end), 'C');

fprintf(...
  '\n Run=''%s'': Continue secondary family of constrained trajectory segments from point %d in run ''%s''.\n', ...
  'brat4', BPlab, 'brat3');

bd4 = coco(prob, 'brat4', [], 1, {'C' 'p'}, [0 5]);

%% Graphical representation of stored solutions

figure(3); clf; hold on; box on; grid on
thm = struct('oid', 'bvp.seg1');
coco_plot_bd(thm, 'brat3', 'C', '||bvp.seg1.x||_{L_2[0,1]}')
thm.lspec   = {'r-', 'LineWidth', 1};
thm.special = {'BP'};
coco_plot_bd(thm, 'brat4', 'C', '||bvp.seg1.x||_{L_2[0,1]}')
hold off
