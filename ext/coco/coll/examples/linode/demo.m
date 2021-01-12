%% Continuation of periodic orbits from Example 10.2 in Recipes for Continuation

% Periodic solutions of the differential equations x''+x'+px=cos(t) are
% given by x(t)=(sin(t)+(p-1)*cos(t))/(p^2-2*p+2). A sample solution and
% the approximant obtained by linear interpolation of the base points found
% using continuation is shown in Figures 1 and 2 for the autonomous and
% non-autonomous encodings, respectively. Figure 3 shows the computed L2
% norms for each of the encodings and the exact value
% sqrt(2*pi)/sqrt(p^2-2*p+2) for the non-autonomous case.

%% Autonomous encoding

% The continuation problem structure encoded below includes a single
% monitor function that evaluates to the problem parameter, and the
% corresponding inactive continuation parameter 'p' representing the
% stiffness. Its dimensional deficit equals 0. A one-dimensional family of
% periodic orbits is obtained by releasing 'p' and allowing it to vary
% during continuation. The implementation is for an autonomous encoding of
% the harmonically excited linear oscillator.

[t0, x0] = ode45(@(t,x) linode_aut(x,1), [0 2*pi], [0; 1; 0]); % Approximate periodic orbit
prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 15, 'var', true);
prob = coco_set(prob, 'cont', 'NAdapt', 1);
coll_args = {@linode_aut, @linode_aut_DFDX, @linode_aut_DFDP, t0, x0, 'p', 1};
prob = ode_isol2coll(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_aut_bc, data, 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx]));

if exist('po_mult_add', 'file')
  prob = po_mult_add(prob, ''); % Store Floquet multipliers with bifurcation data
end

fprintf('\n Run=''%s'': Continue segments along primary branch.\n', ...
  'aut_run1');

coco(prob, 'aut_run1', [], 1, {'p' 'coll.err_TF'}, [0.2 2]);

% To extract the Floquet multipliers from the bifurcation data, use the
% following command:

% bd  = coco_bd_read('aut_run1');
% las = coco_bd_col(bd, {'multipliers(1)', 'multipliers(2)', 'multipliers(3)'});

%% Restart from a previously stored solution

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 1);
prob = ode_coll2coll(prob, '', 'aut_run1', 3);
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_aut_bc, data, 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx]));

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s''.\n', ...
  'aut_run2', 3, 'aut_run1');

coco(prob, 'aut_run2', [], 1, {'p' 'coll.err_TF'}, [0.2 2]);

%% Restart from a stored solution and include variational problem

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 1);
prob = ode_coll2coll(prob, '', 'aut_run1', 3, '-var', eye(3));
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_aut_bc, data, 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx]));

% Hold the initial condition of solution to variational problem fixed
[data, uidx] = coco_get_func_data(prob, 'coll.var', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.coll_var.v0_idx,:), ...
  {'s1' 's2' 's3' 's4' 's5' 's6' 's7' 's8' 's9'});

fprintf(...
  '\n Run=''%s'': Continue segments with variational problem from point %d in run ''%s''.\n', ...
  'aut_run_var', 3, 'aut_run1');

coco(prob, 'aut_run_var', [], 1, {'p' 'coll.err_TF'}, [0.05 3]);

% Use the following commands to extract the monodromy matrix associated
% with a particular solution.

% chart = coco_read_solution('coll.var', 'aut_run_var', 2, 'chart');
% data  = coco_read_solution('coll', 'aut_run_var', 2, 'data');
% M = chart.x(data.coll_var.v1_idx);

% The eigenvalues of the monodromy matrix are the Floquet multipliers.
% These can be obtained explicitly as 1 and exp(-pi+/-sqrt(1-4*p)*pi).

% p = chart.x(data.coll_seg.maps.p_idx);
% sort([eig(M) [1; exp(-pi+sqrt(1-4*p)*pi); exp(-pi-sqrt(1-4*p)*pi)]])

%% Non-autonomous encoding

% The continuation problem structure encoded below includes a single
% monitor function that evaluates to the problem parameter, and the
% corresponding inactive continuation parameter 'p' representing the
% stiffness. Its dimensional deficit equals 0. A one-dimensional family of
% periodic orbits is obtained by releasing 'p' and allowing it to vary
% during continuation. The implementation is for a non-autonomous encoding
% of the harmonically excited linear oscillator.

[t0, x0]   = ode45(@(t,x) linode_het(t,x,1), [0 2*pi], [0; 1]); % Approximate periodic orbit
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'coll', 'NTST', 15, 'var', true);
prob = coco_set(prob, 'cont', 'NAdapt', 1);
coll_args = {@linode_het, @linode_het_DFDX, @linode_het_DFDP, ...
  @linode_het_DFDT, t0, x0, 'p', 1};
prob = ode_isol2coll(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_het_bc, data, 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));

if exist('po_mult_add', 'file')
  prob = po_mult_add(prob, ''); % Store Floquet multipliers with bifurcation data
end

fprintf('\n Run=''%s'': Continue segments along primary branch.\n', ...
  'het_run1');

coco(prob, 'het_run1', [], 1, {'p' 'coll.err_TF'}, [0.2 2]);

% To extract the Floquet multipliers from the bifurcation data, use the
% following command:

% bd  = coco_bd_read('het_run1');
% las = coco_bd_col(bd, {'multipliers(1)', 'multipliers(2)'});

%% Restart from a previously stored solution

prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'cont', 'NAdapt', 1);
prob = ode_coll2coll(prob, '', 'het_run1', 3);
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_het_bc, data, 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s''.\n', ...
  'het_run2', 3, 'het_run1');

coco(prob, 'het_run2', [], 1, {'p' 'coll.err_TF'}, [0.2 2]);

%% Restart from a stored solution and include variational problem

prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'coll', 'NTST', 25);
prob = ode_coll2coll(prob, '', 'het_run1', 3, '-var', eye(2));
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_het_bc, data, 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));

% Hold the initial condition of solution to variational problem fixed
[data, uidx] = coco_get_func_data(prob, 'coll.var', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.coll_var.v0_idx(:)), ...
  {'s1' 's2' 's3' 's4'});

fprintf(...
  '\n Run=''%s'': Continue segments with variational problem from point %d in run ''%s''.\n', ...
  'het_run_var', 3, 'het_run1');

coco(prob, 'het_run_var', [], 1, {'p' 'coll.err_TF'}, [0.05 3]);

% Use the following commands to extract the monodromy matrix associated
% with a particular solution.

% chart = coco_read_solution('coll.var', 'het_run_var', 2, 'chart');
% data  = coco_read_solution('coll', 'het_run_var', 2, 'data');
% M = chart.x(data.coll_var.v1_idx);

% The eigenvalues of the monodromy matrix are the Floquet multipliers.
% These can be obtained explicitly as exp(-pi+/-sqrt(1-4*p)*pi).

% p = chart.x(data.coll_seg.maps.p_idx);
% sort([eig(M) [exp(-pi+sqrt(1-4*p)*pi); exp(-pi-sqrt(1-4*p)*pi)]])

%% Graphical representation of stored solutions

% Figure 1
figure(1); clf; hold on; grid on; box on; axis([0 2*pi -1 1 -1 1])
coco_plot_sol('aut_run_var', 5, '', 't', ...
  {'t', 'p'}, @(t,p) (sin(t)+(p-1)*cos(t))/(p^2-2*p+2), ...
  {'t', 'p'}, @(t,p) (cos(t)-(p-1)*sin(t))/(p^2-2*p+2))
thm = struct();
thm.sol.RO = {'ro', 'MarkerSize', 6};
thm.xlab = 'x_3 mod 2\pi';
coco_plot_sol(thm, 'aut_run_var', 5, '', 'x' , @(x) mod(x(:,3),2*pi), ...
  'x', 'x')
view(-10,50)

% Figure 2
figure(2); hold on; grid on; box on; axis([-1 1 -1 1])
thm = struct();
thm.sol.RO = {'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'w'};
coco_plot_sol(thm, 'het_run_var', 5, '', ...
  {'t', 'p'}, @(t,p) (sin(t)+(p-1)*cos(t))/(p^2-2*p+2), ...
  {'t', 'p'}, @(t,p) (cos(t)-(p-1)*sin(t))/(p^2-2*p+2))
coco_plot_sol('het_run_var', 5, '', 'x', 'x')
hold off

% Figure 3
figure(3); clf; hold on; grid on; box on
coco_plot_bd('aut_run_var', 'p', '||x||_{L_2[0,T]}')
coco_plot_bd('het_run_var', 'p', '||x||_{L_2[0,T]}')
thm = struct();
thm.lspec = {'ro', 'MarkerSize', 4};
thm.ylab  = '||x||_{L_2[0,T]}';
coco_plot_bd(thm, 'aut_run_var', 'p', 'p', ...
  @(p) sqrt(2*pi)./sqrt(p.^2-2*p+2));
hold off
