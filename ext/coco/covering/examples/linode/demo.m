%% Continuation of periodic orbits from Example 10.2 in Recipes for Continuation

% Periodic solutions of the differential equations x''+x'+px=cos(t) are
% given by x(t)=(sin(t)+(p-1)*cos(t))/(p^2-2*p+2).

%% Autonomous encoding

coco_func_data.pointers('set', []);
p0   = 1;
[t0, x0] = ode45(@(t,x) linode_aut(x,p0), [0 2*pi], [0; 1/sqrt(p0^2-2*p0+2); atan(1-p0)]); % Approximate periodic orbit
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'MXCL', false);
prob = coco_set(prob, 'cont', 'atlas', 'kd', 'NPR', 1, ...
  'NAdapt', 4, 'R', .1, 'R_max', .1, 'R_min', .1);
coll_args = {@linode_aut, @linode_aut_DFDX, @linode_aut_DFDP, t0, x0, 'p', p0};
prob = ode_isol2coll(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_aut_bc, [], 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx]));

fprintf('\n Run=''%s'': Continue periodic orbits for autonomous encoding.\n', ...
  'auto');

coco(prob, 'auto', [], 1, {'p' 'coll.err_TF'}, [0.2 2]);

%% Non-autonomous encoding

coco_func_data.pointers('set', []); % only necessary for atlas_kd

prob = coco_prob;
prob = coco_set(prob, 'all', 'CleanData', true);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'coll', 'NTST', 50, 'MXCL', false);
p0   = [1; 1];
t0   = (0:2*pi/99:2*pi)';
x0   = [p0(2)*(sin(t0)+(p0(1)-1)*cos(t0))/(p0(1)^2-2*p0(1)+2) ...
  p0(2)*(cos(t0)-(p0(1)-1)*sin(t0))/(p0(1)^2-2*p0(1)+2) ];
prob = ode_isol2coll(prob, '', @linode, t0, x0, {'p1' 'p2'}, p0);
[data, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'po', @linode_bc, [], 'zero', ...
  'uidx', uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));
prob = coco_add_pars(prob, 'vel', uidx(maps.x0_idx(2)), 'vel', 'active');

%% Continue one-dimensional family of periodic orbits with atlas_1d
prob = coco_set(prob, 'cont', 'atlas', '1d', 'PtMX', 200,  ...
  'NAdapt', 1, 'h0', .2, 'h_max', .2, 'h_min', .2);
coco(prob, 'dim_1_atlas_1d', [], 1, 'p1', [0.2 2]);

%% Continue one-dimensional family of periodic orbits with atlas_kd
prob = coco_set(prob, 'cont', 'atlas', 'kd',  'PtMX', 200,  ...
  'NAdapt', 1, 'R', .2, 'R_max', 2, 'R_min', .02);
coco(prob, 'dim_1_atlas_kd', [], 1, {'p1' 'vel'}, [0.2 2]);

%% Continue two-dimensional family of periodic orbits with atlas_kd
prob = coco_set(prob, 'cont', 'atlas', 'kd',  'PtMX', 1000, ...
  'NAdapt', 1, 'R', .1, 'R_max', 10, 'R_min', .01);
coco(prob, 'dim_2_atlas_kd', [], 2, {'p1' 'p2' 'vel'}, {[0.2 2], [0 1]});

%% Graphical representation of stored solutions

% Figure 1
figure(1); clf; hold on; grid on; box on
thm = struct();
thm.lspec = {'bx', 'MarkerSize', 6};
thm.ylab  = '||x||_{L_2[0,T]}';
coco_plot_bd(thm, 'auto', 'p', '||x||_{L_2[0,T]}')
coco_plot_bd(thm, 'dim_1_atlas_1d', 'p1', '||x||_{L_2[0,T]}')
thm = struct();
thm.lspec = {'ro', 'MarkerSize', 8};
thm.ylab  = '||x||_{L_2[0,T]}';
coco_plot_bd(thm, 'dim_1_atlas_1d', 'p1', 'p1', ...
  @(p) sqrt(2*pi)./sqrt(p.^2-2*p+2));
coco_plot_bd(thm, 'auto', 'p', 'p', ...
  @(p) sqrt(2*pi./(p.^2-2*p+2)+8*pi^3/3+4*pi^2*atan(1-p)+3*atan((1-p).^2)));
hold off

% Figure 2
figure(2); clf; grid on
atlas = coco_bd_read('dim_1_atlas_kd', 'atlas');
plot_atlas_kd(atlas.charts,2)

% Figure 3
figure(3); clf; grid on
atlas = coco_bd_read('dim_2_atlas_kd', 'atlas');
plot_atlas_kd(atlas.charts,3)

