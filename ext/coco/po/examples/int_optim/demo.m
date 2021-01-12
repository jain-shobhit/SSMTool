%% Constrained optimization along a family of periodic orbits
%
% We use the method of successive continuation to look for stationary
% values of an integral functional along a family of periodic orbits of a
% modified Van-der-Pol oscillator.

%% Initial encoding

% Construct equilibrium point zero problem
funcs  = {@mvdP, @mvdP_dx, @mvdP_dp, @mvdP_dxdx, @mvdP_dxdp, @mvdP_dpdp};
x0     = [0;0;0];
pnames = {'p1', 'p2', 'p3', 'p4'};
p0     =  [0.5; 4; 0; 2];

prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, x0, pnames, p0);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');

bd1  = coco(prob, 'ep_run', [], 1, 'p3', [0 1]);

%% Start continuation of periodic orbits from the Hopf bifurcation

HBlab = coco_bd_labs(bd1, 'HB');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [0 100]);
prob = coco_set(prob, 'po', 'bifus', false);

% periodic orbit problem with adjoint
prob = ode_HB2po(prob, '', 'ep_run', HBlab);
prob = adjt_isol2po(prob, '');

% integral functional with adjoint
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data = int_init_data(fdata, 'po.orb');
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);
prob = coco_add_func(prob, 'po.orb.int', @int, @int_du, data, ...
  'inactive', 'po.orb.int', 'uidx', uidx, 'remesh', @int_remesh);

[fdata, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
data = adjt_int_init_data(fdata, 'po.orb');
opt  = data.coll_opt;
aidx = axidx([opt.xcn_idx; opt.T_idx]);
prob = coco_add_adjt(prob, 'po.orb.int', @adjt_int, @adjt_int_du, data, ...
    'd.po.orb.int', 'aidx', aidx, 'remesh', @adjt_int_remesh);

fprintf(...
  '\n Run=''%s'': Continue periodic orbits from point %d in run ''%s''.\n', ...
  'po_run', HBlab, 'ep_run');

cont_args = {1, {'po.orb.int', 'p3', 'd.p1', 'd.p2', 'd.p4', ...
  'd.po.orb.int'}};
bd2  = coco(prob, 'po_run', [], cont_args{:});

%% Switch to branch of nontrivial multipliers

BPlab = coco_bd_labs(bd2, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 1);
prob = coco_set(prob, 'po', 'bifus', false);

% periodic orbit problem with adjoint
prob = ode_BP2po(prob, '', 'po_run', BPlab);
prob = adjt_BP2po(prob, '', 'po_run', BPlab);

% integral functional with adjoint
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data = int_init_data(fdata, 'po.orb');
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);
prob = coco_add_func(prob, 'po.orb.int', @int, @int_du, data, ...
  'inactive', 'po.orb.int', 'uidx', uidx, 'remesh', @int_remesh);

chart = coco_read_solution('po_run', BPlab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
[chart, lidx] = coco_read_adjoint('po.orb.int', 'po_run', BPlab, ...
  'chart', 'lidx');

[fdata, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
data = adjt_int_init_data(fdata, 'po.orb');
opt  = data.coll_opt;
aidx = axidx([opt.xcn_idx; opt.T_idx]);
prob = coco_add_adjt(prob, 'po.orb.int', @adjt_int, @adjt_int_du, data, ...
  'd.po.orb.int', 'aidx', aidx, 'remesh', @adjt_int_remesh, ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

cont_args = {1, {'d.po.orb.int', 'po.orb.int', 'p3', 'd.p1', 'd.p2', ...
  'd.p4'}, [0 1]};
bd3  = coco(prob, 'po_run_lagrange1', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters

EPlab = coco_bd_labs(bd3, 'EP');

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 10);
prob = coco_set(prob, 'po', 'bifus', false);

% periodic orbit problem with adjoint
prob = ode_po2po(prob, '', 'po_run_lagrange1', EPlab(2));
prob = adjt_po2po(prob, '', 'po_run_lagrange1', EPlab(2));

% integral functional with adjoint
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data = int_init_data(fdata, 'po.orb');
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);
prob = coco_add_func(prob, 'po.orb.int', @int, @int_du, data, ...
  'inactive', 'po.orb.int', 'uidx', uidx, 'remesh', @int_remesh);

chart = coco_read_adjoint('po.orb.int', 'po_run_lagrange1', EPlab(2), ...
  'chart');
[fdata, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
data = adjt_int_init_data(fdata, 'po.orb');
opt  = data.coll_opt;
aidx = axidx([opt.xcn_idx; opt.T_idx]);
prob = coco_add_adjt(prob, 'po.orb.int', @adjt_int, @adjt_int_du, data, ...
  'd.po.orb.int', 'aidx', aidx, 'remesh', @adjt_int_remesh, ...
  'l0', chart.x);

prob = coco_add_event(prob, 'OPT', 'BP', 'd.p2', '==', 0);

cont_args = {1, {'d.p2', 'po.orb.int', 'p3', 'd.p1', 'p2', 'd.p4'}};
bd4  = coco(prob, 'po_run_lagrange2', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters

OPTlab = coco_bd_labs(bd4, 'OPT');

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 10);
prob = coco_set(prob, 'po', 'bifus', false);

% periodic orbit problem with adjoint
prob = ode_po2po(prob, '', 'po_run_lagrange2', OPTlab);
prob = adjt_po2po(prob, '', 'po_run_lagrange2', OPTlab);

% integral functional with adjoint
[fdata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
data = int_init_data(fdata, 'po.orb');
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);
prob = coco_add_func(prob, 'po.orb.int', @int, @int_du, data, ...
  'inactive', 'po.orb.int', 'uidx', uidx, 'remesh', @int_remesh);

chart = coco_read_adjoint('po.orb.int', 'po_run_lagrange2', OPTlab, ...
  'chart');
[fdata, axidx] = coco_get_adjt_data(prob, 'po.orb.coll', 'data', 'axidx');
data = adjt_int_init_data(fdata, 'po.orb');
opt  = data.coll_opt;
aidx = axidx([opt.xcn_idx; opt.T_idx]);
prob = coco_add_adjt(prob, 'po.orb.int', @adjt_int, @adjt_int_du, data, ...
  'd.po.orb.int', 'aidx', aidx, 'remesh', @adjt_int_remesh, ...
  'l0', chart.x);

prob = coco_add_event(prob, 'OPT', 'BP', 'd.p1', '==', 0);
cont_args = {1, {'d.p1', 'po.orb.int', 'p3', 'p1', 'p2', 'd.p4'}};
bd5  = coco(prob, 'po_run_lagrange3', [], cont_args{:});
