%% Inflection points along the frequency-response curve for a linear oscillator
%
% This demo illustrates the successive continuation approach to searching
% for stationary points that satisfy the necessary conditions for a local
% extremum in the difference between the maximal values of x1 along two
% periodic solutions to the harmonically excited linear oscillator
%
%   x1' = x2, x2' = -x1-x2+cos x3, x3' = omega
%
% with values of omega that differ by a small fixed amount. In the first
% stage of continuation, local extrema in 'ampdiff' are located along a
% one-dimensional solution manifold with trivial Lagrange multipliers. Each
% of these is a branch point, from which emanates a secondary
% one-dimensional submanifold along which the Lagrange multipliers take on
% nontrivial values. We terminate continuation along such a manifold when
% 'd.ampdiff' equals 1.

%% Initial encoding

prob = coco_prob;
prob = coco_set(prob, 'cont', 'NAdapt', 1);

fcns = {@linode, @linode_dx, @linode_dp, @linode_dxdx, ...
  @linode_dxdp, @linode_dpdp};
bc_funcs = {@linode_bc, @linode_bc_du, @linode_bc_dudu};

%% first run to find initial fold
prob1 = coco_set(prob, 'coll', 'NTST', 30);

% zero problems
[t0, x0]  = ode45(@(t,x) linode(x, 0.98), [0 2*pi/0.98], ...
  [1.01958; 0; 1.56164*0.98]);
coll_args = [fcns, {t0, x0, 0.98}];
bvp_args  = [coll_args, [{'om'}, bc_funcs]];
prob1 = ode_isol2bvp(prob1, 'orb1', bvp_args{:});

[t0, x0]  = ode45(@(t,x) linode(x, 0.88), [0 2*pi/0.88], ...
  [1.10077; 0; 1.49982*0.88]);
coll_args = [fcns, {t0, x0, 0.88}];
bvp_args  = [coll_args, bc_funcs];
prob1 = ode_isol2bvp(prob1, 'orb2', bvp_args{:});

[data1, uidx1] = coco_get_func_data(prob1, 'orb1.bvp.seg1.coll', ...
  'data', 'uidx');
maps1 = data1.coll_seg.maps;
[data2, uidx2] = coco_get_func_data(prob1, 'orb2.bvp.seg1.coll', ...
'data', 'uidx');
maps2 = data2.coll_seg.maps;

prob1 = coco_add_glue(prob1, 'glue', uidx1(maps1.p_idx), ...
  uidx2(maps2.p_idx), -0.1);
prob1 = coco_add_glue(prob1, 'amp', uidx1(maps1.x0_idx(1)), ...
  uidx2(maps2.x0_idx(1)), 'ampdiff', 'inactive');

% adjoints
prob1 = adjt_isol2bvp(prob1, 'orb1');
prob1 = adjt_isol2bvp(prob1, 'orb2');

[data1, axidx1] = coco_get_adjt_data(prob1, 'orb1.bvp.seg1.coll', ...
  'data', 'axidx');
opt1 = data1.coll_opt;
[data2, axidx2] = coco_get_adjt_data(prob1, 'orb2.bvp.seg1.coll', ...
  'data', 'axidx');
opt2 = data2.coll_opt;

prob1 = coco_add_adjt(prob1, 'glue', 'aidx', ...
  [axidx1(opt1.p_idx); axidx2(opt2.p_idx)]);
prob1 = coco_add_adjt(prob1, 'amp', 'd.ampdiff', 'aidx', ...
  [axidx1(opt1.x0_idx(1)); axidx2(opt2.x0_idx(1))]);

% continuation
cont_pars = {'ampdiff' 'om' 'd.ampdiff'};
coco(prob1, 'linode1', [], 1, cont_pars, {[-1 1] [0.4 1.4]});

%% branch switch from fold to grow nontrivial adjoint
bd1   = coco_bd_read('linode1');
BPlab = coco_bd_labs(bd1, 'BP');

% zero problems
prob2 = ode_BP2bvp(prob, 'orb1', 'linode1', BPlab(1));
prob2 = ode_BP2bvp(prob2, 'orb2', 'linode1', BPlab(1));

[data1, uidx1] = coco_get_func_data(prob2, 'orb1.bvp.seg1.coll', ...
  'data', 'uidx');
maps1 = data1.coll_seg.maps;
[data2, uidx2] = coco_get_func_data(prob2, 'orb2.bvp.seg1.coll', ...
'data', 'uidx');
maps2 = data2.coll_seg.maps;

prob2 = coco_add_glue(prob2, 'glue', uidx1(maps1.p_idx), ...
  uidx2(maps2.p_idx), -0.1);
prob2 = coco_add_glue(prob2, 'amp', uidx1(maps1.x0_idx(1)), ...
  uidx2(maps2.x0_idx(1)), 'ampdiff', 'inactive');

% branch switch data
chart = coco_read_solution('linode1', BPlab(1), 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

% adjoints
prob2 = adjt_BP2bvp(prob2, 'orb1', 'linode1', BPlab(1));
prob2 = adjt_BP2bvp(prob2, 'orb2', 'linode1', BPlab(1));

[data1, axidx1] = coco_get_adjt_data(prob2, 'orb1.bvp.seg1.coll', ...
  'data', 'axidx');
opt1 = data1.coll_opt;
[data2, axidx2] = coco_get_adjt_data(prob2, 'orb2.bvp.seg1.coll', ...
  'data', 'axidx');
opt2 = data2.coll_opt;

[chart, lidx] = coco_read_adjoint('glue', 'linode1', BPlab(1), ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'glue', 'aidx', ...
  [axidx1(opt1.p_idx); axidx2(opt2.p_idx)], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('amp', 'linode1', BPlab(1), ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'amp', 'd.ampdiff', 'aidx', ...
  [axidx1(opt1.x0_idx(1)); axidx2(opt2.x0_idx(1))], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

% continuation
cont_pars = {'d.ampdiff' 'ampdiff' 'om'};
coco(prob2, 'linode2', [], 1, cont_pars, {[0 1] [-1 1] [0.4 1.4]});
