coco_use_recipes_toolbox bvp_v1 coll_v1 % Add atlas2d_v6 atlas algorithm and bvp_v1 and coll_v1 toolboxes to search path

addpath ../atlas2d_v6

%% Resonance manifold

p0 = [3.5; 0.35; 0.1];
[t, x0] = ode45(@(t,x) lang(x, p0), [0 5.3], [0.3; 0; 0.4]);

prob = coco_prob();
prob = bvp_isol2seg(prob, '', @lang, @lang_DFDX, @lang_DFDP,  ...
  t, x0, {'om' 'ro' 'eps'}, p0, @po_bc, @po_bc_DFDX);
[data, uidx] = coco_get_func_data(prob, 'bvp.seg.coll', 'data', 'uidx');
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'R', .01, 'PtMX', 5000, 'NPR', 10, ...
  'indcs', uidx([data.x0_idx; data.p_idx]));
coco(prob, 'run1', [], 2, {'ro' 'eps'}, {[] [0.1 0.5]});

coco_use_recipes_toolbox % Remove atlas2d_v6 atlas algorithm and bvp_v1 and coll_v1 toolboxes from search path

%% Graphical representation of stored solutions

atlas = coco_bd_read('run1', 'atlas');

figure(1); clf; hold on
plot_trisurf(atlas.charts, 153, 154, 1)
grid on; view([30 20]); 
alpha 0.4; hold off

figure(2); clf; hold on
plot_trisurf(atlas.charts, 3, 154, 1)
grid on; view([30 20]); 
alpha 0.4; hold off

rmpath ../atlas2d_v6
