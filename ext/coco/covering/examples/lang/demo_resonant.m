% Analysis of a two-dimensional manifold of resonant periodic orbits using
% the atlas_kd algorithm. Examples illustrate the use of mesh adaptation
% and projections onto invariant geometries.

%% Continuation without adaptation, projection onto 'x01', 'x02', 'rho', and 'eps'
coco_func_data.pointers('set', []);
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

p0 = [3.5; 0.35; 0.1];
[t, x0] = ode45(@(t,x) lang(x, p0), [0 5.3], [0.3; 0; 0.7]); % Approximate periodic orbit

prob = coco_set(prob, 'coll', 'NTST', 50);
prob = ode_isol2bvp(prob, '', @lang, @lang_DFDX, @lang_DFDP,  ...
  t, x0, p0, {'om' 'rho' 'eps'}, @po_bc, @po_bc_DFDX);
[data, uidx] = coco_get_func_data(prob, 'bvp.seg1.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'project', ...
  uidx(maps.xbp_idx(1:2)), {'x01' 'x02'});

prob = coco_set(prob, 'cont', 'PtMX', 2000, 'NPR', 100, ...
  'almax', 20, 'NAdapt', 0, 'R_max', 0.4);
coco(prob, 'run_no_adapt_proj', [], 2, {'x01' 'x02' 'rho' 'eps'}, ...
  {[] [] [] [0 0.5]});

atlas = coco_bd_read('run_no_adapt_proj', 'atlas');
figure(1)
clf; hold on
plot_atlas_kd(atlas.charts,1,2,4)
axis equal; axis tight;
axis([-inf inf -inf inf 0 0.5])
view(60,30);
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.9
light('Position', [-0.1 0 0], 'Style', 'local')

figure(2)
clf; hold on
plot_atlas_kd(atlas.charts,1,3,4)
axis equal; axis tight;
axis([-inf inf -inf inf 0 0.5])
view(-60,20)
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.9
light('Position', [-0.1 0 0], 'Style', 'local')

%% Continuation with adaptation, projection onto 'x01', 'x02', 'rho', and 'eps'

prob = coco_set(prob, 'cont', 'PtMX', 2000, 'NPR', 100, ...
  'almax', 20, 'NAdapt', 1, 'R_max', 0.4);
coco(prob, 'run_with_adapt_proj', [], 2, {'x01' 'x02' 'rho' 'eps'}, ...
  {[] [] [] [0 0.5]});

atlas = coco_bd_read('run_with_adapt_proj', 'atlas');
figure(5)
clf; hold on
plot_atlas_kd(atlas.charts,1,2,4)
axis equal; axis tight;
axis([-inf inf -inf inf 0 0.5])
view(60,30);
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.9
light('Position', [-0.1 0 0], 'Style', 'local')

figure(6)
clf; hold on
plot_atlas_kd(atlas.charts,1,3,4)
axis equal; axis tight;
axis([-inf inf -inf inf 0 0.5])
view(-60,20)
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.9
light('Position', [-0.1 0 0], 'Style', 'local')

%% Continuation without adaptation, projection onto all unknowns

coco_func_data.pointers('set', []);
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

p0 = [3.5; 0.35; .1];
[t, x0] = ode45(@(t,x) lang(x, p0), [0 5.3], [0.3; 0; 0.7]); % Approximate periodic orbit

prob = coco_set(prob, 'coll', 'NTST', 50);
prob = ode_isol2bvp(prob, '', @lang, @lang_DFDX, @lang_DFDP,  ...
  t, x0, p0, {'om' 'rho' 'eps'}, @po_bc, @po_bc_DFDX);
[data, uidx] = coco_get_func_data(prob, 'bvp.seg1.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
xbp  = maps.xbp_idx;
labels = cell(1,numel(xbp));
ranges = cell(1,numel(xbp));
for i=1:numel(xbp)
  labels{i} = sprintf('x%d', xbp(i));
end

prob = coco_add_pars(prob, 'project', ...
  uidx([xbp; maps.T_idx]), [labels {'T'}]);

prob = coco_set(prob, 'cont', 'PtMX', 2000, 'NPR', 100, ...
  'almax', 20, 'NAdapt', 0, 'R_max', 0.4);
coco(prob, 'run_no_adapt_no_proj', [], 2, [labels {'T' 'rho' 'eps'}], ...
  [ranges {[] [] [0 0.5]}]);

atlas = coco_bd_read('run_no_adapt_no_proj', 'atlas');
figure(5)
clf; hold on
plot_atlas_kd(atlas.charts,1,2,753)
axis equal; axis tight;
axis([-inf inf -inf inf 0 0.5])
view(60,30); alpha 0.5
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.9
light('Position', [-0.1 0 0], 'Style', 'local')

figure(4)
clf; hold on
plot_atlas_kd(atlas.charts,2,752,753)
axis equal; axis tight;
axis([-inf inf -inf inf 0 0.5])
view(60,30)
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.9
light('Position', [-0.1 0 0], 'Style', 'local')
