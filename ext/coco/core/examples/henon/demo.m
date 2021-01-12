% Demo for periodic orbits of Henon map
% 
% For additional detail, see Henon demo in Recipes for Continuation and
% coco/recipes/demo_henon.m
%
% For functions that use no function data, the anonymous conversion of a
% plain function to a coco function (courtesy of Jan Sieber) eliminates the
% need for an external function definition.

clear
fcn = @(f) @(p,d,u) deal(d, f(u)); % convert plain fcn to coco fcn

%% Initial guess for period-4 orbit of Henon map
p0 = [ 1.0  0.3 ];
x0 = [ 1.3  0.0; -0.7  0.4; 1.0 -0.2; -0.1  0.3 ];

%% Henon map and its residual
henon = @(x,a,b) [ x(2)+1-a*x(1)^2; b*x(1) ];
henon_res = @(x,p,y) y-henon(x, p(1), p(2));
ip = 1:2;
ix = 3:4;
iy = 5:6;
f  = @(u) henon_res(u(ix), u(ip), u(iy));
period = size(x0,1);

%% Append 4 copies of the Henon map to the continuation problem
prob = coco_prob;
for i=1:period
  inxt = mod(i, period)+1;
  prob = coco_add_func(prob, ['henon' num2str(i)], fcn(f), [], 'zero', ...
    'u0', [p0 x0(i,:) x0(inxt,:)]');
end
prob = coco_add_pars(prob, 'x1', ix(1), 'x1');
prob = coco_add_pars(prob, 'pars', ip, {'a' 'b'});

%% Add glue between parameters, as well as x(i+1) and y(i)
for i=2:period
  uidx = coco_get_func_data(prob, ['henon' num2str(i)], 'uidx');
  prob = coco_add_glue(prob, ['pglue' num2str(i)], ip, uidx(ip));
end

for i=1:period
  uidx = coco_get_func_data(prob, ['henon' num2str(i)], 'uidx');
  inxt = mod(i, period)+1;
  uidxnxt = coco_get_func_data(prob, ['henon' num2str(inxt)], 'uidx');
  prob = coco_add_glue(prob, ['yglue' num2str(i)], uidx(iy), uidxnxt(ix));
end

%% Continue along family of period-4 orbits under variations in a.
bd1 = coco(prob, 'run1', [], 1, {'a' 'x1'}, [0.8 1.2]);

%% Switch branch at branch point
lab = coco_bd_labs(bd1, 'BP');

% reconstruct continuation problem structure
chart = coco_read_solution('run1', lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
    
prob = coco_prob;
for i=1:period
  fid = ['henon' num2str(i)];
  [chart, uidx] = coco_read_solution(fid, 'run1', lab, 'chart', 'uidx');
  inxt = mod(i, period)+1;
  prob = coco_add_func(prob, fid, fcn(f), [], 'zero', ...
    'u0', chart.x, 't0', cdata.v(uidx));
end
prob = coco_add_pars(prob, 'x1', ix(1), 'x1');
prob = coco_add_pars(prob, 'pars', ip, {'a' 'b'});

for i=2:period
  uidx = coco_get_func_data(prob, ['henon' num2str(i)], 'uidx');
  prob = coco_add_glue(prob, ['pglue' num2str(i)], ip, uidx(ip));
end

for i=1:period
  uidx = coco_get_func_data(prob, ['henon' num2str(i)], 'uidx');
  inxt = mod(i, period)+1;
  uidxnxt = coco_get_func_data(prob, ['henon' num2str(inxt)], 'uidx');
  prob = coco_add_glue(prob, ['yglue' num2str(i)], uidx(iy), uidxnxt(ix));
end

% continue along secondary branch
bd2 = coco(prob, 'run2', [], 1, {'a' 'x1'}, [0.8 1.2]);

%% Visualize the result
figure(1); clf; hold on
thm = struct();
thm.special = {'FP'};
thm.FP      = {'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 10};
coco_plot_bd(thm, 'run1', 'a', 'x1')
thm = struct();
thm.ylab  = 'x_1';
thm.lspec =  {'k:'  'LineWidth'  2};
coco_plot_bd(thm, 'run2', 'a', 'x1')
grid on
