function figure_13_3
% Figure 13.3: Initial steps of covering the manifold of Example 13.2 with
% the 2-dimensional atlas algorithm from Sect. 13.3. After four steps we
% obtain an atlas with one interior and three boundary charts (a), where
% each chart has all other charts as neighbors. The piecewise-linear
% approximation of the atlas boundary is the union of the bold line
% segments. The algorithm continues to grow the atlas and expands the
% boundary outwards (b).

% Generate data
if ~coco_exist({'pillow' 'run1'}, 'run')
  coco_use_recipes_toolbox atlas2d_v4
  prob = coco_prob();
  prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
  prob = coco_set(prob, 'cont', 'PtMX', 3);
  prob = coco_set(prob, 'cont', 'almax', 30);
  prob = coco_set(prob, 'cont', 'h', 0.15);
  prob = coco_set(prob, 'cont', 'NPR', 100);
  prob = coco_add_func(prob, 'pillow', @pillow, [], 'zero', ...
    'u0', [1; 0; 0]);
  prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
  coco(prob, {'pillow' 'run1'}, [], 2, {'x' 'y' 'z'});
  coco_use_recipes_toolbox
end
if ~coco_exist({'pillow' 'run2'}, 'run')
  coco_use_recipes_toolbox atlas2d_v4
  prob = coco_prob();
  prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
  prob = coco_set(prob, 'cont', 'PtMX', 20);
  prob = coco_set(prob, 'cont', 'almax', 30);
  prob = coco_set(prob, 'cont', 'h', 0.15);
  prob = coco_set(prob, 'cont', 'NPR', 100);
  prob = coco_add_func(prob, 'pillow', @pillow, [], 'zero', ...
    'u0', [1; 0; 0]);
  prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
  coco(prob, {'pillow' 'run2'}, [], 2, {'x' 'y' 'z'});
  coco_use_recipes_toolbox
end

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 2 -0.22 0.17 -0.2 0.2])
view([90 0])

atlas = coco_bd_read({'pillow' 'run1'}, 'atlas');
plot_charts(atlas.charts, true)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 2 -0.4 0.4 -0.4 0.4])
view([90 0])

atlas = coco_bd_read({'pillow' 'run2'}, 'atlas');
plot_charts(atlas.charts, true)

hold off

end

function plot_charts(charts, id_flag)
f   = @(x) x([1 2 3])';
is_boundary = @(x) any(x.nb==0);
scale = 1;
boundary = [];
tri = [];
X   = [];
C   = [];
N   = numel(charts);
for k=1:N
  chart = charts{k};
  xx = chart.x;
  x  = f(xx);
  X  = [ X ; x ];
  xo = size(X,1);
  NN = size(chart.s,2);
  if isfield(chart, 'v')
    chart.v = scale*chart.v;
  else
    chart.v = (scale*chart.R)*ones(size(chart.s,2),1);
  end
  nb = circshift(chart.nb,[0 -1]);
  for l=1:NN-1;
    x   = f(xx+min(chart.R,chart.v(l))*(chart.TS*chart.s(:,l)));
    X   = [ X ; x ];
    tri = [tri ; xo xo+l xo+l+1];
    C   = [ C ; is_boundary(chart)+1 ];
    if nb(l)==0
      boundary = [boundary ; [xo+l xo+l+1]];
    end
  end
  x   = f(xx+min(chart.R,chart.v(NN))*(chart.TS*chart.s(:,NN)));
  X   = [ X ; x ];
  tri = [tri ; xo xo+NN xo+1];
  C   = [ C ; is_boundary(chart)+1 ];
  if nb(NN)==0
    boundary = [boundary ; [xo+NN xo+1]];
  end
end

cmap = [0.9 0.9 0.9 ; 0.8 0.8 0.8];
colormap(cmap);
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'EdgeColor', 0.6*[1 1 1], 'LineWidth', 0.5);
if id_flag
  xoff = 0.005;
  for k=1:size(boundary,1)
    x = X(boundary(k,:),1);
    y = X(boundary(k,:),2);
    z = X(boundary(k,:),3);
    plot3(x+xoff,y,z, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black');
  end
  for k=1:numel(charts)
    text(charts{k}.x(1)+xoff, charts{k}.x(2), charts{k}.x(3), ...
      sprintf('%d', charts{k}.id), 'FontSize', 22, ...
      'FontName', 'Helvetica', 'FontWeight', 'bold', ...
      'Interpreter', 'none', 'Color', 'black', ...
      'BackgroundColor', 'none', 'EdgeColor', 'none', ...
      'LineStyle', '-', 'LineWidth', 1, 'Margin', 0.5, ...
      'HorizontalAlignment', 'center', ...
      'VerticalAlignment', 'middle', ...
      'BackgroundColor', cmap(is_boundary(charts{k})+1,:));
  end
end

end
