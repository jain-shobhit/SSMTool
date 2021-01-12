function cover_figure_1

% Generate data
if ~coco_exist('pillow', 'run')
  coco_use_recipes_toolbox atlas2d_v4
  
  prob = coco_prob();
  prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
  prob = coco_set(prob, 'cont', 'PtMX', 1000);
  prob = coco_set(prob, 'cont', 'almax', 30);
  prob = coco_set(prob, 'cont', 'h', 0.1);
  prob = coco_set(prob, 'cont', 'NPR', 100);
  prob = coco_add_func(prob, 'pillow', @pillow, [], 'zero', ...
    'u0', [1; 0; 0]);
  prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
  coco(prob, 'pillow', [], 2, {'x' 'y' 'z'}, ...
    {[], [-0.5 0.5], [-0.35 0.35]});
  
  coco_use_recipes_toolbox
end

% Extract data
atlas = coco_bd_read('pillow', 'atlas');

% Plot data
figure(1)
clf
hold on
axis([0 2 -0.5 0.5 -0.35 0.35]*1.04)
axis off
view([90 0])

cmap = [0.2 0.2 0.9; 0.15 0.15 0.85];
ecol = [0.1 0.1 0.8];
plot_charts(atlas.charts, cmap, ecol, 0.75)

hold off

end

function plot_charts(charts, cmap, ecol, lw)
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
  X  = [ X ; x ]; %#ok<*AGROW>
  xo = size(X,1);
  NN = size(chart.s,2);
  if isfield(chart, 'v')
    chart.v = chart.v;
  else
    chart.v = chart.R*ones(size(chart.s,2),1);
  end
  nb = circshift(chart.nb,[0 -1]);
  for l=1:NN-1;
    if chart.R<=chart.v(l)
      rr = scale*chart.R;
    else
      rr = chart.v(l);
    end
    x   = f(xx+rr*(chart.TS*chart.s(:,l)));
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

colormap(cmap);
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'EdgeColor', ecol, ...
'LineWidth', lw);

end
