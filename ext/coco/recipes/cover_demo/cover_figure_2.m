function cover_figure_2

p = fileparts(mfilename('fullpath'));
f = fullfile(p, 'cache_klein.mat');
if exist(f, 'file')==2
  load(f);
else
  if ~coco_exist('klein', 'run')    
    prob = coco_prob();
    prob = coco_set(prob, 'cont', 'atlas', @atlas2_x.create);
    prob = coco_set(prob, 'cont', 'PtMX', 20000);
    prob = coco_set(prob, 'cont', 'almax', 15);
    prob = coco_set(prob, 'cont', 'h', 0.1);
    prob = coco_set(prob, 'cont', 'NPR', 100);
    prob = coco_add_func(prob, 'klein', @klein, [], 'zero', ...
      'u0', [0; -1-sqrt(2); 0]);
    prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
    coco(prob, 'klein', [], 2, {'x' 'y' 'z'});
    
    coco_use_recipes_toolbox
  end
  atlas = coco_bd_read('klein', 'atlas');
  [tri, X, C] = get_charts(atlas, 1, 2, 3);
  N = 12000;
  M = max(C);
  I = N+find(rand(1,M-N)<=(1-linspace(0,nthroot(0.95,3),M-N).^3));
  r = arrayfun(@(x) (x<N) || any(x==I), C);
  save(f, 'tri', 'X', 'C', 'r');
  
  % save(f, 'tri', 'X', 'C');
end

figure(1)
clf
hold on
view([-175 50])

cmap = autumn(256);
colormap(cmap(32:end,:));
trisurf(tri(r,:), -X(:,1), X(:,3), X(:,2), X(:,2), 'LineStyle', 'none');
axis equal
axis off

camproj perspective
campos([-2.8029   33.5060   35.8342])
camva(7.5)
camtarget([0.0919   -0.7667   -0.2964])

lighting flat
set(gca, 'AmbientLightColor', 0.75*[1 1 1])
material dull

light('Position', [ 2  2 4], 'Color', 0.15*[1 1 1]);
light('Position', [ 1  2 4], 'Color', 0.25*[1 1 1]);
light('Position', [ 0  2 4], 'Color', 0.33*[1 1 1]);
light('Position', [-1  2 4], 'Color', 0.33*[1 1 1]);
light('Position', [-2  2 4], 'Color', 0.25*[1 1 1]);
light('Position', [-2  2 2], 'Color', 0.15*[1 1 1]);
light('Position', [-2  0 2], 'Color', 0.10*[1 1 1]);
light('Position', [-2 -2 2], 'Color', 0.07*[1 1 1]);

hold off

end

function [tri, X, C] = get_charts(charts, ix, iy, iz)
if nargin>=4
  f = @(x) x([ix iy iz])';
else
  f = @(x) x';
end
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
  for l=1:NN-1
    x   = f(xx+min(chart.R,chart.v(l))*(chart.TS*chart.s(:,l)));
    X   = [ X ; x ];
    tri = [tri ; xo xo+l xo+l+1];
    if chart.ep_flag<=1
      C = [ C ; k   ];
    else
      C = [ C ; 0 ];
    end
  end
  x   = f(xx+min(chart.R,chart.v(NN))*(chart.TS*chart.s(:,NN)));
  X   = [ X ; x ];
  tri = [tri ; xo xo+NN xo+1];
  if chart.ep_flag<=1
    C = [ C ; k   ];
  else
    C = [ C ; 0 ];
  end
end
end
