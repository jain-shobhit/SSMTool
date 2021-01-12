function figure_14_1
% Figure 14.1: Applying the modified 2-dimensional covering algorithm from
% Sect. 14.2 to the constrained cylinder in Example 14.4 results in the
% covering shown in different views in panels (a) and (b). A triangulation
% of the surface can be constructed using the neighbor information stored
% with each chart. A resulting triangulated surface is shown in panel (c)
% in the same view as in panel (b).

% Generate data
if ~coco_exist('cylinder', 'run')
  run demo_atlas2d_v5
end

% Extract data
atlas = coco_bd_read('cylinder', 'atlas');

% Plot data: panel (a)
figure(1)
clf
grid on
hold on

plot_charts(atlas.charts, false);

camproj('perspective')
view([135 10])
axis([-0.03 2.03 -1.03 1.03 -1.5 1.5], 'equal')
hold off

% Plot data: panel (b)
figure(2)
clf
grid on
hold on

plot_charts(atlas.charts, false);

camproj('perspective')
view([-45 5])
axis([-0.03 2.03 -1.03 1.03 -1.5 1.5], 'equal')
hold off

% Plot data: panel (c)
figure(3)
clf
grid on
hold on

x0 = [0 0 0];
n  = [1;1;-2];
cfunc = @(x) -n'*(x-x0)';
plot_trisurf(atlas.charts, cfunc)

camproj('perspective')
view([-45 5])
axis([-0.03 2.03 -1.03 1.03 -1.5 1.5], 'equal')
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

function plot_trisurf(charts, cfunc)
f   = @(x) x([1 2 3])';
tri = [];
X   = [];
C   = [];
N   = numel(charts);
for k=1:N
  chart = charts{k};
  X     = [ X ; f(chart.x) ];
  C     = [ C  ; cfunc(f(chart.x)) ];
  ic    = [chart.nb chart.nb(1)];
  ix    = chart.id;
  for l=1:numel(ic)-1
    face = sort([ix ic(l) ic(l+1)]);
    if all(face>0) && (isempty(tri) ||  ~ismember(face, tri, 'rows'))
      tri = [tri ; face];
    end
  end
end

cmap = repmat(linspace(0.6,0.9,100)', 1, 3);
colormap(cmap);
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'FaceColor', 'interp', ...
  'EdgeColor', 0.5*[1 1 1]);
end
