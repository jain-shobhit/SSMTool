function figure_14_2
% Figure 14.2: The 2-dimensional atlas algorithm discussed at the end of
% Sect. 14.2 supports starting on the boundary of the computational domain.
% This is illustrated with the manifold of Example 14.5, where we start on
% a boundary defined by one active constraint in panel (a) and at a corner
% defined by two active constraints in panel (b). In both cases, the
% initial solution guess is given by (x,y,z)=(1,-1,0) and we show the atlas
% obtained after 25 continuation steps.

% Generate data
if ~coco_exist('cylinder1', 'run') || ~coco_exist('cylinder2', 'run')
  run demo_cylinder
end

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0.75 1.35 -1.5 -0.5 -0.1 0.35])
view([0.01 0])

atlas = coco_bd_read('cylinder1', 'atlas');
plot_charts(atlas.charts, true);

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0.6 1.1 -1.5 -0.5 -0.1 0.4])
view([0.01 0])

atlas = coco_bd_read('cylinder2', 'atlas');
plot_charts(atlas.charts, true);

hold off

end

function plot_charts(charts, id_flag)
f   = @(x) x([1 2 3])';
is_boundary = @(x) any(x.nb==0) && (x.ep_flag==0 || x.id==26);
scale = 1;
boundary = [];
tri = [];
X   = [];
CC  = [];
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
    if nb(l)==0 && (chart.ep_flag==0 || chart.id==26)
      boundary = [boundary ; [xo+l xo+l+1]];
    end
  end
  x   = f(xx+min(chart.R,chart.v(NN))*(chart.TS*chart.s(:,NN)));
  X   = [ X ; x ];
  tri = [tri ; xo xo+NN xo+1];
  C   = [ C ; is_boundary(chart)+1 ];
  if nb(NN)==0 && (chart.ep_flag==0 || chart.id==26)
    boundary = [boundary ; [xo+NN xo+1]];
  end
  CC = [CC; mean(X(xo:end,:),1)];
end

cmap = [0.9 0.9 0.9 ; 0.8 0.8 0.8];
colormap(cmap);
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'EdgeColor', 0.6*[1 1 1], ...
  'LineWidth', 0.5, 'FaceAlpha', 0.99);
if id_flag
  yoff = 0.1;
  for k=1:size(boundary,1)
    x = X(boundary(k,:),1);
    y = X(boundary(k,:),2);
    z = X(boundary(k,:),3);
    plot3(x,y-yoff,z, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black');
  end
  for k=1:numel(charts)
    text(CC(k,1), CC(k,2)-yoff, CC(k,3), ...
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
