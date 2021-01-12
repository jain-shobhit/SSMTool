function figure_13_4
% Figure 13.4: After executing 1,000 steps of the algorithm from Sect.
% 13.3, we obtain the partial cover of the manifold of Example 13.2 shown
% in panel (a). The local geometry of this partial cover is visualized in
% the close-up of the top-left corner, shown in panel (b).

% Generate data
if ~coco_exist('pillow', 'run')
  run demo_pillow
end

% Plot data: panel (a)
figure(1)
clf
hold on
grid on

atlas = coco_bd_read('pillow', 'atlas');
plot_charts(atlas.charts, false);

view([50 24])
camproj('perspective')
axis([0 1.25 -1.25 1.25 -1.25 1.25], 'equal')
hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on

plot_charts(atlas.charts, false)

camproj('perspective')
campos([10 -15 5])
camva(2.5)
camtarget([0.4 0.15 0.8])
axis([0 1.25 -1.25 1.25 -1.25 1.25], 'equal')
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
    plot3(x+xoff,y,z);
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
