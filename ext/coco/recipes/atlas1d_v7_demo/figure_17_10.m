function figure_17_10
% Figure 17.10: Continuation of the analysis from Fig. 17.8. Using the
% algorithm for branch switching at period-doubling points developed in
% Sect. 17.3.3, we can further complete the bifurcation diagram (a). Panel
% (b) shows a zoom into the closed family, along which period-doubling
% points were detected previously. A continuation of the emerging
% period-doubled orbits results in another closed family, along which
% period-doubling points are again detected (c). Repetition of this
% procedure results in evidence of a sequence of period-doubling
% bifurcations. The onset of this sequence is shown in the enlargement (d),
% which includes local families of period-2, -4, and -8 orbits; see also
% Fig. 17.11. Note that the families of period-doubled solutions were
% computed with the atlas algorithm without automatic branch switching to
% prevent redundant coverage.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run') && coco_exist('run4', 'run') ...
    && coco_exist('run5', 'run'))
  run demo_atlas1d_v7
end

% Extract data
bd  = coco_bd_read('run2');
bd1 = coco_bd_read('run3');
bd2 = coco_bd_read('run4');
bd3 = coco_bd_read('run5');

% Plot data: panel (a)-(c)
alims = {[-1 81 -2 6], [20 45 -1.9 4.1], [26 42 1.4 4]};

for i=1:3
  figure(i)
  clf
  hold on
  box on
  grid on
  axis(alims{i})
  
  if i==1
    plot_bd(bd)
  else
    plot_sub_bd(bd)
  end
  plot_bd(bd1)
  plot_bd(bd2)
  plot_bd(bd3)
  
  hold off
end

% Plot data: panel (d)
figure(4)
clf
hold on
grid on
box on
axis([41.15 41.45 2.15 2.23])

plot_sub_bd(bd, {'Marker', '.', 'MarkerSize', 15})
plot_bd(bd1, {'LineStyle', 'none', 'LineWidth', sqrt(2), ...
  'Color', 'black'}, {'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerFaceColor', 'white'})
plot_bd(bd2, {'LineStyle', 'none', 'LineWidth', sqrt(2), ...
  'Color', 'black'}, {'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerFaceColor', 'white'})
plot_bd(bd3, {'LineStyle', 'none', 'LineWidth', sqrt(2), ...
  'Color', 'black'}, {'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerFaceColor', 'white'})
plot_bd(bd1, {'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black'}, {'Marker', '.', 'MarkerSize', 15})
plot_bd(bd2, {'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black'}, {'Marker', '.', 'MarkerSize', 15})
plot_bd(bd3, {'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black'}, {'Marker', '.', 'MarkerSize', 15})

hold off

end


function plot_bd(bd, style, marker)
A  = coco_bd_col(bd, 'A');
x0 = coco_bd_col(bd, 'X0');
st = coco_bd_col(bd,'hspo.test.stab');

hold on
if nargin>=2
  plot(A,x0(1,:),style{:})
end
if nargin<3
  marker = {'Marker', '.', 'MarkerSize', 12};
end

plot(A(st==0),x0(1,st==0), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12, marker{:});
plot(A(st>0),x0(1,st>0), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12, marker{:})
idx = coco_bd_idxs(bd, 'EP');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white');
idx = coco_bd_idxs(bd, 'FP');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white');
idx = coco_bd_idxs(bd, 'BP');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white');
idx = coco_bd_idxs(bd, 'PD');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white');
end

function plot_sub_bd(bd, marker)
if nargin<2
  marker = {'Marker', '.', 'MarkerSize', 12};
end

labs = coco_bd_col(bd, 'LAB');
pts  = coco_bd_col(bd, 'PT');

i1=find(cellfun(@(x) ~isempty(x)&&x==41, labs));
i2=find(cellfun(@(x) ~isempty(x)&&x==72, labs));
idx = [1 1+(i1:i2)];
plot_bd(bd(idx,:), {'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black'}, marker);

i1=find(cellfun(@(x) ~isempty(x)&&x==73, labs));
i2=find(cellfun(@(x) ~isempty(x)&&x==140, labs));
idx = [1 1+(i1:i2)];
plot_bd(bd(idx,:), {'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black'}, marker);

i1=find(cellfun(@(x) ~isempty(x)&&x==140, labs));
i2=find(cellfun(@(x) ~isempty(x)&&x==153, labs));
idx = [1 1+(i1:i2)];
plot_bd(bd(idx,:), {'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black'}, marker);


i1=find(cellfun(@(x) ~isempty(x)&&x==388, labs));
i2=find(cellfun(@(x) ~isempty(x)&&x==404, labs));
while (pts(i1+1)-pts(i1))==1
  i1 = i1+1;
end
while (pts(i2+1)-pts(i2))==1
  i2 = i2+1;
end
idx = [1 1+(i1+1:i2)];
plot_bd(bd(idx,:), {'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black'}, marker);

i1=find(cellfun(@(x) ~isempty(x)&&x==404, labs));
i2=find(cellfun(@(x) ~isempty(x)&&x==469, labs));
while (pts(i1+1)-pts(i1))==1
  i1 = i1+1;
end
while (pts(i2+1)-pts(i2))==1
  i2 = i2+1;
end
idx = [1 1+(i1+1:i2)];
plot_bd(bd(idx,:), {'LineStyle', '-', 'LineWidth', 1, ...
  'Color', 'black'}, marker);

end
