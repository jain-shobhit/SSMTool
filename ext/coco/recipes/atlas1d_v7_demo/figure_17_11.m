function figure_17_11
% Figure 17.11: A sequence of period-doubled orbits along the families
% shown in Fig. 17.10(d).

if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run') && coco_exist('run4', 'run') ...
    && coco_exist('run5', 'run'))
  run demo_atlas1d_v7
end

coco_use_recipes_toolbox coll_v1 msbvp_v1 hspo_v2

runs  = {'run3' 'run4' 'run5'};
labs  = [18 16 11];

for i=1:3
  figure(i)
  clf
  hold on
  box on
  grid on
  axis([-7 6 -25 25])
  
  sol = msbvp_read_solution('', runs{i}, labs(i)); % Extract trajectory segments
  for k=1:4:numel(sol)
    plot(sol{k}.x(:,1),sol{k}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
    plot(sol{k+1}.x(:,1),sol{k+1}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.3 0.3 0.3], 'Marker', '.', 'MarkerSize', 12)
    plot(sol{k+2}.x(:,1),sol{k+2}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 12)
    plot(sol{k+3}.x(:,1),sol{k+3}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7], 'Marker', '.', 'MarkerSize', 12)
  end
end

coco_use_recipes_toolbox

end
