function figure_17_7
% Figure 17.7: After each passage through a resonance peak along the
% amplitude response curve shown in Fig. 17.6, we observe the addition of
% an oscillation during one half cycle of the forcing.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run') && coco_exist('run4', 'run') ...
    && coco_exist('run5', 'run'))
  run demo_atlas1d_v7
end

coco_use_recipes_toolbox coll_v1 msbvp_v1 hspo_v2

labs  = [1 60 230];
limits = {
  [-1.5 1.5 -1.5 1.5]
  [-5 5 -10 10]
  [-7 7 -25 25]
  };

for i=1:3
  figure(i)
  clf
  hold on
  box on
  grid on
  
  axis(limits{i})
  
  sol = msbvp_read_solution('', 'run1', labs(i)); % Extract trajectory segments
  plot(sol{1}.x(:,1),sol{1}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  plot(sol{2}.x(:,1),sol{2}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)

  hold off
end

coco_use_recipes_toolbox

end
