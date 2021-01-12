function figure_15_5
% Figure 15.5: Some orbits representative of the continuation illustrated
% in Fig. 15.4. Continuation starts with the impacting orbit shown in (a)
% and passes through a grazing point (b), after which nonphysical orbits
% result (c).

% Generate data
if ~(coco_exist('impact1', 'run') && coco_exist('impact2', 'run'))
  run demo_impact
end

coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v1

% Plot data: panels (a)-(c)
labs = [1, 6, 11];

for i=1:3
  figure(i)
  clf
  hold on
  grid on
  box on
  axis([-1.5 2 -2 2])
  
  sol = msbvp_read_solution('', 'impact1', labs(i)); % Extract trajectory segments
  plot(sol{1}.x(:,1), sol{1}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  plot(sol{2}.x(:,1), sol{2}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
  plot([1 1], [-2 2], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black')
  
  hold off
end

coco_use_recipes_toolbox

end
