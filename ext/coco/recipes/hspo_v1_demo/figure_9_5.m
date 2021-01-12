function figure_9_5
% Figure 9.5: Selected members of the family of impacting three-segment
% periodic orbits of the stick-slip oscillator considered in Example 9.4 in
% the coordinates defiend on page 237. We initialize a computation of this
% family by resegmenting the grazing orbit from Example 9.3 and performing
% an initial correction under the additional constraint that the new
% segment has length 0; that is, it corresponds to the point of grazing
% contact. We again observe nonphysical orbit penetrating the impact
% surface. The labels correspond to the session output included in the
% text.

% Generate data
if ~(coco_exist('stickslip1', 'run') && coco_exist('stickslip2', 'run'))
  run demo_stickslip
end

coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v1

% Plot data: panels (a)-(c)
labs = [4, 1, 2];

for i=1:3
  figure(i)
  clf
  hold on
  grid on
  box on
  axis([0 0.8 -0.4 0.4])
  
  plot([0.5 0.5], [-0.4 0.4], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black')
  
  sol = msbvp_read_solution('', 'stickslip2', labs(i)); % Extract solution
  plot(sol{1}.x(:,1), sol{1}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 12)
  plot(sol{2}.x(:,1), sol{2}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black')
  plot(sol{3}.x(:,2), sol{3}.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black')
  
  hold off
end

coco_use_recipes_toolbox

end
