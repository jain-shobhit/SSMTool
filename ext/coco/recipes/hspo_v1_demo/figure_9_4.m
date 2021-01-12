function figure_9_4
% Figure 9.4: Selected members of a family of nonimpacting two-segment
% periodic orbits of the stick-slip oscillator considered in Example 9.3 in
% the coordinates defined on page 237. Orbit 1 is close to impact and orbit
% 2 is nonphysical, because it penetrates a surface of impact.
% Consequently, between 1 and 2 there must exist a grazing orbit; see Fig.
% 9.5. The labels correspond to the session output included in the text.

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
  axis([0 0.8 -0.3 0.3])
  
  plot([0.5 0.5], [-0.3 0.3], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black');
  
  sol = msbvp_read_solution('', 'stickslip1', labs(i)); % Extract solution
  plot(sol{1}.x(:,1), sol{1}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 12)
  plot(sol{2}.x(:,1), sol{2}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black')
  
  hold off
end

coco_use_recipes_toolbox

end
