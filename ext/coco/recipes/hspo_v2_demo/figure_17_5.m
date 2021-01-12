function figure_17_5
% Figure 17.5: The Duffing equation with bang-bang forcing in Example 17.3
% is equivariant under rotations of the (y_1,y_2) plane by pi, combined
% with a phase shift by pi. As a consequence, we expect to observe
% symmetric as well as nonsymmetric oscillations, where existence of a
% nonsymmetric oscillation implies existence of the symmetry-conjugate
% oscillation. In the present case, continuation starts at a family of
% nonsymmetric oscillations, passes through a symmetry-increasing
% bifurcation point, and continues to produce symmetry-conjugate
% oscillations, as shown in the sequence of plots in panels (a)-(c). The
% solutions in panels (a) and (c) are obtained for the same value of omega,
% but on different sides of the branch point. Clearly, one can transform
% (a) into (c) by a rotation by pi and exchanging black for gray.

% Generate data
if ~(coco_exist('run1', 'run'))
  run demo_hspo_v2
end

coco_use_recipes_toolbox coll_v1 msbvp_v1 hspo_v2

% Plot data: panels (a)-(c)
labs  = [9 12 16];

for i=1:3
  figure(i)
  clf
  hold on
  box on
  grid on
  axis([-5.5 5.5 -16.5 16.5])
  
  sol = msbvp_read_solution('', 'run1', labs(i)); % Extract trajectory segments
  plot(sol{1}.x(:,1), sol{1}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  plot(sol{2}.x(:,1),sol{2}.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
  
  hold off
end

coco_use_recipes_toolbox

end
