function figure_16_3
% Figure 16.3: For oscillations obtained at both sides of the main
% resonance peak of the frequency response curve shown in Fig. 16.2, we
% observe a characteristic phase shift between the forcing y_4(t) and the
% response y_1(t), here shown over the normalized time tau=t/T. To the
% left, both oscillations are in phase (a), while to the right we observe
% an antiphase response (c). The transition occurs close to the fold point
% at label 35 (b).

% Generate data
if ~coco_exist('duffing', 'run')
  run demo_duffing
end

coco_use_recipes_toolbox po_v1 coll_v1

% Plot data: panels (a)-(c)
labs = [16, 35, 58];
alim = {[0 1 -2.08 2.08], [0 1 -4.16 4.16], [0 1 -1.04 1.04]};

for i=1:3
  figure(i)
  clf
  hold on
  box on
  grid on
  axis(alim{i})
  
  sol = po_read_solution('', 'duffing', labs(i)); % Extract solution
  sol.t = sol.t/sol.t(end);
  plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  plot(sol.t, sol.x(:,4), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  
  hold off
end

coco_use_recipes_toolbox

end