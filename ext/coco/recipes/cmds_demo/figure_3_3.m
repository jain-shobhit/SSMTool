function figure_3_3
% Figure 3.3: A set of profiles representative of the family of solutions
% to the Brusselator boundary-value problem in Example 3.14. The boundary
% conditions are f(0)=f(1)=g(1)=1 and g(0)=g_0, where g_0 is the active
% continuation parameter. We start with the exact solution f(x)==g(x)==1
% for g_0=1 and obtain solutions with nontrivial shape as the parameter g_0
% varies. The labels correspond to the session output included in the text.

% Generate data
if ~coco_exist('brusselator', 'run')
  run demo_brusselator
end

% Extract data
bd = coco_bd_read('brusselator');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 1 0.8 1.05])

for lab=1:8
  [data sol] = coco_read_solution('finitediff', 'brusselator', lab); % Extract solution
  t = (0:numel(data.f_idx)-1)/(numel(data.f_idx)-1);
  plot(t, sol.x(data.f_idx), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 15)
end

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 1 -0.5 7.5])

for lab=1:8
  [data sol] = coco_read_solution('finitediff', 'brusselator', lab); % Extract solution
  t = (0:numel(data.g_idx)-1)/(numel(data.g_idx)-1);
  plot(t, sol.x(data.g_idx), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 15)
end

hold off

end
