function figure_8_2
% Figure 8.2: Solutions to the Bratu boundary-value problem given by Eq.
% (8.11) and the given boundary conditions at t=0 and t=1. We start with
% the exact solution z(t)==0 at p=0 and obtain a family of solutions with
% increasing amplitudes. The labels correspond to the session output
% included in the text.

% Generate data
if coco_exist('brat1', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_bratu

coco_use_recipes_toolbox coll_v1 bvp_v1 

% Extract data
bd   = coco_bd_read('brat1');
labs = coco_bd_labs(bd);

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 1 -0.1 1.5])

for lab=labs
  sol = bvp_read_solution('', 'brat1', lab); % Extract solution
  plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 1 -5 5])

for lab=labs
  sol = bvp_read_solution('', 'brat1', lab); % Extract solution
  plot(sol.t, sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

coco_use_recipes_toolbox

end
