function figure_9_3
% Figure 9.3: Families of multisegment periodic orbits from Example 9.2.
% The initial guess is shown as the dashed closed curve in panel (a). The
% initial correction step converges to orbit 1. Selected members of the
% family of orbits resulting from a continuation in the parameter beta are
% shown in (a). We restart a continuation in alpha from orbit 9,
% highlighted in gray. Some members of this family are shown in (b). The
% gray orbit is identical in both panels. The labels correspond to the
% session output included in the text.

% Generate data
if coco_exist('pwlin1', 'run') && coco_exist('pwlin2', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_pwlin

coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v1

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-1.5 4 -3 5])

plot([0 0], [-3 5], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.5 0.5 0.5])

bd   = coco_bd_read('pwlin1'); % Extract bifurcation data
labs = coco_bd_labs(bd);       % Extract labels
for lab=labs
  [sol data] = msbvp_read_solution('', 'pwlin1', lab); % Extract solution
  for i=1:data.nsegs
    if lab==9
      plot(sol{i}.x(:,1), sol{i}.x(:,2), 'LineStyle', '-', ...
        'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'Marker', '.', ...
        'MarkerSize', 12)
    elseif lab==6
    else
      plot(sol{i}.x(:,1), sol{i}.x(:,2),'LineStyle', '-', ...
        'LineWidth', 2, 'Color', 'black', 'Marker', '.', ...
        'MarkerSize', 12)
    end
  end
end
plot(x1(:,1), x1(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.9 0.9 0.9]) %#ok<NODEF>
plot(x2(:,1), x2(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.9 0.9 0.9]) %#ok<NODEF>
plot(x1(:,1), x1(:,2), 'LineStyle', '--', 'LineWidth', 2, ...
  'Color', 'black')
plot(x2(:,1), x2(:,2), 'LineStyle', '--', 'LineWidth', 2, ...
  'Color', 'black')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-1.5 4 -3 5])

plot([0 0], [-3 5], 'LineStyle', '-', 'LineWidth', 1, ...
  'Color', [0.3 0.3 0.3])

bd   = coco_bd_read('pwlin2'); % Extract bifurcation data
labs = coco_bd_labs(bd);       % Extract labels
for lab=labs
  [sol data] = msbvp_read_solution('', 'pwlin2', lab); % Extract solution
  for i=1:data.nsegs
    if abs(sol{i}.p(1)-1)<10*eps
      plot(sol{i}.x(:,1), sol{i}.x(:,2), 'LineStyle', '-', ...
        'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'Marker', '.', ...
        'MarkerSize', 12)
    else
      plot(sol{i}.x(:,1), sol{i}.x(:,2),'LineStyle', '-', ...
        'LineWidth', 2, 'Color', 'black', 'Marker', '.', ...
        'MarkerSize', 12)
    end
  end
end

hold off

coco_use_recipes_toolbox

end
