function figure_20_19
% Figure 20.19: The logarithm of the condition number of the Jacobian of
% the restricted continuation problem plotted against the point number
% normalized by the maximum point number during continuation of the canard
% family of the dynamical system given by the vector field in Eq. (20.45)
% using the (co)moving-mesh adaptive discretization strategies in Sects.
% 18.3.1, 20.2.1, and 20.2.2, respectively.

% Note that machine-dependent round-off errors may result in differences in
% condition numbers.

% Generate data
if ~(coco_exist('3', 'run') && coco_exist('4', 'run') ...
    && coco_exist('5', 'run'))
  run demo_canard
end

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([0 1 0 20])

bd = coco_bd_read('5');
pt = coco_bd_col(bd, 'PT');
pt = pt/pt(end);
y = log10(coco_bd_col(bd, 'lsol.cond'));
plot(pt, y, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])

bd = coco_bd_read('3');
pt = coco_bd_col(bd, 'PT');
pt = pt/pt(end);
y = log10(coco_bd_col(bd, 'lsol.cond'));
plot(pt, y, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black')

bd = coco_bd_read('4');
pt = coco_bd_col(bd, 'PT');
pt = pt/pt(end);
y = log10(coco_bd_col(bd, 'lsol.cond'));
plot(pt, y, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'black')

hold off

end
