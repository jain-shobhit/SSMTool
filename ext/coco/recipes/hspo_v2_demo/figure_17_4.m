function figure_17_4
% Figure 17.4: Frequency response curve of the Duffing oscillator from
% Example 17.3 for forcing amplitude A=26 under variations in the forcing
% frequency omega, while the remaining parameters are set to lambda=0.2 and
% alpha=epsilon=1. Here, black dots correspond to periodic orbits with all
% Floquet multipliers within the unit circle, whereas gray dots correspond
% to periodic orbits with at least one Floquet multiplier outside of the
% unit circle. A number of toolbox and atlas events are detected along this
% curve. As a special case, the branch point at label 12 is a pitchfork
% bifurcation point that marks a characteristic transition from
% nonsymmetric responses through a symmetric one to symmetry-conjugate
% responses. This transition is illustrated with sample solutions at labels
% 9, 12, and 16 in Fig. 17.5

% Generate data
if ~(coco_exist('run1', 'run'))
  run demo_hspo_v2
end

% Extract data
bd1 = coco_bd_read('run1');

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([0.84 1.1 -2 4])

om = coco_bd_col(bd1, 'om');
x0 = coco_bd_col(bd1, 'X0');
st = coco_bd_col(bd1,'hspo.test.stab');
plot(om, x0(1,:), 'LineStyle', '-', 'LineWidth', 1, 'Color', 'black')
plot(om(st==0), x0(1,st==0), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(om(st>0), x0(1,st>0), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd1, 'EP');
plot(om(idx), x0(1,idx),'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black','Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'FP');
plot(om(idx), x0(1,idx),'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black','Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'BP');
plot(om(idx), x0(1,idx),'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black','Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'PD');
plot(om(idx), x0(1,idx),'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black','Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
