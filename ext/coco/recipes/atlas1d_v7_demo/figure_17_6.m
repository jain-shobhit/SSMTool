function figure_17_6
% Figure 17.6: Amplitude response curve of the Duffing oscillator from
% Example 17.4 for omega=1 computed with the atlas code developed in Sect.
% 16.1. Again, lambda=0.2 and alpha=epsilon=1. Along this curve we indicate
% the location of atlas and toolbox ecents. In addition, we encode
% stability information provided by the toolbox monitor function. Here,
% black corresponds to stable and gray to unstable oscillations. We observe
% two resonance peaks, where passage through each peak corresponds to the
% addition of an oscillation during one forcing half cycle, as illustrated
% with sample solutions at labels 1, 60 and 230 in Fig. 17.7. Note that all
% solutions along this response curve are symmetric. The analysis of this
% example continues in Fig. 17.8.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run') && coco_exist('run4', 'run') ...
    && coco_exist('run5', 'run'))
  run demo_atlas1d_v7
end

% Extract data
bd = coco_bd_read('run1');
A  = coco_bd_col(bd, 'A');
x0 = coco_bd_col(bd, 'X0');
st = coco_bd_col(bd,'hspo.test.stab');

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([-1 81 0 6])

plot(A(st==0),x0(1,st==0), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(A(st>0),x0(1,st>0), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd, 'EP');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd, 'FP');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd, 'BP');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd, 'PD');
plot(A(idx),x0(1,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
