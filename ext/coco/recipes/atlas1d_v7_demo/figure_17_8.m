function figure_17_8
% Figure 17.8: Amplitude response curve of the Duffing oscillator from
% Example 17.4 for omega=1 computed with the atlas code developed in Sect.
% 17.3, which implements automatic branch switching at branch points.
% Again, lambda=0.2 and alpha=epsilon=1. In addition to the results shown
% in Fig. 17.6, we obtain a number of families of nonsymmetric oscillations
% emerging from the branch points located during continuation; see also
% Fig. 17.9. Three of these form closed families, and the others terminate
% at the computational boundary. Along the closed family in the center of
% the figure, a number of toolbox events are located, including
% period-doubling bifurcation points. The analysis of this example
% continues in Fig. 17.10. The overlap of two curves at the period-doubling
% point 409 is an artifact of projection.

% Note that machine-dependent round-off errors may result in different
% solution labels for special points.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run') && coco_exist('run4', 'run') ...
    && coco_exist('run5', 'run'))
  run demo_atlas1d_v7
end

% Extract data
bd = coco_bd_read('run2');
A  = coco_bd_col(bd, 'A');
x0 = coco_bd_col(bd, 'X0');
st = coco_bd_col(bd,'hspo.test.stab');

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([-1 81 -2 6])

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
