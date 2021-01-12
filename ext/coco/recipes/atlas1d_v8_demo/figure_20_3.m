function figure_20_3
% Figure 20.3: Interpolation of the increasing test function
% f(t,p)=tanh(pt)/tanh(p) using an equidistributed moving mesh with seven
% mesh points for increasing values of p, as described in Sect. 20.1.3. The
% circles at the bottom indicate the mesh in t. The quality of
% interpolation of the steep front improves compared with the graphs shown
% in Fig. 20.2. The result of using an equidistributed moving mesh with
% variable number of mesh points is shown in Fig. 20.4.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run'))
  run demo_atlas1d_v8
end

% Extract data
bd = coco_bd_read('run2');
labs = coco_bd_labs(bd, 'UZ');

% Plot data: panels (a)-(d)
for i=1:numel(labs)
  figure(i)
  clf
  hold on
  grid on
  box on
  
  [soldata sol] = coco_read_solution('tanh','run2',labs(i));
  plot(-1:.01:1,tanh(sol.x(soldata.p_idx)*(-1:.01:1))/tanh(sol.x(soldata.p_idx)), ...
    'LineStyle', '-', 'LineWidth', 2, 'Color', [0.4 0.4 0.4])
  plot(soldata.t, sol.x(soldata.x_idx), 'LineStyle', '-', ...
    'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 15)
  plot(soldata.t, 0*soldata.t-1.2, 'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, ...
    'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
  
  axis('tight')
  xlims = get(gca, 'XLim');
  ylims = get(gca, 'YLim');
  axis([xlims(1)-0.005*(xlims(2)-xlims(1)) ...
    xlims(2)+0.005*(xlims(2)-xlims(1)) ...
    ylims(1)-0.01*(ylims(2)-ylims(1)) ylims(2)+0.01*(ylims(2)-ylims(1))])
  
  hold off
end

end
